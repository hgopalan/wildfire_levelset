#!/usr/bin/env python3
"""
two_way_coupling_example.py - Minimal Viable Product for two-way wind-fire coupling

Demonstrates the complete MVP workflow:
1. Fire solver generates surface fluxes
2. Fluxes are passed to wind solver
3. Wind solver recomputes 3D wind with heat forcing
4. Wind is read back into fire solver

Usage:
    PYTHONPATH=build/python python3 src/python/two_way_coupling_example.py inputs.i
"""

import sys
import os
import numpy as np
from typing import Optional

# Import wildfire and coupling modules
try:
    from wildfire_solver import WildfireSolver
except ImportError as e:
    print(f"ERROR: Could not import wildfire_solver: {e}")
    print("Make sure to set PYTHONPATH to include build/python")
    sys.exit(1)

# Import coupling utilities
try:
    from surface_flux_extractor import SurfaceFluxExtractor
    from wind_solver_interface import SyntheticWindSolver
except ImportError as e:
    print(f"ERROR: Could not import coupling modules: {e}")
    sys.exit(1)


class TwoWayCoupledSimulation:
    """
    Two-way coupled wind-fire simulation controller.

    Manages the interaction between fire and wind solvers:
    - Extract fire surface fluxes
    - Pass to wind solver
    - Compute wind with heat forcing
    - Update fire wind field
    """

    def __init__(self, fire_inputs: str, wind_solver: Optional[object] = None):
        """
        Initialize coupled simulation.

        Parameters:
            fire_inputs: Path to fire solver inputs file
            wind_solver: Wind solver instance (uses SyntheticWindSolver if None)
        """
        # Initialize fire solver
        print("Initializing fire solver...")
        self.fire = WildfireSolver(fire_inputs)

        # Initialize wind solver
        if wind_solver is None:
            print("Initializing wind solver (synthetic)...")
            self.wind = SyntheticWindSolver()
            self.wind.initialize()
        else:
            self.wind = wind_solver

        # Check domain compatibility
        wind_info = self.wind.get_domain_info()
        if (self.fire.nx != wind_info['nx'] or
            self.fire.ny != wind_info['ny']):
            print(f"WARNING: Fire grid ({self.fire.ny}×{self.fire.nx}) != "
                  f"Wind grid ({wind_info['ny']}×{wind_info['nx']})")
            print("         Assuming grids are aligned; results may be inaccurate")

        # Initialize flux extractor
        self.flux_extractor = SurfaceFluxExtractor(
            (self.fire.ny, self.fire.nx),
            self.fire.dx, self.fire.dy
        )

        # Coupling parameters
        self.coupling_enabled = True
        self.heat_scaling = 1.0  # Feedback strength: 1.0 = full, 0.0 = none
        self.diagnostics = {
            'step': [],
            'time': [],
            'total_heat_MW': [],
            'max_intensity': [],
            'burned_area': [],
        }

    def step(self) -> bool:
        """
        Execute one coupled iteration.

        Returns:
            True if step succeeded, False otherwise
        """
        try:
            # Step 1: Get current fire state
            state = self.fire.get_state()
            current_time = state['time']

            # Step 2: Extract surface fluxes from fire
            if self.coupling_enabled:
                flux_dict = self.flux_extractor.get_flux_dict(state)

                # Apply heat scaling (for sensitivity studies)
                flux_dict['heat_flux'] *= self.heat_scaling
                flux_dict['buoyancy'] *= self.heat_scaling

                # Step 3: Pass fluxes to wind solver and recompute wind
                success = self.wind.solve(current_time, heat_source=flux_dict)
                if not success:
                    print(f"ERROR: Wind solver failed at t={current_time:.2f}s")
                    return False

                # Step 4: Extract new wind field and update fire solver
                u_3d, v_3d, w_3d = self.wind.get_velocity_arrays()
                wind_info = self.wind.get_domain_info()

                # Update fire with new wind (3D)
                self.fire.update_wind_3d(
                    u_3d, v_3d, w_3d,
                    wind_info['nz'],
                    wind_info['zmin'],
                    wind_info['zmax']
                )

                # Diagnostic: print heat feedback info
                if flux_dict['total_heat_MW'] > 0.01:
                    print(f"    Heat feedback: {flux_dict['total_heat_MW']:.2f} MW, "
                          f"Max intensity: {flux_dict['intensity'].max():.1f} kW/m")

            # Step 5: Advance fire simulation
            result = self.fire.step()
            if not result['success']:
                print(f"ERROR: Fire solver failed at step {self.fire.step_num}")
                return False

            # Step 6: Collect diagnostics
            state = self.fire.get_state()
            phi = state['phi']
            intensity = state['intensity']
            burned_area = np.sum(phi <= 0.0) * self.fire.dx * self.fire.dy

            self.diagnostics['step'].append(self.fire.step_num)
            self.diagnostics['time'].append(current_time)
            self.diagnostics['total_heat_MW'].append(
                self.flux_extractor.compute_total_heat_release(state)
            )
            self.diagnostics['max_intensity'].append(intensity.max())
            self.diagnostics['burned_area'].append(burned_area)

            return True

        except Exception as e:
            print(f"ERROR in coupled step: {e}")
            import traceback
            traceback.print_exc()
            return False

    def run(self, num_steps: int = 20, plot_interval: int = 5) -> bool:
        """
        Run coupled simulation for specified number of steps.

        Parameters:
            num_steps: Number of timesteps to run
            plot_interval: Write plotfile every N steps

        Returns:
            True if simulation completed successfully
        """
        print(f"\n{'='*70}")
        print(f"Two-Way Coupled Wind-Fire Simulation")
        print(f"{'='*70}")
        print(f"Heat coupling: {'ENABLED' if self.coupling_enabled else 'DISABLED'}")
        print(f"Heat scaling factor: {self.heat_scaling:.2f}")
        print(f"Target steps: {num_steps}\n")

        for step in range(num_steps):
            success = self.step()
            if not success:
                print(f"✗ Simulation stopped at step {step+1}")
                return False

            # Print progress
            state = self.fire.get_state()
            diag = self.diagnostics

            if (step + 1) % 3 == 0 or step == 0:
                print(f"  Step {step+1:3d}/{num_steps}: "
                      f"t={state['time']:7.2f} s, "
                      f"burned={diag['burned_area'][-1]:8.1f} m², "
                      f"max_intensity={diag['max_intensity'][-1]:6.1f} kW/m, "
                      f"heat={diag['total_heat_MW'][-1]:6.2f} MW")

            # Write plotfile
            if (step + 1) % plot_interval == 0:
                plotfile = f"plt_coupled_{(step+1)//plot_interval:03d}"
                self.fire.write_plotfile(plotfile)
                print(f"           → Wrote {plotfile}")

        print(f"\n{'='*70}")
        print(f"✓ Simulation complete!")
        print(f"{'='*70}\n")
        return True

    def print_summary(self):
        """Print summary statistics."""
        diag = self.diagnostics

        if len(diag['time']) == 0:
            print("No diagnostic data available")
            return

        print(f"\nSimulation Summary")
        print(f"{'-'*70}")
        print(f"Final time: {diag['time'][-1]:.2f} s")
        print(f"Total steps: {diag['step'][-1]}")
        print(f"Final burned area: {diag['burned_area'][-1]:.1f} m²")
        print(f"Max fire intensity: {max(diag['max_intensity']):.1f} kW/m")
        print(f"Peak heat release: {max(diag['total_heat_MW']):.2f} MW")

        # Compute averages
        avg_intensity = np.mean(diag['max_intensity'])
        avg_heat = np.mean(diag['total_heat_MW'])
        print(f"Avg fire intensity: {avg_intensity:.1f} kW/m")
        print(f"Avg heat release: {avg_heat:.2f} MW")

    def finalize(self):
        """Clean up solvers."""
        self.fire.finalize()
        self.wind.finalize()
        print("\n✓ Solvers finalized")


def main():
    """Main entry point."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Two-way coupled wind-fire simulation (MVP)"
    )
    parser.add_argument("inputs_file", help="Fire solver inputs file")
    parser.add_argument("-n", "--num-steps", type=int, default=20,
                        help="Number of timesteps (default: 20)")
    parser.add_argument("-p", "--plot-interval", type=int, default=5,
                        help="Write plotfile every N steps (default: 5)")
    parser.add_argument("-H", "--heat-scale", type=float, default=1.0,
                        help="Heat feedback scaling factor (default: 1.0)")
    parser.add_argument("--no-coupling", action="store_true",
                        help="Disable heat feedback (one-way coupling)")

    args = parser.parse_args()

    # Check inputs file exists
    if not os.path.exists(args.inputs_file):
        print(f"ERROR: Inputs file not found: {args.inputs_file}")
        return 1

    try:
        # Create coupled simulation
        sim = TwoWayCoupledSimulation(args.inputs_file)

        # Configure coupling
        sim.coupling_enabled = not args.no_coupling
        sim.heat_scaling = args.heat_scale

        # Run simulation
        success = sim.run(num_steps=args.num_steps, plot_interval=args.plot_interval)

        # Print summary
        if success:
            sim.print_summary()
            sim.finalize()
            return 0
        else:
            sim.finalize()
            return 1

    except Exception as e:
        print(f"FATAL ERROR: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
