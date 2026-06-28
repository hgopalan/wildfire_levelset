#!/usr/bin/env python3
"""
massconsistent_fire_coupling.py - Two-way coupled fire-wind simulation

Implements bi-directional coupling between wildfire_levelset fire solver
and massconsistent_amr wind solver for accurate fire-atmosphere interaction.

FEATURES:
=========
- Two-way coupling: Wind affects fire, fire heat affects wind
- Rothermel rate of spread (ROS) calculations with wind interaction
- Surface heat flux extraction and feedback to wind solver
- Domain synchronization and validation
- Comprehensive diagnostics and error handling

USAGE:
======
    from massconsistent_fire_coupling import FireWindCoupler
    
    coupler = FireWindCoupler(
        fire_inputs="fire_inputs.i",
        wind_inputs="wind_inputs.txt",
        max_time=3600.0  # 1 hour
    )
    
    results = coupler.run()
    coupler.finalize()

REQUIREMENTS:
=============
- wildfire_levelset built with: -DLEVELSET_BUILD_PYTHON_BINDINGS=ON
- massconsistent_amr built with: -DBUILD_PYTHON_BINDINGS=ON
- PYTHONPATH set to include both: /path/to/wildfire/build/python:/path/to/wind/build/python
- NumPy

TWO-WAY COUPLING FLOW:
=====================
For each timestep:
    1. Wind Solve: Compute 3D wind field with fire heating from previous step
    2. Extract 3D Wind: u_3d, v_3d, w_3d from massconsistent_amr
    3. Update Fire: Pass 3D wind to fire solver
    4. Fire Step: Advance fire using Rothermel ROS
    5. Extract Heat: Get surface heat flux from fire solver
    6. Add to Wind: Pass heat back to wind solver for next iteration
    7. Output: Write diagnostics, statistics, plotfiles
"""

import numpy as np
import sys
import os
from pathlib import Path
from typing import Dict, Tuple, Optional


class FireWindCoupler:
    """
    Two-way coupled fire-wind simulation controller.
    
    Manages coupled iteration between Rothermel fire spread and 3D wind solver,
    including heat feedback from fire to atmosphere.
    """
    
    def __init__(
        self,
        fire_inputs: str,
        wind_inputs: Optional[str] = None,
        max_time: Optional[float] = None,
        max_steps: Optional[int] = None,
        wind_update_interval: int = 1,
        output_interval: float = 60.0,
        verbose: bool = True,
        use_synthetic_wind: bool = False
    ):
        """
        Initialize coupled fire-wind solver system.
        
        Parameters:
            fire_inputs (str): Path to fire solver inputs file
            wind_inputs (str, optional): Path to wind solver inputs file
            max_time (float, optional): Maximum simulation time (seconds)
            max_steps (int, optional): Maximum number of coupled steps
            wind_update_interval (int): Update wind every N fire steps
            output_interval (float): Write output every N seconds
            verbose (bool): Enable detailed logging
            use_synthetic_wind (bool): Use synthetic wind instead of massconsistent_amr
        
        Raises:
            ImportError: If required modules not available
            RuntimeError: If initialization fails
        """
        self.fire_inputs = fire_inputs
        self.wind_inputs = wind_inputs
        self.max_time = max_time
        self.max_steps = max_steps
        self.wind_update_interval = wind_update_interval
        self.output_interval = output_interval
        self.verbose = verbose
        self.use_synthetic_wind = use_synthetic_wind
        
        # Import fire solver
        try:
            from wildfire_solver import WildfireSolver
            self.WildfireSolver = WildfireSolver
        except ImportError as e:
            raise ImportError(
                "Could not import WildfireSolver. "
                "Build with: cmake -S . -B build -DLEVELSET_BUILD_PYTHON_BINDINGS=ON\n"
                f"Error: {e}"
            )
        
        # Try to import wind solver (optional for synthetic wind mode)
        self.wind = None
        if wind_inputs is not None and not use_synthetic_wind:
            try:
                # Try to import massconsistent_amr wind solver
                try:
                    import pyWindSolver as wind_module
                    self.wind = wind_module.WindSolver(wind_inputs)
                    self._log("✓ Wind solver (massconsistent_amr) initialized")
                except (ImportError, AttributeError):
                    # Fallback: try alternate import name
                    try:
                        from wind_solver import WindSolver
                        self.wind = WindSolver(wind_inputs)
                        self._log("✓ Wind solver initialized")
                    except ImportError:
                        self._log("⚠ Wind solver not available, using synthetic wind")
                        self.use_synthetic_wind = True
            except Exception as e:
                self._log(f"⚠ Warning: Could not initialize wind solver: {e}")
                self.use_synthetic_wind = True
        
        # Initialize fire solver
        self.fire = None
        self._initialize_fire()
        
        # Tracking variables
        self.coupled_step = 0
        self.next_output_time = self.fire.time + output_interval
        self.output_count = 0
        self.heat_accumulation = None
        
        # Statistics tracking
        self.stats_history = {
            'time': [],
            'burned_area': [],
            'max_ros': [],
            'max_intensity': [],
            'heat_released': []
        }
    
    def _log(self, msg: str):
        """Print log message if verbose enabled."""
        if self.verbose:
            print(msg)
    
    def _initialize_fire(self):
        """Initialize fire solver from inputs file."""
        try:
            self.fire = self.WildfireSolver(self.fire_inputs)
            self._log(f"✓ Fire solver initialized from {self.fire_inputs}")
            self._log(f"  Domain: X=[{self.fire.xmin:.1f}, {self.fire.xmax:.1f}], "
                     f"Y=[{self.fire.ymin:.1f}, {self.fire.ymax:.1f}]")
            self._log(f"  Grid: {self.fire.nx} × {self.fire.ny}, "
                     f"dx={self.fire.dx:.2f} m, dy={self.fire.dy:.2f} m")
        except Exception as e:
            raise RuntimeError(f"Failed to initialize fire solver: {e}")
    
    def _check_domain_compatibility(self):
        """Verify fire and wind solver domains are compatible."""
        if self.wind is None:
            return True
        
        tol = 1.0  # meter tolerance
        
        fire_bounds = (self.fire.xmin, self.fire.xmax, self.fire.ymin, self.fire.ymax)
        wind_bounds = (self.wind.xmin, self.wind.xmax, self.wind.ymin, self.wind.ymax)
        
        if any(abs(a - b) > tol for a, b in zip(fire_bounds, wind_bounds)):
            self._log("⚠ Warning: Fire and wind domain bounds differ:")
            self._log(f"  Fire: X=[{self.fire.xmin:.1f}, {self.fire.xmax:.1f}], "
                     f"Y=[{self.fire.ymin:.1f}, {self.fire.ymax:.1f}]")
            self._log(f"  Wind: X=[{self.wind.xmin:.1f}, {self.wind.xmax:.1f}], "
                     f"Y=[{self.wind.ymin:.1f}, {self.wind.ymax:.1f}]")
            return False
        
        self._log("✓ Fire and wind domains are compatible")
        return True
    
    def _generate_synthetic_wind_3d(self, time: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray, int, float, float]:
        """
        Generate synthetic 3D wind field for testing without massconsistent_amr.
        
        In a real coupling, this would be replaced with actual wind solver calls.
        
        Returns:
            (u_3d, v_3d, w_3d, nz, zmin, zmax)
        """
        # Grid dimensions
        nz = 10
        zmin = 0.0
        zmax = 1000.0
        
        # Create coordinate arrays
        x = np.linspace(self.fire.xmin, self.fire.xmax, self.fire.nx)
        y = np.linspace(self.fire.ymin, self.fire.ymax, self.fire.ny)
        z = np.linspace(zmin, zmax, nz)
        
        # Base wind (varies with time, simulating diurnal cycle)
        hour = (time / 3600.0) % 24
        diurnal_factor = 0.7 + 0.3 * np.sin(2 * np.pi * (hour - 6) / 24)
        
        u_ref = 5.0 * diurnal_factor  # m/s easterly
        v_ref = 1.0 * diurnal_factor  # m/s southerly
        z_ref = 10.0  # reference height
        z0 = 0.1  # roughness length
        
        # Initialize 3D arrays (Fortran order for compatibility)
        u_3d = np.zeros((nz, self.fire.ny, self.fire.nx), dtype=np.float64, order='F')
        v_3d = np.zeros((nz, self.fire.ny, self.fire.nx), dtype=np.float64, order='F')
        w_3d = np.zeros((nz, self.fire.ny, self.fire.nx), dtype=np.float64, order='F')
        
        # Create log-law profile with spatial variation
        for k in range(nz):
            z_height = z[k] if z[k] > z0 else z0
            height_factor = np.log(z_height / z0) / np.log(z_ref / z0)
            
            for j in range(self.fire.ny):
                for i in range(self.fire.nx):
                    # Spatial modulation (simulates terrain effects)
                    x_norm = (x[i] - self.fire.xmin) / (self.fire.xmax - self.fire.xmin)
                    y_norm = (y[j] - self.fire.ymin) / (self.fire.ymax - self.fire.ymin)
                    
                    spatial_mod = 1.0 + 0.2 * np.sin(2 * np.pi * x_norm) * np.cos(2 * np.pi * y_norm)
                    
                    u_3d[k, j, i] = u_ref * height_factor * spatial_mod
                    v_3d[k, j, i] = v_ref * height_factor * spatial_mod
                    w_3d[k, j, i] = 0.0  # No vertical velocity initially
        
        return u_3d, v_3d, w_3d, nz, zmin, zmax
    
    def _update_wind_in_fire_solver(self):
        """
        Get 3D wind and update fire solver.
        
        For massconsistent_amr: extract wind from solver
        For synthetic: generate test wind
        """
        if self.use_synthetic_wind:
            u_3d, v_3d, w_3d, nz, zmin, zmax = self._generate_synthetic_wind_3d(self.fire.time)
        else:
            if self.wind is None:
                return False
            
            try:
                # Try to extract 3D velocity from wind solver
                if hasattr(self.wind, 'get_velocity_arrays'):
                    u_3d, v_3d, w_3d = self.wind.get_velocity_arrays()
                    nz = self.wind.nz
                    zmin = self.wind.zmin
                    zmax = self.wind.zmax
                elif hasattr(self.wind, 'get_velocity'):
                    vel_dict = self.wind.get_velocity()
                    u_3d = vel_dict['u']
                    v_3d = vel_dict['v']
                    w_3d = vel_dict.get('w', np.zeros_like(u_3d))
                    nz = self.wind.nz
                    zmin = self.wind.zmin
                    zmax = self.wind.zmax
                else:
                    self._log("⚠ Wind solver doesn't have velocity extraction method")
                    return False
                
            except Exception as e:
                self._log(f"⚠ Failed to extract wind from solver: {e}")
                return False
        
        # Update fire solver with 3D wind
        try:
            self.fire.update_wind_3d(u_3d, v_3d, w_3d, nz, zmin, zmax)
            return True
        except Exception as e:
            self._log(f"✗ Failed to update fire wind: {e}")
            return False
    
    def _add_heat_to_wind_solver(self):
        """
        Extract fire heat flux and pass to wind solver.
        
        This implements two-way coupling feedback.
        """
        if self.wind is None:
            return False
        
        try:
            # Get surface heat flux from fire
            flux_data = self.fire.get_surface_fluxes()
            heat_flux = flux_data['heat_flux']
            
            # Try to add heat source to wind solver
            if hasattr(self.wind, 'add_heat_source'):
                self.wind.add_heat_source(heat_flux, flux_data['grid_info'])
                return True
            else:
                if not hasattr(self, '_heat_warning_shown'):
                    self._log("⚠ Wind solver doesn't support add_heat_source() yet")
                    self._heat_warning_shown = True
                return False
                
        except Exception as e:
            self._log(f"⚠ Failed to add heat to wind solver: {e}")
            return False
    
    def _coupled_step(self) -> Dict:
        """
        Execute one coupled timestep.
        
        Workflow:
        1. Update wind in fire solver (from massconsistent_amr or synthetic)
        2. Advance fire solver one timestep (Rothermel ROS calculation)
        3. Extract fire state and heat release
        4. Add heat to wind solver (for two-way coupling)
        5. Return step statistics
        
        Returns:
            dict: Step results with status and metrics
        """
        step_result = {
            'success': True,
            'coupled_step': self.coupled_step,
            'fire_time': self.fire.time,
            'wind_updated': False,
            'heat_added': False,
            'stats': {}
        }
        
        # Only update wind at specified interval
        if self.coupled_step % self.wind_update_interval == 0:
            step_result['wind_updated'] = self._update_wind_in_fire_solver()
        
        # Advance fire solver
        try:
            fire_step_result = self.fire.step()
            if not fire_step_result['success']:
                self._log(f"✗ Fire step failed at step {self.coupled_step}")
                step_result['success'] = False
                return step_result
        except Exception as e:
            self._log(f"✗ Fire step exception: {e}")
            step_result['success'] = False
            return step_result
        
        # Get fire state and statistics
        try:
            fire_stats = self.fire.get_statistics()
            step_result['stats'] = fire_stats
            
            # Add to history
            self.stats_history['time'].append(fire_stats['time'])
            self.stats_history['burned_area'].append(fire_stats['burned_area'])
            self.stats_history['max_ros'].append(fire_stats['max_ros'])
            self.stats_history['max_intensity'].append(fire_stats['max_intensity'])
            self.stats_history['heat_released'].append(fire_stats['total_heat_release'])
            
        except Exception as e:
            self._log(f"⚠ Failed to get fire statistics: {e}")
        
        # Add heat to wind solver (two-way coupling)
        if self.coupled_step % self.wind_update_interval == 0:
            step_result['heat_added'] = self._add_heat_to_wind_solver()
        
        return step_result
    
    def step(self) -> Dict:
        """
        Execute one coupled timestep and manage output.
        
        Returns:
            dict: Step results
        """
        result = self._coupled_step()
        self.coupled_step += 1
        
        # Handle output at specified intervals
        if self.fire.time >= self.next_output_time:
            self._write_output()
            self.next_output_time += self.output_interval
        
        return result
    
    def _write_output(self):
        """Write output files (plotfiles, diagnostics, etc.)."""
        try:
            plotfile_name = f"plt_{self.output_count:05d}"
            self.fire.write_plotfile(plotfile_name)
            
            if self.verbose:
                stats = self.fire.get_statistics()
                print(f"  Output at t={self.fire.time:.2f}s: {plotfile_name}")
                print(f"    Burned area: {stats['burned_area']/1e6:.2f} km²")
                print(f"    Max ROS: {stats['max_ros']:.2f} m/s")
                print(f"    Max intensity: {stats['max_intensity']:.1f} kW/m")
            
            self.output_count += 1
        except Exception as e:
            self._log(f"⚠ Failed to write output: {e}")
    
    def run(self) -> Dict:
        """
        Run the coupled simulation to completion.
        
        Returns:
            dict: Final statistics and summary
        """
        self._log("\n" + "="*70)
        self._log("FIRE-WIND COUPLED SIMULATION")
        self._log("="*70)
        self._log(f"Fire inputs: {self.fire_inputs}")
        if self.wind_inputs:
            self._log(f"Wind inputs: {self.wind_inputs}")
        self._log(f"Synthetic wind: {self.use_synthetic_wind}")
        
        if self.max_time:
            self._log(f"Maximum time: {self.max_time:.1f} s")
        if self.max_steps:
            self._log(f"Maximum steps: {self.max_steps}")
        self._log("="*70 + "\n")
        
        # Domain compatibility check
        self._check_domain_compatibility()
        
        # Main time loop
        self._log("Starting coupled time integration...\n")
        
        while True:
            # Check stopping conditions
            if self.max_time and self.fire.time >= self.max_time:
                self._log(f"\n✓ Reached max time: {self.fire.time:.2f} s")
                break
            
            if self.max_steps and self.coupled_step >= self.max_steps:
                self._log(f"\n✓ Reached max steps: {self.coupled_step}")
                break
            
            # Execute coupled step
            try:
                result = self.step()
                
                if not result['success']:
                    self._log(f"✗ Simulation failed at step {self.coupled_step}")
                    break
                
                # Print progress
                if self.coupled_step % 10 == 0:
                    stats = result['stats']
                    self._log(f"Step {self.coupled_step:6d}: t={self.fire.time:8.2f}s, "
                             f"ROS={stats['max_ros']:6.2f} m/s, "
                             f"Burned={stats['burned_area']/1e6:6.2f} km²")
                
            except KeyboardInterrupt:
                self._log("\n⚠ Simulation interrupted by user")
                break
            except Exception as e:
                self._log(f"\n✗ Unhandled exception: {e}")
                import traceback
                traceback.print_exc()
                break
        
        # Write final output and statistics
        self._write_output()
        self._log("\n" + "="*70)
        self._log("SIMULATION COMPLETE")
        self._log("="*70)
        
        return self._get_final_statistics()
    
    def _get_final_statistics(self) -> Dict:
        """Compute final statistics summary."""
        stats = self.fire.get_statistics()
        
        return {
            'success': True,
            'final_time': self.fire.time,
            'final_step': self.coupled_step,
            'burned_area': stats['burned_area'],
            'perimeter': stats['perimeter'],
            'max_ros': stats['max_ros'],
            'max_intensity': stats['max_intensity'],
            'num_burning_cells': stats['num_burning_cells'],
            'history': self.stats_history
        }
    
    def finalize(self):
        """Clean up and finalize solvers."""
        if self.fire is not None:
            self.fire.finalize()
        
        if self.wind is not None and hasattr(self.wind, 'finalize'):
            try:
                self.wind.finalize()
            except:
                pass
        
        self._log("✓ Solvers finalized")


def main():
    """Example usage and command-line interface."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Two-way coupled fire-wind simulation"
    )
    parser.add_argument("fire_inputs", help="Fire solver inputs file")
    parser.add_argument("--wind-inputs", default=None, help="Wind solver inputs file")
    parser.add_argument("--max-time", type=float, default=3600.0, help="Max sim time (s)")
    parser.add_argument("--max-steps", type=int, default=None, help="Max timesteps")
    parser.add_argument("--synthetic-wind", action="store_true", help="Use synthetic wind")
    parser.add_argument("--output-interval", type=float, default=60.0, help="Output interval (s)")
    
    args = parser.parse_args()
    
    try:
        coupler = FireWindCoupler(
            fire_inputs=args.fire_inputs,
            wind_inputs=args.wind_inputs,
            max_time=args.max_time,
            max_steps=args.max_steps,
            use_synthetic_wind=args.synthetic_wind,
            output_interval=args.output_interval,
            verbose=True
        )
        
        results = coupler.run()
        coupler.finalize()
        
        print("\nFinal Results:")
        print(f"  Burned area: {results['burned_area']/1e6:.2f} km²")
        print(f"  Max ROS: {results['max_ros']:.2f} m/s")
        print(f"  Perimeter: {results['perimeter']/1000:.2f} km")
        
        return 0
        
    except Exception as e:
        print(f"✗ Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
