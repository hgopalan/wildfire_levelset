#!/usr/bin/env python3
"""
example_two_way_coupling.py - Demonstration of two-way wind-fire coupling

This example shows how to couple the wildfire_levelset fire solver with an external
wind solver using two-way coupling, where wind affects fire AND fire affects wind
(via heat feedback).

This is a more advanced coupling mode suitable for:
- Fire-atmosphere interaction studies
- Understanding fire-induced wind effects (updrafts, flow deflection)
- Scenarios where fire heating creates significant atmospheric changes
- Research simulations (more computationally expensive)

The workflow is:
1. Solve the wind field with current heat source
2. Extract 3D wind velocity
3. Update fire solver with wind field
4. Advance fire simulation
5. Extract heat release from fire
6. Pass heat back to wind solver for next solve
7. Repeat steps 1-6

Fire heating DOES affect the wind field in the next iteration.

Usage:
    python3 example_two_way_coupling.py fire_inputs.i wind_inputs.i
    
    Or use the levelset_coupling module:
    
    from levelset_coupling import CoupledWindFireSimulation
    
    coupled = CoupledWindFireSimulation(
        wind_inputs="wind_inputs.i",
        fire_inputs="fire_inputs.i",
        coupling_mode='two_way'
    )
    results = coupled.run(final_time=3600.0, plot_interval=600.0)
    coupled.finalize()
"""

import numpy as np
import sys
import argparse


def run_two_way_coupling(wind_inputs, fire_inputs):
    """
    Run two-way wind-fire coupling simulation.
    
    This example demonstrates the full two-way coupling workflow where:
    1. Wind solver is called with heat source from fire
    2. Wind field is passed to fire solver
    3. Fire is advanced
    4. Heat release is extracted and stored
    5. Heat is added to wind solver for next iteration
    
    Parameters:
        wind_inputs (str): Path to wind solver inputs file
        fire_inputs (str): Path to fire solver inputs file
    
    Returns:
        dict: Simulation results
    """
    import os
    
    # Check if input files exist
    if not os.path.exists(fire_inputs):
        raise FileNotFoundError(f"Fire inputs file not found: {fire_inputs}")
    if not os.path.exists(wind_inputs):
        raise FileNotFoundError(f"Wind inputs file not found: {wind_inputs}")
    
    # Import solvers
    try:
        from wildfire_solver import WildfireSolver
    except ImportError as e:
        raise ImportError(
            "Could not import WildfireSolver. "
            "Make sure wildfire_levelset is built with Python bindings enabled.\n"
            f"Error: {e}"
        )
    
    try:
        from wind_solver import WindSolver
    except ImportError as e:
        raise ImportError(
            "Could not import WindSolver. "
            "Make sure massconsistent_amr is built with Python bindings enabled.\n"
            f"Error: {e}"
        )
    
    # Print header
    print("="*70)
    print("TWO-WAY COUPLING: Wind ↔ Fire (Fire affects wind via heat feedback)")
    print("="*70)
    print()
    
    # Initialize solvers
    print("Initializing solvers...")
    fire = WildfireSolver(fire_inputs)
    print(f"✓ Fire solver initialized")
    print(f"  Grid: {fire.nx} × {fire.ny}")
    print(f"  Domain: X=[{fire.xmin:.1f}, {fire.xmax:.1f}], "
          f"Y=[{fire.ymin:.1f}, {fire.ymax:.1f}]")
    print(f"  Cell size: {fire.dx:.2f} × {fire.dy:.2f} m")
    
    wind = WindSolver(wind_inputs)
    print(f"✓ Wind solver initialized")
    print(f"  Grid: {wind.nx} × {wind.ny} × {wind.nz}")
    print(f"  Domain: X=[{wind.xmin:.1f}, {wind.xmax:.1f}], "
          f"Y=[{wind.ymin:.1f}, {wind.ymax:.1f}]")
    print(f"  Cell size: {wind.dx:.2f} × {wind.dy:.2f} m")
    print()
    
    # Check domain compatibility
    tol = 1.0
    bounds_match = (
        abs(wind.xmin - fire.xmin) <= tol and
        abs(wind.xmax - fire.xmax) <= tol and
        abs(wind.ymin - fire.ymin) <= tol and
        abs(wind.ymax - fire.ymax) <= tol
    )
    spacing_match = (
        abs(wind.dx - fire.dx) <= tol and
        abs(wind.dy - fire.dy) <= tol
    )
    
    if bounds_match and spacing_match:
        print("✓ Wind and fire domains are compatible")
    else:
        if not bounds_match:
            print("⚠️  Warning: Domain bounds don't match perfectly")
        if not spacing_match:
            print("⚠️  Warning: Grid spacing doesn't match perfectly")
    print()
    
    # Simulation parameters
    final_time = 3600.0  # 1 hour in seconds
    plot_interval = 600.0  # Write plots every 10 minutes
    next_plot_time = fire.time + plot_interval
    plot_num = 0
    
    print(f"Parameters:")
    print(f"  Final time: {final_time:.0f} s ({final_time/60:.1f} minutes)")
    print(f"  Plot interval: {plot_interval:.0f} s")
    print(f"  Wind-fire coupling: two-way with heat feedback")
    print()
    
    # Time-stepping loop
    print(f"{'Step':>5} | {'Fire Time (s)':>12} | {'Wind Iters':>10} | "
          f"{'Burned Cells':>12} | {'Heat Release (kW)':>16}")
    print("-" * 75)
    
    step = 0
    total_wind_iters = 0
    total_heat_release = 0.0
    
    try:
        while fire.time < final_time:
            # TWO-WAY COUPLING STEP 1: Solve wind field with current heat source
            # (heat was added in previous step)
            wind.solve()
            wind_iters = wind.iters
            total_wind_iters += wind_iters
            
            # TWO-WAY COUPLING STEP 2: Extract 3D wind velocity
            vel = wind.get_velocity()
            u_3d = vel['u']  # shape (nz, ny, nx)
            v_3d = vel['v']
            w_3d = vel['w']
            
            # TWO-WAY COUPLING STEP 3: Update fire with wind
            fire.update_wind_3d(
                u_3d, v_3d, w_3d,
                wind.nz,
                wind.zmin,
                wind.zmax
            )
            
            # TWO-WAY COUPLING STEP 4: Advance fire simulation
            fire.step()
            
            # TWO-WAY COUPLING STEP 5: Extract heat release from fire
            fire_state = fire.get_state()
            heat_data = fire.compute_heat_release(fire_state)
            surface_flux = heat_data['surface_flux']  # kW/m²
            total_heat = np.sum(surface_flux) * fire.dx * fire.dy / 1000.0  # kW
            total_heat_release += total_heat
            
            # TWO-WAY COUPLING STEP 6: Pass heat back to wind solver
            # This heat will be used in the next wind.solve() call
            try:
                flux_data = fire.get_surface_fluxes()
                grid_info = {
                    'xmin': wind.xmin,
                    'xmax': wind.xmax,
                    'ymin': wind.ymin,
                    'ymax': wind.ymax,
                    'dx': wind.dx,
                    'dy': wind.dy
                }
                wind.add_heat_source(flux_data['heat_flux'], grid_info)
                heat_added = True
            except (AttributeError, Exception) as e:
                heat_added = False
                if step == 0:
                    print(f"⚠️  Note: Wind solver doesn't support heat feedback")
                    print(f"   Continuing with one-way coupling instead")
            
            # Get current state
            phi = fire_state['phi']
            burned_cells = np.sum(phi <= 0)
            
            # Write plotfile if needed
            if fire.time >= next_plot_time:
                plotfile = f"plt_two_way_{plot_num:05d}"
                fire.write_plotfile(plotfile)
                next_plot_time += plot_interval
                plot_num += 1
            
            # Print progress
            if step % 10 == 0:
                heat_indicator = "✓" if heat_added else "✗"
                print(f"{step:5d} | {fire.time:12.1f} | {wind_iters:10d} | "
                      f"{burned_cells:12d} | {total_heat:16.1f}")
            
            step += 1
    
    except Exception as e:
        print(f"\n✗ Simulation failed: {e}")
        raise
    
    print("-" * 75)
    print(f"{step:5d} | {fire.time:12.1f} | {total_wind_iters//step:10d} | "
          f"{burned_cells:12d} | {total_heat_release:16.1f}")
    print()
    print(f"✓ Simulation completed")
    print(f"  Total steps: {step}")
    print(f"  Final time: {fire.time:.1f} s ({fire.time/60:.1f} minutes)")
    print(f"  Burned cells: {burned_cells}")
    print(f"  Total heat released: {total_heat_release:.1f} kW")
    print(f"  Average wind iterations per step: {total_wind_iters/step:.1f}")
    print(f"  Plotfiles written: {plot_num}")
    
    # Cleanup
    fire.finalize()
    wind.finalize()
    
    return {
        'final_time': fire.time,
        'final_step': step,
        'burned_cells': burned_cells,
        'total_heat_release': total_heat_release,
        'total_wind_iters': total_wind_iters,
        'plotfiles': plot_num
    }


def main():
    """Command-line interface for two-way coupling example."""
    parser = argparse.ArgumentParser(
        description='Run two-way wind-fire coupling example',
        epilog='Example: python3 example_two_way_coupling.py wind_inputs.i fire_inputs.i'
    )
    parser.add_argument('wind_inputs', 
                       help='Path to wind solver inputs file')
    parser.add_argument('fire_inputs', 
                       help='Path to fire solver inputs file')
    
    args = parser.parse_args()
    
    # Run simulation
    try:
        results = run_two_way_coupling(args.wind_inputs, args.fire_inputs)
        print(f"\n✓ Results:")
        for key, value in results.items():
            print(f"  {key}: {value}")
        return 0
    except Exception as e:
        print(f"\n✗ Error: {e}", file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
