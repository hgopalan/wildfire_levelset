#!/usr/bin/env python3
"""
example_one_way_coupling.py - Demonstration of one-way wind-to-fire coupling

This example shows how to couple the wildfire_levelset fire solver with an external
wind solver using one-way coupling, where wind affects fire but fire does not
affect wind.

This is the simplest coupling mode and is suitable for:
- Pre-computed wind fields from weather models
- Coupling with atmospheric simulators
- Scenarios where fire-induced wind effects are negligible
- Fast runtime simulations

The workflow is:
1. Solve the wind field (external wind solver)
2. Extract 3D wind velocity
3. Update fire solver with wind field
4. Advance fire simulation
5. Repeat steps 1-4

Fire heating does NOT affect the wind field in this mode.

Usage:
    python3 example_one_way_coupling.py fire_inputs.i wind_inputs.i
    
    Or create a simulation script directly:
    
    from wildfire_solver import WildfireSolver
    import numpy as np
    
    # Initialize fire solver
    fire = WildfireSolver("fire_inputs.i")
    
    # Time loop
    final_time = 3600.0  # 1 hour
    while fire.time < final_time:
        # 1. Get wind field from external source (or synthetic)
        u_3d, v_3d, w_3d = get_wind_field(fire.time)
        
        # 2. Update fire with wind
        fire.update_wind_3d(u_3d, v_3d, w_3d, nz, zmin, zmax)
        
        # 3. Advance fire
        fire.step()
    
    fire.finalize()
"""

import numpy as np
import sys
import argparse


def generate_synthetic_wind_field(nx, ny, nz, 
                                  xmin, xmax, ymin, ymax, zmin, zmax,
                                  time=0.0,
                                  base_u=5.0, base_v=2.0):
    """
    Generate a synthetic 3D wind field for testing.
    
    This creates a simple wind field with:
    - Constant horizontal wind in x and y directions
    - Optional vertical component for convection
    - Height-dependent wind profile (log-law)
    
    Parameters:
        nx, ny, nz (int): Grid dimensions
        xmin, xmax, ymin, ymax, zmin, zmax (float): Domain bounds (meters)
        time (float): Current simulation time (for wind variations)
        base_u (float): Base u-wind component (m/s)
        base_v (float): Base v-wind component (m/s)
    
    Returns:
        tuple: (u_3d, v_3d, w_3d) with shape (nz, ny, nx) each
    """
    # Create meshgrid
    z_levels = np.linspace(zmin, zmax, nz)
    dz = (zmax - zmin) / nz
    
    # Initialize 3D arrays
    u_3d = np.zeros((nz, ny, nx), dtype=np.float64)
    v_3d = np.zeros((nz, ny, nx), dtype=np.float64)
    w_3d = np.zeros((nz, ny, nx), dtype=np.float64)
    
    # Log-law wind profile: u(z) = (u*/κ) * ln(z/z0)
    u_star = 0.4  # friction velocity
    kappa = 0.4   # von Kármán constant
    z0 = 0.1      # roughness length
    
    for k in range(nz):
        z = zmin + (k + 0.5) * dz
        
        # Log-law profile
        height_factor = max(np.log(max(z, z0) / z0), 0.1)
        u_profile = (u_star / kappa) * height_factor
        v_profile = (u_star / kappa) * height_factor * 0.4  # Reduced v-wind
        
        # Set wind field
        u_3d[k, :, :] = base_u * u_profile / np.log(100.0 / z0)  # Normalize to 10m
        v_3d[k, :, :] = base_v * v_profile / np.log(100.0 / z0)
        
        # Vertical wind (weak, for convection)
        w_3d[k, :, :] = 0.0
    
    return u_3d, v_3d, w_3d


def run_one_way_coupling(fire_inputs, fire_only=True, wind_inputs=None):
    """
    Run one-way wind-fire coupling simulation.
    
    Parameters:
        fire_inputs (str): Path to fire solver inputs file
        fire_only (bool): If True, use synthetic wind; if False, couple with wind solver
        wind_inputs (str): Path to wind solver inputs file (if fire_only=False)
    
    Returns:
        dict: Simulation results
    """
    import os
    
    # Check if fire inputs exist
    if not os.path.exists(fire_inputs):
        raise FileNotFoundError(f"Fire inputs file not found: {fire_inputs}")
    
    try:
        from wildfire_solver import WildfireSolver
    except ImportError as e:
        raise ImportError(
            "Could not import WildfireSolver. "
            "Make sure wildfire_levelset is built with Python bindings enabled.\n"
            f"Error: {e}"
        )
    
    # Initialize fire solver
    print("="*70)
    print("ONE-WAY COUPLING: Wind → Fire (Fire does not affect wind)")
    print("="*70)
    print()
    
    fire = WildfireSolver(fire_inputs)
    print(f"✓ Fire solver initialized")
    print(f"  Grid: {fire.nx} × {fire.ny}")
    print(f"  Domain: X=[{fire.xmin:.1f}, {fire.xmax:.1f}], "
          f"Y=[{fire.ymin:.1f}, {fire.ymax:.1f}]")
    print(f"  Cell size: {fire.dx:.2f} × {fire.dy:.2f} m")
    print()
    
    # Initialize wind solver if needed
    wind = None
    if not fire_only and wind_inputs:
        if not os.path.exists(wind_inputs):
            print(f"⚠️  Wind inputs not found: {wind_inputs}")
            print("   Using synthetic wind instead")
        else:
            try:
                from wind_solver import WindSolver
                wind = WindSolver(wind_inputs)
                print(f"✓ Wind solver initialized")
                print(f"  Grid: {wind.nx} × {wind.ny} × {wind.nz}")
                print()
            except Exception as e:
                print(f"⚠️  Could not initialize wind solver: {e}")
                print("   Using synthetic wind instead")
    
    # Simulation parameters
    final_time = 3600.0  # 1 hour in seconds
    plot_interval = 600.0  # Write plots every 10 minutes
    next_plot_time = fire.time + plot_interval
    plot_num = 0
    
    print(f"Parameters:")
    print(f"  Final time: {final_time:.0f} s ({final_time/60:.1f} minutes)")
    print(f"  Plot interval: {plot_interval:.0f} s")
    if wind is not None:
        print(f"  Wind mode: External wind solver")
    else:
        print(f"  Wind mode: Synthetic wind field")
    print()
    
    # Time-stepping loop
    print(f"{'Step':>5} | {'Time (s)':>10} | {'Time (min)':>10} | "
          f"{'Burned Cells':>12} | {'Burned Area (km²)':>16}")
    print("-" * 70)
    
    step = 0
    try:
        while fire.time < final_time:
            # Get wind field
            if wind is not None:
                # Solve wind and extract velocity
                wind.solve()
                vel = wind.get_velocity()
                u_3d = vel['u']
                v_3d = vel['v']
                w_3d = vel['w']
            else:
                # Use synthetic wind
                u_3d, v_3d, w_3d = generate_synthetic_wind_field(
                    fire.nx, fire.ny, 8,  # Assume 8 vertical levels
                    fire.xmin, fire.xmax, fire.ymin, fire.ymax, 0.0, 100.0,
                    time=fire.time,
                    base_u=5.0, base_v=1.0
                )
            
            # Update fire with wind
            fire.update_wind_3d(u_3d, v_3d, w_3d, nz=8, zmin=0.0, zmax=100.0)
            
            # Advance fire
            fire.step()
            
            # Get current state
            state = fire.get_state()
            phi = state['phi']
            burned_cells = np.sum(phi <= 0)
            domain_area = (fire.xmax - fire.xmin) * (fire.ymax - fire.ymin)
            burned_area = burned_cells * fire.dx * fire.dy / 1e6  # km²
            
            # Write plotfile if needed
            if fire.time >= next_plot_time:
                plotfile = f"plt_one_way_{plot_num:05d}"
                fire.write_plotfile(plotfile)
                next_plot_time += plot_interval
                plot_num += 1
            
            # Print progress
            if step % 50 == 0:
                print(f"{step:5d} | {fire.time:10.1f} | {fire.time/60:10.1f} | "
                      f"{burned_cells:12d} | {burned_area:16.3f}")
            
            step += 1
    
    except Exception as e:
        print(f"\n✗ Simulation failed: {e}")
        raise
    
    print("-" * 70)
    print(f"{step:5d} | {fire.time:10.1f} | {fire.time/60:10.1f} | "
          f"{burned_cells:12d} | {burned_area:16.3f}")
    print()
    print(f"✓ Simulation completed")
    print(f"  Total steps: {step}")
    print(f"  Final time: {fire.time:.1f} s ({fire.time/60:.1f} minutes)")
    print(f"  Burned area: {burned_area:.3f} km²")
    print(f"  Plotfiles written: {plot_num}")
    
    # Cleanup
    fire.finalize()
    if wind is not None:
        wind.finalize()
    
    return {
        'final_time': fire.time,
        'final_step': step,
        'burned_cells': burned_cells,
        'burned_area': burned_area,
        'plotfiles': plot_num
    }


def main():
    """Command-line interface for one-way coupling example."""
    parser = argparse.ArgumentParser(
        description='Run one-way wind-fire coupling example',
        epilog='Example: python3 example_one_way_coupling.py regtest/surface_spread/farsite_ellipse/inputs.i'
    )
    parser.add_argument('fire_inputs', 
                       help='Path to fire solver inputs file')
    parser.add_argument('--wind-inputs', 
                       help='Path to wind solver inputs file (optional)')
    parser.add_argument('--synthetic-wind', action='store_true',
                       help='Use synthetic wind field instead of wind solver')
    
    args = parser.parse_args()
    
    # Determine if we should use wind solver or synthetic wind
    fire_only = args.synthetic_wind or (args.wind_inputs is None)
    
    # Run simulation
    try:
        results = run_one_way_coupling(
            args.fire_inputs,
            fire_only=fire_only,
            wind_inputs=args.wind_inputs
        )
        print(f"\n✓ Results: {results}")
        return 0
    except Exception as e:
        print(f"\n✗ Error: {e}", file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
