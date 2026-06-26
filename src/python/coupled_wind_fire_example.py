#!/usr/bin/env python3
"""
coupled_wind_fire_example.py - Demonstration of coupled wind-fire simulation

This script shows how to run a coupled wind-fire simulation where:
1. Wind solver runs for one timestep (simulated with synthetic data)
2. Wind data is passed to the fire solver
3. Fire solver runs for one timestep
4. Process repeats

This demonstrates the workflow that will be used when integrating
with massconsistent_amr or other external wind solvers.

Usage:
    PYTHONPATH=build/python python3 src/python/coupled_wind_fire_example.py inputs.i
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Import the wildfire solver wrapper
from wildfire_solver import WildfireSolver


def generate_synthetic_wind_3d(nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax, time):
    """
    Generate synthetic 3D wind field that varies with time.
    
    This simulates the output from massconsistent_amr or another 3D wind solver.
    In a real coupled simulation, this would be replaced with actual wind solver calls.
    
    Parameters:
        nx, ny, nz: Grid dimensions
        xmin, xmax, ymin, ymax, zmin, zmax: Domain bounds
        time: Current simulation time (seconds)
    
    Returns:
        u_3d, v_3d, w_3d: 3D wind arrays, shape (nz, ny, nx)
    """
    # Create coordinate arrays
    x = np.linspace(xmin, xmax, nx)
    y = np.linspace(ymin, ymax, ny)
    z = np.linspace(zmin, zmax, nz)
    
    # Create meshgrid
    X, Y = np.meshgrid(x, y, indexing='ij')
    
    # Initialize 3D arrays
    u_3d = np.zeros((nz, ny, nx))
    v_3d = np.zeros((nz, ny, nx))
    w_3d = np.zeros((nz, ny, nx))
    
    # Base wind speed (varies with time)
    # Simulates diurnal wind variation
    hour = (time / 3600.0) % 24
    diurnal_factor = 0.7 + 0.3 * np.sin(2 * np.pi * (hour - 6) / 24)
    
    u_ref = 5.0 * diurnal_factor  # m/s westerly
    v_ref = 1.0 * diurnal_factor  # m/s southerly
    z_ref = 10.0  # reference height
    z0 = 0.1  # roughness length
    
    # Create log-law profile with spatial variation
    for k in range(nz):
        z_height = z[k] if z[k] > z0 else z0
        
        # Log-law scaling
        height_factor = np.log(z_height / z0) / np.log(z_ref / z0)
        
        for j in range(ny):
            for i in range(nx):
                # Spatial modulation (simulates terrain effects)
                x_norm = (x[i] - xmin) / (xmax - xmin)
                y_norm = (y[j] - ymin) / (ymax - ymin)
                
                # Create channeling effect
                spatial_mod = 1.0 + 0.3 * np.sin(3 * np.pi * x_norm) * np.cos(2 * np.pi * y_norm)
                
                # Add time-varying vortex
                vortex_x = 0.5 + 0.1 * np.sin(time / 300.0)
                vortex_y = 0.5 + 0.1 * np.cos(time / 300.0)
                dx_vortex = x_norm - vortex_x
                dy_vortex = y_norm - vortex_y
                r_vortex = np.sqrt(dx_vortex**2 + dy_vortex**2)
                vortex_strength = 2.0 * np.exp(-20.0 * r_vortex**2)
                
                u_vortex = -vortex_strength * dy_vortex
                v_vortex = vortex_strength * dx_vortex
                
                u_3d[k, j, i] = u_ref * height_factor * spatial_mod + u_vortex
                v_3d[k, j, i] = v_ref * height_factor * (1.0 + 0.1 * np.cos(2 * np.pi * x_norm)) + v_vortex
                w_3d[k, j, i] = 0.1 * np.sin(np.pi * z_height / zmax)
    
    return u_3d, v_3d, w_3d


def visualize_coupled_simulation(states, wind_states, output_file="coupled_simulation.png"):
    """
    Create visualization of coupled wind-fire simulation results.
    
    Parameters:
        states: List of fire state dictionaries from each timestep
        wind_states: List of wind field dictionaries from each timestep
        output_file: Path to save the visualization
    """
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    
    # Plot initial, middle, and final states
    plot_times = [0, len(states) // 2, len(states) - 1]
    
    for idx, t_idx in enumerate(plot_times):
        state = states[t_idx]
        wind = wind_states[t_idx]
        
        time = state['time']
        phi = state['phi']
        ros = state['ros']
        u_wind = wind['u']
        v_wind = wind['v']
        
        # Top row: Fire front and ROS
        ax_fire = axes[0, idx]
        
        # Plot fire front as contour
        burned = (phi <= 0.0).astype(float)
        im = ax_fire.contourf(burned, levels=[0, 0.5, 1], colors=['white', 'orange', 'red'],
                              alpha=0.7, extent=[state['xmin'], state['xmax'],
                                                state['ymin'], state['ymax']])
        
        # Overlay wind vectors
        skip = max(1, state['nx'] // 20)
        x = np.linspace(state['xmin'], state['xmax'], state['nx'])
        y = np.linspace(state['ymin'], state['ymax'], state['ny'])
        X, Y = np.meshgrid(x[::skip], y[::skip])
        ax_fire.quiver(X, Y, u_wind[::skip, ::skip], v_wind[::skip, ::skip],
                      alpha=0.6, scale=100, color='blue')
        
        ax_fire.set_xlabel('X (m)')
        ax_fire.set_ylabel('Y (m)')
        ax_fire.set_title(f't = {time:.1f} s')
        ax_fire.set_aspect('equal')
        
        # Bottom row: Rate of spread
        ax_ros = axes[1, idx]
        ros_plot = ax_ros.contourf(ros, levels=20, cmap='hot',
                                   extent=[state['xmin'], state['xmax'],
                                          state['ymin'], state['ymax']])
        plt.colorbar(ros_plot, ax=ax_ros, label='ROS (m/s)')
        ax_ros.set_xlabel('X (m)')
        ax_ros.set_ylabel('Y (m)')
        ax_ros.set_title(f'Rate of Spread')
        ax_ros.set_aspect('equal')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"✓ Saved visualization to {output_file}")
    plt.close()


def demonstrate_two_way_coupling_concept():
    """
    Demonstrate the two-way coupling pattern from the problem statement.
    
    This shows the conceptual design for bidirectional coupling:
    
    while fire.time < final_time:
        # Wind → Fire
        wind.solve(fire.time)
        u_3d, v_3d, w_3d = wind.get_velocity_arrays()
        fire.update_wind_3d(u_3d, v_3d, w_3d, wind.nz, wind.zmin, wind.zmax)

        # Advance fire
        fire.step()

        # Fire → Wind (future feature)
        state = fire.get_state()
        heat_release = compute_heat_release(state)
        wind.add_heat_source(heat_release)  # Affects next wind solve
    """
    print("\n" + "="*70)
    print("Two-Way Wind-Fire Coupling Pattern")
    print("="*70)
    print("\nConceptual design for bidirectional coupling:")
    print()
    print("  while fire.time < final_time:")
    print("      # Wind → Fire (already implemented)")
    print("      wind.solve(fire.time)")
    print("      u_3d, v_3d, w_3d = wind.get_velocity_arrays()")
    print("      fire.update_wind_3d(u_3d, v_3d, w_3d, wind.nz, wind.zmin, wind.zmax)")
    print()
    print("      # Advance fire")
    print("      fire.step()")
    print()
    print("      # Fire → Wind (now available via compute_heat_release)")
    print("      state = fire.get_state()")
    print("      heat_release = fire.compute_heat_release(state)")
    print("      # surface_fluxes = fire.get_surface_fluxes()")
    print("      # wind.add_heat_source(heat_release)  # Future: when wind solver supports it")
    print()
    print("="*70)
    print()


def main():
    """Run coupled wind-fire simulation example."""
    
    if len(sys.argv) < 2:
        print("Usage: python coupled_wind_fire_example.py <inputs_file>")
        print("\nExample:")
        print("  PYTHONPATH=build/python python3 src/python/coupled_wind_fire_example.py \\")
        print("    regtest/surface_spread/farsite_ellipse/inputs.i")
        return 1
    
    inputs_file = sys.argv[1]
    
    print("=" * 70)
    print("Coupled Wind-Fire Simulation Example")
    print("=" * 70)
    
    # Show the two-way coupling concept
    demonstrate_two_way_coupling_concept()
    print()
    print("This demonstrates the workflow for coupling an external wind solver")
    print("(e.g., massconsistent_amr) with wildfire_levelset:")
    print()
    print("  1. Wind solver generates 3D wind field")
    print("  2. Wind data passed to fire solver")
    print("  3. Fire solver advances one timestep")
    print("  4. Repeat")
    print()
    print("=" * 70)
    print()
    
    # Initialize fire solver
    print("Initializing fire solver...")
    fire = WildfireSolver(inputs_file)
    
    # 3D wind grid parameters (for spotting calculations)
    # In a real coupled simulation, these would come from massconsistent_amr
    nz = 8
    zmin = 0.0
    zmax = 100.0  # 100 meters above ground level
    
    # Storage for visualization
    states = []
    wind_states = []
    
    # Run coupled simulation
    num_steps = 20
    print(f"\nRunning coupled simulation for {num_steps} timesteps...")
    print()
    
    for step in range(num_steps):
        current_time = fire.time
        
        # Step 1: Generate 3D wind field (simulates massconsistent_amr solver)
        # In a real application, this would be:
        #   u_3d, v_3d, w_3d = wind_solver.solve_and_extract(current_time)
        u_3d, v_3d, w_3d = generate_synthetic_wind_3d(
            fire.nx, fire.ny, nz,
            fire.xmin, fire.xmax,
            fire.ymin, fire.ymax,
            zmin, zmax,
            current_time
        )
        
        # Step 2: Pass 3D wind to fire solver
        # This updates both the 3D wind (for spotting) and 2D column-averaged wind
        # (for surface fire spread)
        fire.update_wind_3d(u_3d, v_3d, w_3d, nz, zmin, zmax)
        
        # Step 3: Advance fire simulation one timestep
        result = fire.step()
        
        # Get current state for visualization
        state = fire.get_state()
        u_2d = state['u_wind']
        v_2d = state['v_wind']
        
        # Step 4: Extract heat release (Fire → Wind coupling)
        # This is now available for passing to wind solver
        heat_data = fire.compute_heat_release(state)
        surface_flux = fire.get_surface_fluxes()
        
        # Store for visualization
        states.append(state)
        wind_states.append({'u': u_2d, 'v': v_2d})
        
        # Print progress with heat release information
        phi = state['phi']
        burned_area = np.sum(phi <= 0.0) * fire.dx * fire.dy
        wind_speed = np.sqrt(u_2d**2 + v_2d**2).mean()
        
        print(f"  Step {step+1:3d}: t={current_time:7.2f} s, "
              f"dt={result['dt']:6.3f} s, "
              f"burned={burned_area:8.1f} m², "
              f"wind={wind_speed:5.2f} m/s, "
              f"heat={heat_data['total_heat_release']:8.0f} kW")
        
        # Write plotfile every 5 steps
        if (step + 1) % 5 == 0:
            plotfile = f"coupled_plt{(step+1)//5:03d}"
            fire.write_plotfile(plotfile)
            print(f"           → Wrote {plotfile}")
    
    # Create visualization
    print()
    print("Creating visualization...")
    visualize_coupled_simulation(states, wind_states, "coupled_wind_fire.png")
    
    # Print summary statistics
    print()
    print("=" * 70)
    print("Simulation Summary")
    print("=" * 70)
    final_state = states[-1]
    phi_final = final_state['phi']
    arrival_final = final_state['arrival_time']
    
    total_burned = np.sum(phi_final <= 0.0) * fire.dx * fire.dy
    total_area = fire.nx * fire.ny * fire.dx * fire.dy
    burned_fraction = total_burned / total_area * 100
    
    # Calculate spread rate
    burned_mask = (phi_final <= 0.0)
    if np.any(burned_mask):
        arrival_burned = arrival_final[burned_mask]
        arrival_burned = arrival_burned[arrival_burned >= 0]
        if len(arrival_burned) > 0:
            mean_arrival_time = arrival_burned.mean()
        else:
            mean_arrival_time = 0.0
    else:
        mean_arrival_time = 0.0
    
    # Get final heat release data
    final_heat = fire.compute_heat_release(final_state)
    
    print(f"Final time: {final_state['time']:.2f} s ({final_state['step']} steps)")
    print(f"Total burned area: {total_burned:.1f} m² ({burned_fraction:.1f}% of domain)")
    print(f"Mean arrival time: {mean_arrival_time:.1f} s")
    print(f"Max ROS: {final_state['ros'].max():.3f} m/s")
    print(f"Max intensity: {final_state['intensity'].max():.1f} kW/m")
    print(f"Max flame length: {final_state['flame_length'].max():.2f} m")
    print()
    print("Heat Release Summary:")
    print(f"  Fire perimeter: {final_heat['perimeter']:.1f} m")
    print(f"  Mean intensity: {final_heat['mean_intensity']:.1f} kW/m")
    print(f"  Total heat release: {final_heat['total_heat_release']:.0f} kW")
    print(f"  Max surface flux: {final_heat['surface_flux'].max():.1f} kW/m²")
    print()
    
    # Finalize solver
    fire.finalize()
    
    print("=" * 70)
    print("✓ Coupled simulation complete!")
    print("=" * 70)
    print()
    print("Heat Release Coupling Implementation Status:")
    print("  ✓ Fire → Wind: compute_heat_release() and get_surface_fluxes() available")
    print("  ✓ Surface flux array ready for wind solver coupling")
    print("  ⧐ Wind → Fire: update_wind_3d() already supports 3D wind input")
    print()
    print("Integration with massconsistent_amr:")
    print("  1. Use CoupledWindFireSolver class for full workflow")
    print("  2. Wind solver needs add_heat_source() method for fire feedback")
    print("  3. Surface fluxes from fire.get_surface_fluxes()['heat_flux']")
    print()
    print("Next steps:")
    print("  1. Implement wind.add_heat_source() in massconsistent_amr")
    print("  2. Use CoupledWindFireSolver for two-way coupling")
    print("  3. Add time-varying fuel moisture from atmospheric model")
    print()
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
