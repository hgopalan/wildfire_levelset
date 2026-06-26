#!/usr/bin/env python3
"""
coupling_modes_example.py - Demonstration of one-way vs two-way wind-fire coupling

This script demonstrates both coupling modes:

1. One-way coupling (wind → fire):
   - Wind solver provides velocity field to fire solver
   - Fire spread is affected by wind
   - Fire does NOT affect wind dynamics
   - Use when wind field is pre-computed or from atmospheric model

2. Two-way coupling (wind ↔ fire):
   - Wind solver provides velocity field to fire solver
   - Fire spread is affected by wind
   - Fire heating affects wind dynamics (when wind solver supports add_heat_source)
   - More realistic but computationally expensive
   - Use for detailed fire-atmosphere interaction studies

Usage:
    PYTHONPATH=build/python python3 src/python/coupling_modes_example.py <inputs_file>

Example:
    PYTHONPATH=build/python python3 src/python/coupling_modes_example.py \\
        regtest/surface_spread/farsite_ellipse/inputs.i
"""

import sys
import numpy as np
from wildfire_solver import WildfireSolver


def generate_synthetic_wind_3d(nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax, time):
    """
    Generate synthetic 3D wind field that varies with time.
    
    This simulates the output from massconsistent_amr or another 3D wind solver.
    In a real application, this would be replaced with actual wind solver calls.
    """
    x = np.linspace(xmin, xmax, nx)
    y = np.linspace(ymin, ymax, ny)
    z = np.linspace(zmin, zmax, nz)
    
    u_3d = np.zeros((nz, ny, nx))
    v_3d = np.zeros((nz, ny, nx))
    w_3d = np.zeros((nz, ny, nx))
    
    # Time-varying wind
    hour = (time / 3600.0) % 24
    diurnal_factor = 0.7 + 0.3 * np.sin(2 * np.pi * (hour - 6) / 24)
    
    u_ref = 5.0 * diurnal_factor  # m/s westerly
    v_ref = 1.0 * diurnal_factor  # m/s southerly
    z_ref = 10.0
    z0 = 0.1
    
    for k in range(nz):
        z_height = z[k] if z[k] > z0 else z0
        height_factor = np.log(z_height / z0) / np.log(z_ref / z0)
        
        for j in range(ny):
            for i in range(nx):
                x_norm = (x[i] - xmin) / (xmax - xmin)
                y_norm = (y[j] - ymin) / (ymax - ymin)
                
                spatial_mod = 1.0 + 0.3 * np.sin(3 * np.pi * x_norm) * np.cos(2 * np.pi * y_norm)
                
                u_3d[k, j, i] = u_ref * height_factor * spatial_mod
                v_3d[k, j, i] = v_ref * height_factor * (1.0 + 0.1 * np.cos(2 * np.pi * x_norm))
                w_3d[k, j, i] = 0.1 * np.sin(np.pi * z_height / zmax)
    
    return u_3d, v_3d, w_3d


def run_one_way_coupling(inputs_file, num_steps=20):
    """
    Run one-way coupling: Wind → Fire (Fire does not affect wind)
    
    In this mode:
    - Wind field is computed independently
    - Wind field is passed to fire solver
    - Fire spread responds to wind
    - Fire does NOT change wind field
    - Suitable for pre-computed wind fields or offline coupling
    """
    print("\n" + "="*80)
    print("ONE-WAY COUPLING DEMONSTRATION: Wind → Fire")
    print("="*80)
    print("\nMode: One-way coupling (wind → fire)")
    print("  • Wind field computed independently (not affected by fire)")
    print("  • Fire spread responds to wind dynamics")
    print("  • Fire heating does NOT modify wind")
    print("  • Suitable for: pre-computed wind, atmospheric model output")
    print("\nWorkflow:")
    print("  while fire.time < final_time:")
    print("      wind.solve()  # Solve wind independently")
    print("      u_3d, v_3d, w_3d = wind.get_velocity()  # Extract 3D wind")
    print("      fire.update_wind_3d(u_3d, v_3d, w_3d, ...)  # Update fire wind")
    print("      fire.step()  # Advance fire")
    print("      # Note: fire heating is NOT fed back to wind solver")
    print()
    
    # Initialize fire solver
    print("Initializing fire solver...")
    fire = WildfireSolver(inputs_file)
    
    # 3D wind parameters
    nz = 8
    zmin = 0.0
    zmax = 100.0
    
    print(f"\nRunning one-way coupled simulation for {num_steps} timesteps...")
    print()
    
    # Run simulation
    for step in range(num_steps):
        current_time = fire.time
        
        # Generate 3D wind field (in real app: wind_solver.solve())
        u_3d, v_3d, w_3d = generate_synthetic_wind_3d(
            fire.nx, fire.ny, nz,
            fire.xmin, fire.xmax,
            fire.ymin, fire.ymax,
            zmin, zmax,
            current_time
        )
        
        # Pass 3D wind to fire solver
        fire.update_wind_3d(u_3d, v_3d, w_3d, nz, zmin, zmax)
        
        # Advance fire
        result = fire.step()
        state = fire.get_state()
        
        # Get heat release (available but NOT fed back to wind in one-way mode)
        heat_data = fire.compute_heat_release(state)
        
        burned_area = heat_data['burned_area']
        wind_speed = np.sqrt(state['u_wind']**2 + state['v_wind']**2).mean()
        
        print(f"  Step {step+1:3d}: t={current_time:7.2f} s, "
              f"dt={result['dt']:6.3f} s, "
              f"burned={burned_area:9.1f} m², "
              f"wind={wind_speed:5.2f} m/s, "
              f"heat={heat_data['total_heat_release']:8.0f} kW "
              f"(not fed back to wind)")
    
    print()
    final_state = fire.get_state()
    final_heat = fire.compute_heat_release(final_state)
    
    print("One-way Coupling Results:")
    print(f"  Final time: {final_state['time']:.2f} s")
    print(f"  Total burned area: {final_heat['burned_area']:.1f} m²")
    print(f"  Fire perimeter: {final_heat['perimeter']:.1f} m")
    print(f"  Mean intensity: {final_heat['mean_intensity']:.1f} kW/m")
    print(f"  Total heat release: {final_heat['total_heat_release']:.0f} kW")
    print(f"  Max ROS: {final_heat['max_ros']:.3f} m/s")
    print()
    
    fire.finalize()
    return final_heat


def run_two_way_coupling(inputs_file, num_steps=20):
    """
    Run two-way coupling: Wind ↔ Fire (Fire affects wind heating)
    
    In this mode:
    - Wind field is computed, possibly affected by fire heating
    - Wind field is passed to fire solver
    - Fire spread responds to wind
    - Fire heating is extracted and available for wind solver
    - Note: Currently awaiting wind.add_heat_source() implementation
    - Suitable for detailed fire-atmosphere interaction studies
    """
    print("\n" + "="*80)
    print("TWO-WAY COUPLING DEMONSTRATION: Wind ↔ Fire")
    print("="*80)
    print("\nMode: Two-way coupling (wind ↔ fire)")
    print("  • Wind field computed with fire heating effects")
    print("  • Fire spread responds to wind dynamics")
    print("  • Fire heating modifies wind field (when wind solver supports add_heat_source)")
    print("  • More realistic fire-atmosphere interaction")
    print("  • Suitable for: detailed studies, fire-induced wind changes")
    print("\nWorkflow:")
    print("  while fire.time < final_time:")
    print("      wind.solve()  # Solve wind (with fire heat effects from previous step)")
    print("      u_3d, v_3d, w_3d = wind.get_velocity()  # Extract 3D wind")
    print("      fire.update_wind_3d(u_3d, v_3d, w_3d, ...)  # Update fire wind")
    print("      fire.step()  # Advance fire")
    print("      heat_flux = fire.get_surface_fluxes()  # Extract heat")
    print("      wind.add_heat_source(heat_flux)  # Feed heat back to wind")
    print()
    
    # Initialize fire solver
    print("Initializing fire solver...")
    fire = WildfireSolver(inputs_file)
    
    # 3D wind parameters
    nz = 8
    zmin = 0.0
    zmax = 100.0
    
    print(f"\nRunning two-way coupled simulation for {num_steps} timesteps...")
    print()
    
    # Run simulation
    for step in range(num_steps):
        current_time = fire.time
        
        # Generate 3D wind field (in real app: wind_solver.solve())
        u_3d, v_3d, w_3d = generate_synthetic_wind_3d(
            fire.nx, fire.ny, nz,
            fire.xmin, fire.xmax,
            fire.ymin, fire.ymax,
            zmin, zmax,
            current_time
        )
        
        # Pass 3D wind to fire solver
        fire.update_wind_3d(u_3d, v_3d, w_3d, nz, zmin, zmax)
        
        # Advance fire
        result = fire.step()
        state = fire.get_state()
        
        # Get heat release and surface fluxes for wind solver
        heat_data = fire.compute_heat_release(state)
        flux_data = fire.get_surface_fluxes()
        
        burned_area = heat_data['burned_area']
        wind_speed = np.sqrt(state['u_wind']**2 + state['v_wind']**2).mean()
        max_flux = flux_data['heat_flux'].max()
        
        print(f"  Step {step+1:3d}: t={current_time:7.2f} s, "
              f"dt={result['dt']:6.3f} s, "
              f"burned={burned_area:9.1f} m², "
              f"wind={wind_speed:5.2f} m/s, "
              f"heat={heat_data['total_heat_release']:8.0f} kW, "
              f"flux_max={max_flux:6.1f} kW/m²")
        
        # In two-way mode, surface fluxes are available for wind solver
        # When wind.add_heat_source() is implemented:
        #   wind.add_heat_source(flux_data['heat_flux'], flux_data['grid_info'])
    
    print()
    final_state = fire.get_state()
    final_heat = fire.compute_heat_release(final_state)
    
    print("Two-way Coupling Results:")
    print(f"  Final time: {final_state['time']:.2f} s")
    print(f"  Total burned area: {final_heat['burned_area']:.1f} m²")
    print(f"  Fire perimeter: {final_heat['perimeter']:.1f} m")
    print(f"  Mean intensity: {final_heat['mean_intensity']:.1f} kW/m")
    print(f"  Total heat release: {final_heat['total_heat_release']:.0f} kW")
    print(f"  Max ROS: {final_heat['max_ros']:.3f} m/s")
    print()
    
    fire.finalize()
    return final_heat


def compare_coupling_modes(one_way_results, two_way_results):
    """
    Compare results between one-way and two-way coupling modes.
    
    Note: In this example, both use synthetic wind, so differences will be
    minimal. In real application with actual wind solver, differences would
    be significant where fire heating affects wind.
    """
    print("\n" + "="*80)
    print("COUPLING MODE COMPARISON")
    print("="*80)
    print()
    
    print("One-way vs Two-way Coupling Results:")
    print()
    print(f"{'Metric':<30} {'One-way':<20} {'Two-way':<20} {'Difference'}")
    print("-" * 80)
    
    metrics = [
        ('Burned Area (m²)', 'burned_area'),
        ('Fire Perimeter (m)', 'perimeter'),
        ('Mean Intensity (kW/m)', 'mean_intensity'),
        ('Total Heat Release (kW)', 'total_heat_release'),
        ('Max ROS (m/s)', 'max_ros')
    ]
    
    for label, key in metrics:
        one_val = one_way_results[key]
        two_val = two_way_results[key]
        diff = two_val - one_val
        pct = (diff / one_val * 100) if one_val != 0 else 0
        
        print(f"{label:<30} {one_val:<20.1f} {two_val:<20.1f} "
              f"{diff:+.1f} ({pct:+.1f}%)")
    
    print()
    print("Notes:")
    print("  • In this example, wind is synthetic and fixed, so differences are minimal")
    print("  • In real application with wind.add_heat_source(), fire heating")
    print("    would cause wind to change, leading to significant differences")
    print("  • Two-way coupling shows more realistic fire-atmosphere interaction")
    print()


def main():
    """Demonstrate both coupling modes."""
    
    if len(sys.argv) < 2:
        print("Usage: python coupling_modes_example.py <inputs_file>")
        print("\nExample:")
        print("  PYTHONPATH=build/python python3 src/python/coupling_modes_example.py \\")
        print("    regtest/surface_spread/farsite_ellipse/inputs.i")
        return 1
    
    inputs_file = sys.argv[1]
    num_steps = 20
    
    print("="*80)
    print("Wind-Fire Coupling Modes Comparison")
    print("="*80)
    print()
    print("This demonstrates the two coupling modes available:")
    print()
    print("1. ONE-WAY COUPLING (wind → fire):")
    print("   - Wind field computed independently")
    print("   - Wind affects fire spread")
    print("   - Fire does NOT affect wind")
    print("   - Use case: Pre-computed wind, offline coupling with weather model")
    print()
    print("2. TWO-WAY COUPLING (wind ↔ fire):")
    print("   - Wind field computed with fire heating effects")
    print("   - Wind affects fire spread")
    print("   - Fire heating affects wind (when add_heat_source is available)")
    print("   - Use case: Detailed fire-atmosphere interaction studies")
    print()
    print("="*80)
    print()
    
    # Run one-way coupling
    one_way_results = run_one_way_coupling(inputs_file, num_steps)
    
    # Run two-way coupling
    two_way_results = run_two_way_coupling(inputs_file, num_steps)
    
    # Compare results
    compare_coupling_modes(one_way_results, two_way_results)
    
    print("="*80)
    print("Summary:")
    print("="*80)
    print()
    print("✓ One-way coupling: IMPLEMENTED")
    print("  - Use WildfireSolver.update_wind_3d() and step()")
    print("  - Use CoupledWindFireSolver(coupling_mode='one_way')")
    print()
    print("✓ Two-way coupling: FRAMEWORK READY")
    print("  - WildfireSolver.compute_heat_release() and get_surface_fluxes() available")
    print("  - Use CoupledWindFireSolver(coupling_mode='two_way')")
    print("  - Awaiting wind.add_heat_source() implementation in massconsistent_amr")
    print()
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
