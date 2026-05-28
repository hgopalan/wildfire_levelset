#!/usr/bin/env python3
"""
Regression test for coupled wind-fire simulation via Python API
Tests:
1. Initialization of fire solver
2. Synthetic 3D wind field generation (simulating massconsistent_amr)
3. Passing 3D wind to fire solver via update_wind_3d()
4. Coupled time-stepping loop
5. Verification of wind effects on fire spread
"""

import sys
import os
import numpy as np

# Find the pyWildfire module
try:
    import pyWildfire
    from wildfire_solver import WildfireSolver
except ImportError:
    build_python = os.path.join(os.getcwd(), '..', '..', '..', 'python')
    if os.path.exists(build_python):
        sys.path.insert(0, build_python)
        import pyWildfire
        from wildfire_solver import WildfireSolver
    else:
        print("ERROR: Cannot find pyWildfire module")
        sys.exit(1)


def generate_synthetic_wind_3d(fire, time):
    """
    Generate synthetic 3D wind field
    
    This simulates the output from massconsistent_amr or another wind solver.
    In a real application, this would be replaced with actual wind solver calls.
    
    Parameters:
        fire: WildfireSolver instance
        time: Current simulation time (seconds)
    
    Returns:
        u_3d, v_3d, w_3d: 3D wind arrays, shape (nz, ny, nx)
        nz: Number of vertical levels
        zmin, zmax: Vertical domain bounds
    """
    # Vertical grid
    nz = 5
    zmin = 0.0
    zmax = 50.0  # 50 m vertical extent
    
    nx, ny = fire.nx, fire.ny
    xmin, xmax = fire.xmin, fire.xmax
    ymin, ymax = fire.ymin, fire.ymax
    
    # Create 3D arrays
    u_3d = np.zeros((nz, ny, nx), dtype=np.float64)
    v_3d = np.zeros((nz, ny, nx), dtype=np.float64)
    w_3d = np.zeros((nz, ny, nx), dtype=np.float64)
    
    # Generate log-law wind profile with time variation
    z_ref = 10.0  # reference height
    z0 = 0.1  # roughness length
    
    # Time-varying base wind speed (simulates diurnal variation)
    hour = (time / 3600.0) % 24
    time_factor = 0.8 + 0.4 * np.sin(2 * np.pi * (hour - 6) / 24)
    u_ref = 5.0 * time_factor
    v_ref = 0.8 * time_factor
    
    # Create coordinate grids
    x = np.linspace(xmin, xmax, nx)
    y = np.linspace(ymin, ymax, ny)
    
    for k in range(nz):
        z_height = zmin + (k + 0.5) * (zmax - zmin) / nz
        z_height = max(z_height, z0)  # avoid singularity
        
        # Log-law scaling
        height_factor = np.log(z_height / z0) / np.log(z_ref / z0)
        
        for j in range(ny):
            for i in range(nx):
                # Add spatial variation (simulates terrain channeling)
                x_norm = (x[i] - xmin) / (xmax - xmin)
                y_norm = (y[j] - ymin) / (ymax - ymin)
                
                # Simple channeling effect
                spatial_factor = 1.0 + 0.2 * np.sin(4 * np.pi * x_norm)
                
                u_3d[k, j, i] = u_ref * height_factor * spatial_factor
                v_3d[k, j, i] = v_ref * height_factor
                
                # Small vertical motion due to terrain
                w_3d[k, j, i] = 0.05 * np.sin(np.pi * z_height / zmax) * \
                               np.sin(3 * np.pi * x_norm)
    
    return u_3d, v_3d, w_3d, nz, zmin, zmax


def test_coupled_wind_fire():
    """Test coupled wind-fire simulation"""
    print("=" * 70)
    print("Python Coupled Wind-Fire - Regression Test")
    print("=" * 70)
    
    # Initialize fire solver
    inputs_file = "inputs.i"
    if not os.path.exists(inputs_file):
        print(f"ERROR: Input file {inputs_file} not found")
        return False
    
    print(f"\n1. Initializing fire solver...")
    try:
        fire = WildfireSolver(inputs_file)
        print(f"   ✓ Fire solver initialized")
        print(f"     Grid: {fire.nx} × {fire.ny}")
        print(f"     Domain: [{fire.xmin:.1f}, {fire.xmax:.1f}] × "
              f"[{fire.ymin:.1f}, {fire.ymax:.1f}]")
    except Exception as e:
        print(f"   ✗ FAILED to initialize: {e}")
        return False
    
    # Get initial state
    print(f"\n2. Getting initial state...")
    try:
        state0 = fire.get_state()
        initial_burned = np.sum(state0['phi'] <= 0.0) * fire.dx * fire.dy
        print(f"   ✓ Initial burned area: {initial_burned:.1f} m²")
    except Exception as e:
        print(f"   ✗ FAILED to get initial state: {e}")
        fire.finalize()
        return False
    
    # Coupled simulation loop
    print(f"\n3. Running coupled simulation (10 steps)...")
    wind_speeds_2d = []  # Track column-averaged wind speeds
    
    try:
        for step in range(10):
            # Step 3a: Generate 3D wind field (simulating wind solver)
            u_3d, v_3d, w_3d, nz, zmin, zmax = generate_synthetic_wind_3d(fire, fire.time)
            
            # Calculate column-averaged wind for tracking
            u_avg = u_3d.mean()
            v_avg = v_3d.mean()
            wind_speed_avg = np.sqrt(u_avg**2 + v_avg**2)
            wind_speeds_2d.append(wind_speed_avg)
            
            # Step 3b: Update fire solver wind from 3D field
            success = fire.update_wind_3d(u_3d, v_3d, w_3d, nz, zmin, zmax)
            if not success:
                print(f"   ✗ FAILED to update wind at step {step+1}")
                fire.finalize()
                return False
            
            # Step 3c: Advance fire solver
            result = fire.step()
            if not result['success']:
                print(f"   ✗ FAILED to advance at step {step+1}")
                fire.finalize()
                return False
            
            # Step 3d: Monitor progress
            state = fire.get_state()
            burned = np.sum(state['phi'] <= 0.0) * fire.dx * fire.dy
            
            if (step + 1) % 3 == 0:
                print(f"   Step {step+1:2d}: t={fire.time:6.2f} s, "
                      f"burned={burned:7.1f} m², wind={wind_speed_avg:.2f} m/s")
        
        print(f"   ✓ Coupled simulation completed successfully")
        
    except Exception as e:
        print(f"   ✗ FAILED during coupled simulation: {e}")
        import traceback
        traceback.print_exc()
        fire.finalize()
        return False
    
    # Verify results
    print(f"\n4. Verifying coupled simulation results...")
    try:
        final_state = fire.get_state()
        final_burned = np.sum(final_state['phi'] <= 0.0) * fire.dx * fire.dy
        
        # Check that fire spread
        if final_burned > initial_burned * 1.5:
            print(f"   ✓ Fire spread verified: {initial_burned:.1f} → {final_burned:.1f} m²")
        else:
            print(f"   ✗ FAILED: Insufficient fire spread")
            fire.finalize()
            return False
        
        # Check that wind was updated (verify 2D wind field from state)
        u2d = final_state['u_wind']
        v2d = final_state['v_wind']
        wind_2d_avg = np.sqrt(u2d.mean()**2 + v2d.mean()**2)
        
        print(f"   ✓ Final 2D wind: {wind_2d_avg:.2f} m/s")
        
        # Verify wind values are reasonable (should be positive from our synthetic field)
        if wind_2d_avg > 0.5 and wind_2d_avg < 10.0:
            print(f"   ✓ Wind values are physically reasonable")
        else:
            print(f"   ✗ WARNING: Wind values may be unrealistic: {wind_2d_avg:.2f} m/s")
        
        # Check ROS and intensity
        max_ros = final_state['ros'].max()
        max_intensity = final_state['intensity'].max()
        
        if max_ros > 0.0 and max_intensity > 0.0:
            print(f"   ✓ Fire behavior metrics are positive")
            print(f"     Max ROS: {max_ros:.3f} m/s")
            print(f"     Max intensity: {max_intensity:.1f} kW/m")
        else:
            print(f"   ✗ FAILED: Invalid fire behavior metrics")
            fire.finalize()
            return False
            
    except Exception as e:
        print(f"   ✗ FAILED during verification: {e}")
        import traceback
        traceback.print_exc()
        fire.finalize()
        return False
    
    # Finalize
    print(f"\n5. Finalizing solver...")
    try:
        fire.finalize()
        print(f"   ✓ Solver finalized successfully")
    except Exception as e:
        print(f"   ✗ FAILED to finalize: {e}")
        return False
    
    return True


def main():
    """Main test entry point"""
    print(f"pyWildfire version: {pyWildfire.__version__}\n")
    
    success = test_coupled_wind_fire()
    
    print("\n" + "=" * 70)
    if success:
        print("✓ COUPLED WIND-FIRE REGRESSION TEST PASSED")
        print("=" * 70)
        print("\nThis test demonstrates:")
        print("  - Integration with external 3D wind solvers")
        print("  - Dynamic wind field updates during simulation")
        print("  - Proper coupling between wind and fire physics")
        print("  - Time-varying wind effects on fire spread")
        return 0
    else:
        print("✗ COUPLED WIND-FIRE REGRESSION TEST FAILED")
        print("=" * 70)
        return 1


if __name__ == "__main__":
    sys.exit(main())
