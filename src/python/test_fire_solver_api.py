#!/usr/bin/env python3
"""
test_fire_solver_api.py - Test script for pyWildfire fire solver Python API

Tests the new fire solver control functions:
- Initialization from inputs file
- Time-stepping
- State extraction
- Wind field updates
- Plotfile writing

Run after building with Python bindings enabled:
    cmake -S . -B build -DLEVELSET_BUILD_PYTHON_BINDINGS=ON -DLEVELSET_DIM_2D=ON
    cmake --build build
    PYTHONPATH=build/python python3 src/python/test_fire_solver_api.py
"""

import sys
import os
import numpy as np

try:
    import pyWildfire
    from wildfire_solver import WildfireSolver
except ImportError as e:
    print(f"Error: Could not import pyWildfire or wildfire_solver")
    print(f"  {e}")
    print("\nMake sure to:")
    print("  1. Build with -DLEVELSET_BUILD_PYTHON_BINDINGS=ON")
    print("  2. Set PYTHONPATH to point to build/python directory")
    print("  Example: PYTHONPATH=build/python python3 src/python/test_fire_solver_api.py")
    sys.exit(1)


def test_basic_initialization():
    """Test 1: Basic fire solver initialization"""
    print("=" * 70)
    print("Test 1: Basic fire solver initialization")
    print("=" * 70)
    
    # Create a minimal inputs file
    inputs_content = """
# Minimal inputs for testing
n_cell_x = 32
n_cell_y = 32
max_grid = 16

plo_x = 0.0
plo_y = 0.0
phi_x = 1000.0
phi_y = 1000.0

# Constant wind
ux = 5.0
uy = 0.0

# Simple circular ignition
ignition.type = circle
ignition.x0 = 500.0
ignition.y0 = 500.0
ignition.radius = 50.0

# CFL and time control
cfl = 0.5
nsteps = 10

# Propagation method
propagation_method = levelset
"""
    
    # Write test inputs file
    test_inputs = "/tmp/test_fire_solver_inputs.i"
    with open(test_inputs, 'w') as f:
        f.write(inputs_content)
    
    try:
        # Test low-level API
        result = pyWildfire.initialize(test_inputs)
        
        print(f"\nInitialization result:")
        print(f"  Success: {result['success']}")
        print(f"  Grid: {result['nx']} × {result['ny']}")
        print(f"  Domain: X=[{result['xmin']}, {result['xmax']}], "
              f"Y=[{result['ymin']}, {result['ymax']}]")
        print(f"  Cell size: dx={result['dx']:.2f} m, dy={result['dy']:.2f} m")
        
        if result['success']:
            print("\n✓ Test PASSED: Initialization successful")
            pyWildfire.finalize()
            return True
        else:
            print("\n✗ Test FAILED: Initialization failed")
            return False
            
    except Exception as e:
        print(f"\n✗ Test FAILED with exception: {e}")
        import traceback
        traceback.print_exc()
        return False
    finally:
        if os.path.exists(test_inputs):
            os.remove(test_inputs)


def test_time_stepping():
    """Test 2: Time-stepping and state extraction"""
    print("\n" + "=" * 70)
    print("Test 2: Time-stepping and state extraction")
    print("=" * 70)
    
    # Create test inputs
    inputs_content = """
n_cell_x = 16
n_cell_y = 16
max_grid = 16

plo_x = 0.0
plo_y = 0.0
phi_x = 500.0
phi_y = 500.0

ux = 3.0
uy = 0.5

ignition.type = circle
ignition.x0 = 250.0
ignition.y0 = 250.0
ignition.radius = 30.0

cfl = 0.5
nsteps = 5
propagation_method = levelset
"""
    
    test_inputs = "/tmp/test_timestep_inputs.i"
    with open(test_inputs, 'w') as f:
        f.write(inputs_content)
    
    try:
        # Initialize
        result = pyWildfire.initialize(test_inputs)
        if not result['success']:
            print("✗ Test FAILED: Could not initialize")
            return False
        
        nx, ny = result['nx'], result['ny']
        dx, dy = result['dx'], result['dy']
        
        # Run a few timesteps
        print(f"\nRunning 3 timesteps...")
        for step in range(3):
            # Advance
            step_result = pyWildfire.advance()
            
            if not step_result['success']:
                print(f"✗ Test FAILED: Timestep {step+1} failed")
                pyWildfire.finalize()
                return False
            
            # Get state
            state = pyWildfire.get_state()
            
            # Check state fields
            phi = state['phi']
            ros = state['ros']
            intensity = state['intensity']
            flame_length = state['flame_length']
            u_wind = state['u_wind']
            v_wind = state['v_wind']
            
            # Verify shapes
            expected_shape = (ny, nx)
            if phi.shape != expected_shape:
                print(f"✗ Test FAILED: phi shape {phi.shape} != {expected_shape}")
                pyWildfire.finalize()
                return False
            
            # Calculate burned area
            burned_area = np.sum(phi <= 0.0) * dx * dy
            
            print(f"  Step {step+1}: t={state['time']:.3f} s, "
                  f"dt={step_result['dt']:.4f} s, "
                  f"burned={burned_area:.1f} m²")
            print(f"           ROS: min={ros.min():.3f}, max={ros.max():.3f} m/s")
            print(f"           Intensity: max={intensity.max():.1f} kW/m")
        
        # Verify fire has spread
        final_state = pyWildfire.get_state()
        final_burned = np.sum(final_state['phi'] <= 0.0) * dx * dy
        
        # Initial burned area (circle of radius 30 m)
        initial_burned = np.pi * 30.0**2
        
        if final_burned > initial_burned * 1.5:
            print(f"\n✓ Test PASSED: Fire spread from {initial_burned:.1f} to "
                  f"{final_burned:.1f} m²")
            pyWildfire.finalize()
            return True
        else:
            print(f"\n✗ Test FAILED: Fire did not spread sufficiently")
            print(f"  Initial: {initial_burned:.1f} m², Final: {final_burned:.1f} m²")
            pyWildfire.finalize()
            return False
            
    except Exception as e:
        print(f"\n✗ Test FAILED with exception: {e}")
        import traceback
        traceback.print_exc()
        return False
    finally:
        if os.path.exists(test_inputs):
            os.remove(test_inputs)


def test_wind_update():
    """Test 3: Wind field update"""
    print("\n" + "=" * 70)
    print("Test 3: Wind field update")
    print("=" * 70)
    
    # Create test inputs
    inputs_content = """
n_cell_x = 16
n_cell_y = 16
max_grid = 16
plo_x = 0.0
plo_y = 0.0
phi_x = 500.0
phi_y = 500.0
ux = 1.0
uy = 0.0
ignition.type = circle
ignition.x0 = 250.0
ignition.y0 = 250.0
ignition.radius = 30.0
cfl = 0.5
propagation_method = levelset
"""
    
    test_inputs = "/tmp/test_wind_update_inputs.i"
    with open(test_inputs, 'w') as f:
        f.write(inputs_content)
    
    try:
        # Initialize
        result = pyWildfire.initialize(test_inputs)
        if not result['success']:
            print("✗ Test FAILED: Could not initialize")
            return False
        
        nx, ny = result['nx'], result['ny']
        
        # Create new wind field (uniform 10 m/s easterly)
        u_new = np.full((ny, nx), 10.0)
        v_new = np.zeros((ny, nx))
        
        print(f"\nUpdating wind field to uniform 10 m/s easterly...")
        success = pyWildfire.update_wind(u_new, v_new)
        
        if not success:
            print("✗ Test FAILED: Wind update failed")
            pyWildfire.finalize()
            return False
        
        # Verify wind was updated
        state = pyWildfire.get_state()
        u_wind = state['u_wind']
        v_wind = state['v_wind']
        
        u_mean = u_wind.mean()
        v_mean = v_wind.mean()
        
        print(f"  Retrieved wind: u_mean={u_mean:.2f} m/s, v_mean={v_mean:.2f} m/s")
        
        if abs(u_mean - 10.0) < 0.1 and abs(v_mean - 0.0) < 0.1:
            print("\n✓ Test PASSED: Wind field updated successfully")
            pyWildfire.finalize()
            return True
        else:
            print(f"\n✗ Test FAILED: Wind field not updated correctly")
            print(f"  Expected: u=10.0, v=0.0")
            print(f"  Got: u={u_mean:.2f}, v={v_mean:.2f}")
            pyWildfire.finalize()
            return False
            
    except Exception as e:
        print(f"\n✗ Test FAILED with exception: {e}")
        import traceback
        traceback.print_exc()
        return False
    finally:
        if os.path.exists(test_inputs):
            os.remove(test_inputs)


def test_wildfire_solver_class():
    """Test 4: High-level WildfireSolver class"""
    print("\n" + "=" * 70)
    print("Test 4: High-level WildfireSolver class")
    print("=" * 70)
    
    # Create test inputs
    inputs_content = """
n_cell_x = 20
n_cell_y = 20
max_grid = 16
plo_x = 0.0
plo_y = 0.0
phi_x = 600.0
phi_y = 600.0
ux = 4.0
uy = 1.0
ignition.type = circle
ignition.x0 = 300.0
ignition.y0 = 300.0
ignition.radius = 40.0
cfl = 0.5
propagation_method = levelset
"""
    
    test_inputs = "/tmp/test_solver_class_inputs.i"
    with open(test_inputs, 'w') as f:
        f.write(inputs_content)
    
    try:
        # Initialize using high-level class
        fire = WildfireSolver(test_inputs)
        
        print(f"\nSolver initialized:")
        print(f"  Grid: {fire.nx} × {fire.ny}")
        print(f"  Domain: X=[{fire.xmin:.1f}, {fire.xmax:.1f}], "
              f"Y=[{fire.ymin:.1f}, {fire.ymax:.1f}]")
        
        # Run a few steps
        print(f"\nRunning 5 steps...")
        for i in range(5):
            result = fire.step()
            state = fire.get_state()
            burned = np.sum(state['phi'] <= 0.0) * fire.dx * fire.dy
            print(f"  Step {i+1}: t={fire.time:.2f} s, burned={burned:.1f} m²")
        
        # Get final state
        final_state = fire.get_state()
        final_burned = np.sum(final_state['phi'] <= 0.0) * fire.dx * fire.dy
        
        fire.finalize()
        
        if final_burned > 0:
            print(f"\n✓ Test PASSED: WildfireSolver class works correctly")
            return True
        else:
            print(f"\n✗ Test FAILED: No fire spread detected")
            return False
            
    except Exception as e:
        print(f"\n✗ Test FAILED with exception: {e}")
        import traceback
        traceback.print_exc()
        return False
    finally:
        if os.path.exists(test_inputs):
            os.remove(test_inputs)


def main():
    """Run all tests"""
    print(f"pyWildfire version: {pyWildfire.__version__}\n")
    
    results = []
    results.append(test_basic_initialization())
    results.append(test_time_stepping())
    results.append(test_wind_update())
    results.append(test_wildfire_solver_class())
    
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Tests run: {len(results)}")
    print(f"Tests passed: {sum(results)}")
    print(f"Tests failed: {len(results) - sum(results)}")
    
    if all(results):
        print("\n✓ ALL TESTS PASSED")
        return 0
    else:
        print("\n✗ SOME TESTS FAILED")
        return 1


if __name__ == "__main__":
    sys.exit(main())
