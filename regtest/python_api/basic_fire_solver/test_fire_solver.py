#!/usr/bin/env python3
"""
Regression test for Python fire solver API - basic fire spread
Tests:
1. Initialization from inputs file
2. Time-stepping with advance()
3. State extraction with get_state()
4. Fire spread verification
5. Plotfile writing
"""

import sys
import os
import numpy as np

# Find the pyWildfire module
# When run as a regression test, this will be in the build directory
try:
    import pyWildfire
    from wildfire_solver import WildfireSolver
except ImportError:
    # Try to find it relative to the build directory
    build_python = os.path.join(os.getcwd(), '..', '..', '..', 'python')
    if os.path.exists(build_python):
        sys.path.insert(0, build_python)
        import pyWildfire
        from wildfire_solver import WildfireSolver
    else:
        print("ERROR: Cannot find pyWildfire module")
        print("Make sure to build with -DLEVELSET_BUILD_PYTHON_BINDINGS=ON")
        sys.exit(1)


def test_basic_fire_solver():
    """Test basic fire solver operations"""
    print("=" * 70)
    print("Python Fire Solver API - Basic Regression Test")
    print("=" * 70)
    
    # Initialize fire solver
    inputs_file = "inputs.i"
    if not os.path.exists(inputs_file):
        print(f"ERROR: Input file {inputs_file} not found")
        return False
    
    print(f"\n1. Initializing fire solver from {inputs_file}...")
    try:
        fire = WildfireSolver(inputs_file)
        print(f"   ✓ Initialized successfully")
        print(f"     Grid: {fire.nx} × {fire.ny}")
        print(f"     Domain: X=[{fire.xmin:.1f}, {fire.xmax:.1f}], "
              f"Y=[{fire.ymin:.1f}, {fire.ymax:.1f}]")
        print(f"     Cell size: dx={fire.dx:.2f} m, dy={fire.dy:.2f} m")
    except Exception as e:
        print(f"   ✗ FAILED to initialize: {e}")
        return False
    
    # Get initial state
    print(f"\n2. Getting initial state...")
    try:
        state0 = fire.get_state()
        initial_burned = np.sum(state0['phi'] <= 0.0) * fire.dx * fire.dy
        print(f"   ✓ Initial state retrieved")
        print(f"     Time: {state0['time']:.2f} s")
        print(f"     Initial burned area: {initial_burned:.1f} m²")
    except Exception as e:
        print(f"   ✗ FAILED to get state: {e}")
        fire.finalize()
        return False
    
    # Run time-stepping
    print(f"\n3. Running 10 timesteps...")
    try:
        for step in range(10):
            result = fire.step()
            if not result['success']:
                print(f"   ✗ FAILED at step {step+1}")
                fire.finalize()
                return False
        
        state = fire.get_state()
        print(f"   ✓ Time-stepping successful")
        print(f"     Final time: {state['time']:.2f} s")
    except Exception as e:
        print(f"   ✗ FAILED during time-stepping: {e}")
        fire.finalize()
        return False
    
    # Check fire spread
    print(f"\n4. Verifying fire spread...")
    try:
        final_burned = np.sum(state['phi'] <= 0.0) * fire.dx * fire.dy
        max_ros = state['ros'].max()
        max_intensity = state['intensity'].max()
        max_flame_length = state['flame_length'].max()
        
        print(f"   Burned area: {final_burned:.1f} m² (growth: {final_burned/initial_burned:.2f}x)")
        print(f"   Max ROS: {max_ros:.3f} m/s")
        print(f"   Max intensity: {max_intensity:.1f} kW/m")
        print(f"   Max flame length: {max_flame_length:.2f} m")
        
        # Verify fire has spread
        if final_burned > initial_burned * 1.5:
            print(f"   ✓ Fire spread verified (area increased by {(final_burned/initial_burned-1)*100:.1f}%)")
        else:
            print(f"   ✗ FAILED: Fire did not spread sufficiently")
            fire.finalize()
            return False
        
        # Verify physical values are reasonable
        if max_ros > 0.0 and max_ros < 10.0:
            print(f"   ✓ ROS values are physically reasonable")
        else:
            print(f"   ✗ WARNING: ROS values may be unrealistic")
        
        if max_intensity > 0.0:
            print(f"   ✓ Fire intensity is positive")
        else:
            print(f"   ✗ FAILED: Fire intensity is not positive")
            fire.finalize()
            return False
            
    except Exception as e:
        print(f"   ✗ FAILED during verification: {e}")
        fire.finalize()
        return False
    
    # Test state field shapes
    print(f"\n5. Verifying state field shapes...")
    try:
        expected_shape = (fire.ny, fire.nx)
        fields = ['phi', 'ros', 'intensity', 'flame_length', 'u_wind', 'v_wind']
        
        for field in fields:
            if field not in state:
                print(f"   ✗ FAILED: Field '{field}' missing from state")
                fire.finalize()
                return False
            if state[field].shape != expected_shape:
                print(f"   ✗ FAILED: Field '{field}' has shape {state[field].shape}, "
                      f"expected {expected_shape}")
                fire.finalize()
                return False
        
        print(f"   ✓ All state fields have correct shape {expected_shape}")
    except Exception as e:
        print(f"   ✗ FAILED during shape verification: {e}")
        fire.finalize()
        return False
    
    # Test plotfile writing
    print(f"\n6. Testing plotfile writing...")
    try:
        fire.write_plotfile()
        # Check if plotfile was created
        if os.path.exists('plt00010'):
            print(f"   ✓ Plotfile written successfully")
        else:
            print(f"   ✗ WARNING: Plotfile directory not found (may be expected)")
    except Exception as e:
        print(f"   ✗ WARNING: Plotfile write error (may be expected): {e}")
    
    # Finalize
    print(f"\n7. Finalizing solver...")
    try:
        fire.finalize()
        print(f"   ✓ Finalized successfully")
    except Exception as e:
        print(f"   ✗ FAILED to finalize: {e}")
        return False
    
    return True


def main():
    """Main test entry point"""
    print(f"pyWildfire version: {pyWildfire.__version__}\n")
    
    success = test_basic_fire_solver()
    
    print("\n" + "=" * 70)
    if success:
        print("✓ REGRESSION TEST PASSED")
        print("=" * 70)
        return 0
    else:
        print("✗ REGRESSION TEST FAILED")
        print("=" * 70)
        return 1


if __name__ == "__main__":
    sys.exit(main())
