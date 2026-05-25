#!/usr/bin/env python3
"""
test_pywildfire.py - Test script for pyWildfire Python bindings

Tests the load_wind_from_arrays function with synthetic 3D wind data.
This demonstrates how to use the Python interface to load wind data
that can be used by wildfire_levelset.

Run after building with Python bindings enabled:
    cmake -S . -B build -DLEVELSET_BUILD_PYTHON_BINDINGS=ON -DLEVELSET_DIM_2D=ON
    cmake --build build
    PYTHONPATH=build/python python3 src/python/test_pywildfire.py
"""

import sys
import numpy as np

try:
    import pyWildfire
except ImportError as e:
    print(f"Error: Could not import pyWildfire module")
    print(f"  {e}")
    print("\nMake sure to:")
    print("  1. Build with -DLEVELSET_BUILD_PYTHON_BINDINGS=ON")
    print("  2. Set PYTHONPATH to point to build/python directory")
    print("  Example: PYTHONPATH=build/python python3 src/python/test_pywildfire.py")
    sys.exit(1)

def test_basic_wind_loading():
    """Test basic wind data loading from numpy arrays"""
    print("=" * 70)
    print("Test 1: Basic wind loading from numpy arrays")
    print("=" * 70)
    
    # Grid dimensions (same as generate_plt_wind.py in regression tests)
    nx, ny, nz = 8, 8, 4
    
    # Domain bounds (UTM Zone 11N, Southern California)
    xmin, xmax = 329900.0, 330500.0
    ymin, ymax = 3774900.0, 3775500.0
    zmin, zmax = 0.0, 40.0
    
    print(f"\nGrid: {nx} × {ny} × {nz}")
    print(f"Domain: X=[{xmin}, {xmax}], Y=[{ymin}, {ymax}], Z=[{zmin}, {zmax}]")
    
    # Create uniform 3D wind field (Fortran order: k varies fastest)
    # Shape should be (nz, ny, nx) but flatten to 1D array
    u_wind = 5.0  # m/s westerly
    v_wind = 0.5  # m/s southerly
    w_wind = 0.0  # m/s (no vertical motion)
    
    # Create 3D arrays
    u_3d = np.full((nz, ny, nx), u_wind, dtype=np.float64)
    v_3d = np.full((nz, ny, nx), v_wind, dtype=np.float64)
    w_3d = np.full((nz, ny, nx), w_wind, dtype=np.float64)
    
    print(f"\nUniform wind: u={u_wind} m/s, v={v_wind} m/s, w={w_wind} m/s")
    
    # Flatten to 1D in Fortran (column-major) order
    u_flat = u_3d.flatten('F')
    v_flat = v_3d.flatten('F')
    w_flat = w_3d.flatten('F')
    
    print(f"Array sizes: {len(u_flat)} points (expected: {nx*ny*nz})")
    
    # Load wind data
    try:
        result = pyWildfire.load_wind_from_arrays(
            nx, ny, nz,
            xmin, xmax, ymin, ymax, zmin, zmax,
            u_flat, v_flat, w_flat
        )
        
        print("\n✓ Successfully loaded wind data")
        print(f"  Valid: {result['valid']}")
        print(f"  2D grid: {result['nx_2d']} × {result['ny_2d']}")
        print(f"  Total 3D points: {result['n_points']}")
        print(f"  2D domain: X=[{result['xmin']:.1f}, {result['xmax']:.1f}], "
              f"Y=[{result['ymin']:.1f}, {result['ymax']:.1f}]")
        
        # Check column-averaged wind
        u2d = result['u2d']
        v2d = result['v2d']
        
        print(f"\n2D column-averaged wind:")
        print(f"  u2d shape: {u2d.shape}")
        print(f"  u2d mean: {u2d.mean():.4f} m/s (expected: {u_wind:.4f})")
        print(f"  v2d mean: {v2d.mean():.4f} m/s (expected: {v_wind:.4f})")
        
        # Verify accuracy
        u_error = abs(u2d.mean() - u_wind)
        v_error = abs(v2d.mean() - v_wind)
        
        if u_error < 1e-6 and v_error < 1e-6:
            print("\n✓ Test PASSED: Column-averaged wind matches input")
            return True
        else:
            print(f"\n✗ Test FAILED: Errors: u={u_error:.2e}, v={v_error:.2e}")
            return False
            
    except Exception as e:
        print(f"\n✗ Test FAILED with exception: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_variable_wind_field():
    """Test with a variable wind field (height-dependent)"""
    print("\n" + "=" * 70)
    print("Test 2: Variable wind field (height-dependent)")
    print("=" * 70)
    
    nx, ny, nz = 8, 8, 4
    xmin, xmax = 329900.0, 330500.0
    ymin, ymax = 3774900.0, 3775500.0
    zmin, zmax = 0.0, 40.0
    
    # Create height-dependent wind (log profile approximation)
    # u increases with height
    u_3d = np.zeros((nz, ny, nx), dtype=np.float64)
    v_3d = np.zeros((nz, ny, nx), dtype=np.float64)
    
    for k in range(nz):
        z_height = zmin + (k + 0.5) * (zmax - zmin) / nz
        # Simple linear increase with height
        u_k = 3.0 + 2.0 * (z_height / zmax)  # 3 m/s at surface, 5 m/s at top
        v_k = 0.2 + 0.3 * (z_height / zmax)  # 0.2 m/s at surface, 0.5 m/s at top
        u_3d[k, :, :] = u_k
        v_3d[k, :, :] = v_k
    
    # Expected averages
    u_expected = u_3d.mean()
    v_expected = v_3d.mean()
    
    print(f"\nVariable wind field:")
    print(f"  u range: [{u_3d.min():.2f}, {u_3d.max():.2f}] m/s")
    print(f"  v range: [{v_3d.min():.2f}, {v_3d.max():.2f}] m/s")
    print(f"  Expected averages: u={u_expected:.3f}, v={v_expected:.3f}")
    
    # Load without w component (test optional argument)
    try:
        result = pyWildfire.load_wind_from_arrays(
            nx, ny, nz,
            xmin, xmax, ymin, ymax, zmin, zmax,
            u_3d.flatten('F'), v_3d.flatten('F')  # No w_array
        )
        
        print("\n✓ Successfully loaded variable wind data (without w component)")
        
        u2d_mean = result['u2d'].mean()
        v2d_mean = result['v2d'].mean()
        
        print(f"  2D averages: u={u2d_mean:.3f}, v={v2d_mean:.3f}")
        
        u_error = abs(u2d_mean - u_expected)
        v_error = abs(v2d_mean - v_expected)
        
        if u_error < 0.01 and v_error < 0.01:
            print("\n✓ Test PASSED: Column-averaged wind correct for variable field")
            return True
        else:
            print(f"\n✗ Test FAILED: Errors: u={u_error:.3f}, v={v_error:.3f}")
            return False
            
    except Exception as e:
        print(f"\n✗ Test FAILED with exception: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    """Run all tests"""
    print(f"pyWildfire version: {pyWildfire.__version__}\n")
    
    results = []
    results.append(test_basic_wind_loading())
    results.append(test_variable_wind_field())
    
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
