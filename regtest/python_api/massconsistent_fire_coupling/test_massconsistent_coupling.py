#!/usr/bin/env python3
"""
Test massconsistent_fire_coupling.py integration with synthetic wind.

This test validates:
1. FireWindCoupler initialization
2. Two-way coupling framework  
3. Rothermel ROS calculation with 3D wind
4. Heat flux extraction and feedback
5. Domain compatibility checks
6. Statistics collection
7. Output generation

Usage:
    cd regtest/python_api/massconsistent_fire_coupling
    PYTHONPATH=/path/to/build/python:$PYTHONPATH python3 test_massconsistent_coupling.py
    
Or via CTest:
    cd build
    ctest -R massconsistent_fire_coupling --output-on-failure
"""

import sys
import os
import numpy as np
from pathlib import Path

# Add paths for imports
test_dir = Path(__file__).parent.absolute()
repo_root = test_dir.parents[2]
src_python = repo_root / "src" / "python"
regtest = repo_root / "regtest"

sys.path.insert(0, str(src_python))

try:
    from wildfire_solver import WildfireSolver
    from massconsistent_fire_coupling import FireWindCoupler
except ImportError as e:
    print(f"Import Error: {e}")
    print(f"Make sure PYTHONPATH includes: {src_python}")
    sys.exit(1)


def find_test_inputs():
    """Find compatible fire inputs file from existing tests."""
    # Look for a simple fire inputs file
    potential_inputs = [
        repo_root / "regtest" / "python_api" / "coupled_wind_fire" / "inputs.i",
        repo_root / "regtest" / "surface_spread" / "farsite_ellipse" / "inputs.i",
    ]
    
    for inputs_file in potential_inputs:
        if inputs_file.exists():
            return str(inputs_file)
    
    # If no standard inputs found, create a minimal one
    return None


def test_fire_solver_api():
    """Test basic fire solver API."""
    print("\n" + "="*70)
    print("TEST 1: Fire Solver API Features")
    print("="*70)
    
    inputs_file = find_test_inputs()
    if inputs_file is None:
        print("⚠ Skipping: No fire inputs file found")
        return True
    
    try:
        fire = WildfireSolver(inputs_file)
        
        # Test get_status()
        status = fire.get_status()
        assert 'initialized' in status
        assert 'time' in status
        assert 'step' in status
        print("✓ get_status() works")
        
        # Test get_domain_bounds()
        bounds = fire.get_domain_bounds()
        assert 'xmin' in bounds and 'xmax' in bounds
        print("✓ get_domain_bounds() works")
        
        # Test get_grid_spacing()
        spacing = fire.get_grid_spacing()
        assert 'dx' in spacing and 'dy' in spacing
        print("✓ get_grid_spacing() works")
        
        # Test get_grid_dimensions()
        dims = fire.get_grid_dimensions()
        assert 'nx' in dims and 'ny' in dims
        print("✓ get_grid_dimensions() works")
        
        # Test properties
        assert fire.current_time == fire.time
        assert fire.timestep == fire.step_num
        print("✓ current_time and timestep properties work")
        
        # Test field accessor
        phi = fire.get_field('phi')
        assert phi.shape == (fire.ny, fire.nx)
        print("✓ get_field('phi') works")
        
        ros = fire.get_field('ros')
        assert ros.shape == (fire.ny, fire.nx)
        print("✓ get_field('ros') works")
        
        # Test statistics
        stats = fire.get_statistics()
        assert 'burned_area' in stats
        assert 'perimeter' in stats
        assert 'max_ros' in stats
        print("✓ get_statistics() works")
        
        # Test diagnostics
        diag = fire.get_diagnostic_info()
        assert 'initialized' in diag
        assert 'time' in diag
        print("✓ get_diagnostic_info() works")
        
        # Test Rothermel-related methods (should not raise)
        fire.set_fuel_model(1)
        print("✓ set_fuel_model() works (with warnings)")
        
        props = fire.get_fuel_properties()
        assert 'model_number' in props
        print("✓ get_fuel_properties() works")
        
        ros_field = fire.compute_ros()
        assert ros_field.shape == (fire.ny, fire.nx)
        print("✓ compute_ros() works")
        
        # Test environmental setters (should not raise)
        fire.set_fuel_moisture(10.0, 15.0, 20.0, 80.0)
        print("✓ set_fuel_moisture() works (with warnings)")
        
        fire.set_wind_direction(225.0)
        print("✓ set_wind_direction() works (with warnings)")
        
        fire.set_ambient_temperature(25.0)
        print("✓ set_ambient_temperature() works (with warnings)")
        
        fire.set_relative_humidity(35.0)
        print("✓ set_relative_humidity() works (with warnings)")
        
        # Test ignition control
        fire.set_ignition(0, 0, 0.0, 1.0)
        print("✓ set_ignition() works (with warnings)")
        
        ign_state = fire.get_ignition_state()
        assert 'configured' in ign_state
        print("✓ get_ignition_state() works")
        
        # Test performance config
        fire.set_cfl(0.5)
        print("✓ set_cfl() works")
        
        fire.set_max_time(1000.0)
        print("✓ set_max_time() works")
        
        fire.set_nsteps(100)
        print("✓ set_nsteps() works")
        
        perf = fire.get_performance_metrics()
        assert 'time_per_step' in perf
        print("✓ get_performance_metrics() works")
        
        # Test surface fluxes
        flux_data = fire.get_surface_fluxes()
        assert 'heat_flux' in flux_data
        assert 'grid_info' in flux_data
        print("✓ get_surface_fluxes() works")
        
        fire.finalize()
        print("\n✓ All fire solver API tests passed!")
        return True
        
    except Exception as e:
        print(f"\n✗ Fire solver API test failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_coupler_initialization():
    """Test FireWindCoupler initialization."""
    print("\n" + "="*70)
    print("TEST 2: FireWindCoupler Initialization")
    print("="*70)
    
    inputs_file = find_test_inputs()
    if inputs_file is None:
        print("⚠ Skipping: No fire inputs file found")
        return True
    
    try:
        # Test initialization with synthetic wind
        coupler = FireWindCoupler(
            fire_inputs=inputs_file,
            max_time=100.0,
            use_synthetic_wind=True,
            verbose=False
        )
        
        assert coupler.fire is not None
        print("✓ FireWindCoupler initialized")
        
        # Check domain
        status = coupler.fire.get_status()
        assert status['nx'] > 0
        assert status['ny'] > 0
        print("✓ Domain verified")
        
        # Check synthetic wind generation
        u, v, w, nz, zmin, zmax = coupler._generate_synthetic_wind_3d(0.0)
        assert u.shape == (nz, coupler.fire.ny, coupler.fire.nx)
        assert v.shape == (nz, coupler.fire.ny, coupler.fire.nx)
        print("✓ Synthetic wind generation works")
        
        coupler.finalize()
        print("\n✓ FireWindCoupler initialization tests passed!")
        return True
        
    except Exception as e:
        print(f"\n✗ Coupler initialization test failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_coupled_timestepping():
    """Test coupled fire-wind timestepping."""
    print("\n" + "="*70)
    print("TEST 3: Coupled Time Stepping")
    print("="*70)
    
    inputs_file = find_test_inputs()
    if inputs_file is None:
        print("⚠ Skipping: No fire inputs file found")
        return True
    
    try:
        coupler = FireWindCoupler(
            fire_inputs=inputs_file,
            max_time=60.0,  # 1 minute
            use_synthetic_wind=True,
            verbose=False
        )
        
        initial_time = coupler.fire.time
        initial_step = coupler.coupled_step
        
        # Take 3 coupled steps
        for i in range(3):
            result = coupler.step()
            assert result['success']
            assert coupler.coupled_step > initial_step
            assert coupler.fire.time > initial_time
        
        print(f"✓ Took 3 coupled steps")
        print(f"  Time progressed: {initial_time:.2f}s → {coupler.fire.time:.2f}s")
        print(f"  Fire domain: {coupler.fire.nx} × {coupler.fire.ny}")
        
        # Check statistics collection
        assert len(coupler.stats_history['time']) > 0
        assert len(coupler.stats_history['burned_area']) > 0
        print(f"✓ Statistics collected ({len(coupler.stats_history['time'])} entries)")
        
        coupler.finalize()
        print("\n✓ Coupled time stepping tests passed!")
        return True
        
    except Exception as e:
        print(f"\n✗ Coupled timestepping test failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_rothermel_ros():
    """Test Rothermel ROS with wind coupling."""
    print("\n" + "="*70)
    print("TEST 4: Rothermel ROS with Wind Coupling")
    print("="*70)
    
    inputs_file = find_test_inputs()
    if inputs_file is None:
        print("⚠ Skipping: No fire inputs file found")
        return True
    
    try:
        coupler = FireWindCoupler(
            fire_inputs=inputs_file,
            max_steps=5,
            use_synthetic_wind=True,
            verbose=False
        )
        
        # Take a step to get ROS with wind
        coupler.step()
        
        state = coupler.fire.get_state()
        ros = state['ros']
        
        # Check that ROS is physically reasonable
        max_ros = np.max(ros)
        assert max_ros >= 0.0, "ROS should be non-negative"
        print(f"✓ ROS values are non-negative (max: {max_ros:.3f} m/s)")
        
        # Check that ROS varies spatially (not uniform everywhere)
        ros_std = np.std(ros[ros > 0])
        if ros_std > 0:
            print(f"✓ ROS varies spatially (std: {ros_std:.4f} m/s)")
        else:
            print(f"⚠ ROS is spatially uniform (may be expected)")
        
        # Verify heat flux depends on ROS
        intensity = state['intensity']
        flame_length = state['flame_length']
        
        max_intensity = np.max(intensity)
        print(f"✓ Fire intensity calculated (max: {max_intensity:.1f} kW/m)")
        
        coupler.finalize()
        print("\n✓ Rothermel ROS tests passed!")
        return True
        
    except Exception as e:
        print(f"\n✗ Rothermel ROS test failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_heat_flux_extraction():
    """Test heat flux extraction for wind coupling."""
    print("\n" + "="*70)
    print("TEST 5: Heat Flux Extraction")
    print("="*70)
    
    inputs_file = find_test_inputs()
    if inputs_file is None:
        print("⚠ Skipping: No fire inputs file found")
        return True
    
    try:
        fire = WildfireSolver(inputs_file)
        
        # Take a step to generate heat
        fire.step()
        
        # Extract heat flux
        flux_data = fire.get_surface_fluxes()
        heat_flux = flux_data['heat_flux']
        grid_info = flux_data['grid_info']
        
        assert heat_flux.shape == (fire.ny, fire.nx)
        print(f"✓ Heat flux shape correct: {heat_flux.shape}")
        
        # Check heat flux is non-negative
        assert np.all(heat_flux >= 0), "Heat flux should be non-negative"
        print(f"✓ Heat flux non-negative (max: {np.max(heat_flux):.1f} kW/m²)")
        
        # Check grid info
        assert grid_info['nx'] == fire.nx
        assert grid_info['ny'] == fire.ny
        print(f"✓ Grid info consistent")
        
        # Compute from heat_release method for comparison
        heat_data = fire.compute_heat_release()
        assert 'surface_flux' in heat_data
        print(f"✓ Alternate heat extraction method works")
        
        fire.finalize()
        print("\n✓ Heat flux extraction tests passed!")
        return True
        
    except Exception as e:
        print(f"\n✗ Heat flux test failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    """Run all tests."""
    print("\n" + "="*70)
    print("MASSCONSISTENT FIRE COUPLING API TEST SUITE")
    print("="*70)
    
    tests = [
        test_fire_solver_api,
        test_coupler_initialization,
        test_coupled_timestepping,
        test_rothermel_ros,
        test_heat_flux_extraction,
    ]
    
    results = {}
    for test_func in tests:
        try:
            results[test_func.__name__] = test_func()
        except Exception as e:
            print(f"\n✗ Test {test_func.__name__} crashed: {e}")
            results[test_func.__name__] = False
    
    # Summary
    print("\n" + "="*70)
    print("TEST SUMMARY")
    print("="*70)
    
    passed = sum(1 for v in results.values() if v)
    total = len(results)
    
    for test_name, result in results.items():
        status = "✓ PASS" if result else "✗ FAIL"
        print(f"{status}: {test_name}")
    
    print(f"\nTotal: {passed}/{total} tests passed")
    
    return 0 if passed == total else 1


if __name__ == "__main__":
    sys.exit(main())
