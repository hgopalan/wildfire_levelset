#!/usr/bin/env python3
"""
test_two_way_coupling_mvp.py - Regression test for two-way coupling MVP

Tests:
1. Surface flux extraction from fire state
2. Wind solver interface implementation
3. Two-way coupled time loop
4. Heat feedback effects (wind modification)
"""

import sys
import os
import numpy as np

# Find and import pyWildfire
try:
    import pyWildfire
    from wildfire_solver import WildfireSolver
except ImportError:
    build_python = os.path.join(os.getcwd(), '..', '..', '..', 'build', 'python')
    if os.path.exists(build_python):
        sys.path.insert(0, build_python)
        import pyWildfire
        from wildfire_solver import WildfireSolver
    else:
        print("ERROR: Cannot find pyWildfire module")
        sys.exit(1)

# Import coupling modules
from surface_flux_extractor import SurfaceFluxExtractor
from wind_solver_interface import SyntheticWindSolver


def test_surface_flux_extraction():
    """Test 1: Surface flux extraction."""
    print("=" * 70)
    print("Test 1: Surface Flux Extraction")
    print("=" * 70)

    # Create mock fire state
    ny, nx = 32, 32
    dx, dy = 30.0, 30.0

    mock_state = {
        'intensity': np.random.exponential(100.0, (ny, nx)) * (np.random.rand(ny, nx) > 0.8),
        'flame_length': np.ones((ny, nx)) * 5.0,
        'phi': np.random.randn(ny, nx) * 50.0 - 20.0,
    }

    # Extract fluxes
    extractor = SurfaceFluxExtractor((ny, nx), dx, dy)
    fluxes = extractor.get_flux_dict(mock_state)

    # Verify outputs
    assert 'heat_flux' in fluxes, "Missing heat_flux key"
    assert 'buoyancy' in fluxes, "Missing buoyancy key"
    assert 'fire_footprint' in fluxes, "Missing fire_footprint key"
    assert 'total_heat_MW' in fluxes, "Missing total_heat_MW key"

    assert fluxes['heat_flux'].shape == (ny, nx), f"Wrong heat_flux shape: {fluxes['heat_flux'].shape}"
    assert fluxes['buoyancy'].shape == (ny, nx), f"Wrong buoyancy shape: {fluxes['buoyancy'].shape}"

    # Check values are reasonable
    assert np.all(fluxes['heat_flux'] >= 0.0), "Heat flux has negative values"
    assert np.all(fluxes['fire_footprint'] >= 0.0), "Fire footprint has negative values"
    assert np.all(fluxes['fire_footprint'] <= 1.0), "Fire footprint has values > 1"

    total_heat = fluxes['total_heat_MW']
    assert total_heat >= 0.0, f"Total heat is negative: {total_heat}"

    print(f"  ✓ Flux shapes correct")
    print(f"  ✓ Heat flux range: {fluxes['heat_flux'].min():.1f} - {fluxes['heat_flux'].max():.1f} W/m²")
    print(f"  ✓ Total heat: {total_heat:.2f} MW")
    print(f"  ✓ Active fire fraction: {fluxes['fire_footprint'].mean()*100:.1f}%")
    print()
    return True


def test_wind_solver_interface():
    """Test 2: Wind solver interface."""
    print("=" * 70)
    print("Test 2: Wind Solver Interface")
    print("=" * 70)

    # Create wind solver
    wind = SyntheticWindSolver()
    assert wind.initialize(), "Wind solver initialization failed"
    print("  ✓ Wind solver initialized")

    # Get domain info
    info = wind.get_domain_info()
    assert 'nx' in info and 'ny' in info and 'nz' in info, "Missing domain info"
    print(f"  ✓ Domain info: {info['nx']}×{info['ny']}×{info['nz']} grid")

    # Solve without heat
    assert wind.solve(0.0), "Wind solver solve failed"
    u, v, w = wind.get_velocity_arrays()
    assert u.shape == (info['nz'], info['ny'], info['nx']), f"Wrong u shape: {u.shape}"
    assert v.shape == (info['nz'], info['ny'], info['nx']), f"Wrong v shape: {v.shape}"
    assert w.shape == (info['nz'], info['ny'], info['nx']), f"Wrong w shape: {w.shape}"
    print(f"  ✓ Wind arrays shape correct: (nz, ny, nx) = {u.shape}")

    # Check wind values
    u_mean_no_heat = np.abs(u).mean()
    print(f"  ✓ Wind without heat: |u|_mean={u_mean_no_heat:.2f} m/s")

    # Solve with heat
    heat_source = {
        'heat_flux': np.ones((info['ny'], info['nx'])) * 50000.0,
        'buoyancy': np.ones((info['ny'], info['nx'])) * 5.0,
        'fire_footprint': np.ones((info['ny'], info['nx'])),
    }
    assert wind.solve(0.0, heat_source=heat_source), "Wind solver solve with heat failed"
    u2, v2, w2 = wind.get_velocity_arrays()
    print(f"  ✓ Wind solver handles heat source")

    # Verify wind was modified by heat
    u_mean_with_heat = np.abs(u2).mean()
    print(f"  ✓ Wind with heat: |u|_mean={u_mean_with_heat:.2f} m/s")
    assert u_mean_with_heat < u_mean_no_heat, "Heat should reduce horizontal wind speed"
    print(f"  ✓ Heat feedback applied: {(1-u_mean_with_heat/u_mean_no_heat)*100:.1f}% reduction")

    assert wind.finalize(), "Wind solver finalization failed"
    print("  ✓ Wind solver finalized")
    print()
    return True


def test_two_way_coupled_loop():
    """Test 3: Two-way coupled time loop."""
    print("=" * 70)
    print("Test 3: Two-Way Coupled Time Loop")
    print("=" * 70)

    # Locate inputs file
    inputs_file = "inputs.i"
    if not os.path.exists(inputs_file):
        print(f"  ⚠ Skipping: {inputs_file} not found in current directory")
        print(f"    Current directory: {os.getcwd()}")
        return True  # Skip but don't fail

    try:
        # Initialize fire solver
        fire = WildfireSolver(inputs_file)
        print(f"  ✓ Fire solver initialized: {fire.nx}×{fire.ny}")

        # Initialize wind solver
        wind = SyntheticWindSolver()
        wind.initialize()
        print(f"  ✓ Wind solver initialized: {wind.nx}×{wind.ny}×{wind.nz}")

        # Initialize flux extractor
        extractor = SurfaceFluxExtractor((fire.ny, fire.nx), fire.dx, fire.dy)
        print(f"  ✓ Flux extractor initialized")

        # Run 5 coupled steps
        num_steps = 5
        print(f"\n  Running {num_steps} coupled steps...")

        for step in range(num_steps):
            # Get fire state
            state = fire.get_state()

            # Extract and pass fluxes to wind
            fluxes = extractor.get_flux_dict(state)
            wind.solve(state['time'], heat_source=fluxes)

            # Get new wind and update fire
            u_3d, v_3d, w_3d = wind.get_velocity_arrays()
            wind_info = wind.get_domain_info()

            fire.update_wind_3d(u_3d, v_3d, w_3d, wind_info['nz'],
                               wind_info['zmin'], wind_info['zmax'])

            # Advance fire
            fire.step()

            # Print progress
            state = fire.get_state()
            phi = state['phi']
            burned = np.sum(phi <= 0.0) * fire.dx * fire.dy
            print(f"    Step {step+1}: t={state['time']:.2f} s, "
                  f"burned={burned:.1f} m², heat={fluxes['total_heat_MW']:.2f} MW")

        fire.finalize()
        wind.finalize()
        print(f"  ✓ Coupled loop completed successfully")
        print()
        return True

    except Exception as e:
        print(f"  ✗ Coupled loop failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_heat_feedback_effects():
    """Test 4: Verify heat feedback modifies wind."""
    print("=" * 70)
    print("Test 4: Heat Feedback Effects")
    print("=" * 70)

    wind = SyntheticWindSolver()
    wind.initialize()

    # Reference wind (no heat)
    wind.solve(0.0, heat_source=None)
    u_ref, v_ref, w_ref = wind.get_velocity_arrays()
    u_magnitude_ref = np.sqrt(u_ref**2 + v_ref**2)

    # Wind with varying heat intensities
    heat_levels = [0.0, 0.25, 0.5, 1.0]  # Fraction of max heat

    print(f"  Testing wind response to varying heat intensity...")
    print(f"  {'Heat Level':<15} {'|U| Mean':<15} {'Change':<15}")
    print(f"  {'-'*45}")

    prev_u_mag = u_magnitude_ref.mean()

    for heat_frac in heat_levels:
        if heat_frac == 0.0:
            heat_source = None
        else:
            max_heat = 100000.0  # 100 kW/m²
            heat_source = {
                'heat_flux': np.ones((wind.ny, wind.nx)) * max_heat * heat_frac,
                'buoyancy': np.ones((wind.ny, wind.nx)) * heat_frac,
                'fire_footprint': np.ones((wind.ny, wind.nx)),
            }

        wind.solve(0.0, heat_source=heat_source)
        u, v, w = wind.get_velocity_arrays()
        u_magnitude = np.sqrt(u**2 + v**2).mean()

        if heat_frac == 0.0:
            change_pct = 0.0
        else:
            change_pct = (1.0 - u_magnitude / u_magnitude_ref.mean()) * 100

        print(f"  {heat_frac*100:>5.0f}%          {u_magnitude:>6.3f} m/s        "
              f"{change_pct:>+6.1f}%")

        # Verify monotonicity (more heat → smaller wind)
        if heat_frac > 0.0:
            assert u_magnitude < prev_u_mag, \
                f"Wind should decrease with heat, but {u_magnitude} >= {prev_u_mag}"
            prev_u_mag = u_magnitude

    wind.finalize()
    print(f"  ✓ Heat feedback effects verified (monotonic wind reduction)")
    print()
    return True


def main():
    """Run all tests."""
    print(f"\nRegression Test Suite: Two-Way Coupling MVP")
    print(f"{'='*70}\n")

    tests = [
        ("Surface Flux Extraction", test_surface_flux_extraction),
        ("Wind Solver Interface", test_wind_solver_interface),
        ("Two-Way Coupled Loop", test_two_way_coupled_loop),
        ("Heat Feedback Effects", test_heat_feedback_effects),
    ]

    passed = 0
    failed = 0

    for name, test_func in tests:
        try:
            if test_func():
                passed += 1
            else:
                failed += 1
                print(f"✗ {name} FAILED\n")
        except Exception as e:
            failed += 1
            print(f"✗ {name} FAILED with exception: {e}\n")
            import traceback
            traceback.print_exc()

    # Summary
    print("=" * 70)
    print(f"Test Summary: {passed} passed, {failed} failed")
    print("=" * 70)

    if failed == 0:
        print("✓ ALL TESTS PASSED")
        return 0
    else:
        print(f"✗ {failed} TEST(S) FAILED")
        return 1


if __name__ == "__main__":
    sys.exit(main())
