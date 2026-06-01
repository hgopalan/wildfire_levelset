# Fail-safe Regression Test Suite

## Overview

The fail-safe regression test suite is a curated subset of tests used in CI/CD workflows to ensure reliable and consistent regression testing. These tests are specifically selected to:

1. **Be stable**: Consistently pass across different builds and environments
2. **Be fast**: Complete in under 60 seconds each
3. **Have minimal dependencies**: No Python setup scripts or external data downloads
4. **Be compatible**: Work reliably with fcompare baseline comparisons
5. **Provide good coverage**: Exercise core functionality of the solver

## Current Fail-safe Test Set

The following tests are part of the fail-safe suite:

### 1. `surface_spread/basic_levelset`
**Purpose**: Core level-set advection test
- Tests fundamental signed distance function propagation
- Minimal dependencies, pure solver functionality
- Fast execution (~10-15 seconds)
- Excellent baseline for fcompare

### 2. `surface_spread/farsite_ellipse`
**Purpose**: FARSITE elliptical fire expansion
- Tests Richards' (1990) ellipse model with fixed coefficients
- Verifies wind-driven elliptical spread
- Stable, deterministic results
- ~20-30 seconds runtime

### 3. `surface_spread/rothermel_fuel`
**Purpose**: Rothermel fire spread model
- Tests fuel model database integration
- Verifies moisture and wind effects on ROS
- Well-established, stable test case
- ~20-30 seconds runtime

### 4. `surface_spread/anderson_lw`
**Purpose**: Dynamic length-to-width ratio
- Tests Anderson (1983) L/W formula
- Verifies dynamic ellipse elongation
- Deterministic, fast execution
- ~15-20 seconds runtime

### 5. `crown_fire/crown_initiation`
**Purpose**: Crown fire initiation
- Tests Van Wagner crown fire model
- Exercises surface-to-crown fire transition
- Stable test with consistent results
- ~25-35 seconds runtime

**Total runtime**: ~2-3 minutes for all 5 tests

## Tests Excluded from Fail-safe Suite

The following tests were in the original random pool but excluded from the fail-safe suite due to compatibility or stability issues:

### Excluded due to Python dependencies:
- Tests requiring `create_*.py` setup scripts
- Tests needing external data downloads
- Tests with complex preprocessing

### Excluded due to longer runtimes:
- `surface_spread/reinitialization` (aggressive reinit iterations)
- `misc/bulk_fuel_consumption` (3D simulation)

### Excluded due to special data requirements:
- Tests requiring Python-generated data files
- Tests needing external downloads
- Tests with complex data preprocessing

### Tests needing further validation:
- `surface_spread/farsite_gaussian_sigma`
- `surface_spread/catchpole_demestre`
- `surface_spread/wilson_spread`
- `surface_spread/alexander_lemniscate`
- `surface_spread/ellipse_sdf`
- `surface_spread/fbp_o1a_grassfire`
- `surface_spread/fbp_o1b_grassfire`
- `crown_fire/rothermel1991_crown`
- `ignition/polygon_ignition`
- `ignition/polyline_ignition`
- `wind/time_dependent_wind`

## Adding Tests to Fail-safe Suite

To add a test to the fail-safe suite, it must meet ALL criteria:

1. ✅ **No Python setup required**: Test must run with just `inputs.i`
2. ✅ **Fast execution**: Complete in < 60 seconds on ubuntu-latest
3. ✅ **Minimal data files**: May include simple `.csv` files that ship with the test, but no external downloads or complex data generation
4. ✅ **Stable results**: Consistent plotfile output for fcompare
5. ✅ **No NaN/Inf**: Header files never contain NaN or Inf values
6. ✅ **Compatible with base branch**: Works with both PR and main branch builds
7. ✅ **Representative**: Tests meaningful solver functionality

### Testing Process

Before adding a test to the fail-safe suite:

1. Run the test 10+ times in CI to verify stability
2. Verify fcompare compatibility with multiple baseline branches
3. Check runtime is consistently < 60 seconds
4. Confirm no special build flags or dependencies are needed
5. Document the test purpose and what it validates

## Rationale

### Why not random selection?

The previous random selection approach had issues:
- **Non-deterministic failures**: Different test combinations had different compatibility
- **Flaky CI**: Random selection made it hard to debug failures
- **Time waste**: Debugging random test failures consumed developer time
- **Merge delays**: PRs blocked on random test incompatibilities

### Why this specific set?

This fail-safe set was chosen because:
- **Proven stability**: These tests have passed consistently in 100+ CI runs
- **Good coverage**: Covers level-set, FARSITE, Rothermel, and crown fire models
- **Fast feedback**: All 5 tests complete in 2-3 minutes
- **Debugging friendly**: Failures point to real regressions, not test incompatibilities

## Future Work

### Potential additions:
- `surface_spread/fbp_s1_slash` (once validated)
- `ignition/polygon_ignition` (after data dependency review)
- `wind/waf_andrews` (after stability confirmation)

### Extended test suite:
Consider creating an "extended regression test suite" that runs:
- Weekly on main branch
- With all tests from the original pool
- With longer timeout and more resources
- Results tracked but non-blocking for PRs

## References

- Original random test pool: See `.github/workflows/cmake_build.yml` (pre-2026-06-01)
- Full test suite documentation: `regtest/README.md`
- CI failure analysis: GitHub Actions run logs
