# Fire Acceleration Model Implementation Summary

## Overview

This implementation adds the full FARSITE temporal fire acceleration model (McAlpine & Wakimoto 1991) to the wildfire_levelset toolkit, including wind-onset time-lag capability for modeling fire intensification after wind changes.

## Research Findings

Extensive research revealed that:
1. **O'Neal/Grala** does NOT appear in FARSITE literature or source code
2. The actual FARSITE model is based on **McAlpine & Wakimoto (1991) / VanWagner's equation**
3. The model uses **R(t) = R_E × (1 - exp(-A × t))** with per-cell temporal tracking
4. Two acceleration constants: A_point = 0.115 1/min (small fires), A_line = 0.300 1/min (large fires)
5. The current implementation used a simplified size-based model (Catchpole et al. 1992)

## Implementation Details

### Models Implemented

1. **Size-Based Model** (Catchpole et al. 1992) - EXISTING
   - Global scaling factor based on fire size
   - α = 1 - exp(-r_fire / L_acc)
   - Simple and efficient

2. **FARSITE Temporal Model** (McAlpine & Wakimoto 1991) - NEW
   - Per-cell temporal tracking
   - R(t) = R_E × (1 - exp(-A × t))
   - Automatic point/line constant switching

3. **Wind-Onset Time-Lag** - NEW EXTENSION
   - Optional wind response delay
   - R_target(t) = R_prev + (R_E - R_prev) × (1 - exp(-dt/tau_wind))
   - Realistic fire response to sudden wind changes

### Code Changes

**New/Modified Files:**
1. `src/fire_acceleration.H` - Complete rewrite with both models
2. `src/parse_inputs.H` - Extended AccelerationParams
3. `src/parse_inputs.cpp` - Parameter initialization and validation
4. `src/multifab_setup.H` - Added accel_state_mf (3 components per cell)
5. `src/main.cpp` - Initialize state and pass dt to acceleration
6. `docs/usage.rst` - Comprehensive documentation
7. `FARSITE_FEATURES_IMPLEMENTATION.md` - Feature summary

**New Test Suite:**
- `regtest/surface_spread/fire_acceleration/`
  - inputs.size_based - Size-based model test
  - inputs.temporal_point - FARSITE temporal test
  - inputs.wind_lag - Wind-lag test
  - run_regtest.py - Automated test runner

### Parameters Added

```
acceleration.enable = 1                # Enable model
acceleration.use_temporal = 1          # 0=size-based, 1=FARSITE
acceleration.L_acc = 50.0              # Size-based length scale [m]
acceleration.A_point = 0.115           # Point constant [1/min]
acceleration.A_line = 0.300            # Line constant [1/min]
acceleration.perim_limit = 402.3       # Perimeter threshold [m]
acceleration.enable_wind_lag = 1       # Enable wind-lag
acceleration.tau_wind = 180.0          # Wind time constant [s]
```

## Validation Plan

1. **Build Test**: Verify code compiles with AMReX
2. **Size-Based Test**: Confirm backward compatibility
3. **Temporal Point Test**: Verify ~20 min to 90% equilibrium
4. **Temporal Line Test**: Verify ~8 min to 90% equilibrium (faster than point)
5. **Wind-Lag Test**: Confirm ROS lags behind equilibrium after wind step

## Expected Behavior

- **Point fires**: Slow acceleration (A = 0.115 1/min)
- **Line fires**: Fast acceleration (A = 0.300 1/min)
- **Wind increases**: Exponential ramp-up with tau_wind time constant
- **Wind decreases**: Immediate adjustment to new equilibrium
- **Large fires**: Minimal impact (already at steady-state)

## References

1. McAlpine, R.S. & Wakimoto, R.H. (1991). The acceleration of fire from point source to equilibrium spread. Forest Science, 37(5), 1314–1337.
2. Alexander, M.E., Stocks, B.J. & Lawson, B.D. (1992). Fire behavior in Black Spruce-lichen woodland. Info. Rep. NOR-X-310. USDA/CFS.
3. Catchpole, E.A., de Mestre, N.J. & Gill, A.M. (1992). Intensity of fire at its perimeter. Australian Journal of Ecology, 17(1), 1–4.
4. Finney, M.A. (1998/2004). FARSITE: Fire Area Simulator. USDA Forest Service RMRS-RP-4.

## Status

✅ Implementation complete
✅ Documentation updated
✅ Regression tests created
⏳ Build and runtime validation pending (requires AMReX)
