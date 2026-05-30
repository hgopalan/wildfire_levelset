# Implementation Summary: 10 New Wildfire Features

This document summarizes the implementation of 10 new wildfire features based on peer-reviewed literature and operational fire codes.

## Overview

All features have been implemented with minor code changes (50-200 lines each), following the established code patterns. The core algorithms are complete and ready for integration into the main simulation loop.

## Completed Features

### 1. Fuel Continuity Factor (Heterogeneous Fuel Distributions)
**Files Modified:**
- `src/fuel_database.H` - Added `continuity` field to FuelModel struct
- `src/parse_inputs.H` - Added `fuel_continuity` to RothermelParams

**Implementation:**
- Continuity factor (0-1) multiplies computed ROS
- 1.0 = continuous fuel bed (default)
- 0.5 = 50% fuel coverage (patchy fuels)
- 0.0 = no fuel (gaps)

**Literature:**
- Finney, M.A. (2006). FlamMap capabilities. USDA Forest Service RMRS-P-41.
- Parsons, R.A., et al. (2017). Ecological Modelling, 222(3):679-691.

**Integration Needed:** Apply continuity factor in ROS computation in main loop

---

### 2. NFDRS Fire Danger Class (Operational Categories)
**Files Modified:**
- `src/fire_intensity_class.H` - Added NFDRS classification functions

**Implementation:**
- `classify_nfdrs_danger()` - 5 classes: Low/Moderate/High/Very High/Extreme
- `compute_nfdrs_danger_class()` - MultiFab version
- Based on fireline intensity thresholds

**Literature:**
- Deeming, J.E., et al. (1977). NFDRS-1978. USDA GTR INT-39.
- Andrews, P.L. & Bradshaw, L.S. (1997). GTR INT-367.

**Integration Needed:** Add to plotfile output

---

### 3. Crown Fraction Burned (CFB Diagnostic)
**Files Modified:**
- `src/crown_initiation.H` - Added CFB calculation function

**Implementation:**
- `compute_crown_fraction_burned()` function
- CFB = (R - R_surface) / (R_active_crown - R_surface)
- Distinguishes passive vs active crown fire

**Literature:**
- Scott, J.H. & Reinhardt, E.D. (2001). USDA RMRS-RP-29.
- Cruz, M.G., et al. (2005). Can. J. Forest Res., 35(7):1626-1639.

**Integration Needed:** Add to output variables

---

### 4. Effective Wind Speed (Wind + Slope Combined)
**Files Created:**
- `src/effective_wind_speed.H` - Complete implementation

**Implementation:**
- Vector combination: U_eff = sqrt(U_wind² + U_slope_equiv²)
- Slope converted to equivalent wind speed
- `compute_effective_wind_speed_field()` for MultiFab

**Literature:**
- Rothermel, R.C. (1983). GTR INT-143.
- Sharples, J.J. (2008). Int. J. Wildland Fire, 17(1):1-17.

**Integration Needed:** Add to plotfile output

---

### 5. Thomas Flame Length Model (Alternative to Byram)
**Files Modified:**
- `src/compute_fire_behavior.H` - Added Thomas formula option

**Implementation:**
- Thomas (1963): L = 0.0266 * I^0.667
- Byram (1959): L = 0.0775 * I^0.46 (default)
- Selectable via `flame_length_model` parameter

**Literature:**
- Thomas, P.H. (1963). Combustion Symposium, 9(1):844-859.
- Nelson, R.M. Jr. (1980). USDA Research Paper SE-205.

**Integration Needed:** Add parameter to parse_inputs.cpp

---

### 6. Fuel Boundary Smoothing (Eliminate Discontinuities)
**Files Created:**
- `src/fuel_boundary_smoothing.H` - Complete implementation

**Implementation:**
- Distance-weighted ROS blending at fuel boundaries
- 2-3 cell transition zone
- Uses existing two_fuel_blending methods
- `apply_fuel_boundary_smoothing()` function

**Literature:**
- Finney, M.A. (1998). USDA RMRS-RP-4.
- Finney, M.A. (2006). USDA RMRS-P-41.

**Integration Needed:** Call in main loop after ROS computation

---

### 7. CSIRO Grassfire Acceleration (Non-Equilibrium Growth)
**Files Modified:**
- `src/cheney_gould_model.H` - Added acceleration functions

**Implementation:**
- `csiro_grassfire_acceleration()` - a(t) = 1 - exp(-t/t_accel)
- `cheney_gould_ros_with_acceleration()` - combined function
- Default t_accel = 600s (10 minutes)

**Literature:**
- Cheney, N.P. & Gould, J.S. (1995). Int. J. Wildland Fire, 5(4):237-247.
- Gould, J.S., et al. (2007). Project Vesta. Ensis-CSIRO.

**Integration Needed:** Add parameters and use in grassfire simulations

---

### 8. Burnout Time Separation (Flaming vs Smoldering)
**Files Modified:**
- `src/burnout_time.H` - Added separation functions

**Implementation:**
- `compute_flaming_smoldering_times()` - splits total residence time
- `get_flaming_fraction_by_fuel_type()` - fuel-dependent fractions
- Fine fuels: 70% flaming, Medium: 40%, Heavy: 20%, Duff: 10%

**Literature:**
- Anderson, H.E. (1969). USDA Research Paper INT-69.
- Byram, G.M. (1959). Forest Fire: Control and Use, pp. 61-89.
- Frandsen, W.H. (1997). Can. J. Forest Res., 27(9):1471-1477.

**Integration Needed:** Expose as user-configurable parameters

---

### 9. Simard Moisture Model (Exponential Time-Lag)
**Files Created:**
- `src/simard_moisture.H` - Complete implementation

**Implementation:**
- `simard_moisture_update()` - M(t+dt) = M_eq + (M(t) - M_eq) * exp(-dt/tau)
- Standard time-lags: 1-hr, 10-hr, 100-hr
- `compute_equilibrium_moisture_from_rh()` - equilibrium from RH
- Simpler alternative to full differential equations

**Literature:**
- Simard, A.J. (1968). Forestry Branch Info. Rep. FF-X-14.
- Fosberg, M.A. (1975). USDA Research Paper RM-152.

**Integration Needed:** Add as moisture model option

---

### 10. Post-Frontal Smoldering (Residual Heat Release)
**Files Modified:**
- `src/duff_moisture_smoldering.H` - Extended with post-frontal features

**Implementation:**
- `compute_residual_smoldering_heat_release()` - exponential decay
- Fuel-type-dependent time constants (30 min - 6 hours)
- `compute_total_post_frontal_heat_release()` - combines all sources
- Important for smoke and air quality modeling

**Literature:**
- Frandsen, W.H. (1997). Can. J. Forest Res., 27(9):1471-1477.
- Urbanski, S.P. (2014). Forest Ecol. Mgmt., 317:1-8.

**Integration Needed:** Track time since front passage, add to heat release output

---

## Code Quality

All implementations follow established patterns:
- GPU-compatible (AMREX_GPU_HOST_DEVICE macros)
- Properly documented with references
- Physically sound (clamped ranges, safe defaults)
- Consistent with existing code style
- 50-200 lines per feature (minor additions)

## Next Steps for Full Integration

1. **parse_inputs.cpp**: Add parameters for all new features
2. **main.cpp**: Call new functions in simulation loop
3. **plot_results.H**: Add new diagnostics to plotfile output
4. **Regression tests**: Create tests for each feature
5. **Documentation**: Update usage guide with new parameters

## Summary

All 10 features from wildfire literature are now implemented and ready for integration. Each feature adds significant scientific value with minimal code complexity, following the design goal of "easy features with high value."
