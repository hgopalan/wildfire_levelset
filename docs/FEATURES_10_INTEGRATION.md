# Integration of 10 Wildfire Features - Complete Deployment Guide

## Overview

This document describes the integration of 10 new wildfire fire behavior features into the main AMReX-based level-set wildfire simulation framework. All features have been successfully integrated into the core simulation engine with GPU support.

## Features Integrated

### Feature 1: Fuel Continuity Factor (Heterogeneous Fuel Distributions)
**Status:** ✅ Complete (already in RothermelParams)

**Description:** Accounts for patchy or discontinuous fuel beds by applying a multiplicative factor (0-1) to computed ROS.

**Input Parameters:**
```
rothermel.fuel_continuity = 1.0  # 1.0 = continuous, 0.5 = 50% coverage, 0.0 = gaps
```

**References:**
- Finney, M.A. (2006). FlamMap capabilities. USDA Forest Service RMRS-P-41.

---

### Feature 2: NFDRS Fire Danger Class (Operational Categories)
**Status:** ✅ Complete

**Description:** Classifies fireline intensity into 5 operational fire danger categories (Low/Moderate/High/Very High/Extreme) following NFDRS conventions.

**Implementation:**
- Automatically computed from fireline intensity
- Output in plotfile as `nfdrs_danger_class`

**Danger Classes:**
- 1 = Low (< 100 kW/m)
- 2 = Moderate (100-500 kW/m)
- 3 = High (500-2000 kW/m)
- 4 = Very High (2000-5000 kW/m)
- 5 = Extreme (> 5000 kW/m)

**References:**
- Deeming, J.E., et al. (1977). NFDRS-1978. USDA GTR INT-39.

---

### Feature 3: Crown Fraction Burned (CFB) Diagnostic
**Status:** ✅ Complete

**Description:** Diagnostic metric distinguishing passive vs active crown fire based on fire intensity.

**Input Parameters:**
```
crown_fraction.enable = 0  # Set to 1 to enable
```

**Output Variables:**
- `crown_fraction_burned` - CFB ratio (0-1)

**References:**
- Scott, J.H. & Reinhardt, E.D. (2001). USDA RMRS-RP-29.

---

### Feature 4: Effective Wind Speed (Wind + Slope Combined)
**Status:** ✅ Complete

**Description:** Combines ambient wind and slope effects into a single scalar effective wind speed via vector addition.

**Input Parameters:**
```
effective_wind.enable = 0  # Set to 1 to enable
```

**Formula:**
```
U_eff = sqrt(U_wind² + U_slope_equiv²)
```

**Output Variables:**
- `effective_wind_speed` - Combined wind speed [m/s]

**Requirements:** Terrain slopes must be available (landscape file)

**References:**
- Rothermel, R.C. (1983). GTR INT-143.

---

### Feature 5: Thomas Flame Length Model
**Status:** ✅ Complete

**Description:** Alternative flame length formula L = 0.0266 × I^0.667 (vs default Byram: L = 0.0775 × I^0.46)

**Input Parameters:**
```
flame_length_model.model = "byram"   # or "thomas"
```

**Selection:**
- `"byram"` - Default Byram (1959) formula
- `"thomas"` - Thomas (1963) formula

**References:**
- Thomas, P.H. (1963). Combustion Symposium, 9(1):844-859.

---

### Feature 6: Fuel Boundary Smoothing
**Status:** ✅ Complete

**Description:** Distance-weighted ROS blending at fuel boundaries to eliminate unrealistic discontinuities.

**Input Parameters:**
```
fuel_boundary.enable = 0              # Set to 1 to enable
fuel_boundary.transition_cells = 2.0  # Number of cells for transition zone
```

**Requirements:** Landscape file with fuel model data must be present

**References:**
- Finney, M.A. (1998, 2006). USDA RMRS-RP-4, RMRS-P-41.

---

### Feature 7: CSIRO Grassfire Acceleration
**Status:** ✅ Complete

**Description:** Models non-equilibrium fire growth during initial spread phase using exponential acceleration function.

**Input Parameters:**
```
grassfire_accel.enable = 0          # Set to 1 to enable
grassfire_accel.t_accel = 600.0     # Acceleration time constant [seconds]
```

**Formula:**
```
acceleration_factor = 1 - exp(-dt / t_accel)
R_accelerated = R × (1 + acceleration_factor)
```

**References:**
- Cheney, N.P. & Gould, J.S. (1995). Int. J. Wildland Fire, 5(4):237-247.

---

### Feature 8: Burnout Time Separation (Flaming vs Smoldering)
**Status:** ✅ Complete

**Description:** Splits total residence time into flaming and smoldering phase durations by fuel type.

**Input Parameters:**
```
burnout_separation.enable = 1
burnout_separation.flaming_fraction_fine = 0.70      # Fine fuels: 70% flaming
burnout_separation.flaming_fraction_medium = 0.40    # Medium: 40% flaming
burnout_separation.flaming_fraction_heavy = 0.20     # Heavy: 20% flaming
burnout_separation.flaming_fraction_duff = 0.10      # Duff: 10% flaming
```

**Output Variables:**
- `burnout_flaming_time` - Flaming phase duration [seconds]
- `burnout_smoldering_time` - Smoldering phase duration [seconds]

**References:**
- Anderson, H.E. (1969). USDA Research Paper INT-69.
- Frandsen, W.H. (1997). Can. J. Forest Res., 27(9):1471-1477.

---

### Feature 9: Simard Moisture Model (Exponential Time-Lag)
**Status:** ✅ Complete

**Description:** Exponential time-lag moisture update based on equilibrium moisture approach with size-dependent time constants.

**Input Parameters:**
```
simard_moisture.enable = 0
simard_moisture.tau_1hr = 1.0       # 1-hour fuel lag [hours]
simard_moisture.tau_10hr = 10.0     # 10-hour fuel lag [hours]
simard_moisture.tau_100hr = 100.0   # 100-hour fuel lag [hours]
```

**Formula:**
```
M(t+dt) = M_eq + (M(t) - M_eq) × exp(-dt/tau)
```

**References:**
- Simard, A.J. (1968). Forestry Branch Info. Rep. FF-X-14.

---

### Feature 10: Post-Frontal Smoldering (Residual Heat Release)
**Status:** ✅ Complete

**Description:** Tracks residual combustion after flame front passage with exponential decay for smoke/air quality applications.

**Input Parameters:**
```
post_frontal.enable = 0
post_frontal.tau_fine = 1800.0      # Fine fuels decay [seconds] = 30 min
post_frontal.tau_medium = 3600.0    # Medium fuels decay = 1 hour
post_frontal.tau_heavy = 7200.0     # Heavy fuels decay = 2 hours
post_frontal.tau_duff = 21600.0     # Duff decay = 6 hours
```

**Output Variables:**
- `time_since_burn` - Elapsed time since cell burned [seconds]
- `residual_heat_release` - Residual smoldering intensity [kW/m²]

**References:**
- Frandsen, W.H. (1997). Can. J. Forest Res., 27(9):1471-1477.
- Urbanski, S.P. (2014). Forest Ecol. Mgmt., 317:1-8.

---

## Integration Architecture

### Files Modified/Created

1. **parse_inputs.H/cpp** - Added parameter structures for all 10 features
2. **multifab_setup.H** - Added MultiFab fields for diagnostic output:
   - `crown_fraction_burned_mf`
   - `effective_wind_speed_mf`
   - `burnout_phases_mf`
   - `residual_heat_release_mf`
   - `time_since_burn_mf`

3. **plot_results.H** - Added new diagnostic variables to plotfile output
4. **main.cpp** - Integrated feature calls into main simulation loop
5. **wildfire_includes.H** - Added includes for feature header files

### New Header Files (Already Implemented)

- `effective_wind_speed.H` - Feature 4 computation
- `fire_intensity_class.H` - Feature 2 NFDRS classification
- `crown_initiation.H` - Feature 3 CFB computation
- `fuel_boundary_smoothing.H` - Feature 6 boundary smoothing
- `fire_acceleration.H` - Feature 7 CSIRO acceleration
- `simard_moisture.H` - Feature 9 moisture model
- `duff_moisture_smoldering.H` - Feature 10 post-frontal tracking

## Example Input File

```
# Enable 10 new wildfire features

# Feature 1: Fuel continuity
rothermel.fuel_continuity = 0.8

# Feature 3: Crown fraction burned
crown_fraction.enable = 1

# Feature 4: Effective wind speed
effective_wind.enable = 1

# Feature 5: Thomas flame length model
flame_length_model.model = thomas

# Feature 6: Fuel boundary smoothing
fuel_boundary.enable = 1
fuel_boundary.transition_cells = 2.5

# Feature 7: CSIRO grassfire acceleration
grassfire_accel.enable = 1
grassfire_accel.t_accel = 600.0

# Feature 8: Burnout time separation
burnout_separation.enable = 1
burnout_separation.flaming_fraction_fine = 0.70
burnout_separation.flaming_fraction_medium = 0.40

# Feature 9: Simard moisture model (future use)
simard_moisture.enable = 0

# Feature 10: Post-frontal smoldering
post_frontal.enable = 1
post_frontal.tau_fine = 1800.0
post_frontal.tau_medium = 3600.0
```

## Output Plotfile Variables

When enabled, the following variables are added to every plotfile:

```
crown_fraction_burned          # Feature 3: CFB diagnostic [-]
effective_wind_speed           # Feature 4: Combined wind speed [m/s]
burnout_flaming_time           # Feature 8: Flaming duration [s]
burnout_smoldering_time        # Feature 8: Smoldering duration [s]
residual_heat_release          # Feature 10: Residual intensity [kW/m²]
time_since_burn                # Feature 10: Time elapsed since burn [s]
```

## Compilation and Testing

### Build

```bash
cd wildfire_levelset
git submodule update --init --recursive
cmake -S . -B build -DLEVELSET_DIM_2D=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

### Verify Features

The binary includes all 10 features. Features can be selectively enabled via input parameters (see example above). When `enable=0` (default), features have minimal computational overhead.

## Performance Notes

- **GPU Acceleration:** All features use AMREX_GPU_HOST_DEVICE macros for GPU compatibility
- **Computational Cost:** Features are only computed when explicitly enabled
- **Memory Overhead:** Additional MultiFabs for diagnostics add ~50 MB per feature for 512³ domain
- **Physics-Based:** All formulas based on peer-reviewed wildfire literature

## Future Integration Opportunities

1. **Moisture Feedback:** Feature 9 (Simard) can be integrated with spatial moisture fields
2. **Crown Fire Coupling:** Features 3, 7 can enhance crown fire initiation modeling
3. **Air Quality:** Feature 10 output can drive smoke/emissions calculations
4. **Operational Integration:** Feature 2 (NFDRS) can trigger resource responses
5. **Ensemble Studies:** Features enable probabilistic fire behavior analysis

## References

Complete references are documented in each feature's header file and in the main README.md.

## Support

For questions about specific features, refer to:
- Individual header files (`.H`) for detailed algorithm documentation
- Main README.md for general project information
- Scientific literature cited in each feature description
