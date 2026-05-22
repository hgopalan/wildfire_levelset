# FARSITE-Compatible Features Implementation Summary

This document summarizes FARSITE-compatible features added to the wildfire_levelset toolkit.

## Feature 1: Fire Intensity Class Raster

### Overview
Classifies Byram (1959) fireline intensity into 6 discrete classes following the FARSITE / NFFL convention, matching the classification used in the HTML fire report.

### Implementation
- **New file**: `src/fire_intensity_class.H`
  - `classify_fire_intensity()`: GPU-compatible function for per-cell classification
  - `compute_fire_intensity_class()`: MultiFab-level computation
- **Modified files**:
  - `src/multifab_setup.H`: Added `fire_intensity_class_mf` MultiFab
  - `src/plot_results.H`: Added computation and output
  - `src/parse_inputs.H`: Added to plot_vars documentation

### Classes
| Class | Intensity (kW/m) | Flame Length (m) | Description |
|-------|------------------|------------------|-------------|
| I     | < 10             | < 0.5            | Creeping / smouldering |
| II    | 10-175           | 0.5-1.5          | Low intensity surface fire |
| III   | 175-500          | 1.5-2.5          | Moderate intensity |
| IV    | 500-2000         | 2.5-5            | High intensity; spotting begins |
| V     | 2000-10000       | 5-11             | Very high; active crown fire |
| VI    | > 10000          | > 11             | Extreme; firestorm |

### Usage
The field is automatically computed and written to every plotfile as `fire_intensity_class`.

To selectively output only this field:
```
plot_vars = fire_intensity_class fireline_intensity
```

---

## Feature 2: Spot Fire Ignition Delay

### Overview
Models the delay between firebrand landing and actual ignition based on fuel moisture content, following FARSITE's approach to more realistically represent fire spread timing.

### Implementation
- **New file**: `src/spot_ignition_delay.H`
  - `compute_spot_ignition_delay()`: Delay calculation based on moisture deficit
  - `compute_spatial_spot_ignition_delay()`: Per-cell delay with spatial moisture
- **Modified files**:
  - `src/parse_inputs.H`: Added delay parameters to `SpottingParams`
  - `src/parse_inputs.cpp`: Added parameter initialization and validation

### Model
The delay formula is:
```
τ_delay [s] = τ_base × (M_x / (M_x - M_f))
```

Where:
- `τ_base` = base delay time at zero moisture [s] (default: 120 s)
- `M_x` = extinction moisture content [fraction]
- `M_f` = current fuel moisture content [fraction]

At extinction moisture (M_f → M_x), delay approaches infinity (no ignition).
At very dry conditions (M_f → 0), delay approaches τ_base.

### Usage
Enable ignition delay in your inputs file:
```
spotting.enable = 1
spotting.enable_delay = 1
spotting.tau_base = 120.0  # base delay in seconds
```

### Status
✅ Parameters and delay model implemented
⏳ Integration with firebrand_spotting.H pending (requires MultiFab for tracking landing times)

---

## Feature 3: Pasquill-Gifford Atmospheric Stability Class

### Overview
Integrates Pasquill-Gifford atmospheric stability classification into the Briggs plume rise model, allowing FARSITE-compatible adjustment of smoke plume rise based on atmospheric conditions.

### Implementation
- **Modified files**:
  - `src/smoke_plume_rise.H`: 
    - Added `get_stability_correction_factor()` function
    - Updated header documentation with Pasquill-Gifford references
    - Modified `compute_smoke_plume_rise()` to apply stability correction
  - `src/parse_inputs.H`: Added stability class parameters to `SmokePlumeParams`
  - `src/parse_inputs.cpp`: Added parameter initialization and validation

### Stability Classes
| Class | Conditions | Correction Factor | Effect on Plume Rise |
|-------|------------|-------------------|---------------------|
| A | Extremely unstable (strong solar radiation, light winds) | 1.20 | +20% enhanced rise |
| B | Moderately unstable | 1.10 | +10% enhanced rise |
| C | Slightly unstable | 1.05 | +5% enhanced rise |
| D | Neutral (overcast or moderate wind) | 1.00 | No correction |
| E | Slightly stable (clear night, light winds) | 0.85 | -15% suppressed rise |
| F | Moderately stable (clear night, very light winds) | 0.70 | -30% suppressed rise |

### Physics
The base Briggs (1965/1969) formula computes:
```
Δh = 1.6 × F_B^(1/3) × x_f^(2/3) / u
```

With Pasquill-Gifford correction:
```
Δh_corrected = Δh × stability_factor
```

This accounts for:
- **Unstable conditions (A-C)**: Enhanced vertical mixing → greater plume rise
- **Stable conditions (E-F)**: Suppressed vertical motion → reduced plume rise

### Usage
Enable stability correction in your inputs file:
```
smoke_plume.enable = 1
smoke_plume.use_stability_correction = 1
smoke_plume.stability_class = A  # Options: A, B, C, D, E, or F
```

Default (without correction) is equivalent to:
```
smoke_plume.use_stability_correction = 0  # or omit
smoke_plume.stability_class = D
```

---

## References

### Fire Intensity Classification
- Byram, G.M. (1959). Combustion of forest fuels. In: Davis, K.P. (ed.) Forest Fire: Control and Use. McGraw-Hill.
- Andrews, P.L. et al. (2011). BehavePlus fire modeling system, version 5.0: Design and features. USDA Forest Service RMRS-GTR-249.

### Spot Fire Ignition Delay
- Albini, F.A. (1979). Spot fire distance from burning trees. USDA Forest Service INT-56.
- Finney, M.A. (1998). FARSITE: Fire Area Simulator. USDA Forest Service RMRS-RP-4.
- Anderson, H.E. (1970). Forest fuel ignitibility. Fire Technology 6(4):312-319.

### Pasquill-Gifford Stability
- Briggs, G.A. (1965). A plume rise model compared with observations. JAPCA 15(9):433-438.
- Briggs, G.A. (1969). Plume rise. USAEC Critical Review Series, TID-25075.
- Pasquill, F. (1961). The estimation of the dispersion of windborne material. Meteorological Magazine 90:33-49.
- Gifford, F.A. (1961). Use of routine meteorological observations for estimating atmospheric dispersion. Nuclear Safety 2:47-51.
- Turner, D.B. (1970). Workbook of atmospheric dispersion estimates. US EPA AP-26.

---

## Feature 4: FARSITE Temporal Fire Acceleration Model

### Overview
Implements the full FARSITE temporal acceleration model (McAlpine & Wakimoto 1991 / VanWagner's equation) with per-cell temporal tracking and wind-onset time-lag capability. Models the delayed approach to quasi-steady-state spread rate in small fires or after wind changes.

### Implementation
- **Modified file**: `src/fire_acceleration.H`
  - `apply_fire_acceleration_size_based()`: Original Catchpole et al. 1992 size-based model
  - `apply_fire_acceleration_temporal()`: New FARSITE temporal model with VanWagner's equation
  - `apply_fire_acceleration()`: Main dispatch function
- **Modified files**:
  - `src/parse_inputs.H`: Extended `AccelerationParams` with temporal model parameters
  - `src/parse_inputs.cpp`: Added parameter initialization and validation
  - `src/multifab_setup.H`: Added `accel_state_mf` for per-cell state tracking
  - `src/main.cpp`: Initialize temporal state and pass dt to acceleration function
  - `docs/usage.rst`: Comprehensive documentation of both models

### Models

#### Size-Based Model (Catchpole et al. 1992)
```
α = 1 - exp(-r_fire / L_acc)
R_eff = α × R_QSS

where r_fire = √(A_burned / π)
```

- Global scaling based on current fire size
- Simple and computationally efficient
- Fire achieves ~63% of QSS at r_fire = L_acc, ~95% at r_fire = 3×L_acc

#### FARSITE Temporal Model (McAlpine & Wakimoto 1991)
```
R(t) = R_E × (1 - exp(-A × t))

where:
  A = A_point (0.115 1/min) for small fires (perim < 402.3 m)
  A = A_line (0.300 1/min) for large fires (perim >= 402.3 m)
  t = elapsed time since entering current acceleration state [s]
```

- Per-cell temporal tracking of acceleration state
- Automatic switching between point and line acceleration constants
- Each cell tracks: current ROS, elapsed time, equilibrium ROS

#### Wind-Onset Time-Lag Extension
```
When R_E changes (wind change):
  R_target(t) = R_prev + (R_E - R_prev) × (1 - exp(-dt/tau_wind))
  R(t) = R_target × (1 - exp(-A × t_accel))
```

- Optional wind-lag model for realistic wind response
- Fire ROS ramps up exponentially after wind changes
- Controlled by time constant tau_wind (default: 180 s = 3 min)

### Usage

Enable size-based model (backward compatible):
```
acceleration.enable = 1
acceleration.L_acc = 50.0
```

Enable FARSITE temporal model:
```
acceleration.enable = 1
acceleration.use_temporal = 1
acceleration.A_point = 0.115
acceleration.A_line = 0.300
acceleration.perim_limit = 402.3
```

Enable wind-onset time-lag:
```
acceleration.enable = 1
acceleration.use_temporal = 1
acceleration.enable_wind_lag = 1
acceleration.tau_wind = 180.0
```

### Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `acceleration.enable` | 0 | Enable acceleration model (1=yes, 0=no) |
| `acceleration.use_temporal` | 0 | 0=size-based, 1=FARSITE temporal |
| `acceleration.L_acc` | 50.0 m | Length scale for size-based model |
| `acceleration.A_point` | 0.115 1/min | Point ignition acceleration constant |
| `acceleration.A_line` | 0.300 1/min | Line ignition acceleration constant |
| `acceleration.perim_limit` | 402.3 m | Perimeter threshold (20 chains) |
| `acceleration.enable_wind_lag` | 0 | Enable wind-onset time-lag (1=yes, 0=no) |
| `acceleration.tau_wind` | 180.0 s | Wind response time constant (3 min) |

### Validation

Regression tests created in `regtest/surface_spread/fire_acceleration/`:
- **inputs.size_based**: Validates size-based model (Catchpole et al. 1992)
- **inputs.temporal_point**: Validates FARSITE temporal model with point ignition
- **inputs.wind_lag**: Validates wind-onset time-lag capability

Expected behavior:
- Point ignition: ~20 min to reach 90% equilibrium (A = 0.115 1/min)
- Line ignition: ~8 min to reach 90% equilibrium (A = 0.300 1/min)
- Wind-lag: ROS lags behind equilibrium by ~3×tau_wind after sudden wind increase

### References
- McAlpine, R.S. & Wakimoto, R.H. (1991). "The acceleration of fire from point source to equilibrium spread." *Forest Science*, 37(5), 1314–1337.
- Alexander, M.E., Stocks, B.J. & Lawson, B.D. (1992). "Fire behavior in Black Spruce-lichen woodland." Info. Rep. NOR-X-310. USDA/CFS.
- Catchpole, E.A., de Mestre, N.J. & Gill, A.M. (1992). "Intensity of fire at its perimeter." *Australian Journal of Ecology*, 17(1), 1–4.
- Finney, M.A. (1998/2004). "FARSITE: Fire Area Simulator." USDA Forest Service RMRS-RP-4.

### Status
✅ Size-based model implemented and working
✅ FARSITE temporal model implemented with per-cell state tracking
✅ Wind-onset time-lag capability implemented
✅ Regression tests created
✅ Documentation updated

---

## Files Modified

### New Files
1. `src/fire_intensity_class.H` - Fire intensity classification functions
2. `src/spot_ignition_delay.H` - Ignition delay model functions
3. `regtest/surface_spread/fire_acceleration/` - Fire acceleration regression tests

### Modified Files
1. `src/multifab_setup.H` - Added fire_intensity_class_mf MultiFab and accel_state_mf
2. `src/plot_results.H` - Added fire intensity class computation and output
3. `src/parse_inputs.H` - Added parameters for all features including acceleration
4. `src/parse_inputs.cpp` - Added parameter initialization and validation
5. `src/smoke_plume_rise.H` - Added Pasquill-Gifford stability correction
6. `src/fire_acceleration.H` - Complete rewrite with FARSITE temporal model
7. `src/main.cpp` - Initialize acceleration state and pass dt to acceleration function
8. `docs/usage.rst` - Comprehensive documentation for all features

---

## Testing Recommendations

### Fire Intensity Class
1. Run existing regression test with fire_intensity_class in plot_vars
2. Verify classification matches expected intensity ranges
3. Compare with HTML report intensity table

### Spot Fire Ignition Delay
1. Create test with spotting enabled and varying moisture conditions
2. Verify delay increases as M_f approaches M_x
3. Test spatial moisture variation

### Pasquill-Gifford Stability
1. Run smoke plume test with different stability classes
2. Verify:
   - Class A gives ~20% higher plume rise than Class D
   - Class F gives ~30% lower plume rise than Class D
3. Compare with/without stability correction enabled
