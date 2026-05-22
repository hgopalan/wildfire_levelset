# Cell Size Effects Correction Feature

## Overview

This feature implements empirical correction factors for fire spread rates based on grid cell resolution. It accounts for the physical reality that finer computational grids better resolve fire-front geometry, leading to higher effective spread rates when compared to coarser grids.

## Physical Basis

Fire spread models like Rothermel assume a representative fuel bed at a given scale. When implementing these models on computational grids:
- **Coarser grids**: Numerical diffusion and averaging effects cause under-prediction of fire spread
- **Finer grids**: Better resolution of fire-front topology leads to higher effective spread rates

FARSITE addresses this through empirical cell size correction factors, which this implementation provides.

## Implementation Details

### Parameters

Add the following to your inputs file to enable cell size correction:

```
cellsize.enable = 1
cellsize.dx_ref = 30.0
cellsize.correction_exponent = 0.1
```

- **cellsize.enable** (default: 0)
  - 0: Disable cell size correction (backward compatible, default)
  - 1: Enable empirical correction factors

- **cellsize.dx_ref** (default: 30.0) [meters]
  - Reference cell size for correction formula
  - Typical value: 30 m (common for FARSITE simulations)

- **cellsize.correction_exponent** (default: 0.1)
  - Power-law exponent in correction formula
  - Range: [0.0, 1.0]
  - Larger values = stronger resolution effect
  - Empirical value: 0.1 (balances realism and stability)

### Correction Formula

The rate of spread is adjusted as:

```
R_corrected = R_base × (dx_ref / dx_min)^exponent
```

where:
- `R_base` = rate of spread without correction (m/s)
- `dx_ref` = reference cell size (default: 30 m)
- `dx_min` = actual minimum cell size in domain = min(dx, dy, dz)
- `exponent` = power-law exponent (default: 0.1)

### Examples

With `dx_ref = 30.0` and `exponent = 0.1`:

| Cell Size (m) | dx/dx_ref | Correction Factor | Effect |
|---|---|---|---|
| 10 | 3.0 | (3.0)^0.1 = 1.116 | +11.6% faster |
| 15 | 2.0 | (2.0)^0.1 = 1.072 | +7.2% faster |
| 20 | 1.5 | (1.5)^0.1 = 1.041 | +4.1% faster |
| 30 | 1.0 | (1.0)^0.1 = 1.000 | No change |
| 60 | 0.5 | (0.5)^0.1 = 0.933 | -6.7% slower |
| 100 | 0.3 | (0.3)^0.1 = 0.878 | -12.2% slower |

## Usage

### FARSITE Propagation (Recommended)

Cell size correction is designed for FARSITE elliptical propagation and applies only when `propagation_method = farsite`:

```
propagation_method = farsite
cellsize.enable = 1
cellsize.dx_ref = 30.0
cellsize.correction_exponent = 0.1
```

The correction is applied to the rate of spread computed by the selected fire spread model (Rothermel, Balbi, FBP, Cheney-Gould, or Lautenberger).

### Level-Set Propagation

Cell size correction is **not applied** to level-set propagation, even if `cellsize.enable = 1`. The level-set method naturally incorporates resolution effects through its CFL-based timestep constraint. Parameters are silently ignored.

```
propagation_method = levelset
cellsize.enable = 1  # Ignored; will not affect levelset spread
```

### MTT Propagation

Cell size correction is **not applied** to minimum travel time (MTT) propagation. Parameters are silently ignored.

```
propagation_method = mtt
cellsize.enable = 1  # Ignored; will not affect MTT spread
```

## Backward Compatibility

- **Default behavior preserved**: Cell size correction is disabled by default (`cellsize.enable = 0`)
- **Existing simulations unaffected**: All existing input files run without change
- **Optional feature**: Only applies when explicitly enabled
- **Non-destructive**: Original ROS computation remains unchanged when disabled

## Theory References

### FARSITE Implementation
- Finney, M.A. (1998). FARSITE: Fire Area Simulator. USDA Forest Service RMRS-RP-4.
  - Section 3.4: Discussion of cell size effects on fire spread
  - Empirical observations showing resolution-dependent behavior

### Fire Behavior Theory
- Rothermel, R.C. (1972). A mathematical model for predicting fire spread in wildland fuels. USDA Forest Service Res. Pap. INT-115.
  - Foundation for spread models used with this correction

- Anderson, H.E. (1982). Aids to determining fuel models for estimating fire behavior. USDA Forest Service Gen. Tech. Rep. INT-122.
  - Scale considerations in fuel model selection

## Files Modified

1. **src/parse_inputs.H**
   - Added `CellSizeParams` struct with parameters

2. **src/parse_inputs.cpp**
   - Parameter parsing and validation for `cellsize.*` parameters

3. **src/cell_size_correction.H** (new)
   - GPU-compatible correction factor computation

4. **src/compute_rothermel_R.H**
   - Modified to extract cell sizes and apply correction factor to final ROS

5. **src/main.cpp**
   - Updated all `compute_rothermel_R` calls to pass correction parameters
   - Applied only for FARSITE initialization and main loop

6. **regtest/CMakeLists.txt**
   - Added 3 regression tests for the feature

7. **regtest/surface_spread/cell_size_correction_***/**
   - Test cases for baseline, with correction, and levelset validation

## Testing

Three regression tests validate the feature:

1. **cell_size_correction_baseline**: FARSITE without correction (control)
2. **cell_size_correction_with_correction**: FARSITE with correction enabled (experimental)
3. **cell_size_correction_levelset**: Level-set with correction parameters (validates they're ignored)

Run tests with:
```bash
ctest -L regtest -R cell_size_correction
```

## Notes for Users

- **Validation recommended**: Test the correction factor with your specific fuel and terrain before using in production simulations
- **Sensitivity analysis**: Consider running simulations at multiple resolutions to validate the correction exponent
- **Documentation**: Include the correction parameters in simulation metadata and reports for reproducibility
- **Coupling with other features**: Works with all fire spread models (Rothermel, Balbi, FBP, etc.) and all FARSITE fire shape options

## Future Enhancements

Potential extensions:
- Fuel-type-dependent exponents
- Direction-dependent corrections (finer resolution parallel to wind)
- Adaptive exponents based on ROS magnitude
- Integration with machine learning calibration
