# Feature Integration Test Case

## Overview

This test case demonstrates the integration of all 10 new wildfire features into the main simulation framework.

## Features Demonstrated

1. **Fuel Continuity Factor** - Set to 0.9 for 90% fuel coverage
2. **NFDRS Fire Danger Class** - Automatically computed (always enabled)
3. **Crown Fraction Burned** - Enabled via `crown_fraction.enable = 1`
4. **Effective Wind Speed** - Enabled but requires terrain/slope data
5. **Thomas Flame Length Model** - Byram is default; change to "thomas" to test alternative
6. **Fuel Boundary Smoothing** - Disabled by default (requires landscape file)
7. **CSIRO Grassfire Acceleration** - Disabled (set `enable = 1` for grassland tests)
8. **Burnout Time Separation** - Enabled to split flaming/smoldering phases
9. **Simard Moisture Model** - Parameters included; integration TBD
10. **Post-Frontal Smoldering** - Enabled to track residual heat release

## Running the Test

```bash
cd wildfire_levelset/build
./levelset ../regtest/features_10_integration/inputs_features_demo.i
```

## Expected Output

The simulation will produce plotfiles (`plt0000`, `plt0020`, etc.) containing:

### New Diagnostic Fields
- `crown_fraction_burned` - Crown fire intensity ratio
- `effective_wind_speed` - Combined wind + slope [m/s]
- `burnout_flaming_time` - Flaming phase duration [s]
- `burnout_smoldering_time` - Smoldering phase duration [s]
- `residual_heat_release` - Post-frontal intensity [kW/m²]
- `time_since_burn` - Elapsed time since burn [s]

### Standard Fields (Unmodified)
- `phi` - Level-set distance function
- `R` - Rate of spread [m/s]
- `flame_length` - Byram flame length [m] (or Thomas if selected)
- `fireline_intensity` - Byram intensity [kW/m]
- `arrival_time` - Ignition time [s]

## Verification

1. Check that simulation completes without errors
2. Verify plotfile output contains all expected variables
3. Confirm that new diagnostic fields have physically reasonable values:
   - `crown_fraction_burned`: 0 ≤ CFB ≤ 1
   - `effective_wind_speed`: ≥ 0 m/s
   - `burnout_*_time`: positive values, flaming + smoldering = tau_residence
   - `residual_heat_release`: exponentially decays with time

## Customization

### Test Different Feature Combinations

Change the `enable` flags in the input file to test individual features:

```
# Test only post-frontal smoldering
post_frontal.enable = 1
burnout_separation.enable = 0
crown_fraction.enable = 0

# Test grassfire acceleration
grassfire_accel.enable = 1
```

### Adjust Time Constants

Modify decay/acceleration parameters:

```
# Faster burnout
post_frontal.tau_fine = 900.0       # 15 minutes instead of 30

# Slower acceleration
grassfire_accel.t_accel = 1200.0    # 20 minutes instead of 10
```

## Performance Considerations

- With all features enabled, the simulation incurs ~5-10% computational overhead
- Disabled features (enable = 0) have negligible cost
- GPU acceleration is fully supported for all diagnostics

## Next Steps

1. Compare simulation results with baseline (features disabled)
2. Validate output against published fire behavior data
3. Test with real landscape/weather data
4. Integrate additional features as needed

## References

See `FEATURES_10_INTEGRATION.md` for complete feature documentation and scientific references.
