# WindNinja Ridge/Canyon Speed-up Regression Test

## Purpose

Tests the WindNinja-style ridge speed-up and canyon channeling wind-terrain model
(`wind_terrain.model = windninja_ridge_canyon`, Option 7).

This model accounts for the empirical observation that:
- **Ridges**: Wind accelerates when it climbs upslope (due to convergence at the crest)
- **Canyons/Valleys**: Wind accelerates when it flows down a channel (drainage flow)

## Physical Model

The effective wind speed is modified based on the alignment between the wind direction
and the terrain gradient:

```
alignment = dot(wind_xy, upslope_unit) / |wind|
```

Where `upslope_unit = grad(z) / |grad(z)|` points in the direction of steepest ascent.

- **Ridge (alignment > 0)**: wind is blowing upslope
  ```
  f = 1 + k_ridge * tan_phi * alignment
  ```

- **Canyon (alignment < 0)**: wind is blowing downslope
  ```
  f = 1 + k_canyon_wn * tan_phi * |alignment|
  ```

The effective wind is `U_eff = U * f` (always ≥ ambient wind speed).

## Configuration

- **Domain**: 1000 m × 1000 m (same Gaussian hill as `terrain_wind` test)
- **Grid**: 100 × 100 cells (10 m resolution)
- **Fire spread model**: Rothermel (1972)
- **Wind-terrain model**: Option 7 `windninja_ridge_canyon`
- **k_ridge**: 1.5 (stronger ridge speed-up)
- **k_canyon_wn**: 0.8 (moderate canyon channeling)
- **Wind**: Approximately eastward (from west), from the Gaussian hill wind file

## Expected Behavior

On the Gaussian hill with eastward wind:
- **Western face (upslope)**: Wind aligns with upslope direction → `alignment > 0` → ridge
  speed-up increases Rothermel wind factor → higher ROS on the approach to the crest
- **Eastern face (downslope)**: Wind aligns downslope (past crest) → `alignment < 0` → 
  canyon channeling → additional ROS increase on lee side

The fire should spread faster both on the upslope (windward) face due to ridge speed-up
and on the lee side due to canyon channeling, creating an asymmetric spread pattern
compared to the plain Rothermel model (Option 1) or Pimont exponential correction
(Option 6).

## Comparison with Other Options

| Option | Model | Effect on Gaussian Hill |
|--------|-------|------------------------|
| 1 | `none` | Symmetric, base ROS only |
| 4 | `canyon_wind` | Isotropic amplification ∝ tan(phi), no wind direction dependence |
| 6 | `pimont` | Exponential amplification, no wind direction dependence |
| 7 | `windninja_ridge_canyon` | Directional: ridge on windward face, canyon on lee |

## Files

- `inputs.i`: Input parameters
- `gaussian_hill_terrain.csv`: Terrain elevation data (reused from `terrain_wind`)
- `gaussian_hill_wind.csv`: Wind velocity data (reused from `terrain_wind`)

## Run Command

```bash
cd regtest/windninja_ridge_canyon
../../build/levelset inputs.i
```

> **Note**: This test requires a 2D build (`-DLEVELSET_DIM_2D=ON`).

## References

- Forthofer, J.M. (2007). *Modeling Wind in Complex Terrain for Use in Fire Spread Prediction*. Colorado State University MS thesis.
- Rothermel, R.C. (1983). *How to Predict the Spread and Intensity of Forest and Range Fires*. USDA Forest Service GTR INT-143.
