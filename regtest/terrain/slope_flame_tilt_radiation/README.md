# Slope-Dependent Flame Tilt for Radiation Test

## Purpose

Tests the **slope-dependent flame tilt** feature that enhances radiation preheating on upslope fires by accounting for flame angle influenced by terrain slope.

This demonstrates that:
1. `radiation_preheating.use_slope_tilt = 1` enables slope-flame tilt coupling
2. Upslope fires have flames tilted forward by terrain angle
3. Forward-tilted flames increase preheating distance downwind
4. Fire spread accelerates more on steep slopes compared to flat terrain
5. Flame tilt angle is computed from wind and slope effects

## Test Configuration

| Parameter | Value | Description |
|-----------|-------|-------------|
| Domain | 1000 m × 1000 m | Simple upslope terrain |
| Grid | 64 × 64 cells | ~15.6 m resolution |
| Terrain slope | 30° (0.577) | Steep slope in +x direction |
| Elevation range | 0–577 m | Linear rise across domain |
| Wind | 4 m/s eastward | Aligned with slope |
| Fuel model | FM4 chaparral | Southern California fuel |
| Fuel moisture | 8% | Moderate conditions |
| Slope tilt coupling | Enabled | Key test feature |

## Expected Behavior

1. **Fire ignition**: Circular fire starts at base of slope
2. **Upslope spread**: Fire accelerates in +x direction (upslope)
   - Standard Rothermel slope factor + radiation preheating
   - Enhanced by forward-tilted flames
3. **Cross-slope spread**: Slower spread in ±y directions (no slope effect)
4. **Perimeter asymmetry**: Elongated ellipse biased toward upslope

## Files

- `inputs.i`: Test configuration
- `slope_terrain.csv`: 30° upslope terrain data
- `README.md`: This file

## Run Command

```bash
cd regtest/terrain/slope_flame_tilt_radiation
../../build/levelset inputs.i
```

> **Note**: Requires a 2D build (`-DLEVELSET_DIM_2D=ON`).

## References

- Butler, B.W., et al. (2004). A radiation-driven model for crown fire spread. *Canadian Journal of Forest Research*, 34(8), 1588-1599.
- Viegas, D.X. (2004). Slope and wind effects on fire propagation. *International Journal of Wildland Fire*, 13(2), 143-156.
