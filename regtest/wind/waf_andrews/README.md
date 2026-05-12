# WAF Andrews Test

## Purpose

Tests the Andrews (2018) logarithmic Wind Adjustment Factor (WAF) together with
the Maximum Effective Wind Speed (MEWS) cap.

The Andrews WAF converts a 20-ft (6.1 m) open wind speed to the midflame height
wind speed required by the Rothermel (1972) model:

```
WAF = 1.83 / ln((20 + 0.36 × h) / (0.13 × h))    [h = fuel bed depth in ft]
```

For the FM4 chaparral fuel (h = 6 ft) the WAF ≈ 0.36.  A 5 m/s ambient wind is
reduced to approximately 1.8 m/s at midflame height.

The MEWS cap (`rothermel.use_wind_limit = 1`) prevents unrealistically large ROS
by limiting the effective wind speed to the Rothermel calibration range.

## Configuration

| Parameter | Value | Notes |
|-----------|-------|-------|
| Domain | 500 m × 500 m | UTM Zone 11N, S. California |
| Grid | 64 × 64 | 2-D build required |
| Fuel | FM4 chaparral | h = 6 ft, M_f = 8 % |
| Wind | 5 m/s eastward | 20-ft reference wind |
| `waf_formula` | `andrews` | logarithmic Albini & Baughman (1979) |
| `use_waf` | 1 | WAF enabled |
| `use_wind_limit` | 1 | MEWS cap enabled |

## Expected Behaviour

- Fire spreads eastward with a reduced (WAF-adjusted) effective wind.
- ROS is lower than without the WAF, consistent with midflame wind sheltering.
- The fire perimeter is elongated in the downwind direction.

## Build

```bash
cmake -S . -B build -DLEVELSET_DIM_2D=ON
cmake --build build -j
```

## Run Command

```bash
./build/levelset regtest/wind/waf_andrews/inputs.i
```

## References

- Andrews, P.L. (2018). *The Rothermel Surface Fire Spread Model and Associated
  Developments: A Comprehensive Explanation.* USDA For. Serv. Gen. Tech. Rep.
  RMRS-GTR-371.
- Albini, F.A. & Baughman, R.G. (1979). *Estimating Windspeeds for Predicting
  Wildland Fire Behavior.* USDA For. Serv. Res. Pap. INT-221.
