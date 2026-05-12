# WAF BehavePlus Test

## Purpose

Tests the BehavePlus simplified linear Wind Adjustment Factor (WAF) for open and
shrub fuel models.

BehavePlus computes the midflame WAF using a linear approximation:

```
WAF = 0.36 + 0.004 × h_in    [h_in = fuel bed depth in inches]
```

For the FM4 chaparral fuel (h = 6 ft = 72 in): WAF = 0.36 + 0.288 = **0.648**.
Compare with the Andrews logarithmic formula for the same fuel: WAF ≈ 0.36.

When a landscape file provides per-cell canopy cover and canopy height, the
closed-canopy path uses exponential Beer–Lambert attenuation:

```
WAF_canopy = WAF_open(h_c) × exp(−α_c × f_c)
```

where `α_c` is the canopy attenuation coefficient (`rothermel.waf_canopy_alpha`,
default 1.5).  This test exercises the open-fuel formula only (no landscape file).

## Configuration

| Parameter | Value | Notes |
|-----------|-------|-------|
| Domain | 500 m × 500 m | UTM Zone 11N, S. California |
| Grid | 64 × 64 | 2-D build required |
| Fuel | FM4 chaparral | h = 6 ft (72 in), M_f = 8 % |
| Wind | 5 m/s eastward | 20-ft reference wind |
| `waf_formula` | `behaviorplus` | linear 0.36 + 0.004×h_in |
| `waf_canopy_alpha` | 1.5 | exponential sheltering coefficient |
| `use_waf` | 1 | WAF enabled |
| `use_wind_limit` | 1 | MEWS cap enabled |

## WAF Formula Comparison (FM4, h = 6 ft)

| Formula | WAF | Midflame wind (5 m/s input) |
|---------|-----|-----------------------------|
| Andrews logarithmic | 0.36 | 1.80 m/s |
| BehavePlus linear | 0.648 | 3.24 m/s |

The BehavePlus formula gives a higher midflame wind for deep fuel beds, producing
a faster spread rate than the Andrews formula at the same ambient wind speed.

## Expected Behaviour

- Fire spreads eastward at a higher rate than the `waf_andrews` test (for FM4).
- ROS is higher because the BehavePlus WAF is larger for deep chaparral fuel.
- The fire perimeter is more elongated in the downwind direction.

## Build

```bash
cmake -S . -B build -DLEVELSET_DIM_2D=ON
cmake --build build -j
```

## Run Command

```bash
./build/levelset regtest/wind/waf_behaviorplus/inputs.i
```

## References

- Andrews, P.L. (2018). *The Rothermel Surface Fire Spread Model and Associated
  Developments: A Comprehensive Explanation.* USDA For. Serv. Gen. Tech. Rep.
  RMRS-GTR-371.
- BehavePlus Fire Modeling System. USDA Forest Service.
  https://www.frames.gov/behaviorplus/home
- Massman, W.J. (1987). A comparative study of some mathematical models of the
  mean wind structure and aerodynamic drag of plant canopies.
  *Boundary-Layer Meteorology*, 40(1–2), 179–197.
