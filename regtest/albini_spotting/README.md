# Albini Spotting Regtest

Tests the Albini (1983) firebrand spotting model with 2-D trajectory integration.

## What is tested

- `albini_spotting.enable = 1` activates the new model.
- Byram's fire line intensity is derived from the Rothermel rate-of-spread field.
- Albini (1983) lofting height formula: `H_z = 12.2 * I_B^(1/3)` [m, kW/m].
- Stochastic firebrand launch with intensity-weighted probability `P_base`.
- 2-D forward-Euler particle trajectory integrates horizontal wind `(u, v)` for
  flight time `t_f = H_z / v_t`.
- Spot ignitions are placed as circles of radius `spot_radius` in the level-set
  field `phi`.
- Diagnostic fields `albini_Hz`, `albini_count`, `albini_dist`, `albini_active`
  are written to every plot file.

## Domain

1 km × 1 km with 100×100 cells (dx = 10 m).  Wind is 5 m/s in x + 1 m/s in y.
FARSITE ellipse spread is used (`skip_levelset = 1`).

## Expected result

Spot fires appear downwind of the main fire front.  With `random_seed = 42` the
result is deterministic.
