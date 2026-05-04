# Cheney & Gould (1995 / 1998) Grassland Fire Spread Test

## Purpose

Tests fire spread using the Cheney & Gould empirical grassland fire spread model.
This model was calibrated against experimental grassland fires in Australia and
provides better agreement with observed grassland spread rates than Rothermel in
open grass fuels.

## Model Overview

The head-fire rate of spread is a piecewise-linear function of 10-m open wind
speed with moisture and curing corrections:

- **For U₁₀ ≤ 5 km/h**: `R = (0.165 + 0.534 × U₁₀) × exp(−0.108 × MC) × CF`
- **For U₁₀ > 5 km/h**: `R = (−0.020 + 0.640 × U₁₀) × exp(−0.108 × MC) × CF`

where `MC` is dead fine fuel moisture content [%] and `CF` is the curing fraction.

## Configuration

- **Domain**: 1000 m × 1000 m flat grassland patch
- **Grid**: 64 × 64 × 64 cells
- **Initial condition**: Sphere at domain centre (500 m, 500 m), radius 20 m
- **Wind**: 10 m/s from the west (= 36 km/h 10-m open wind)
- **Fuel moisture** (`cheney_gould.moisture`): 8% dead fine fuel moisture
- **Curing** (`cheney_gould.curing`): 1.0 (fully cured grass)
- **Simulation time**: 600 s
- **Output interval**: Every 20 steps
- **Propagation method**: `levelset`

## Expected Behavior

With U₁₀ = 36 km/h (well into the high-wind regime), MC = 8%, and CF = 1.0:

```
f_MC = exp(−0.108 × 8) ≈ 0.420
R_kmh = (−0.020 + 0.640 × 36) × 0.420 ≈ 9.67 km/h
R_ms  ≈ 2.69 m/s
```

The fire should expand asymmetrically, spreading rapidly downwind (east) and
slowly upwind (west), forming an elongated elliptical perimeter.

## Variants

Modify these parameters to explore model sensitivity:

| Parameter | Effect |
|-----------|--------|
| `cheney_gould.moisture` | Higher moisture → lower ROS (exponential dampening) |
| `cheney_gould.curing`   | Lower curing → proportionally lower ROS |
| `u_x`                   | Wind speed controls regime (low ≤ 5 km/h vs high > 5 km/h) |

Low-wind example (U₁₀ < 5 km/h, i.e. u_x < 1.39 m/s):
```
cheney_gould.moisture = 15.0
cheney_gould.curing   = 0.80
u_x = 1.0
```

## Run Command

```bash
./build/levelset regtest/cheney_gould_grassfire/inputs.i
```

## References

- Cheney, N.P. & Gould, J.S. (1995). "Fire growth in grassland fuels."
  *Int. J. Wildland Fire*, 5(4), 237–247.

- Cheney, N.P., Gould, J.S. & Catchpole, W.R. (1998). "Prediction of fire
  spread in grasslands." *Int. J. Wildland Fire*, 8(1), 1–13.
