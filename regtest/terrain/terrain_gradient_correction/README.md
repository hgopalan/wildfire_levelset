# Terrain Gradient Correction Test

## Purpose

Verifies that the level-set gradient magnitude `|∇φ|` is computed on the
terrain surface rather than on the flat Cartesian plane.

When a terrain file is present, `godunov_norm_grad_phi` replaces the flat
grid spacings with terrain arc-length spacings:

```
eff_dx = dx * sqrt(1 + (dz/dx)^2)
eff_dy = dy * sqrt(1 + (dz/dy)^2)
```

This test uses a **steep Gaussian hill** (H = 200 m, σ = 100 m) where the
slope magnitude at the inflection ring (r = σ) reaches ≈ 1.21 (about 50°),
giving an effective spacing factor of `sqrt(1 + 1.21²) ≈ 1.57`.  Without the
correction the gradient would be over-estimated by 57 % on the steepest flanks,
causing the fire to spread too slowly upslope.

## Configuration

| Parameter      | Value                        |
|----------------|------------------------------|
| Domain         | 1000 m × 1000 m              |
| Grid           | 64 × 64 cells (≈ 15.6 m res) |
| Terrain        | Gaussian hill H=200 m, σ=100 m |
| Max slope      | ≈ 1.21 (50°) at r = σ        |
| Gradient correction factor | up to 1.57× |
| Wind           | 5 m/s westerly               |
| Fuel           | FM4 Chaparral, M_f = 8 %     |
| Ignition       | Sphere r = 50 m at western edge |
| Propagation    | `levelset` (exercises terrain-corrected WENO3 gradient) |

## Terrain Model

The Gaussian hill is:

```
z(x,y) = H * exp(-r² / (2σ²))
```

where `r² = (x - x_c)² + (y - y_c)²`, `H = 200 m`, `σ = 100 m`.

The terrain slope components are:

```
dz/dx = -H (x - x_c) / σ² * exp(-r² / (2σ²))
dz/dy = -H (y - y_c) / σ² * exp(-r² / (2σ²))
```

Maximum slope magnitude ≈ H/σ · e^(-1/2) ≈ 1.21.

## Expected Behaviour

- Fire starts as a sphere on the western edge.
- Wind drives it eastward toward and up the steep hill.
- With terrain-corrected gradients, `|∇φ|` on the steep flanks is
  smaller (≈ 1/1.57 of the flat value), correctly reflecting the longer
  surface arc between grid nodes.
- The solver should complete 100 steps and write 4 plotfiles without
  NaN/Inf values or negative time steps.

## Files

| File                        | Description                        |
|-----------------------------|------------------------------------|
| `inputs.i`                  | Solver input parameters            |
| `steep_gaussian_terrain.csv`| 1600-point XYZ terrain elevation   |
| `README.md`                 | This file                          |

## Run Command

```bash
cd regtest/terrain/terrain_gradient_correction
../../../build/levelset inputs.i
```

## Notes

- Requires a 2-D build: `cmake -DLEVELSET_DIM_2D=ON`
- The `reinitialization` path in `numerical_schemes.H` still uses flat
  Cartesian spacings (default slope = 0) — this is intentional because
  the signed-distance property is defined in the horizontal plane.
