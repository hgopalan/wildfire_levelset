# Alexander et al. Lemniscate (Limaçon) Fire Shape Test

## Purpose

Tests the Alexander et al. lemniscate / Limaçon fire spread shape model
(`farsite.fire_shape_model = lemniscate`).

## Background

The lemniscate (Limaçon) model describes the fire spread wavelet as a
cardioid-family curve in polar coordinates:

```
r(θ) = coeff_b + (coeff_a − coeff_c) / 2 · cosθ
```

where θ is the angle between the fire-front normal and the wind direction.

Properties:

| Angle | Spread rate |
|-------|-------------|
| θ = 0  (head fire)   | `(coeff_b + (coeff_a − coeff_c)/2) × R_base × dt` |
| θ = π/2 (flank fire) | `coeff_b × R_base × dt` |
| θ = π  (backing fire)| `(coeff_b − (coeff_a − coeff_c)/2) × R_base × dt` |

When `coeff_b = (coeff_a + coeff_c) / 2` (as in this regtest), the formula is
a true Limaçon: the head rate becomes exactly `coeff_a × R_base × dt` and the
backing rate becomes exactly `coeff_c × R_base × dt`.  For other values of
`coeff_b` the flank rate still equals `coeff_b` exactly, while the head and
backing rates are shifted by `coeff_b − (coeff_a + coeff_c)/2`.  Values below
zero are clamped to zero (non-convex shapes are suppressed).

## Configuration

- **Domain**: 100 m × 100 m (UTM Zone 11N, Southern California)
- **Grid**: 64² cells
- **Initial condition**: Box source centred in domain
- **Wind**: constant 0.3 m/s in x-direction
- **FARSITE coefficients**: a = 1.0 (head), b = 0.6 (flank = arithmetic mean of a and c), c = 0.2 (backing)

## Expected Behavior

The fire spreads in a lemniscate / cardioid shape.  Because `coeff_b =
(coeff_a + coeff_c) / 2 = 0.6`, the head fire spreads at `coeff_a × R_base =
1.0 × R_base`, the flanks at `0.6 × R_base`, and the backing fire at `0.2 ×
R_base`.

## Run Command

```bash
./build/levelset regtest/alexander_lemniscate/inputs.i
```

## Reference

Alexander, M.E., Stocks, B.J., & Lawson, B.D. (1991). "Fire behaviour in
black spruce-lichen woodland: the Porter Lake project." *Information Report
NOR-X-310*, Canadian Forest Service, Edmonton, Alberta.
