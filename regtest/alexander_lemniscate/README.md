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
| θ = 0  (head fire)   | `coeff_a × R_base × dt` |
| θ = π/2 (flank fire) | `coeff_b × R_base × dt` |
| θ = π  (backing fire)| `coeff_c × R_base × dt` |

When `coeff_b = (coeff_a + coeff_c) / 2`, the formula is a true Limaçon and
the head/backing rates match `coeff_a` and `coeff_c` exactly.  If `coeff_b`
is set independently, the head and backing rates may differ from `coeff_a` and
`coeff_c` but the flank spread matches `coeff_b` exactly.  Values below zero
are clamped to zero (non-convex shapes are suppressed).

## Configuration

- **Domain**: 100 m × 100 m (UTM Zone 11N, Southern California)
- **Grid**: 64² cells
- **Initial condition**: Box source centred in domain
- **Wind**: constant 0.3 m/s in x-direction
- **FARSITE coefficients**: a = 1.0 (head), b = 0.6 (flank = arithmetic mean of a and c), c = 0.2 (backing)

## Expected Behavior

The fire spreads in a lemniscate / cardioid shape that is more elongated than a
circle but less sharply anisotropic than an ellipse.  Head fire spreads fastest
(rate `coeff_a`), flanks at rate `coeff_b`, and backing fire slowest (`coeff_c`).

## Run Command

```bash
./build/levelset regtest/alexander_lemniscate/inputs.i
```

## Reference

Alexander, M.E., Stocks, B.J., & Lawson, B.D. (1991). "Fire behaviour in
black spruce-lichen woodland: the Porter Lake project." *Information Report
NOR-X-310*, Canadian Forest Service, Edmonton, Alberta.
