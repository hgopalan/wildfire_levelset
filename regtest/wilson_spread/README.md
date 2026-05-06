# Wilson (1988) Ellipse-from-Focus Fire Shape Test

## Purpose

Tests the Wilson (1988) single-ellipse fire spread shape model
(`farsite.fire_shape_model = wilson`).

## Background

Wilson's model places the fire's point of origin at the **rear focus** of a
single ellipse (rather than at the centre, as in the double-ellipse model).
This is physically motivated: from a single ignition point, the fire perimeter
at time *t* is an ellipse whose far end (head fire) is displaced by `R_H·t`
and whose near end (backing fire) is displaced by `R_B·t`.

The spread rate at angle θ from the wind direction is given by the conic-section
polar equation from the rear focus:

```
r(θ) = ℓ / (1 − e · cosθ)
```

where:

- `R_H = coeff_a × R_base × dt` (head spread)
- `R_B = coeff_c × R_base × dt` (backing spread)
- `a_ell = (R_H + R_B) / 2`
- `e = (R_H − R_B) / (R_H + R_B)`  (eccentricity)
- `ℓ = a_ell (1 − e²) = 2 R_H R_B / (R_H + R_B)` (semi-latus rectum)

The flank spread (at θ = π/2) equals the harmonic mean of the head and backing
rates: `ℓ = 2·R_H·R_B / (R_H + R_B)`.  The `coeff_b` parameter is not used by
this model.

## Configuration

- **Domain**: 100 m × 100 m (UTM Zone 11N, Southern California)
- **Grid**: 64² cells
- **Initial condition**: Box source centred in domain
- **Wind**: constant 0.3 m/s in x-direction
- **FARSITE coefficients**: a = 1.0 (head), c = 0.2 (backing)

## Expected Behavior

The fire spreads as a single ellipse.  Compared to the Richards (1990) model,
the fire shape is physically more consistent: the ignition point sits exactly
at the rear focus and the flank spread is derived rather than independently
specified.  Head and backing fire rates match `coeff_a` and `coeff_c` exactly.

## Run Command

```bash
./build/levelset regtest/wilson_spread/inputs.i
```

## Reference

Wilson, A.A.G. (1988). "Width of firebreak that is necessary to stop grass
fires – some field experiments and a simulation model." *Canadian Journal of
Forest Research*, 18(6), 682–687.
