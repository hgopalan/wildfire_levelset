# Catchpole & de Mestre (1986) Double-Ellipse Fire Shape Test

## Purpose

Tests the Catchpole & de Mestre (1986) double-ellipse fire spread shape model
(`farsite.fire_shape_model = catchpole_demestre`).

## Background

The double-ellipse model represents the Huygens wavelet emitted by each
fire-front point as two half-ellipses:

- **Head fire half** (downwind, cos θ ≥ 0): ellipse with semi-axes `a` (along
  wind) and `b` (crosswind). The polar spread rate is:

  `r(θ) = a·b / sqrt(b²·cos²θ + a²·sin²θ) × R_head × dt`

- **Backing fire half** (upwind, cos θ < 0): ellipse with semi-axes `c`
  (upwind) and `b` (crosswind).

  `r(θ) = c·b / sqrt(b²·cos²θ + c²·sin²θ) × R_head × dt`

Unlike the default Richards (1990) model (which uses a linear blend), this
model is geometrically exact for the double-ellipse shape, giving a sharper
transition at the flank.

## Configuration

- **Domain**: 100 m × 100 m (UTM Zone 11N, Southern California)
- **Grid**: 64² cells
- **Initial condition**: Box source centred in domain
- **Wind**: constant 0.3 m/s in x-direction
- **FARSITE coefficients**: a = 1.0 (head), b = 0.4 (flank), c = 0.2 (backing)

## Expected Behavior

The fire spreads as a double-ellipse.  Head spread is fastest; backing slowest.
Compared to the Richards (1990) model, the flanks are slightly more rounded and
the transition between head and backing is smoother.

## Run Command

```bash
./build/levelset regtest/catchpole_demestre/inputs.i
```

## Reference

Catchpole, E.A. & de Mestre, N.J. (1986). "Physical models for a spreading
line fire." *Australian Forestry*, 49(2), 102–111.
