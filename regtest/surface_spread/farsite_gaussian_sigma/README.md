# FARSITE Gaussian-Smoothed Spread-Point Stamping Test

## Purpose

Tests the `farsite.gaussian_sigma` option, which replaces the default
single-cell phi stamping with a signed-distance-function (SDF) disk of
radius *sigma* centred on each propagated fire-front position.  This
produces a smooth, differentiable phi field at each FARSITE time step —
analogous to the `fire_gaussian_sigma` / `fire_points_file` initialisation
approach — reducing numerical noise along the advancing fire perimeter.

## Configuration

- **Domain**: 1000 m × 1000 m (50 × 50 cells, Δx = Δy = 20 m)
- **Grid**: 50 × 50 cells
- **Initial condition**: Box ignition (x: 50–100 m, y: 300–700 m)
- **Wind**: Constant (u_x = 5 m/s)
- **Time steps**: 20
- **Fuel model**: FM4 (chaparral), M_f = 0.08
- **FARSITE parameters**: Richards (1990), a = 1.0, b = 0.4, c = 0.2, L/W = 3.0
- **Gaussian sigma**: `farsite.gaussian_sigma = 40.0` m (2 cells)

## Expected Behaviour

The fire advances eastward in an elongated ellipse.  With smoothing enabled,
the stamped phi field shows a continuous negative region surrounded by a
smooth zero contour, rather than isolated single-cell markers.  The run
should complete without errors.

## Run Command

```bash
./build/levelset regtest/surface_spread/farsite_gaussian_sigma/inputs.i
```
