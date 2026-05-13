# ember_cascade_flux — Regression Test

Tests the flux-based ember cascade model (`ember_cascade_flux.H`) with
3-D wind data from a synthetic massconsistent_amr plotfile.

## Purpose

Verifies that:

1. The Albini (1983) plume-lofting height `H_z = 12.2 × I_B^(1/3)` is
   correctly computed per fire-front cell.
2. The Gaussian transport kernel accumulates a landing-flux density
   (`ember_cascade_flux` plotfile field) over the downwind domain.
3. The Poisson ignition threshold (`N_min_density`) correctly gates which
   cells receive spot-fire ignitions (`ember_cascade_ignition` field).
4. Height-averaged wind from the 3-D plt file is read and bilinearly
   interpolated onto the fire-model grid for the mean landing displacement.

## Test scenario

| Parameter | Value |
|---|---|
| Domain | 600 m × 600 m (UTM Zone 11N) |
| Grid | 64 × 64, max_grid=32 |
| Fuel | FM4 chaparral (M_f = 7 %) |
| Wind | u = 6 m/s, v = 1 m/s |
| Ignition | Sphere, radius 25 m |
| `I_B_min` | 10 kW/m |
| `k_flux` | 1.0 |
| `sigma_base` | 30 m |
| `k_sigma` | 0.1 |
| `N_min_density` | 5×10⁻⁴ embers/m²/s |

## Running manually

```bash
cd regtest/spotting/ember_cascade_flux
python3 generate_plt_wind.py          # creates plt_wind_3d/ directory
../../../build/levelset inputs.i
```

## Running via CTest

The Python setup step is registered as a CTest fixture, so it runs
automatically before the solver step:

```bash
cd build
ctest -R ember_cascade_flux --output-on-failure
```

## plt file format

`generate_plt_wind.py` writes a minimal single-level AMReX plotfile (`plt_wind_3d/`)
on an 8 × 8 × 4 grid with variables `u`, `v`, `w`, `vel_magnitude`.
The reader in `src/plt_wind_reader.H` column-averages `u` and `v` over the
vertical levels to produce a 2-D height-averaged wind field, which is then
bilinearly interpolated onto the fire-model cell centres.

## Output fields

| Field | Description |
|---|---|
| `ember_cascade_flux` | Accumulated ember landing-flux density [embers/m²/s] |
| `ember_cascade_ignition` | 1.0 where a spot fire was placed this check step |
