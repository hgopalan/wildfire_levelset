# Balbi + Viegas (Six Options) + Heat Flux Regression Test

## Purpose

Tests the following new capabilities:

1. **Balbi (2009) + Viegas diagnostics**: When `fire_spread_model = balbi`, the
   Viegas diagnostic uses the Balbi amplitude coefficient `A` and buoyancy velocity
   `v_b` instead of Rothermel's `R_0` and wind factor `phi_w`. The Viegas ROS is:
   ```
   R_V = A * (1 + sin(alpha_w) - cos(alpha_w)) * Phi_s_V
   ```
   where `alpha_w` is the wind-only flame tilt angle and `Phi_s_V = exp(a_V * tan_phi)`
   is the Viegas exponential slope enhancement factor.

2. **Option 2 (viegas_ros) with Balbi**: `wind_terrain.model = viegas_ros` overrides
   the Balbi ROS with `max(R_balbi, R_viegas_balbi)` in eruptive cells
   (`tan_phi > viegas.tan_phi_c`).

3. **Heat flux MultiFab with Balbi v_b augmentation**: When `heat_flux.enable_upward = 1`,
   the Balbi buoyancy velocity is augmented by the fire-induced buoyancy:
   ```
   v_b_Q  = k_upward * sqrt(g * Q * ref_height / (rho_air * Cp_air * T_a))
   v_b_eff = sqrt(v_b_fuel^2 + v_b_Q^2)
   ```
   A larger `v_b` makes the flame more vertical (less forward tilt), reducing
   radiant heat to unburned fuel and thereby reducing ROS for strongly buoyant fires.

## Configuration

- **Domain**: 1000 m × 1000 m (Gaussian hill, same as `terrain_wind` test)
- **Grid**: 100 × 100 cells (10 m resolution)
- **Fire spread model**: Balbi (2009)
- **Propagation**: level-set
- **Wind-terrain**: Option 2 (viegas_ros, Balbi baseline)
- **Heat flux**: 5000 W/m² uniform (moderate crown fire)
- **Fuel**: FM4 Chaparral
- **Viegas**: enabled (automatically by `wind_terrain.model = viegas_ros`)

## How to Switch Wind-Terrain Options

Edit `inputs.i` and change `wind_terrain.model` to any of:

| Value | Option | Description |
|-------|--------|-------------|
| `viegas_ros` | 2 | Replace Balbi ROS with Viegas-Balbi ROS in eruptive cells |
| `viegas_wind` | 3 | Add Viegas buoyancy wind in eruptive cells only |
| `canyon_wind` | 4 | Rothermel (1983) canyon wind amplification |
| `viegas_neto` | 5 | Viegas & Neto (1994) buoyancy wind at all cells |
| `pimont` | 6 | Pimont et al. (2009) exponential slope correction |
| `windninja_ridge_canyon` | 7 | WindNinja ridge speed-up and canyon channeling |

## Expected Behavior

- **Without heat flux**: Fire spreads according to Balbi ROS with Viegas slope enhancement.
  Cells on the Gaussian hill crest (high slope) are flagged as eruptive and R is increased.
- **With heat flux (5000 W/m²)**: Balbi `v_b` is augmented by ~0.9 m/s (at 300 K ambient),
  making the flame slightly more vertical. This reduces forward radiation and slightly
  reduces ROS on low-slope areas. The effect diminishes as slope increases (slope tilt
  dominates over wind tilt on steep terrain).
- **Diagnostic fields**: `viegas_ROS` shows Balbi-based Viegas ROS; `viegas_eruptive_flag`
  marks slope > 22°; `viegas_ROS_excess` is positive where Viegas exceeds Balbi.

## Files

- `inputs.i`: Input parameters
- `gaussian_hill_terrain.csv`: Terrain elevation data (2500 points, reused from terrain_wind test)
- `gaussian_hill_wind.csv`: Wind velocity data (2500 points, reused from terrain_wind test)

## Run Command

```bash
cd regtest/balbi_viegas_heatflux
../../build/levelset inputs.i
```

> **Note**: This test requires a 2D build (`-DLEVELSET_DIM_2D=ON`).

## References

- Balbi, J.H. et al. (2009). A physical model for wildland fires. *Combustion and Flame*, 156(12), 2217-2230.
- Viegas, D.X. (2004). Slope and wind effects on fire propagation. *Int. J. Wildland Fire*, 13(2), 143-156.
- WindNinja: Forthofer, J.M. (2007). Modeling Wind in Complex Terrain for Use in Fire Spread Prediction. Colorado State University MS thesis.
