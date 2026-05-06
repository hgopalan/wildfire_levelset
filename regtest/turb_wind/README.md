# Turbulent Wind Perturbation Regression Test

## Purpose

Tests the Ornstein-Uhlenbeck (OU) turbulent wind perturbation model with a
Gaussian spatial correlation kernel (`turb_wind.model = ou_process`).

This demonstrates that:
1. The OU process produces temporally correlated wind gusts with decorrelation
   time `1/theta = 1/0.05 = 20 s`.
2. The Gaussian spatial kernel (`L_c = 100 m`) produces spatially coherent
   perturbations — neighbouring cells receive similar gusts while cells
   separated by more than `L_c` receive nearly independent perturbations.
3. The perturbed wind is correctly combined with the base wind and fed into
   the Rothermel ROS computation each timestep.
4. The fire perimeter evolves realistically with stochastic wind variability
   around the mean 5 m/s eastward base flow.

## Physical Model

### Ornstein-Uhlenbeck process

The wind perturbation `(u', v')` follows the discrete-time OU update:

```
alpha      = exp(-theta * dt)
sigma_step = sigma * sqrt(1 - alpha^2)
u'(t+dt)   = alpha * u'(t) + sigma_step * eta_u(t)
v'(t+dt)   = alpha * v'(t) + sigma_step * eta_v(t)
```

- **Stationary std dev** of each cell's perturbation = `sigma` = 1.0 m/s
- **Temporal decorrelation time** = `1/theta` = 20 s
- The perturbed wind is `u_eff = u_base + u'`, `v_eff = v_base + v'`

### Gaussian spatial correlation kernel

When `L_c > 0` the OU noise term `eta` is spatially correlated:

1. Generate per-cell standard normal white noise `xi(x,y)`.
2. Apply separable 1-D Gaussian smoothing in `x` then `y` with kernel
   standard deviation `sigma_k = L_c / dx = 100 / 10 = 10 cells`.
3. Rescale the smoothed field to unit RMS variance so `sigma` remains the
   stationary standard deviation of each cell's perturbation.
4. Use the rescaled noise as `eta` in the OU update above.

Spatial autocorrelation function of the perturbation field:
```
C(r) ≈ sigma^2 * exp(-r^2 / (2 * L_c^2))
```
(Gaussian envelope due to the separable Gaussian kernel.)

## Configuration

| Parameter | Value | Description |
|-----------|-------|-------------|
| Domain | 1000 m × 1000 m | Flat terrain |
| Grid | 100 × 100 cells | 10 m resolution |
| Base wind | 5 m/s east | Constant background |
| `turb_wind.theta` | 0.05 s⁻¹ | Reversion rate → 20 s decorrelation time |
| `turb_wind.sigma` | 1.0 m/s | Stationary perturbation std dev |
| `turb_wind.L_c` | 100 m | Spatial correlation length (= 10 cells) |
| `turb_wind.random_seed` | 42 | Fixed seed for reproducibility |
| Fire spread | Rothermel FM4 | Southern California chaparral |
| Propagation | Level-set | `propagation_method = levelset` |
| Simulation time | 1200 s | |

## Expected Behaviour

- Fire starts as a sphere of radius 30 m at (150, 500) m.
- Mean wind drives fire eastward; ROS ≈ Rothermel prediction for FM4 at 5 m/s.
- Wind gusts (OU fluctuations up to ±1–2 σ) cause the fire perimeter to
  deviate slightly from the perfectly-elliptic shape seen with constant wind.
- Perimeter is spatially irregular at scales below L_c = 100 m; smooth at
  larger scales (consistent with the spatial kernel).
- The `velx` / `vely` plotfile fields show the perturbed wind (base + u', v').

## Files

- `inputs.i`: Input parameters

## Run Command

```bash
cd regtest/turb_wind
../../build/levelset inputs.i
```

> **Note**: Requires a 2D build (`-DLEVELSET_DIM_2D=ON`).

## Comparison with Constant Wind

Run with `turb_wind.model = none` (or remove the `turb_wind.*` lines) to
see the perfectly smooth elliptical fire perimeter from the constant-wind
base case.  The turbulent run should produce a visibly rougher perimeter
with the same overall eastward spread rate.

## References

- Uhlenbeck, G.E. & Ornstein, L.S. (1930). On the theory of Brownian motion.
  *Phys. Rev.* 36(5):823–841.
- Finney, M.A. et al. (2011). Role of wind, fuel moisture, and terrain in
  controlling fire movement. *Ecosphere* 2(3):art17.
- Cruz, M.G. et al. (2019). Uncertainty in wildfire behaviour research.
  *Current Forestry Reports* 5:155–172.
