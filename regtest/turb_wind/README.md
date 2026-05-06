# Turbulent Wind Perturbation Regression Test

## Purpose

Tests the **Random Fourier Feature (RFF) spectral noise** turbulent wind
perturbation model with Ornstein-Uhlenbeck temporal evolution
(`turb_wind.model = spectral_noise`).

This demonstrates that:
1. `N_modes = 32` wavenumber pairs are drawn at initialisation from the 2-D
   isotropic Gaussian power spectrum with length scale `L_c = 100 m`.
2. Each mode has a scalar OU amplitude that evolves with decorrelation time
   `1/theta = 20 s`; modes are independent so the composite field is both
   spatially structured and temporally correlated.
3. The perturbation field is reconstructed on the GPU at every timestep as a
   cosine superposition weighted by the current OU amplitudes.
4. The perturbed wind is fed into the Rothermel ROS computation, causing
   the fire perimeter to evolve with realistic stochastic variability around
   the mean 5 m/s eastward base flow.

## Available turbulent wind models

| `turb_wind.model` | Description |
|---|---|
| `none` | No perturbation (deterministic baseline) |
| `ou_process` (L_c=0) | Domain-uniform OU: all cells get the same gust |
| `ou_process` (L_c>0) | Per-cell OU driven by Gaussian-smoothed white noise |
| `spectral_noise` | RFF spectral noise with OU amplitudes (this test) |
| `direction_walk` | Bounded direction random walk; speed preserved |

## Physical Model

### Random Fourier Feature (RFF) spectral noise

At initialisation `N_modes` wavenumber pairs and phases are sampled once:

```
kx_n, ky_n  ~  N(0, (1/L_c)^2)   (isotropic Gaussian PSD)
phi_u_n, phi_v_n  ~  U[0, 2π)
A_u_n = A_v_n = 0
```

At every timestep the OU amplitudes are updated on the CPU:

```
alpha = exp(-theta * dt)
A_u_n(t+dt) = alpha * A_u_n(t) + sqrt(1 - alpha^2) * N(0,1)
A_v_n(t+dt) = alpha * A_v_n(t) + sqrt(1 - alpha^2) * N(0,1)
```

The perturbation field is then reconstructed on the GPU:

```
u'(x,y) = sigma * sqrt(2/N) * Σ_n  A_u_n * cos(kx_n*x + ky_n*y + phi_u_n)
v'(x,y) = sigma * sqrt(2/N) * Σ_n  A_v_n * cos(kx_n*x + ky_n*y + phi_v_n)
```

This gives:
- **Temporal decorrelation time** = `1/theta = 20 s`
- **Stationary std dev** of each cell's perturbation ≈ `sigma = 1.0 m/s`
- **Approximate 2-D autocorrelation**: `C(r) ≈ sigma^2 * exp(-r^2 / (2*L_c^2))`
- **Gaussian power spectrum**: energy concentrated at wavenumbers `k ~ 1/L_c`

## Configuration

| Parameter | Value | Description |
|-----------|-------|-------------|
| Domain | 1000 m × 1000 m | Flat terrain |
| Grid | 100 × 100 cells | 10 m resolution |
| Base wind | 5 m/s east | Constant background |
| `turb_wind.theta` | 0.05 s⁻¹ | OU reversion → 20 s decorrelation time |
| `turb_wind.sigma` | 1.0 m/s | Perturbation std dev per cell |
| `turb_wind.L_c` | 100 m | Spatial correlation length (= 10 cells) |
| `turb_wind.N_modes` | 32 | Number of Fourier modes |
| `turb_wind.random_seed` | 42 | Fixed seed for reproducibility |
| Fire spread | Rothermel FM4 | Southern California chaparral |
| Propagation | Level-set | `propagation_method = levelset` |
| Simulation time | 1200 s | |

## Expected Behaviour

- Fire starts as a sphere of radius 30 m at (150, 500) m.
- Mean wind drives fire eastward; mean ROS ≈ Rothermel prediction for FM4 at 5 m/s.
- Wind gusts (OU fluctuations ±1–2 σ) cause the fire perimeter to deviate
  from the perfectly-elliptic shape seen with constant wind.
- Perimeter roughness reflects the spatial structure of the spectral modes:
  smooth at scales ≫ L_c, irregular below L_c.
- The `velx`/`vely` plotfile fields show the perturbed wind (base + u', v').

## Files

- `inputs.i`: Input parameters

## Run Command

```bash
cd regtest/turb_wind
../../build/levelset inputs.i
```

> **Note**: Requires a 2D build (`-DLEVELSET_DIM_2D=ON`).

## Comparison Runs

Replace the `turb_wind.*` block to compare models:

```
# Baseline (no turbulence)
# [remove all turb_wind lines]

# Domain-uniform OU
turb_wind.model = ou_process
turb_wind.theta = 0.05
turb_wind.sigma = 1.0

# OU + Gaussian spatial kernel
turb_wind.model = ou_process
turb_wind.theta = 0.05
turb_wind.sigma = 1.0
turb_wind.L_c = 100.0

# Direction walk
turb_wind.model = direction_walk
turb_wind.sigma_theta = 0.1
turb_wind.theta_max = 0.5236
```

## References

- Kraichnan, R.H. (1970). Diffusion by a random velocity field.
  *Phys. Fluids* 13(1):22–31.
- Rahimi, A. & Recht, B. (2007). Random features for large-scale kernel
  machines. *NIPS 2007*.
- Uhlenbeck, G.E. & Ornstein, L.S. (1930). On the theory of Brownian motion.
  *Phys. Rev.* 36(5):823–841.
- Finney, M.A. et al. (2011). Role of wind, fuel moisture, and terrain in
  controlling fire movement. *Ecosphere* 2(3):art17.
