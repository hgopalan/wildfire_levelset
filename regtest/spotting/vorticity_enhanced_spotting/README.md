# Vorticity-Enhanced Spotting Regtest

Tests the integration of the Weise & Biging (1996) fire whirl model with the Albini (1983) 
firebrand spotting model to produce vorticity-enhanced spotting distances.

## What is tested

- **Fire whirl model**: `weise_biging.enable = 1` computes fire whirl characteristics:
  - Whirl height, radius, angular velocity, tangential velocity
  - Based on fireline intensity, flame length, and wind speed
  
- **Vorticity-enhanced spotting**: `weise_biging.enhance_spotting = 1` couples whirl effects to Albini spotting:
  - Enhancement factor: `1 + α × (v_θ/U) × (Ω×r_w/U)`
  - Multiplies base Albini spotting distance
  - Clipped to range [1.0, 5.0] to avoid extreme values
  
- **Albini spotting model**: Same as base `albini_spotting` regtest:
  - Lofting height: `H_z = 12.2 × I_B^(1/3)` [m, kW/m]
  - 2-D trajectory integration through wind field
  - Stochastic firebrand launch

- **Diagnostic output**: All Weise & Biging fields are written to plotfiles:
  - `weise_flame_height`, `weise_flame_tilt`
  - `weise_whirl_height`, `weise_whirl_radius`
  - `weise_angular_velocity`, `weise_max_tang_vel`

## Domain

500 m × 500 m with 64×64 cells (dx ≈ 7.8 m).  
Wind: 8 m/s eastward (higher than base Albini test to generate stronger whirls).  
Propagation: Level-set method.

## Expected result

- Spot fires appear downwind of the main fire front
- **Enhanced spotting distance** compared to base Albini test (without whirl enhancement)
- With `random_seed = 42` the result is deterministic
- Weise & Biging diagnostic fields show non-zero whirl characteristics at fire front

## Physics verification (FM4, u = 8 m/s, v_t = 1 m/s, α = 1.0)

Assuming typical values at the fire front:
- R ≈ 120 ft/min (higher due to stronger wind)
- I_B ≈ 12000 kW/m
- H_z ≈ 269 m (from Albini formula)
- Flight time ≈ 269 s
- Base spotting distance ≈ 8 m/s × 269 s ≈ 2152 m

Fire whirl (Weise & Biging):
- Whirl height H_w ≈ 269 m (matches vertical flame height)
- Whirl radius r_w ≈ 0.1 × 269 m ≈ 27 m
- Tangential velocity v_θ ≈ 8 m/s (same order as wind)
- Enhancement factor ≈ 1 + 1.0 × (8/8) × (1×27/8) ≈ 1 + 3.375 ≈ 4.375 (clipped to 5.0 max)

**Enhanced spotting distance ≈ 2152 m × min(4.375, 5.0) ≈ 9400 m** (far exceeds domain, will hit boundary)

In practice, the enhancement will vary across the fire front based on local intensity and wind conditions.
