# Ember Accumulation and Probabilistic Ignition Test

## Purpose

Tests the **ember accumulation and probabilistic ignition** model that tracks landed ember density with exponential decay (burnout) and computes ignition probability based on accumulated coverage.

This demonstrates that:
1. `ember_accumulation.enable = 1` activates ember density tracking
2. Ember density accumulates from Albini spotting landing flux
3. Landed embers decay exponentially with burnout time constant
4. Probabilistic ignition occurs when accumulated density exceeds threshold
5. Fuel moisture dampens ignition probability
6. Diagnostic fields track ember density and ignition events

## Physical Model

### Ember Density Evolution

Ember density at each cell evolves according to:

```
ρ_ember(t+dt) = ρ_ember(t) + F_landing × dt - ρ_ember(t) × k_decay × dt
```

Where:
- `F_landing` = ember landing flux [embers/m²/s] from Albini spotting
- `k_decay` = burnout decay rate [1/s] (default: 1/600 = 10-min burnout)

### Probabilistic Ignition

At each check interval, ignition probability is computed:

```
P_ignition = 1 - exp(-λ)
λ = k_ignition × ρ_ember × dt_check / ρ_threshold
```

Where:
- `ρ_threshold` = minimum ember density for ignition [embers/m²]
- `k_ignition` = ignition rate coefficient
- `dt_check` = time between ignition checks [s]

### Moisture Damping

When fuel moisture M_f is available, ignition probability is reduced:

```
f_moisture = max(0, (M_x - M_f) / M_x)
P_ignition_eff = P_ignition × f_moisture
```

Higher moisture → lower ignition probability.

## Test Configuration

| Parameter | Value | Description |
|-----------|-------|-------------|
| Domain | 1000 m × 1000 m | Flat terrain |
| Grid | 64 × 64 cells | ~15.6 m resolution |
| Fuel model | FM4 chaparral | High-intensity fuel |
| Fuel moisture | 6% | Dry conditions |
| Wind | 8 m/s eastward | Strong wind for spotting |
| Simulation time | ~2 hours | Long duration for accumulation |
| **Albini Spotting** | | |
| Launch probability | 0.02 | Higher rate for testing |
| Terminal velocity | 1.2 m/s | Firebrand fall speed |
| Min Byram intensity | 50 kW/m | Threshold for lofting |
| Spot radius | 10 m | Ignition zone size |
| **Ember Accumulation** | | |
| Decay rate k_decay | 0.00167 s⁻¹ | 10-min burnout (1/600) |
| Density threshold | 8 embers/m² | Min for ignition |
| Ignition coefficient | 1.2 | Rate multiplier |
| Check interval | 30 s | Ignition evaluation period |
| Moisture damping | Enabled | M_f effect on P_ignition |

## Expected Behavior

1. **Main fire**: Advances eastward with level-set propagation
2. **Ember launch**: Albini model generates firebrands from high-intensity front cells
3. **Ember landing**: Firebrands land downwind, creating landing flux F_landing
4. **Density accumulation**: Repeated landings build up ember density ρ_ember
5. **Exponential decay**: Landed embers burn out with 10-min time constant
6. **Probabilistic ignition**: High-density spots ignite stochastically
7. **Secondary fires**: Ignited spots grow and merge with main fire

## Plotfile Diagnostics

The following fields are written to plotfiles:

- `ember_density`: Accumulated ember density [embers/m²]
- `ember_flux`: Instantaneous landing flux [embers/m²/s] from Albini
- `albini_Hz`: Lofting height [m]
- `albini_count`: Cumulative spot fire count
- `albini_dist`: Landing distance [m]
- `phi`: Level-set function showing fire perimeter

## Comparison Tests

**Without ember accumulation** (`ember_accumulation.enable = 0`):
- Albini spots ignite immediately upon landing
- No density tracking or probabilistic ignition
- All spots with P_catch > random threshold ignite

**With ember accumulation** (`ember_accumulation.enable = 1`):
- Ember density accumulates over time
- Single embers may not ignite immediately
- Repeated landings at same location increase ignition probability
- Embers burn out if not enough accumulate

## Files

- `inputs.i`: Test configuration
- `README.md`: This file

## Run Command

```bash
cd regtest/spotting/ember_accumulation
../../build/levelset inputs.i
```

> **Note**: Requires a 2D build (`-DLEVELSET_DIM_2D=ON`).

## References

- Sardoy, N., et al. (2007). Modeling transport and combustion of firebrands. *Combustion and Flame*, 150(3), 151-169.
- Koo, E., et al. (2010). Modelling firebrand transport in wildfires. *Proceedings of the Combustion Institute*, 33(2), 2449-2456.
- Albini, F.A. (1983). Transport of firebrands by line thermals. *Combustion Science and Technology*, 32(5-6), 277-288.
