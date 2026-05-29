# McArthur Moisture Scaling Regression Test

## Purpose

Tests the **McArthur-style temperature/RH-dependent moisture response time scaling** feature that adjusts fuel moisture drying/wetting rates based on ambient conditions.

This demonstrates that:
1. `precipitation_moisture.use_mcarthur_scaling = 1` activates the McArthur temperature/RH adjustment
2. Moisture response times are scaled based on temperature and relative humidity
3. Higher temperatures → faster drying (reduced τ)
4. Higher relative humidity → slower drying (increased τ)
5. The scaling affects fire spread by modifying fuel moisture evolution after precipitation

## Test Configuration

| Parameter | Value | Description |
|-----------|-------|-------------|
| Domain | 1000 m × 1000 m | Flat terrain, UTM Zone 11N |
| Grid | 64 × 64 cells | ~15.6 m resolution |
| Initial precipitation | 10 mm | Recent rainfall raises moisture |
| Temperature | 25°C | Moderately warm |
| Relative humidity | 30% | Moderately dry |
| Base τ (1-hr fuel) | 3600 s (1 hour) | Standard drying time |
| Base τ (10-hr fuel) | 14400 s (4 hours) | |
| Base τ (100-hr fuel) | 86400 s (24 hours) | |
| Fuel model | FM4 chaparral | Southern California fuel type |
| Wind | 5 m/s eastward | Moderate wind |
| Simulation time | ~4 hours | Sufficient to observe moisture decay |

## Expected Behavior

1. **Initial state**: Fire starts in fuel with elevated moisture from recent rain
2. **Moisture evolution**: Fuel moisture decays exponentially toward equilibrium
3. **McArthur acceleration**: Decay rate is modified by T/RH conditions
4. **Fire spread**: ROS increases as fuel moisture decreases during simulation
5. **Plotfile output**: `moisture_d1`, `moisture_d10`, `moisture_d100` fields show spatial moisture

## Comparison Tests

To verify McArthur scaling, compare runs with different settings:

**Without McArthur scaling** (`use_mcarthur_scaling = 0`):
- Fuel moisture decays with unmodified time constants

**With McArthur scaling** (`use_mcarthur_scaling = 1`):
- Fuel moisture decay rate adjusted for temperature and humidity

**Different T/RH scenarios**:
- Hot, dry (T = 35°C, RH = 15%): Faster drying
- Cool, humid (T = 5°C, RH = 80%): Slower drying

## Files

- `inputs.i`: Test configuration
- `README.md`: This file

## Run Command

```bash
cd regtest/moisture/mcarthur_scaling
../../build/levelset inputs.i
```

> **Note**: Requires a 2D build (`-DLEVELSET_DIM_2D=ON`).

## References

- McArthur, A.G. (1967). Fire behaviour in eucalypt forests. Commonwealth of Australia Forestry and Timber Bureau Leaflet 107.
- Nelson, R.M. (2000). Prediction of diurnal change in 10-h fuel stick moisture content. *Canadian Journal of Forest Research*, 30(7), 1071-1087.
