# Periodic Wind Gust Factor Test

## Purpose

Tests the **periodic wind gust factor** feature that modulates the base wind field with a sinusoidal pattern representing thermal turbulence or atmospheric oscillations.

This demonstrates that:
1. `turb_wind.enable = 1` with periodic gust settings activates wind modulation
2. Wind speed varies sinusoidally: `V(t) = V_base × (1 + A × sin(2πt/T))`
3. Fire spread rate responds to gust cycles (accelerating/decelerating)
4. Fire perimeter exhibits asymmetry from time-varying wind
5. Multiple gust periods can be simulated

## Physical Model

### Periodic Gust Modulation

The base wind velocity is multiplied by a time-varying factor:

```
gust_factor(t) = 1 + A × sin(2π × t / T + φ)
V(t) = V_base × gust_factor(t)
```

Where:
- `A` = gust amplitude (dimensionless, typically 0.2–0.5)
- `T` = gust period [s] (e.g., 600 s = 10 min for thermal gusts)
- `φ` = initial phase [radians]

**Example**: With V_base = 5 m/s, A = 0.4, T = 600 s:
- Maximum wind: V_max = 5 × 1.4 = 7.0 m/s (at t = T/4)
- Minimum wind: V_min = 5 × 0.6 = 3.0 m/s (at t = 3T/4)
- Mean wind: V_mean = 5.0 m/s

## Test Configuration

| Parameter | Value | Description |
|-----------|-------|-------------|
| Domain | 1000 m × 1000 m | Flat terrain |
| Grid | 64 × 64 cells | ~15.6 m resolution |
| Base wind | 5 m/s eastward | Mean wind speed |
| Fuel model | FM4 chaparral | Southern California fuel |
| Fuel moisture | 8% | Moderate conditions |
| Simulation time | ~2 hours | Covers multiple gust cycles |
| **Periodic Gust** | | |
| Amplitude A | 0.4 (40%) | Wind varies ±40% around mean |
| Period T | 600 s (10 min) | Thermal turbulence timescale |
| Phase φ | 0 radians | Start at mean wind |
| Min wind | 3.0 m/s | Occurs periodically |
| Max wind | 7.0 m/s | Occurs periodically |

## Expected Behavior

1. **Fire ignition**: Circular fire starts at domain center
2. **Gust cycles**: Wind oscillates between 3.0 and 7.0 m/s
3. **ROS variation**: Fire spread accelerates and decelerates with wind
4. **Perimeter asymmetry**: Fire shape reflects gust phase
5. **Net effect**: Mean spread rate affected by non-linear ROS-wind relationship

## Files

- `inputs.i`: Test configuration
- `README.md`: This file

## Run Command

```bash
cd regtest/wind/periodic_gust_factor
../../build/levelset inputs.i
```

> **Note**: Requires a 2D build (`-DLEVELSET_DIM_2D=ON`).

## References

- Finney, M.A., et al. (2011). Role of wind, fuel moisture, and terrain in controlling fire movement. *Ecosphere*, 2(3), art17.
- Clements, C.B. (2010). Thermodynamic structure of a grass fire plume. *International Journal of Wildland Fire*, 19(7), 895-902.
