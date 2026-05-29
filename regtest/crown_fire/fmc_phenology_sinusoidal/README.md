# FMC Sinusoidal Phenology Test

## Purpose

Tests the **sinusoidal FMC phenology** model that varies foliar moisture content seasonally based on a simple sinusoidal curve.

## Test Configuration

| Parameter | Value | Description |
|-----------|-------|-------------|
| Domain | 1000 m × 1000 m | Flat terrain |
| Grid | 64 × 64 cells | ~15.6 m resolution |
| Wind | 6 m/s eastward | Moderate wind |
| Fuel model | FM4 chaparral | Surface fuel |
| **FMC Phenology** | | |
| Model | Sinusoidal | Simple seasonal curve |
| Mean FMC | 100% | Annual average |
| Amplitude | ±40% | Seasonal range: 60–140% |
| Peak DOY | 150 | Late May (spring greenup) |
| Start DOY | 200 | Mid-July (post-peak) |
| **Expected FMC** | ~70% | Mid-summer value |

## Expected Behavior

1. **Fire ignition**: Surface fire starts in dry chaparral
2. **FMC variation**: FMC varies seasonally per sinusoidal curve
3. **Crown fire transition**: Depends on FMC value and surface intensity
4. **Seasonal effect**: Lower summer FMC → easier crown initiation

## Files

- `inputs.i`: Test configuration
- `README.md`: This file

## Run Command

```bash
cd regtest/crown_fire/fmc_phenology_sinusoidal
../../build/levelset inputs.i
```

> **Note**: Requires a 2D build (`-DLEVELSET_DIM_2D=ON`).

## References

- Van Wagner, C.E. (1977). Conditions for the start and spread of crown fire. *Canadian Journal of Forest Research*, 7(1), 23-34.
- Jolly, W.M., et al. (2005). A generalized, bioclimatic index to predict foliar phenology. *Global Change Biology*, 11(4), 619-632.
