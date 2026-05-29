# FMC Growing Degree Day (GDD) Phenology Test

## Purpose

Tests the **GDD-based FMC phenology** model that varies foliar moisture content based on accumulated growing degree days, representing temperature-driven spring greenup.

## Physical Model

### Growing Degree Day Accumulation

Daily GDD accumulation:
```
GDD_daily = max(0, T_mean - T_base)
GDD_total = Σ GDD_daily
```

### FMC Progression

```
if GDD < GDD_start:
    FMC = FMC_dormant
elif GDD < GDD_peak:
    FMC = FMC_dormant + (FMC_peak - FMC_dormant) × (GDD - GDD_start) / (GDD_peak - GDD_start)
else:
    FMC = FMC_peak
```

## Test Configuration

| Parameter | Value | Description |
|-----------|-------|-------------|
| Domain | 1000 m × 1000 m | Flat terrain |
| Grid | 64 × 64 cells | ~15.6 m resolution |
| Wind | 6 m/s eastward | Moderate wind |
| **GDD Phenology** | | |
| Model | GDD | Temperature-driven |
| Dormant FMC | 80% | Winter/cured state |
| Peak FMC | 140% | Full greenup |
| GDD start | 50 | Greenup begins |
| GDD peak | 400 | Full greenup threshold |
| Base temperature | 5°C | Minimum for growth |
| Mean temperature | 15°C | Daily average |
| Start DOY | 100 | Early April |
| **Expected GDD** | ~1000 | 100 days × 10°/day |
| **Expected FMC** | 140% | Full greenup |

## Expected Behavior

1. **Fire ignition**: Surface fire in dry fuel
2. **High FMC**: Spring greenup conditions (140%)
3. **Crown fire resistance**: High FMC → high I_o threshold
4. **Limited crown spread**: Crown fire unlikely or slow

## Files

- `inputs.i`: Test configuration
- `README.md`: This file

## Run Command

```bash
cd regtest/crown_fire/fmc_phenology_gdd
../../build/levelset inputs.i
```

> **Note**: Requires a 2D build (`-DLEVELSET_DIM_2D=ON`).

## References

- Jolly, W.M., et al. (2005). A generalized, bioclimatic index to predict foliar phenology. *Global Change Biology*, 11(4), 619-632.
- Richardson, A.D., et al. (2006). Phenology of a northern hardwood forest canopy. *Global Change Biology*, 12(7), 1174-1188.
