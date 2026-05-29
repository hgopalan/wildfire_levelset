# Integration Scenario B: Crown Fire with Phenology

## Purpose

Demonstrates **active crown fire spread** combined with **seasonal foliar moisture phenology** in a Sierra Nevada mixed-conifer forest during the critical spring greenup period.

## Integrated Features

### 2026 Enhancement Features
- ✅ **GDD-based FMC phenology**: Growing degree day model for spring greenup
- ✅ **McArthur moisture scaling**: Temperature/RH-dependent surface fuel drying

### Crown Fire Features
- **Cruz crown fire model**: Empirical active crown fire rate of spread
- **Van Wagner crown fire initiation**: CFB (crown fraction burned) criteria
- **Rothermel1991 crown fire spread**: Surface-to-crown fire transition

### Existing Features
- **Terrain effects**: Moderate eastward slope (10% grade)
- **Radiation preheating**: View factor-based fuel preheating
- **Bulk fuel consumption**: Post-frontal burnout simulation
- **Fire behavior diagnostics**: Intensity, flame length tracking

## Scenario Description

**Location**: Sierra Nevada mixed-conifer forest, California  
**Terrain**: Moderate eastward slope (300 m rise over 3 km)  
**Vegetation**: Dense mixed-conifer (ponderosa pine, white fir, incense cedar)  
**Canopy**: CBD = 0.18 kg/m³, CBH = 3 m  
**Ignition**: Single point at western edge  
**Weather**: Moderate spring afternoon wind (8 m/s easterly)  
**Season**: Late May (spring greenup transition period)  
**Duration**: 2 hours

## Test Configuration

| Category | Parameter | Value |
|----------|-----------|-------|
| **Domain** | Size | 3 km × 3 km |
| | Grid | 96 × 96 cells |
| **Weather** | Wind | 8 m/s easterly |
| | Temperature | 22°C |
| | Humidity | 30% |
| **Crown** | CBD | 0.18 kg/m³ |
| | CBH | 3 m |
| | FMC | 95% (from phenology) |
| **Phenology** | Model | Growing Degree Days (GDD) |
| | GDD current | 250 (mid-greenup) |
| | FMC range | 85-130% |

## Physical Interpretation

This scenario captures the critical transition period when:

1. **Spring greenup increases foliar moisture**: GDD accumulation drives FMC from dormant (85%) toward peak (130%)
2. **Crown fire vulnerability changes**: Higher FMC reduces crown fire spread rate and probability
3. **Surface-to-crown transition**: Van Wagner criteria determine when surface fire transitions to crown fire
4. **Seasonal fire behavior**: Demonstrates how phenology affects fire behavior throughout the growing season

## Expected Behavior

- Active crown fire initiation where surface fire intensity exceeds threshold
- Crown fire spread rate modulated by foliar moisture (slower than late-summer conditions)
- Terrain enhancement of fire spread on upslope (eastward) direction
- Transition between surface and crown fire regimes
- Bulk fuel consumption in fire wake

## Comparison Tests

Compare with:
- **Late summer scenario**: Same domain/weather but with dormant FMC (85%) → faster crown fire
- **Peak greenup scenario**: GDD = 800, FMC = 130% → much slower/no crown fire
- **No phenology**: Fixed FMC → loses seasonal variation

## Files

- `inputs.i`: Complete scenario configuration
- `slope_terrain.csv`: Eastward slope elevation data (10% grade)
- `README.md`: This file

## Run Command

```bash
cd regtest/integration/crown_fire_phenology
../../build/levelset inputs.i
```

> **Note**: Requires a 2D build (`-DLEVELSET_DIM_2D=ON`).

## Scientific References

- Cruz, M.G., Alexander, M.E., & Wakimoto, R.H. (2005). "Development and testing of models for predicting crown fire rate of spread in conifer forest stands." *Canadian Journal of Forest Research*, 35(7), 1626-1639.
- Van Wagner, C.E. (1977). "Conditions for the start and spread of crown fire." *Canadian Journal of Forest Research*, 7(1), 23-34.
- Rothermel, R.C. (1991). "Predicting behavior and size of crown fires in the northern Rocky Mountains." *USDA Forest Service Research Paper INT-438*.
