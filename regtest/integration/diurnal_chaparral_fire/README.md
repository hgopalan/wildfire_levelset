# Integration Scenario A: Complete Diurnal Wildfire

## Purpose

Demonstrates a **realistic operational forecast scenario** combining multiple 2026 enhancement features in a single 24-hour Southern California chaparral wildfire simulation.

## Integrated Features

### 2026 Enhancement Features
- ✅ **McArthur moisture scaling**: Temperature/RH-dependent drying rates
- ✅ **FMC sinusoidal phenology**: Seasonal foliar moisture (mid-summer)
- ✅ **Ember accumulation**: Density tracking with probabilistic ignition
- ✅ **Periodic wind gusts**: Afternoon thermal turbulence (15-min cycles)
- ✅ **Slope-dependent flame tilt**: Enhanced upslope radiation preheating

### Existing Features
- **FARSITE elliptical expansion**: Anderson dynamic L/W ratio
- **Terrain effects**: Gaussian hill with slope/aspect
- **Albini spotting**: Physics-based firebrand transport
- **Radiation preheating**: View factor-based fuel preheating
- **Burn period gate**: 10:00–20:00 local time burning window

## Scenario Description

**Location**: Southern California chaparral landscape  
**Terrain**: Gaussian hill (150 m height, 400 m characteristic width)  
**Vegetation**: FM4 mature chaparral (6-8 ft tall)  
**Ignition**: Single point at western edge at 10:00 local time  
**Weather**: Hot, dry afternoon with periodic gusts  
**Duration**: 24 hours (July 19–20, 2024)

## Test Configuration

| Category | Parameter | Value |
|----------|-----------|-------|
| **Domain** | Size | 2 km × 2 km |
| | Grid | 128 × 128 cells |
| **Weather** | Base wind | 4 m/s westerly |
| | Temperature | 28°C |
| | Humidity | 25% |
| | Gust amplitude | 35% |

## Files

- `inputs.i`: Complete scenario configuration
- `terrain_hill.csv`: Gaussian hill elevation/slope data
- `README.md`: This file

## Run Command

```bash
cd regtest/integration/diurnal_chaparral_fire
../../build/levelset inputs.i
```

> **Note**: Requires a 2D build (`-DLEVELSET_DIM_2D=ON`).
