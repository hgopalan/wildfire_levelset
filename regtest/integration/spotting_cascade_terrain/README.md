# Integration Scenario C: Spotting Cascade with Terrain

## Purpose

Demonstrates **multi-generation spotting cascade** across **complex mountain terrain** during extreme fire weather conditions, showcasing how spot fires can create multiple generations of fire ignitions.

## Integrated Features

### 2026 Enhancement Features
- ✅ **Ember accumulation**: Density tracking with probabilistic ignition
- ✅ **Periodic wind gusts**: Extreme fire weather gusting (45% amplitude)
- ✅ **Slope-dependent flame tilt**: Enhanced upslope spread
- ✅ **McArthur moisture scaling**: Temperature/RH-dependent drying

### Spotting Features
- **Albini physics-based spotting**: Thermal plume lofting with 2-D trajectory integration
- **Multi-generation cascade**: Spot fires can generate their own spot fires
- **Spotting diagnostics**: Track up to 5 generations of spot fire lineage

### Existing Features
- **Complex terrain**: Ridge-valley system with realistic topography
- **FARSITE elliptical expansion**: Anderson dynamic L/W ratio
- **Radiation preheating**: View factor-based fuel preheating
- **Bulk fuel consumption**: Post-frontal burnout

## Scenario Description

**Location**: Mountain terrain with ridge-valley system  
**Terrain**: Two parallel N-S ridges (250m & 200m high) with central valley  
**Vegetation**: Dense chaparral (FM4)  
**Ignition**: Single point on windward (SW) ridge  
**Weather**: Extreme conditions - 12 m/s wind, 35°C, 15% RH  
**Duration**: 4 hours

## Test Configuration

| Category | Parameter | Value |
|----------|-----------|-------|
| **Domain** | Size | 4 km × 4 km |
| | Grid | 128 × 128 cells |
| **Weather** | Wind | 12 m/s (easterly) |
| | Temperature | 35°C |
| | Humidity | 15% |
| | Gust amplitude | 45% |
| **Spotting** | Base probability | 0.025 |
| | Intensity threshold | 80 kW/m |
| | Max loft height | 1200 m |
| | Ember lifetime | ~16 min |
| **Terrain** | Ridge 1 height | 250 m |
| | Ridge 2 height | 200 m |
| | Valley depth | 80 m |

## Physical Interpretation

This scenario captures extreme fire behavior where:

1. **Primary fire generates embers**: High-intensity fire on windward ridge produces firebrands
2. **Embers cross valley**: Strong winds carry embers across topographic barriers
3. **Spot fires ignite downwind**: Accumulated embers ignite new fires on leeward ridge
4. **Cascade effect**: Spot fires become intense enough to generate their own embers
5. **Multi-generation spread**: 3-5 generations of spot fires leapfrog across terrain

## Expected Behavior

- Rapid fire spread on upslope terrain faces
- Multiple spot fire ignitions downwind of main fire
- Spotting across valley to leeward ridge
- Secondary spotting from intense spot fires
- Non-contiguous fire perimeter with islands
- Exponential area growth when spotting dominates

## Spotting Cascade Mechanics

**Generation 0** (primary fire):
- Originates from initial ignition
- Generates embers when intensity > 80 kW/m

**Generation 1** (first-order spots):
- Ignited by embers from Generation 0
- Can generate their own embers once established

**Generation 2-5** (higher-order spots):
- Created by embers from previous generation spots
- Tracked in spotting diagnostics file

## Comparison Tests

Compare with:
- **No spotting**: Same scenario without spotting → much slower spread, no valley crossing
- **Flat terrain**: Remove ridges → different spotting patterns
- **Calm conditions**: Reduce wind to 4 m/s → local spotting only, no cascade
- **Single generation**: Disable cascade → limits spread distance

## Files

- `inputs.i`: Complete scenario configuration
- `ridge_valley_terrain.csv`: Complex terrain elevation data
- `README.md`: This file

## Run Command

```bash
cd regtest/integration/spotting_cascade_terrain
../../build/levelset inputs.i
```

> **Note**: Requires a 2D build (`-DLEVELSET_DIM_2D=ON`).

## Output Files

- `fire_stats.csv`: Fire behavior time series
- `spotting_events.csv`: Log of all spotting events with generation tracking
- `plt*`: Plotfiles showing fire progression and ember density

## Scientific References

- Albini, F.A. (1983). "Transport of firebrands by line thermals." *Combustion Science and Technology*, 32(5-6), 277-288.
- Koo, E., et al. (2010). "Modelling firebrand transport in wildfires using HIGRAD/FIRETEC." *International Journal of Wildland Fire*, 19(2), 200-216.
- Storey, M.A., et al. (2021). "Examining the influence of fire weather and land use on the 2017 Chile fires." *Frontiers in Climate*, 3, 641056.
