# FARSITE Wind Stream Simulation Test

Tests the FARSITE empirical wind stream simulation model (Option 8) which uses terrain curvature to identify ridges and valleys, then applies empirical speed modifications and direction deflections.

## Key Features Tested

1. **Terrain curvature computation**: Second derivatives of elevation field
2. **Ridge acceleration**: Wind speeds up over convex terrain (positive curvature)
3. **Lee-side sheltering**: Wind slows in protected areas downwind of ridges
4. **Valley channeling**: Wind is funneled through concave terrain (negative curvature)
5. **Wind direction deflection**: Wind follows terrain contours

## Test Configuration

- **Terrain**: Gaussian hill (reused from windninja_ridge_canyon test)
- **Wind**: Spatially-varying field from file
- **Fuel**: FM4 (chaparral)
- **Grid**: 64×64 cells, 1000m×1000m domain
- **Propagation**: Level-set method

## Expected Behavior

The fire should:
- Accelerate on the windward (upslope) face of the hill
- Slow in the sheltered lee-side region
- Show direction deflection where wind crosses terrain contours

## Parameters

- `k_ridge_farsite = 1.5`: Ridge speed-up coefficient
- `k_shelter = 0.6`: Lee-side sheltering coefficient
- `k_valley = 0.8`: Valley channeling coefficient
- `k_deflection = 0.3`: Wind direction deflection coefficient
- `min_curvature = 0.0001`: Minimum curvature threshold

## References

- Finney, M.A. (1998). FARSITE: Fire Area Simulator—Model Development and Evaluation. USDA Forest Service Research Paper RMRS-RP-4.
