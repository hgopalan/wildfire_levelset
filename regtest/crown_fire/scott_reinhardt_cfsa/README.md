# Scott & Reinhardt Crown Fire Surface Area (CFSA) Test

## Purpose

Tests the Scott & Reinhardt (2001) Crown Fire Surface Area (CFSA) model implementation.

## Physics

The CFSA represents the effective burning surface area per unit ground area during active crown fire. It depends on:
- Canopy bulk density (CBD)
- Canopy depth (CH - CBH)
- Surface area coefficient (k_sa)

### Model Equation

```
CFSA = min(CBD × (CH - CBH) × k_sa, max_cfsa)
```

where:
- CBD = 0.20 kg/m³ (canopy bulk density)
- CH = 15.0 m (canopy height)
- CBH = 3.0 m (canopy base height)
- k_sa = 6.0 m²/kg (surface area coefficient)
- max_cfsa = 3.0 (maximum CFSA)

### Expected Results

```
Canopy depth = CH - CBH = 15.0 - 3.0 = 12.0 m
CFSA (uncapped) = 0.20 × 12.0 × 6.0 = 14.4
CFSA (capped) = min(14.4, 3.0) = 3.0
```

The test should produce active crown fire with CFSA ≈ 3.0 in the burned region.

## Reference

Scott, J.H. & Reinhardt, E.D. (2001). "Assessing Crown Fire Potential by Linking Models of Surface and Crown Fire Behavior." USDA Forest Service Research Paper RMRS-RP-29.

## Validation

- Verify CFSA values are computed and stored in output
- Check CFSA ≈ 3.0 in active crown fire regions
- Verify CFSA = 0 in unburned regions
- Confirm crown fire heat release scales with CFSA
