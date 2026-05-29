# Multi-Layer Canopy Wind Profile Test

## Purpose

Tests the multi-layer canopy wind profile model based on exponential attenuation (Massman 1997, Inoue 1963).

## Physics

Wind speed decreases exponentially from canopy top to surface. The model divides the vertical domain into layers and computes wind speed at each layer midpoint.

### Model Equations

**Within canopy (z < h):**
```
u(z) = u_h × exp(-α × (1 - z/h))
```

**Above canopy (z > h):**
```
u(z) = u_h × [ln((z-d)/z₀) / ln((h-d)/z₀)]
```

where:
- u_h = wind speed at canopy top [m/s]
- α = attenuation coefficient = 2.5 (from LAI ≈ 4)
- h = canopy height = 20.0 m
- d = displacement height = 0.7 × h = 14.0 m
- z₀ = roughness length = 0.1 × h = 2.0 m

### Expected Results

For input wind at 10m (below canopy top):
```
Layer 1 (z ≈ 2m):   u ≈ 2.5 m/s (strong attenuation near surface)
Layer 2 (z ≈ 6m):   u ≈ 4.0 m/s
Layer 3 (z ≈ 10m):  u = 8.0 m/s (reference height)
Layer 4 (z ≈ 14m):  u ≈ 9.5 m/s
Layer 5 (z ≈ 18m):  u ≈ 10.5 m/s (near canopy top)
```

## References

- Massman, W.J. (1997). "An analytical one-dimensional model of momentum transfer by vegetation of arbitrary structure." Boundary-Layer Meteorology, 83(3), 407-421.
- Inoue, E. (1963). "On the turbulent structure of airflow within crop canopies." Journal of the Meteorological Society of Japan, 41, 317-326.

## Validation

- Verify exponential wind decrease from canopy top to surface
- Check surface wind is significantly reduced (30-40% of canopy top wind)
- Confirm logarithmic profile above canopy
- Validate layer-specific wind speeds affect fire spread rate
