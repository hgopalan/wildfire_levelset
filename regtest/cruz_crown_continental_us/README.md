# Cruz, Alexander & Wakimoto (2005) Crown Fire Spread Test — Continental US

## Purpose

Tests fire spread using the Cruz, Alexander & Wakimoto (2005) simple algebraic
crown fire spread model. This model was calibrated on active crown fire
observations from **western North American conifer forests**, including data
from the Sierra Nevada, Rocky Mountains, and Pacific Northwest of the
**continental United States**.

## Model Overview

The head-fire rate of spread is a single closed-form algebraic expression:

```
R [m/min] = 11.02 × U₁₀^0.90 × CBD^0.19 × exp(−0.17 × MC₁₀)
```

where:
- `U₁₀` — 10-m open wind speed [km/h]
- `CBD` — canopy bulk density [kg/m³]
- `MC₁₀` — 10-h timelag fine fuel moisture content [%]

The model is fully algebraic (one closed-form equation, no iteration) and
simple to evaluate per grid cell.

**Applicable conditions (Cruz et al. 2005):**
- Wind: 5–60 km/h
- CBD: 0.05–0.40 kg/m³ (western US conifer canopies)
- MC₁₀: 4–25%

## Scenario: Sierra Nevada Mixed-Conifer Crown Fire

This test represents a summer crown fire in a Sierra Nevada white fir / lodgepole
pine stand under typical severe fire weather conditions (hot, dry, windy day).

### Configuration

| Parameter | Value | Justification |
|-----------|-------|---------------|
| Domain | 2000 m × 2000 m | Representative forest patch |
| Grid | 64 × 64 × 64 cells | 31.25 m resolution |
| Initial condition | Sphere at (1000, 1000) m, radius 30 m | Point ignition |
| Wind (`u_x`) | 6.0 m/s (= 21.6 km/h) | Moderate Diablo/Santa Ana conditions |
| `cruz_crown.CBD` | 0.15 kg/m³ | Dense Sierra Nevada mixed-conifer stand |
| `cruz_crown.MC10` | 8.0% | Dry summer conditions (Jul–Aug) |
| Simulation time | 1800 s (30 min) | Short-term crown fire spread |
| Output interval | Every 20 steps | ~5–10 min resolution |
| Propagation | `levelset` | Level-set method |

### Expected Behavior

With U₁₀ = 21.6 km/h, CBD = 0.15 kg/m³, MC₁₀ = 8%:

```
R = 11.02 × 21.6^0.90 × 0.15^0.19 × exp(−0.17 × 8)
  ≈ 11.02 × 15.36 × 0.698 × 0.259
  ≈ 30.6 m/min
  ≈ 0.51 m/s
```

The fire should expand asymmetrically, spreading rapidly downwind (east) and
more slowly upwind (west), forming an elongated elliptical perimeter consistent
with observed western US crown fire behavior.

After 30 minutes, the downwind extent should be on the order of 900–1000 m
from the ignition point.

## Variants

Modify these parameters to explore model sensitivity:

| Parameter | Effect |
|-----------|--------|
| `cruz_crown.CBD` | Higher CBD → faster crown fire spread (denser canopy provides more fuel continuity) |
| `cruz_crown.MC10` | Higher moisture → slower spread (exponential dampening) |
| `u_x` | Wind speed is the dominant driver (power-law exponent 0.90) |

**Rocky Mountain ponderosa pine example** (sparser canopy, drier conditions):
```
cruz_crown.CBD  = 0.08   # kg/m³ — sparse ponderosa pine
cruz_crown.MC10 = 5.0    # %     — extreme drought conditions
u_x = 8.0               # m/s   — strong downslope wind
```

**Pacific Northwest Douglas-fir example** (denser canopy, moderate moisture):
```
cruz_crown.CBD  = 0.25   # kg/m³ — dense Douglas-fir stand
cruz_crown.MC10 = 12.0   # %     — moderate moisture
u_x = 5.0               # m/s   — sustained wind
```

## Run Command

```bash
./build/levelset regtest/cruz_crown_continental_us/inputs.i
```

## References

- Cruz, M.G., Alexander, M.E., and Wakimoto, R.H. (2005). "Development and
  testing of models for predicting crown fire rate of spread in conifer forest
  stands." *Canadian Journal of Forest Research*, 35(7), 1626–1639.
  https://doi.org/10.1139/x05-068

- Scott, J.H. and Reinhardt, E.D. (2001). "Assessing crown fire potential by
  linking models of surface and crown fire behavior." USDA Forest Service
  Research Paper RMRS-RP-29. (Background on CBD and crown fire in western US.)
