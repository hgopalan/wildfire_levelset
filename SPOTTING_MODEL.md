# Firebrand Spotting Model Documentation

## Overview

The firebrand spotting model simulates the generation of new fire ignition points ahead of the main fire front due to wind-borne firebrands (burning embers). This model works within the FARSITE framework and uses a probability-based approach to generate spot fires stochastically.

## Physical Background

### Firebrand Transport

Firebrands are burning pieces of vegetation (bark, branches, leaves) that are lofted by convection and transported downwind by ambient winds. When they land in receptive fuel beds, they can create new ignition points (spot fires) ahead of the main fire.

Key factors affecting spotting:
- **Wind speed**: Higher winds increase lofting and transport distance
- **Fire intensity**: Hotter fires produce more and larger firebrands
- **Fuel moisture**: Drier fuels produce more viable firebrands
- **Topography**: Affects wind patterns and firebrand trajectories

### Importance in Wildfire Modeling

Spotting is a critical fire behavior phenomenon that:
- Accelerates fire spread beyond the continuous fire front
- Creates non-contiguous fire perimeters
- Challenges suppression efforts
- Is responsible for many wildfire-related structure losses

## Model Formulation

### Spotting Probability

The probability of spotting from a fire front cell is computed as:

```
P_spot = P_base × f_wind × f_intensity × f_moisture
```

Where:
- `P_base`: Base spotting probability (user-specified, typically 0.01-0.05)
- `f_wind`: Wind factor (0-1)
- `f_intensity`: Fire intensity factor (0-1)
- `f_moisture`: Fuel moisture factor (0-1)

#### Wind Factor

```
f_wind = 1 - exp(-k_wind × U²/100)
```

- `U`: Wind speed in ft/min (from Rothermel model)
- `k_wind`: Wind coefficient (user-specified, typically 0.1-0.5)
- Higher wind speeds increase probability non-linearly

#### Intensity Factor

```
f_intensity = min(1.0, I / I_critical)
```

- `I`: Fire intensity approximated as `R × h_heat` (BTU/ft²/min)
- `R`: Rate of spread from Rothermel model (ft/min)
- `h_heat`: Heat content of fuel (BTU/lb)
- `I_critical`: Critical intensity threshold (user-specified, typically 500-2000 BTU/ft²/min)

#### Moisture Factor

```
f_moisture = max(0, (M_x - M_f) / M_x)
```

- `M_x`: Moisture of extinction (from Rothermel params)
- `M_f`: Fuel moisture content (from Rothermel params)
- Drier fuels (lower M_f) produce more viable firebrands

### Spotting Distance

Two distance distribution models are available:

#### Lognormal Distribution (Recommended)

```
ln(d) ~ N(μ, σ²)
```

Where:
- `μ = ln(d_mean) - 0.5σ²` ensures mean distance equals `d_mean`
- `σ = d_sigma` controls distribution width
- Lognormal is physically realistic for firebrand transport

**Typical values**:
- `d_mean = 0.1-0.3` (10-30% of domain size)
- `d_sigma = 0.4-0.8` (wider distribution for more variability)

#### Exponential Distribution

```
d ~ Exp(λ)
```

Where:
- `λ = d_lambda` is the decay rate
- Mean distance = `1/λ`
- Simpler model, less realistic

### Spotting Direction

Spots are placed primarily downwind with lateral dispersion:

```
direction = wind_direction + lateral_offset
```

Where:
- `wind_direction`: Unit vector in wind direction
- `lateral_offset ~ N(0, σ_lateral²)` 
- `σ_lateral = tan(lateral_spread_angle) × distance`

This creates a cone of possible landing locations with apex at the fire front.

### Spot Fire Ignition

Each generated spot creates a small ignition zone:
- Shape: Circular (2D) or spherical (3D)
- Radius: `spot_radius` (user-specified, typically 2-3 grid cells)
- Implementation: Set `phi` to signed distance from spot center within radius

## Parameters

### Required Parameters

All spotting parameters use the `spotting.` prefix:

| Parameter | Description | Typical Range | Default |
|-----------|-------------|---------------|---------|
| `enable` | Enable spotting (0=off, 1=on) | 0 or 1 | 0 |
| `P_base` | Base spotting probability | 0.01-0.10 | 0.02 |
| `k_wind` | Wind speed coefficient | 0.1-0.5 | 0.3 |
| `I_critical` | Critical fire intensity (BTU/ft²/min) | 500-2000 | 1000.0 |
| `d_mean` | Mean spotting distance | 0.05-0.3 | 0.1 |
| `d_sigma` | Distance std dev (lognormal) | 0.3-1.0 | 0.5 |
| `d_lambda` | Distance decay rate (exponential) | 5.0-20.0 | 10.0 |
| `distance_model` | "lognormal" or "exponential" | - | "lognormal" |
| `lateral_spread_angle` | Lateral dispersion angle (degrees) | 10-30 | 15.0 |
| `spot_radius` | Radius of spot fire ignition zone | 0.01-0.05 | 0.02 |
| `random_seed` | RNG seed (0=time-based) | any int | 0 |
| `check_interval` | Check for spotting every N steps | 1-10 | 5 |

### Parameter Selection Guidelines

#### Spotting Frequency (`P_base`, `check_interval`)
- **Low frequency**: `P_base = 0.01-0.02`, `check_interval = 10`
- **Moderate frequency**: `P_base = 0.03-0.05`, `check_interval = 5`
- **High frequency**: `P_base = 0.08-0.10`, `check_interval = 3`

Start conservative (low frequency) and increase as needed.

#### Distance Scaling
Scale `d_mean` relative to domain size and cell resolution:
- For 1.0 × 1.0 domain with 64 cells: `d_mean = 0.1-0.2`
- For larger domains: `d_mean` can be absolute (e.g., 100m if using real units)

#### Spot Size
- Minimum: `spot_radius ≥ 2 × cell_size`
- Typical: `spot_radius = 2-3 × cell_size`
- Too small: spots may not be resolved
- Too large: unrealistic ignition zones

#### Wind Sensitivity
- **Low sensitivity**: `k_wind = 0.1` (spotting occurs even in light winds)
- **Moderate sensitivity**: `k_wind = 0.3` (default)
- **High sensitivity**: `k_wind = 0.5-1.0` (spotting only in strong winds)

## Integration with FARSITE

The spotting model requires FARSITE to be enabled:

```bash
# FARSITE must be enabled
farsite.enable = 1
skip_levelset = 1  # Recommended for FARSITE mode

# Then enable spotting
spotting.enable = 1
```

**Execution Order** (each timestep):
1. FARSITE computes elliptical fire spread
2. Fire front advances to new positions
3. Spotting probability computed at new fire front
4. Firebrand spots generated stochastically
5. New ignition points added to phi field
6. Process repeats

## Output Fields

The spotting model adds four output fields to plotfiles:

| Field | Description | Range | Units |
|-------|-------------|-------|-------|
| `spot_prob` | Spotting probability | 0.0-1.0 | dimensionless |
| `spot_count` | Number of firebrands generated | ≥ 0 | count |
| `spot_dist` | Spotting distance field | ≥ 0 | sim units |
| `spot_active` | Active spot fire flag | 0 or 1 | binary |

### Visualization Tips

In ParaView or VisIt:
1. **Fire perimeter**: Threshold `phi < 0` to see burned regions
2. **Spotting zones**: Color by `spot_prob` to see high-probability areas
3. **Spot locations**: Track where `phi` becomes negative in unburned regions
4. **Temporal evolution**: Animate to observe spot fire generation and growth

## Example Configurations

### Conservative Spotting (Few Spots)
```bash
spotting.enable = 1
spotting.P_base = 0.01
spotting.k_wind = 0.5
spotting.check_interval = 10
spotting.d_mean = 0.1
spotting.d_sigma = 0.3
```

### Moderate Spotting (Realistic)
```bash
spotting.enable = 1
spotting.P_base = 0.03
spotting.k_wind = 0.3
spotting.check_interval = 5
spotting.d_mean = 0.15
spotting.d_sigma = 0.5
```

### Aggressive Spotting (Many Spots)
```bash
spotting.enable = 1
spotting.P_base = 0.08
spotting.k_wind = 0.2
spotting.check_interval = 3
spotting.d_mean = 0.2
spotting.d_sigma = 0.8
```

### Exponential Distance Model
```bash
spotting.enable = 1
spotting.distance_model = exponential
spotting.d_lambda = 8.0  # mean = 1/λ = 0.125
```

## Validation and Testing

### Expected Behaviors

A properly configured spotting model should exhibit:
1. **Downwind bias**: Most spots in wind direction
2. **Distance distribution**: Variety of distances, centered on `d_mean`
3. **Lateral spread**: Some spots perpendicular to wind
4. **Multiple generations**: Spot fires can themselves generate new spots
5. **Merging**: Spot fires grow and eventually merge with main fire

### Validation Checks

Run diagnostic cases to verify:
- **No wind**: Minimal or no spotting (`u_x = u_y = u_z = 0`)
- **Low intensity**: Reduced spotting with low fire intensity
- **High moisture**: Reduced spotting with high fuel moisture (`M_f ≈ M_x`)
- **Reproducibility**: Same results with fixed `random_seed`

### Common Issues

**No spots generated**:
- Check that `farsite.enable = 1` and `skip_levelset = 1`
- Increase `P_base` or decrease `check_interval`
- Verify wind speed > 0
- Check that fire front exists (phi ≈ 0)

**Too many spots**:
- Decrease `P_base`
- Increase `check_interval`
- Increase `k_wind` for higher wind sensitivity

**Spots too close/far**:
- Adjust `d_mean` for mean distance
- Adjust `d_sigma` for distance variability

## Limitations and Future Work

### Current Limitations

1. **CPU-only random numbers**: RNG uses `std::random`, limiting GPU portability
2. **Simplified physics**: No firebrand lofting height or burnout modeling
3. **No fuel receptivity**: All locations equally receptive to spot ignition
4. **No terrain interaction**: Topography doesn't affect landing probability

### Potential Enhancements

1. **GPU-compatible RNG**: Use AMReX random for device kernels
2. **Firebrand trajectory model**: Physics-based ballistic transport
3. **Fuel moisture/type effects**: Vary receptivity by fuel properties
4. **Terrain effects**: Upslope/downslope influence on landing zones
5. **Spot fire growth model**: Different growth rates for spot fires
6. **Statistical calibration**: Fit parameters to observed fire data

## References

1. Koo, E., et al. (2010). "Modelling firebrand transport in wildfires using HIGRAD/FIRETEC." *International Journal of Wildland Fire*, 19(6), 722-729.

2. Sardoy, N., et al. (2007). "Modeling transport and combustion of firebrands from burning trees." *Combustion and Flame*, 150(3), 151-169.

3. Tarifa, C. S., et al. (1965). "On the flight paths and lifetimes of burning particles of wood." *Symposium (International) on Combustion*, 10(1), 1021-1037.

4. Albini, F. A. (1979). "Spot fire distance from burning trees—a predictive model." USDA Forest Service, Intermountain Forest and Range Experiment Station, Research Paper INT-56.

5. Finney, M. A. (1998). "FARSITE: Fire Area Simulator—model development and evaluation." Research Paper RMRS-RP-4, USDA Forest Service.

## Contact

For questions or issues related to the spotting model, please open an issue on the GitHub repository.
