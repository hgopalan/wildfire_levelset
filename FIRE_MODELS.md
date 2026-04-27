# Fire Spread Models Documentation

This document provides detailed information about the fire spread models implemented in the wildfire level-set solver.

## Overview

The solver now supports four fire spread modeling approaches:

1. **Basic Model**: Simple wind-driven spread (original implementation)
2. **Rothermel Model**: Includes terrain slope corrections
3. **FARSITE Model**: Full model with Anderson L/W ratio and combined wind/terrain effects
4. **Firebrand Spotting Model**: Probability-based spotting integrated with FARSITE
5. **Bulk Fuel Consumption Model**: Computes fraction of fuel consumed (FARSITE only)

## Model Selection

Use runtime parameters to select the desired model:

```bash
# Basic model (default)
./build/levelset u_x=0.5 u_y=0.0

# Rothermel model with terrain
./build/levelset u_x=0.5 u_y=0.0 use_terrain_effects=true terrain_slope=15.0 terrain_aspect=90.0

# FARSITE model with L/W ratio
./build/levelset u_x=0.5 u_y=0.0 use_terrain_effects=true use_farsite_model=true terrain_slope=15.0 terrain_aspect=90.0

# FARSITE with firebrand spotting
./build/levelset skip_levelset=1 farsite.enable=1 spotting.enable=1 u_x=0.4

# FARSITE with bulk fuel consumption
./build/levelset skip_levelset=1 farsite.enable=1 farsite.use_bulk_fuel_consumption=1 u_x=0.4
```

## Anderson L/W Ratio (FARSITE)

### Background
Anderson (1983) developed an empirical relationship for the length-to-width (L/W) ratio of wildfire spread patterns based on wind speed. This relationship is fundamental to the FARSITE fire behavior model.

### Formula
```
L/W = 0.936 × exp(0.2566 × U) + 0.461 × exp(-0.1548 × U) - 0.397
```

where:
- U = wind speed at midflame height (mph)
- L/W = 1.0 for calm conditions (circular fire)
- L/W increases with wind speed (elliptical fire)

### Implementation Details
- Wind speed is converted from m/s to mph (factor: 2.23694)
- L/W ratio is bounded between 1.0 and 8.0 for physical realism
- The ratio affects the directional rate of spread

### Physical Meaning
- L/W = 1.0: Circular fire (no wind)
- L/W = 2.0: Fire spreads twice as fast in wind direction vs. perpendicular
- L/W = 4.0: Highly elongated fire shape (strong winds)

## Rothermel Slope Factor

### Background
Rothermel (1972) developed a mathematical model for wildland fire spread that includes terrain slope effects. Fires spread faster uphill and slower downhill.

### Formula
```
φ_s = 5.275 × tan²(θ)
```

where:
- θ = terrain slope angle (radians)
- φ_s = slope factor (dimensionless multiplier)

### Implementation Details
- Slope input is in degrees, converted to radians internally
- The coefficient 5.275 is for typical fuel bed properties
- Effect is directional: maximum uphill, zero perpendicular, negative downhill

### Physical Meaning
- 10° slope: ~5% increase in uphill spread rate
- 20° slope: ~20% increase in uphill spread rate
- 30° slope: ~50% increase in uphill spread rate

## FARSITE Combined Model

### Background
FARSITE (Fire Area Simulator) combines wind and slope effects vectorially, accounting for their relative directions.

### Algorithm
1. Calculate wind factor: φ_w = exp(0.05 × U)
2. Calculate slope factor: φ_s (Rothermel equation)
3. Determine alignment between wind and slope directions
4. Combine effects: ROS = base × (1 + φ_w + φ_s × alignment)

### Vectorial Combination
- **Aligned** (wind uphill): Effects add constructively
- **Perpendicular**: Only wind effect applies
- **Opposed** (wind downhill): Effects partially cancel

### Implementation Details
- Wind direction computed from velocity components: atan2(v_y, v_x)
- Slope aspect represents uphill direction
- Alignment factor: cos(wind_direction - slope_aspect)
- Only positive alignment contributes to slope effect

## Bulk Fuel Consumption Fraction Model (FARSITE Only)

### Background
The bulk fuel consumption fraction model computes what fraction of available fuel is consumed as the fire passes through. Not all fuel is consumed during fire passage - the fraction depends on fire intensity, residence time, and fuel properties. This model is currently only available in the FARSITE pathway.

### Physical Basis
Fuel consumption fraction is influenced by:
- **Fire Intensity**: Higher intensity fires consume more fuel
- **Residence Time**: Longer residence allows more complete combustion
- **Fuel Moisture**: Lower moisture increases consumption (via intensity)
- **Fuel Type**: Implicitly through Rothermel parameters

### Formula
The model uses a sigmoid-like transition function:

```
I_norm = (I_R / I_ref) × √(τ / τ_ref)
f_c = f_c_min + (f_c_max - f_c_min) × [0.5 × (1 + tanh(I_norm - 1))]
```

where:
- I_R = Reaction intensity from Rothermel model [BTU/ft²/min]
- I_ref = Reference intensity (1000 BTU/ft²/min)
- τ = Residence time [seconds]
- τ_ref = Reference residence time (60 seconds)
- f_c_min = Minimum consumption fraction (default: 0.5)
- f_c_max = Maximum consumption fraction (default: 0.9)
- f_c = Computed fuel consumption fraction (0.0-1.0)

### Parameter Guidelines

| Parameter | Description | Typical Range | Default |
|-----------|-------------|---------------|---------|
| `farsite.use_bulk_fuel_consumption` | Enable model (0=off, 1=on) | 0 or 1 | 0 |
| `farsite.tau_residence` | Residence time [seconds] | 30-120 | 60.0 |
| `farsite.f_consumed_min` | Minimum fraction | 0.3-0.6 | 0.5 |
| `farsite.f_consumed_max` | Maximum fraction | 0.8-1.0 | 0.9 |

### Physical Interpretation

**Low Intensity Fires** (I_R < 500 BTU/ft²/min):
- Fast-moving grass fires
- f_c ≈ f_c_min (typically 0.5)
- 50% of fuel remains unconsumed

**Moderate Intensity Fires** (500-2000 BTU/ft²/min):
- Typical brush and forest fires
- f_c transitions between min and max
- Consumption depends on residence time

**High Intensity Fires** (I_R > 2000 BTU/ft²/min):
- Crown fires, high fuel loads
- f_c ≈ f_c_max (typically 0.9)
- 90% of fuel consumed

### Implementation Details
- Computed at each fire front cell during FARSITE spread calculation
- Uses reaction intensity I_R from Rothermel model
- Output field `fuel_consumption` in plotfiles shows f_c values
- Model is applied locally (spatially varying if I_R varies)
- Does NOT affect current spread rate (for future coupling)

### Example Usage

**Enable with default parameters:**
```bash
./build/levelset skip_levelset=1 farsite.enable=1 \
  farsite.use_bulk_fuel_consumption=1 u_x=0.4
```

**Fast-moving grass fire (lower consumption):**
```bash
./build/levelset skip_levelset=1 farsite.enable=1 \
  farsite.use_bulk_fuel_consumption=1 \
  farsite.tau_residence=30 \
  farsite.f_consumed_min=0.3 \
  farsite.f_consumed_max=0.7 \
  u_x=0.6
```

**Slow-moving crown fire (higher consumption):**
```bash
./build/levelset skip_levelset=1 farsite.enable=1 \
  farsite.use_bulk_fuel_consumption=1 \
  farsite.tau_residence=120 \
  farsite.f_consumed_min=0.6 \
  farsite.f_consumed_max=0.95 \
  rothermel.fuel_model=TU5
```

### Output
The `fuel_consumption` field in plotfiles contains the computed consumption fraction (0.0-1.0) at each grid point. Values are only meaningful at fire front locations where phi ≈ 0.

### Residual Fuel
The residual fuel fraction is simply:
```
f_residual = 1 - f_c
```

This can be important for:
- Multi-pass fire simulations
- Fuel recovery/regrowth modeling
- Post-fire effects analysis

### Future Enhancements
Potential improvements for future versions:
- Couple consumption to spread rate (reduced fuel → reduced ROS)
- Fuel moisture effects on consumption (beyond intensity)
- Time-integrated consumption (track fuel state over time)
- Fuel type-specific consumption models
- Validation against experimental data

### References
1. Albini, F. A. (1976). "Estimating wildfire behavior and effects." USDA Forest Service Research Paper INT-30.
2. Byram, G. M. (1959). "Combustion of forest fuels." In Forest Fire: Control and Use.
3. Andrews, P. L. (2018). "The Rothermel surface fire spread model and associated developments: A comprehensive explanation." USDA Forest Service General Technical Report RMRS-GTR-371.

## Terrain Data

### Slope
- Units: degrees
- Range: 0° (flat) to 90° (vertical)
- Typical values: 0-30° for most terrain

### Aspect
- Units: degrees (input), radians (internal)
- Convention:
  - 0° = East
  - 90° = North
  - 180° = West
  - 270° = South
- Represents direction of steepest uphill slope

### Uniform Terrain
For simple cases, specify constant slope and aspect:
```bash
terrain_slope=15.0    # 15-degree slope
terrain_aspect=45.0   # Northeast facing
```

### Variable Terrain
For complex terrain, use `init_terrain_from_elevation()`:
- Requires elevation data as a MultiFab
- Computes slope and aspect from elevation gradients
- Needs ghost cells for boundary handling

## Rate of Spread Calculation

### Basic Model
```
ROS = wind_magnitude
```

### Rothermel Model
```
ROS = wind_magnitude × (1 + slope_factor × uphill_component)
```

### FARSITE Model
```
ROS = wind_magnitude × combined_factor × (1 + (L/W - 1) × 0.3)
```

where:
- combined_factor accounts for wind and terrain
- 0.3 is an empirical scaling factor to moderate L/W influence

## Level-Set Equation

The governing equation for the fire front propagation is:

```
∂φ/∂t = -ROS × (|∇φ| - ε × Δφ)
```

where:
- φ = level-set function (zero contour = fire front)
- ROS = rate of spread (from fire models above)
- |∇φ| = gradient magnitude (computed with WENO5-Z)
- Δφ = Laplacian (artificial viscosity term)
- ε = viscosity coefficient (default: 0.4)

## References

1. Anderson, H. E. (1983). "Predicting wind-driven wild land fire size and shape." Research Paper INT-305, USDA Forest Service.

2. Rothermel, R. C. (1972). "A mathematical model for predicting fire spread in wildland fuels." Research Paper INT-115, USDA Forest Service.

3. Finney, M. A. (1998). "FARSITE: Fire Area Simulator - Model development and evaluation." Research Paper RMRS-RP-4, USDA Forest Service.

4. Mandel, J., et al. (2014). "A wildland fire model with data assimilation." Mathematics and Computers in Simulation, 79(3), 584-606.

## Example Use Cases

### Flat Terrain, Moderate Wind
```bash
./build/levelset u_x=1.0 u_y=0.0 nsteps=500
```
Expected: Elliptical fire shape elongated in x-direction

### Uphill Fire with Wind Alignment
```bash
./build/levelset u_x=1.0 u_y=0.0 use_terrain_effects=true terrain_slope=20.0 terrain_aspect=0.0
```
Expected: Accelerated spread (wind and slope aligned)

### Uphill Fire with Cross-Wind
```bash
./build/levelset u_x=1.0 u_y=0.0 use_terrain_effects=true terrain_slope=20.0 terrain_aspect=90.0
```
Expected: Asymmetric fire shape (wind and slope perpendicular)

### FARSITE Model with Complex Conditions
```bash
./build/levelset u_x=0.8 u_y=0.3 use_terrain_effects=true use_farsite_model=true \
  terrain_slope=25.0 terrain_aspect=30.0 nsteps=1000 plot_int=50
```
Expected: Complex fire shape reflecting combined wind/terrain effects

## Firebrand Spotting Model

### Background

The firebrand spotting model simulates the generation of new fire ignition points ahead of the main fire front. Firebrands (burning embers) are lofted by fire convection and transported downwind, creating spot fires when they land in receptive fuel beds. This is a critical fire behavior phenomenon that:
- Accelerates fire spread beyond the continuous fire front
- Creates non-contiguous fire perimeters
- Is responsible for many structure losses in wildfire events

### Model Overview

The spotting model uses a probability-based approach:
1. Identifies fire front locations (phi ≈ 0)
2. Computes spotting probability based on wind, fire intensity, and fuel moisture
3. Generates new ignition points stochastically
4. Places spot fires at downwind locations with lateral dispersion

### Integration with FARSITE

**Important**: The spotting model works within the FARSITE framework and requires:
- `farsite.enable = 1` (FARSITE must be enabled)
- `skip_levelset = 1` (recommended for FARSITE mode)
- `spotting.enable = 1` (enable spotting)

The spotting model does NOT work with level-set advection. It is designed specifically for the FARSITE elliptical expansion approach.

### Probability Function

Spotting probability at each fire front cell:
```
P_spot = P_base × f_wind × f_intensity × f_moisture
```

Where:
- `P_base`: Base probability (user-specified, 0.01-0.10)
- `f_wind = 1 - exp(-k_wind × U²/100)`: Wind factor (higher wind → higher probability)
- `f_intensity = min(1, I/I_critical)`: Fire intensity factor
- `f_moisture = (M_x - M_f)/M_x`: Fuel moisture factor

### Distance Distribution

Two models available:
1. **Lognormal** (recommended): `ln(d) ~ N(μ, σ²)` where `μ = ln(d_mean) - 0.5σ²`
2. **Exponential**: `d ~ Exp(λ)` where mean = `1/λ`

Lognormal is more physically realistic for firebrand transport.

### Directional Dispersion

Spots are placed primarily downwind with lateral spread:
- Primary direction: Wind direction
- Lateral offset: Gaussian with std dev = `tan(lateral_angle) × distance`
- Creates a cone of possible landing locations

### Key Parameters

| Parameter | Description | Typical Range |
|-----------|-------------|---------------|
| `spotting.P_base` | Base probability | 0.01-0.10 |
| `spotting.k_wind` | Wind coefficient | 0.1-0.5 |
| `spotting.I_critical` | Critical intensity (BTU/ft²/min) | 500-2000 |
| `spotting.d_mean` | Mean distance | 0.1-0.3 |
| `spotting.d_sigma` | Distance std dev (lognormal) | 0.3-1.0 |
| `spotting.lateral_spread_angle` | Lateral angle (degrees) | 10-30 |
| `spotting.spot_radius` | Spot fire radius | 0.01-0.05 |
| `spotting.check_interval` | Check every N steps | 1-10 |

See `SPOTTING_MODEL.md` for comprehensive parameter documentation.

### Example Usage

**Basic spotting with moderate frequency:**
```bash
./build/levelset skip_levelset=1 farsite.enable=1 spotting.enable=1 \
  u_x=0.4 spotting.P_base=0.03 spotting.d_mean=0.15
```

**Aggressive spotting with high frequency:**
```bash
./build/levelset skip_levelset=1 farsite.enable=1 spotting.enable=1 \
  u_x=0.5 spotting.P_base=0.08 spotting.k_wind=0.2 \
  spotting.check_interval=3 spotting.d_mean=0.2
```

**Using exponential distance model:**
```bash
./build/levelset skip_levelset=1 farsite.enable=1 spotting.enable=1 \
  u_x=0.4 spotting.distance_model=exponential spotting.d_lambda=8.0
```

### Expected Behavior

A properly configured spotting simulation should show:
1. Fire starts at initial ignition point
2. FARSITE propagates fire elliptically in wind direction
3. New spot fires appear ahead of main fire (primarily downwind)
4. Distance distribution follows chosen model (lognormal/exponential)
5. Some lateral dispersion perpendicular to wind
6. Spot fires grow and eventually merge with main fire
7. Multiple generations of spotting possible

### Output Fields

Spotting adds four fields to plotfiles:
- `spot_prob`: Spotting probability (0.0-1.0)
- `spot_count`: Number of firebrands generated
- `spot_dist`: Spotting distance field
- `spot_active`: Active spot fire flag (0 or 1)

### Validation

Validate spotting behavior by checking:
- **No wind case**: Minimal or no spotting
- **Low intensity**: Reduced spotting frequency
- **High moisture**: Reduced spotting frequency
- **Reproducibility**: Fixed random seed gives consistent results

### Limitations

Current implementation:
- CPU-only random number generation (not GPU-compatible)
- Simplified physics (no firebrand lofting height or burnout)
- No fuel receptivity variations
- No terrain interaction in landing probability

See `SPOTTING_MODEL.md` for detailed physics, implementation, and future enhancements.

## Notes and Limitations

### Current Implementation
- Terrain is currently uniform (constant slope/aspect)
- Wind field is constant in space and time
- No fuel or moisture variability
- Firebrand spotting available (see Spotting Model section above)

### Future Enhancements
- Variable terrain from elevation data
- Non-uniform wind fields from file input
- Fuel bed properties and moisture content
- Two-way coupling with atmospheric models
- Enhanced spotting with firebrand physics

### Numerical Considerations
- WENO5-Z scheme maintains accuracy on complex fire fronts
- RK3 time integration provides stability
- Reinitialization preserves signed distance property
- Artificial viscosity (ε) prevents front sharpening

## Validation

The models have been implemented following published formulations, but should be validated against:
- Experimental fire data
- Operational fire model outputs (FARSITE, FlamMap)
- Field observations

Users should verify results are physically reasonable for their specific applications.
