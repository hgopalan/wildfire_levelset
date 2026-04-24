# Fire Spread Models Documentation

This document provides detailed information about the fire spread models implemented in the wildfire level-set solver.

## Overview

The solver now supports three fire spread modeling approaches:

1. **Basic Model**: Simple wind-driven spread (original implementation)
2. **Rothermel Model**: Includes terrain slope corrections
3. **FARSITE Model**: Full model with Anderson L/W ratio and combined wind/terrain effects

## Model Selection

Use runtime parameters to select the desired model:

```bash
# Basic model (default)
./build/levelset u_x=0.5 u_y=0.0

# Rothermel model with terrain
./build/levelset u_x=0.5 u_y=0.0 use_terrain_effects=true terrain_slope=15.0 terrain_aspect=90.0

# FARSITE model with L/W ratio
./build/levelset u_x=0.5 u_y=0.0 use_terrain_effects=true use_farsite_model=true terrain_slope=15.0 terrain_aspect=90.0
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

## Notes and Limitations

### Current Implementation
- Terrain is currently uniform (constant slope/aspect)
- Wind field is constant in space and time
- No fuel or moisture variability
- Simplified fire physics (no spotting, crown fire, etc.)

### Future Enhancements
- Variable terrain from elevation data
- Non-uniform wind fields from file input
- Fuel bed properties and moisture content
- Two-way coupling with atmospheric models
- Spotting and ember transport

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
