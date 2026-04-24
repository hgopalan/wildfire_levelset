# Gaussian Hill Test Case with Anderson L/W Ratio

## Overview

This test case demonstrates wildfire spread on a Gaussian-shaped hill using the Anderson (1983) length-to-width (L/W) ratio formulation for elliptical fire spread, combined with Rothermel (1972) terrain effects.

## Physical Setup

### Terrain
- **Type**: Simplified Gaussian hill representation
- **Slope**: Constant slope approximation (15% grade / 8.5° in x-direction)
- **Actual Gaussian hills** have spatially varying slope: z(x,y) = H·exp(-r²/(2σ²))
- This simplified version uses constant slope to represent the hillside conditions

### Fire Source
- **Location**: Base of hill at (500m, 500m, 10m)
- **Initial shape**: Sphere with 25m radius
- **Fuel model**: Anderson fuel model 4 (Chaparral)

### Wind Conditions
- **Speed**: 5 m/s (~11 mph) in x-direction (upslope)
- **Expected L/W ratio**: ~2.3 (from Anderson 1983 formula)
- **Effect**: Fire elongates in wind direction, enhanced by terrain slope

## Fire Spread Models

### Rothermel (1972) Model
- Computes base rate of spread (ROS) including:
  - Fuel properties (chaparral characteristics)
  - Wind speed effects (φ_w factor)
  - **Terrain slope effects** (φ_s factor): ~5.275·tan²(slope)
  - Moisture content

### FARSITE Ellipse Model with Anderson L/W Ratio
- **Anderson (1983) L/W formula**: 
  ```
  L/W = 0.936·exp(0.2566·U) + 0.461·exp(-0.1548·U) - 0.397
  ```
  where U is wind speed in mph
  
- **Dynamic behavior**: L/W ratio varies with local wind speed
- **Conversion**: L/W ratio → Richards' elliptical coefficients (a, b, c)
- **Terrain integration**: Slope factor from Rothermel enhances upslope spread

## Expected Results

1. **Fire shape**: Elliptical, elongated in x-direction (wind/upslope)
2. **L/W ratio**: ~2.3 for 11 mph wind
3. **Enhanced upslope spread**: Combined wind + terrain effects
4. **Asymmetric spread**: 
   - Head fire (downwind/upslope): maximum spread
   - Flank fire (crosswind): intermediate spread  
   - Backing fire (upwind/downslope): minimum spread

## Comparison with Fixed Coefficients

To compare Anderson L/W vs. fixed Richards' coefficients:

1. **Anderson mode** (current): `farsite.use_anderson_LW = 1`
   - L/W adapts to wind speed dynamically
   
2. **Fixed mode**: `farsite.use_anderson_LW = 0`
   - Uses constant coefficients a=1.0, b=0.5, c=0.2
   - Equivalent to fixed L/W ≈ 1.2

## Running the Test

```bash
cd tests/gaussian_hill_anderson
../../build/main inputs.i
```

## References

- Anderson, H.E. (1983). "Predicting wind-driven wild land fire size and shape." 
  Research Paper INT-305, USDA Forest Service.
  
- Rothermel, R.C. (1972). "A mathematical model for predicting fire spread in 
  wildland fuels." Research Paper INT-115, USDA Forest Service.
  
- Richards, G.D. (1990). "An elliptical growth model of forest fire fronts and 
  its numerical solution." International Journal of Numerical Methods in Engineering.
