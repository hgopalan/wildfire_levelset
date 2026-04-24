# Gaussian Hill with Fixed Richards' Coefficients

## Overview

This is a comparison test case to `gaussian_hill_anderson`. It uses the same physical setup but with **fixed Richards' elliptical coefficients** instead of the dynamic Anderson (1983) L/W ratio.

## Key Difference

| Aspect | gaussian_hill_anderson | gaussian_hill_fixed_coeffs |
|--------|------------------------|----------------------------|
| L/W Mode | Anderson (1983) - Dynamic | Fixed Richards' coefficients |
| L/W at 5 m/s wind | ~2.3 | ~1.2 |
| Parameter | `farsite.use_anderson_LW = 1` | `farsite.use_anderson_LW = 0` |
| Fire elongation | More elongated | Less elongated |

## Expected Results

The fire spread should be:
- **Less elongated** in the wind direction compared to Anderson test
- Still affected by terrain slope (Rothermel φ_s factor)
- More circular overall shape due to lower L/W ratio

## Purpose

Compare these two test cases side-by-side to see:
1. Effect of Anderson's wind-speed-dependent L/W ratio
2. Difference between dynamic and fixed ellipse parameterization
3. How terrain effects (slope factor) work in both models

## Running Both Tests

```bash
# Run Anderson L/W test
cd tests/gaussian_hill_anderson
../../build/main inputs.i

# Run fixed coefficients test
cd ../gaussian_hill_fixed_coeffs
../../build/main inputs.i
```

Then compare the output plotfiles to see the difference in fire shape evolution.
