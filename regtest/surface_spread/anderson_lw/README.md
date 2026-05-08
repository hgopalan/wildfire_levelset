# Anderson L/W Ratio Test

## Purpose
Tests the Anderson (1983) dynamic length-to-width ratio calculation for fire ellipse elongation based on wind speed.

## Configuration
- **Domain**: 1.0 x 1.0 x 1.0 unit cube
- **Grid**: 64³ cells
- **Initial condition**: Box source at center (0.45-0.55)
- **Velocity**: Moderate wind (0.4, 0.0, 0.0) ≈ 10 mph
- **Time steps**: 100
- **Output interval**: Every 10 steps
- **FARSITE**: Enabled with Anderson L/W ratio
- **Fuel**: FM4 Chaparral

## Anderson (1983) Formula
```
L/W = 0.936 * exp(0.2566 * U) + 0.461 * exp(-0.1548 * U) - 0.397
```
where U is wind speed in mph.

## Expected Behavior
- For ~10 mph wind: L/W ≈ 2.5
- Fire ellipse elongates in wind direction
- Coefficients computed dynamically:
  - a = 1.0 (head fire)
  - b ≈ 0.24 (flank fire)
  - c = 0.2 (backing fire)

## Comparison Test
Compare with `farsite_ellipse` test which uses fixed coefficients:
- Fixed L/W = 3.0
- Fixed b = 0.4

The Anderson model produces more realistic ellipse shapes based on actual wind conditions.

## Run Command
```bash
./build/levelset regtest/anderson_lw/inputs.i
```

## References
Anderson, H.E. (1983). "Predicting wind-driven wild land fire size and shape." Research Paper INT-305, USDA Forest Service.
