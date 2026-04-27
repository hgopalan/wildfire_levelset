# Bulk Fuel Consumption Fraction Regression Test

## Purpose
Tests the bulk fuel consumption fraction model that computes what fraction of available fuel is consumed as fire passes through.

## Configuration
- **Domain**: 1.0 x 1.0 x 1.0 unit cube
- **Grid**: 64³ cells
- **Initial condition**: Sphere at center with radius 0.12
- **Velocity**: Moderate wind (0.25, 0.1, 0.0)
- **Time steps**: 80
- **Output interval**: Every 10 steps
- **Fuel model**: FM4 (Chaparral, 6 ft)
- **Fuel moisture**: 7%

## Fuel Consumption Parameters
- **tau_residence**: 60.0 seconds (moderate residence time)
- **f_consumed_min**: 0.5 (50% minimum consumption for fast/low-intensity fires)
- **f_consumed_max**: 0.9 (90% maximum consumption for slow/high-intensity fires)
- **Anderson L/W**: Enabled for realistic fire shape

## Expected Behavior
1. Fire spreads according to FARSITE model
2. `fuel_consumption` field computed at each fire front cell
3. Low intensity areas: f_c ≈ 0.5 (50% consumption)
4. High intensity areas: f_c ≈ 0.9 (90% consumption)
5. Transition between min and max follows sigmoid based on normalized intensity

## Model Physics
- **Normalized intensity**: I_norm = (I_R / I_ref) * √(τ / τ_ref)
- **Reference intensity**: I_ref = 1000 BTU/ft²/min
- **Reference residence**: τ_ref = 60 seconds
- **Transition function**: f_c = f_c_min + (f_c_max - f_c_min) * [0.5 * (1 + tanh(I_norm - 1))]

## Physical Interpretation
- **I_R < 500 BTU/ft²/min**: Fast grass fires → f_c ≈ 0.5
- **500 < I_R < 2000**: Transition zone → 0.5 < f_c < 0.9
- **I_R > 2000 BTU/ft²/min**: Crown fires → f_c ≈ 0.9

## Output Fields
- `fuel_consumption`: Computed consumption fraction (0.0-1.0)
- `phi`: Fire perimeter indicator function
- Other standard FARSITE output fields

## Run Command
```bash
./build/levelset regtest/bulk_fuel_consumption/inputs.i
```

## Success Criteria
- Simulation completes without errors
- `fuel_consumption` values between 0.5 and 0.9
- Higher values at more intense fire front locations
- Spatially varying consumption based on local fire intensity

## Validation
Check that:
- All consumption values in valid range [0.5, 0.9]
- Values increase with fire intensity
- Smooth spatial variation (no discontinuities)
