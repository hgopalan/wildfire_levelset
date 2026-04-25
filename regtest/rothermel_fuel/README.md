# Rothermel Fuel Model Test

## Purpose
Tests fire spread with different fuel types from the NFFL (Northern Forest Fire Laboratory) fuel model database.

## Configuration
- **Domain**: 1.0 x 1.0 x 1.0 unit cube
- **Grid**: 64³ cells
- **Initial condition**: Sphere at center with radius 0.1
- **Velocity**: Moderate wind (0.2, 0.1, 0.0)
- **Time steps**: 100
- **Output interval**: Every 10 steps
- **Fuel model**: FM1 (Short Grass)
- **Fuel moisture**: 6%

## Expected Behavior
Fire spreads according to Rothermel's model for short grass fuel. Higher rate of spread expected due to low fuel moisture and fine fuel characteristics.

## Other Fuel Models to Test
You can modify `rothermel.fuel_model` to test other fuels:
- FM1: Short Grass (1 ft)
- FM2: Timber (Grass and Understory)
- FM3: Tall Grass (2.5 ft)
- FM4: Chaparral (6 ft)
- FM5: Brush (2 ft)
- FM8: Closed Timber Litter
- FM10: Timber (Litter and Understory)

## Run Command
```bash
./build/levelset regtest/rothermel_fuel/inputs.i
```
