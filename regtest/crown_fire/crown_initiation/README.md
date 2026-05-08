# Crown Fire Initiation Regression Test

## Purpose
Tests Van Wagner's (1977) threshold model for crown fire initiation combined with Rothermel surface fire spread.

## Configuration
- **Domain**: 1.0 x 1.0 x 1.0 unit cube
- **Grid**: 64³ cells
- **Initial condition**: Sphere at center with radius 0.1
- **Velocity**: Moderate wind (0.3, 0.1, 0.0)
- **Time steps**: 50
- **Output interval**: Every 10 steps
- **Fuel model**: FM10 (Timber with Litter and Understory)
- **Fuel moisture**: 8%

## Crown Fire Parameters
- **CBH** (Canopy Base Height): 4.0 m
- **CBD** (Canopy Bulk Density): 0.15 kg/m³
- **FMC** (Foliar Moisture Content): 100%
- **Crown fraction weight**: 1.0 (full crown fire effect)
- **Units**: Metric (m, kW/m)

## Expected Behavior
1. Fire initiates as surface fire
2. When surface fire intensity exceeds critical threshold I_o, crown fire initiates
3. Crown fire spreads faster than surface fire
4. `crown_fraction` field shows where crown fire is active (0.0-1.0)
5. Total fire spread rate increases significantly when crown fire is active

## Van Wagner Model Physics
- **Critical intensity**: I_o = 0.010 * CBH * (460 + 25.9 * FMC)
- For CBH=4.0m, FMC=100%: I_o ≈ 122 kW/m
- **Crown ROS**: R_crown = 3.0 / CBD = 3.0 / 0.15 = 20 m/min
- Crown fire dominates when initiated due to much faster spread rate

## Output Fields
- `phi`: Fire perimeter indicator function
- `crown_fraction`: Fraction of spread from crown fire (0.0 = surface only, 1.0 = crown dominated)
- `fuel_consumption`: Fuel consumption fraction
- Other standard FARSITE output fields

## Run Command
```bash
./build/levelset regtest/crown_initiation/inputs.i
```

## Success Criteria
- Simulation completes without errors
- Crown fraction > 0 at fire front for high intensity areas
- Crown fraction = 0 for low intensity or no fire areas
- Spread rate increases when crown fire initiates
