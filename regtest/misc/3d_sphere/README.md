# 3D Fire Spread Test

## Purpose
Tests full 3D fire spread simulation with FARSITE elliptical expansion, Rothermel model, and terrain slope in 3D.

## Configuration
- **Domain**: 1.0 x 1.0 x 1.0 unit cube
- **Grid**: 48³ cells (reduced for computational efficiency)
- **Initial condition**: Sphere at (0.3, 0.3, 0.3) with radius 0.12
- **Velocity**: 3D wind (0.25, 0.15, 0.05)
- **Time steps**: 100
- **Output interval**: Every 20 steps
- **FARSITE**: Enabled with Anderson L/W ratio
- **Fuel**: FM3 Tall Grass
- **Terrain**: 10% slope in x-direction

## Features Tested
1. **3D Level-Set Advection**: Full 3D signed distance function evolution
2. **3D FARSITE**: Elliptical expansion in 3D space
3. **3D Rothermel**: Fire spread model with 3D wind and terrain
4. **3D Reinitialization**: Maintaining |∇φ| = 1 in 3D
5. **Terrain Effects**: Upslope fire acceleration in 3D

## Build Requirements
This test requires a 3D build (default):
```bash
cmake -S . -B build
cmake --build build -j
```

For 2D builds, this test will not run correctly.

## Expected Behavior
- Fire spreads as a 3D ellipsoid
- Elongation in primary wind direction (x)
- Some elongation in secondary wind direction (y)
- Upslope acceleration in x-direction
- Spherical shape preserved during advection via reinitialization

## Run Command
```bash
./build/levelset regtest/3d_sphere/inputs.i
```

## Visualization
Use ParaView or VisIt to visualize the 3D fire spread:
1. Open `plt*` directories
2. Display the `phi` field
3. Create isosurface at phi = 0 (fire boundary)
4. Animate over time steps

## Computational Notes
- 3D simulations are more expensive than 2D
- Grid size reduced to 48³ for faster testing
- For production runs, increase to 64³ or higher
