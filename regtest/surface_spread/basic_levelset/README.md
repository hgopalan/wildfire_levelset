# Basic Level-Set Advection Test

## Purpose
Tests basic level-set advection with a spherical source and constant velocity field.

## Configuration
- **Domain**: 1.0 x 1.0 x 1.0 unit cube
- **Grid**: 64³ cells
- **Initial condition**: Sphere at (0.3, 0.5, 0.5) with radius 0.15
- **Velocity**: Constant (0.25, 0.1, 0.0)
- **Time steps**: 100
- **Output interval**: Every 20 steps

## Expected Behavior
The sphere should advect diagonally across the domain, maintaining its shape through periodic reinitialization.

## Run Command
```bash
./build/levelset regtest/basic_levelset/inputs.i
```
