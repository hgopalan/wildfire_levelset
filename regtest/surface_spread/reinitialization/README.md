# Level-Set Reinitialization Test

## Purpose
Tests periodic reinitialization to maintain the signed distance property of the level-set function during advection.

## Configuration
- **Domain**: 1.0 x 1.0 x 1.0 unit cube
- **Grid**: 64³ cells
- **Initial condition**: Sphere at (0.25, 0.25, 0.25) with radius 0.15
- **Velocity**: Diagonal (0.3, 0.2, 0.1)
- **Time steps**: 200
- **Output interval**: Every 10 steps
- **Reinitialization**:
  - Frequency: Every 5 time steps (aggressive)
  - Iterations: 30
  - Time step: 0.3

## Background
During advection, the level-set function φ can deviate from being a signed distance function (|∇φ| = 1). Reinitialization solves:

```
∂φ/∂τ + sign(φ₀)(|∇φ| - 1) = 0
```

to restore the signed distance property while preserving the zero level-set (the interface).

## Expected Behavior
- The sphere maintains its shape better with frequent reinitialization
- Interface remains sharp and well-defined
- No artificial diffusion or numerical artifacts
- Compare with `basic_levelset` which uses less frequent reinitialization (every 20 steps)

## Reinitialization Parameters
- `reinit_int`: Steps between reinitialization (-1 to disable)
- `reinit_iters`: Number of pseudo-time iterations
- `reinit_dtau`: Pseudo-time step size

## Run Command
```bash
./build/levelset regtest/reinitialization/inputs.i
```

## Analysis
Compare outputs with `basic_levelset`:
- Check |∇φ| - should be closer to 1.0
- Interface sharpness
- Shape preservation during advection
