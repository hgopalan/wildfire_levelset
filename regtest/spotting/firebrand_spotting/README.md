# Firebrand Spotting Model Test

## Purpose
Tests the firebrand spotting model with FARSITE elliptical fire expansion. This test demonstrates the generation of new ignition points ahead of the main fire based on a probability function.

## Configuration
- **Domain**: 1.0 x 1.0 x 1.0 unit cube
- **Grid**: 64³ cells
- **Initial condition**: Small sphere (radius 0.1) at (0.3, 0.5, 0.5) - point ignition
- **Velocity**: Constant wind (0.4, 0.0, 0.0) in x-direction
- **Time steps**: 50
- **Output interval**: Every 5 steps
- **FARSITE**: Enabled with fixed coefficients
  - a = 1.0 (head fire)
  - b = 0.5 (flank fire)
  - c = 0.2 (backing fire)
  - L/W ratio = 3.0
- **Spotting Parameters**:
  - Base probability: 5%
  - Wind coefficient: 0.3
  - Critical intensity: 800 BTU/ft²/min
  - Mean spotting distance: 0.15 (15% of domain)
  - Distance model: Lognormal with σ = 0.4
  - Lateral spread angle: 20 degrees
  - Spot radius: 0.03 (3% of domain)
  - Check interval: Every 3 timesteps
  - Random seed: 12345 (for reproducibility)

## Expected Behavior
1. Fire starts at initial ignition point (0.3, 0.5, 0.5)
2. FARSITE model propagates fire elliptically in wind direction
3. Spotting model generates new ignition points:
   - Primarily downwind (positive x-direction)
   - Distance follows lognormal distribution centered at 0.15
   - Some lateral dispersion perpendicular to wind
4. Multiple generations of spot fires should appear
5. Spot fires grow and merge with main fire over time

## Validation Points
- **Spotting occurs**: Check that new ignitions appear ahead of main fire
- **Downwind bias**: Spots should be primarily in +x direction
- **Distance distribution**: Most spots around mean distance (0.15)
- **Lateral spread**: Some spots offset perpendicular to wind
- **Reproducibility**: Fixed random seed ensures consistent results

## Output Fields
The plotfiles contain the following fields:
- `phi`: Level-set function (fire indicator)
- `velx`, `vely`, `velz`: Velocity components
- `farsite_dx`, `farsite_dy`, `farsite_dz`: FARSITE spread displacements
- `R`: Rothermel rate of spread
- `spot_prob`: Spotting probability field
- `spot_count`: Number of firebrands generated per cell
- `spot_dist`: Spotting distance field
- `spot_active`: Flag for active spot fires

## Run Command
```bash
./build/levelset regtest/firebrand_spotting/inputs.i
```

## Visualization Tips
When visualizing in ParaView or VisIt:
1. Look at `phi < 0` to see burned regions
2. Examine `spot_prob` to see where spotting is likely
3. Track temporal evolution to observe spot fire generation
4. Compare early timesteps (main fire only) vs. later (with spots)
