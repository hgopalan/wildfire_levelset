# FARSITE Ellipse Model Test

## Purpose
Tests the FARSITE elliptical fire expansion model (Richards 1990) with fixed coefficients.

## Configuration
- **Domain**: 1.0 x 1.0 x 1.0 unit cube
- **Grid**: 64³ cells
- **Initial condition**: Box source (0.4-0.6 in all dimensions)
- **Velocity**: Constant wind (0.3, 0.0, 0.0) in x-direction
- **Time steps**: 100
- **Output interval**: Every 10 steps
- **FARSITE**: Enabled with fixed coefficients
  - a = 1.0 (head fire)
  - b = 0.4 (flank fire)
  - c = 0.2 (backing fire)
  - L/W ratio = 3.0

## Expected Behavior
The fire should spread elliptically with elongation in the wind direction. Head fire spreads fastest, backing fire slowest.

## Run Command
```bash
./build/levelset regtest/farsite_ellipse/inputs.i
```
