# External Terrain and Wind Test

## Purpose
Tests fire spread over a Gaussian hill with spatially-varying terrain slopes and wind field. This demonstrates the code's capability to handle external terrain data and non-uniform wind fields.

## Configuration
- **Domain**: 1000 m x 1000 m
- **Grid**: 100 x 100 cells (10 m resolution)
- **Terrain**: Gaussian hill
  - Center: (500, 500) m
  - Height: 100 m
  - Width (σ): 150 m
- **Wind field**: Spatially-varying
  - Base wind: 5 m/s in x-direction
  - Speedup over hill crest: up to 7.5 m/s
  - Includes turbulent fluctuations
- **Initial condition**: Line fire at x = 100-120 m, y = 200-800 m
- **Time steps**: 200
- **Output interval**: Every 20 steps
- **Fuel**: FM4 Chaparral
- **FARSITE**: Enabled with Anderson L/W ratio

## Terrain Model
The Gaussian hill is defined by:
```
z(x,y) = H * exp(-r²/(2σ²))
```
where:
- H = 100 m (peak height)
- σ = 150 m (width parameter)
- r² = (x - x_center)² + (y - y_center)²

## Wind Model
The wind field approximates flow over a hill with:
```
u(x,y) = u_base * (1 + 0.5 * z/H) + turbulence
```
This gives up to 50% wind speedup at the hill crest.

## Expected Behavior
- Fire starts as a line at the western edge
- Wind drives fire eastward toward and over the hill
- Fire accelerates upslope due to terrain effect (φ_s factor)
- Wind speedup at crest further accelerates fire
- Elliptical spread pattern with Anderson L/W ratio
- Fire shape influenced by both wind and terrain

## Files
- `inputs.i`: Input parameters
- `gaussian_hill_terrain.csv`: Terrain elevation data (2500 points)
- `gaussian_hill_wind.csv`: Wind velocity data (2500 points)

## Run Command
```bash
cd regtest/terrain_wind
../../build/levelset inputs.i
```

## Notes
- This test requires 2D build (`-DLEVELSET_DIM_2D=ON`)
- Terrain slopes computed via inverse distance weighting
- Wind field interpolated to grid using IDW
- Both terrain and wind files use unstructured data format
