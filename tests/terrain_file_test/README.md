# Terrain File Test

This test demonstrates the use of spatially-varying terrain slopes read from a file.

## Test Setup

- **Domain**: 1000m x 1000m (64 x 64 cells)
- **Fuel Model**: FM4 (Chaparral)
- **Wind**: 2.0 m/s in +X direction
- **Terrain**: Simple slope rising from 0m to 200m over 1000m horizontal distance
  - Slope angle: ~11.3 degrees (tan ≈ 0.2)
  - Slope direction: uphill in +X direction
- **Initial Fire**: Sphere centered at (500, 500) with radius 50m

## Expected Behavior

The fire should spread faster uphill (+X direction) due to the terrain slope effect in the Rothermel model. The slope factor φ_s is computed from the terrain data at each grid cell, allowing for spatially-varying fire spread rates based on local terrain characteristics.

## Running the Test

From the repository root:

```bash
cd tests/terrain_file_test
../../build/levelset inputs.i
```

## Terrain File Format

The `simple_slope_terrain.txt` file contains elevation data in X Y Z format:
- Lines starting with '#' are comments
- Each data line has three values: X-coordinate, Y-coordinate, Elevation
- The terrain slopes are automatically computed using inverse distance weighting interpolation

## Comparison

To compare with constant slope behavior, you can run with constant slopes instead:

```bash
../../build/levelset inputs.i rothermel.terrain_file="" rothermel.slope_x=0.2 rothermel.slope_y=0.0
```

This should produce similar results since the terrain in this test case has a uniform slope.
