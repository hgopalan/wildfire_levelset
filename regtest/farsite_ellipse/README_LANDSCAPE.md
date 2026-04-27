# FARSITE Landscape File Support

## Overview

This directory contains test cases for the FARSITE ellipse model with landscape file support.

## Landscape File Format

The landscape file format is an ASCII text file with the following columns:

```
X Y ELEVATION SLOPE ASPECT FUEL_MODEL
```

Where:
- **X, Y**: Coordinates in meters
- **ELEVATION**: Elevation above sea level in meters
- **SLOPE**: Slope angle in degrees (0-90)
- **ASPECT**: Slope aspect in degrees (0-360)
  - 0° = North
  - 90° = East
  - 180° = South
  - 270° = West
- **FUEL_MODEL**: NFFL fuel model number (optional, defaults to 0)

Lines starting with `#` are treated as comments and ignored.

## Sample Landscape Files

### socal_chaparral_landscape.lcp

A 100m × 100m sample landscape file representing typical Southern California chaparral terrain:

- **Terrain**: Moderate slopes (10-35 degrees) representing coastal mountain foothills
- **Aspect**: Predominantly southwest-facing (fire-prone orientation)
- **Fuel**: NFFL Fuel Model 4 (chaparral)
- **Elevation Range**: 100-188 meters
- **Grid**: 11×11 points spaced 10m apart

This configuration is representative of fire-prone Southern California terrain where Santa Ana wind events drive rapid wildfire spread upslope.

## Usage

To use a landscape file in your simulation:

1. Create or obtain a landscape file in the format described above
2. In your input file, specify the landscape file path:

```
rothermel.landscape_file = socal_chaparral_landscape.lcp
```

3. When a landscape file is specified, the slope and elevation data from any terrain file will be ignored, and the landscape file data will be used instead.

## Test Cases

### inputs_landscape.i

A test case demonstrating:
- Use of landscape file for terrain data
- FARSITE ellipse model with Anderson L/W ratio
- Western wind (simulating Santa Ana conditions)
- Line fire ignition at western edge
- 100m × 100m domain matching the sample landscape

Run this test with:
```bash
./wildfire_levelset inputs_landscape.i
```

## Notes

- The landscape file takes precedence over terrain files when both are specified
- The code uses inverse distance weighting (IDW) interpolation to map landscape data points to the simulation grid
- Slope and aspect from the landscape file are converted to slope_x and slope_y components for use in the Rothermel fire spread model
