# Velocity File Test

This test demonstrates the velocity field initialization from unstructured 2D data.

## Files
- `velocity_data.txt`: Example unstructured velocity data in X Y U V format
- `inputs.i`: Input file that specifies `velocity_file` parameter

## Data Format
The velocity file should contain unstructured 2D data with four columns:
```
X Y U V
```
where:
- X, Y: Spatial coordinates
- U, V: Velocity components in X and Y directions

Lines starting with `#` are treated as comments and empty lines are ignored.

## Interpolation Method
The code uses Inverse Distance Weighting (IDW) interpolation with power=2 to map the unstructured data points onto the structured AMReX MultiFab grid.

## Usage
In your inputs file, specify:
```
velocity_file = path/to/velocity_data.txt
```

If `velocity_file` is not specified or is empty, the code falls back to using constant velocity values (`u_x`, `u_y`, `u_z`).

## Note
This feature is only available in 2D mode (when `LEVELSET_DIM_2D=ON`). In 3D mode, the code will always use constant velocity initialization.
