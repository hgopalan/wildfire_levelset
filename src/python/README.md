# pyWildfire Python Bindings

Python bindings for wildfire_levelset that enable direct loading of 3D wind field data from Python without file I/O.

## Features

- Load 3D wind velocity fields from numpy arrays
- Automatic computation of 2D column-averaged wind for fire spotting
- Zero-copy data transfer (when using pyAMReX MultiFab in future)
- Compatible with massconsistent_amr output via pyAMReX

## Building

### Prerequisites

- Python 3.11+
- NumPy
- CMake 3.20+
- C++17 compiler
- pybind11 (automatically fetched if not found)

### Build Instructions

```bash
cd wildfire_levelset

# Configure with Python bindings enabled
cmake -S . -B build \
  -DLEVELSET_DIM_2D=ON \
  -DLEVELSET_BUILD_PYTHON_BINDINGS=ON

# Build
cmake --build build -j

# The Python module will be in build/python/pyWildfire.*.so
```

### Installation (Optional)

```bash
cmake --install build --prefix /usr/local
```

Or add to PYTHONPATH:

```bash
export PYTHONPATH=$PWD/build/python:$PYTHONPATH
```

## Usage

### Basic Example: Loading from NumPy Arrays

```python
import numpy as np
import pyWildfire

# Define 3D wind grid
nx, ny, nz = 8, 8, 4  # Grid dimensions

# Physical domain (UTM Zone 11N, Southern California)
xmin, xmax = 329900.0, 330500.0  # meters
ymin, ymax = 3774900.0, 3775500.0
zmin, zmax = 0.0, 40.0  # meters above ground level

# Create uniform 3D wind field
u = np.full((nz, ny, nx), 5.0)  # 5 m/s westerly
v = np.full((nz, ny, nx), 0.5)  # 0.5 m/s southerly
w = np.full((nz, ny, nx), 0.0)  # no vertical motion

# Load and compute column-averaged wind
# Arrays must be flattened in Fortran (column-major) order
result = pyWildfire.load_wind_from_arrays(
    nx, ny, nz,
    xmin, xmax, ymin, ymax, zmin, zmax,
    u.flatten('F'),  # Flatten in Fortran order
    v.flatten('F'),
    w.flatten('F')   # Optional
)

# Access results
print(f"2D grid: {result['nx_2d']} × {result['ny_2d']}")
print(f"Column-averaged u-wind: {result['u2d']}")  # shape: (ny, nx)
print(f"Column-averaged v-wind: {result['v2d']}")  # shape: (ny, nx)
```

### Advanced Example: Height-Dependent Wind Profile

```python
import numpy as np
import pyWildfire

nx, ny, nz = 16, 16, 8
xmin, xmax = 329000.0, 331000.0
ymin, ymax = 3774000.0, 3776000.0
zmin, zmax = 0.0, 100.0

# Create log-law wind profile
u_3d = np.zeros((nz, ny, nx))
v_3d = np.zeros((nz, ny, nx))

for k in range(nz):
    z = zmin + (k + 0.5) * (zmax - zmin) / nz
    # Log-law profile: u(z) = (u*/κ) * ln(z/z0)
    u_star = 0.4  # friction velocity
    kappa = 0.4   # von Kármán constant
    z0 = 0.1      # roughness length
    
    u_k = (u_star / kappa) * np.log(max(z, z0) / z0)
    u_3d[k, :, :] = u_k
    v_3d[k, :, :] = 0.0

result = pyWildfire.load_wind_from_arrays(
    nx, ny, nz,
    xmin, xmax, ymin, ymax, zmin, zmax,
    u_3d.flatten('F'), v_3d.flatten('F')
)

print(f"Mean column-averaged wind: {result['u2d'].mean():.2f} m/s")
```

### Future: Integration with pyAMReX (massconsistent_amr)

Once pyAMReX integration is complete, you will be able to:

```python
import pyamrex.space3d as amr
import pyWildfire

# Run massconsistent_amr solver (pseudocode)
wind_mf = solve_wind_field(...)  # Returns pyAMReX MultiFab

# Extract data and load into wildfire_levelset format
u_data = wind_mf.to_numpy('u')
v_data = wind_mf.to_numpy('v')
w_data = wind_mf.to_numpy('w')
geom = wind_mf.geometry()

result = pyWildfire.load_wind_from_multifab(
    wind_mf, geom, var_names=['u', 'v', 'w']
)
```

## API Reference

### `load_wind_from_arrays()`

Load 3D wind field from numpy arrays and compute 2D column-averaged wind.

**Parameters:**
- `nx, ny, nz` (int): Grid dimensions in x, y, z directions
- `xmin, xmax, ymin, ymax, zmin, zmax` (float): Physical domain bounds in meters
- `u_array` (ndarray): 3D u-velocity component (m/s), flattened in Fortran order
- `v_array` (ndarray): 3D v-velocity component (m/s), flattened in Fortran order
- `w_array` (ndarray, optional): 3D w-velocity component (m/s), flattened in Fortran order

**Returns:**
Dictionary with keys:
- `valid` (bool): Whether data loaded successfully
- `nx_2d, ny_2d` (int): 2D grid dimensions
- `n_points` (int): Total number of 3D points
- `xmin, xmax, ymin, ymax` (float): 2D domain bounds
- `u2d` (ndarray): Column-averaged u-component, shape (ny, nx)
- `v2d` (ndarray): Column-averaged v-component, shape (ny, nx)

**Notes:**
- Arrays must be in Fortran (column-major) order with indexing [k, j, i]
- The 2D wind is computed by averaging over all vertical levels at each (x, y) column
- This matches the format expected by wildfire_levelset for spotting calculations

## Testing

Run the test suite:

```bash
# Set PYTHONPATH to find the module
export PYTHONPATH=$PWD/build/python:$PYTHONPATH

# Run tests
python3 src/python/test_pywildfire.py
```

Expected output:
```
pyWildfire version: 0.1.0

======================================================================
Test 1: Basic wind loading from numpy arrays
======================================================================
...
✓ Test PASSED: Column-averaged wind matches input

======================================================================
Test 2: Variable wind field (height-dependent)
======================================================================
...
✓ Test PASSED: Column-averaged wind correct for variable field

======================================================================
SUMMARY
======================================================================
Tests run: 2
Tests passed: 2
Tests failed: 0

✓ ALL TESTS PASSED
```

## Integration with massconsistent_amr

The primary use case for these bindings is to enable direct data transfer from
massconsistent_amr (3D mass-consistent wind solver) to wildfire_levelset without
writing intermediate plotfiles to disk.

### Current Workflow (File-Based)
```
massconsistent_amr → plotfile (disk) → plt_wind_reader.H → wildfire_levelset
```

### New Workflow (Memory-Based)
```
massconsistent_amr → pyAMReX MultiFab → pyWildfire → wildfire_levelset
```

Benefits:
- Eliminates disk I/O overhead
- Faster data transfer (zero-copy when possible)
- Enables coupled simulations in single Python script
- Easier to integrate into workflows

## Troubleshooting

### Import Error

If you get `ImportError: No module named 'pyWildfire'`:
1. Make sure you built with `-DLEVELSET_BUILD_PYTHON_BINDINGS=ON`
2. Set PYTHONPATH: `export PYTHONPATH=/path/to/build/python:$PYTHONPATH`
3. Check that `pyWildfire.*.so` exists in `build/python/`

### Array Size Mismatch

If you get array size mismatch errors:
- Verify your arrays are the correct size: `nx * ny * nz`
- Make sure you flatten with Fortran order: `array.flatten('F')`
- Check that your array dtype is `float64` (double precision)

### AMReX Initialization

The module automatically initializes AMReX on first use. If you see AMReX-related
errors, make sure you're not calling `amrex::Initialize()` elsewhere in your code.

## License

Same as wildfire_levelset (see LICENSE file in repository root).
