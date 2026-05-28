# Python API Regression Tests

This directory contains regression tests for the Python bindings (`pyWildfire`) that enable programmatic control of the wildfire solver from Python.

## Overview

| Test | Description | Purpose |
|------|-------------|---------|
| `basic_fire_solver/` | Basic fire solver API | Tests initialization, time-stepping, state extraction, and finalization |
| `coupled_wind_fire/` | Coupled wind-fire simulation | Demonstrates integration with external 3D wind solvers |

## Running the Tests

### Prerequisites

1. Build with Python bindings enabled:
   ```bash
   cmake -S . -B build -DLEVELSET_DIM_2D=ON -DLEVELSET_BUILD_PYTHON_BINDINGS=ON
   cmake --build build -j
   ```

2. Ensure Python 3.6+ and NumPy are installed

### Run via CTest (Recommended)

From the build directory:
```bash
cd build

# Run all Python API tests
ctest -L python_api --output-on-failure

# Run individual tests
ctest -R python_api_basic_fire_solver --output-on-failure
ctest -R python_api_coupled_wind_fire --output-on-failure
```

### Run Directly

Set `PYTHONPATH` to include the build directory:
```bash
export PYTHONPATH=/path/to/build/python:$PYTHONPATH

# Run basic test
cd regtest/python_api/basic_fire_solver
python3 test_fire_solver.py

# Run coupled test
cd ../coupled_wind_fire
python3 test_coupled_wind_fire.py
```

## Integrating with Mass Consistent Wind Solver

The `coupled_wind_fire` test demonstrates how to integrate with external wind solvers. Currently, it uses synthetic wind data, but it can be easily adapted for real wind solvers.

### Using massconsistent_amr

[massconsistent_amr](https://github.com/hgopalan/massconsistent_amr) provides a mass-consistent diagnostic wind solver with Python bindings.

#### 1. Build massconsistent_amr with Python Bindings

```bash
# Clone the repository
git clone --recurse-submodules https://github.com/hgopalan/massconsistent_amr.git
cd massconsistent_amr

# Build with Python bindings
cmake -S . -B build -DBUILD_PYTHON_BINDINGS=ON
cmake --build build -j
```

#### 2. Set PYTHONPATH

```bash
export PYTHONPATH=/path/to/wildfire_levelset/build/python:\
/path/to/massconsistent_amr/build/python:$PYTHONPATH
```

#### 3. Create Coupled Simulation Script

```python
from wildfire_solver import WildfireSolver
from pyWindSolver import WindSolver  # massconsistent_amr Python module

# Initialize both solvers
fire = WildfireSolver("fire_inputs.i")
wind = WindSolver("wind_inputs.txt")

# Coupled time loop
while fire.time < final_time:
    # Solve wind field
    wind.solve(fire.time)
    u_3d, v_3d, w_3d = wind.get_velocity_arrays()
    
    # Update fire wind
    fire.update_wind_3d(u_3d, v_3d, w_3d, wind.nz, wind.zmin, wind.zmax)
    
    # Advance fire
    fire.step()
    
    # Optional: Extract fire state for two-way coupling
    state = fire.get_state()

fire.finalize()
wind.finalize()
```

### Replacing with Other Wind Solvers

The Python API is designed to work with any wind solver that can provide 3D velocity fields. Here's how to integrate different wind solvers:

#### Option 1: AMReX-based Wind Solver (like massconsistent_amr)

If your wind solver is also AMReX-based:
1. Add Python bindings using pybind11 (follow massconsistent_amr pattern)
2. Implement `get_velocity_arrays()` to return NumPy arrays
3. Follow the coupling example above

#### Option 2: External Wind Solver (WindNinja, WRF, etc.)

For external solvers without native Python bindings:

```python
import subprocess
import numpy as np
from wildfire_solver import WildfireSolver

def run_wind_solver_and_load(time, fire_domain):
    """
    Run external wind solver and load results
    """
    # Example: Call WindNinja
    subprocess.run([
        'WindNinja_cli',
        '--domain', f'{fire_domain.xmin},{fire_domain.xmax},...',
        '--time', str(time),
        '--output', 'wind_output.asc'
    ])
    
    # Read output and convert to 3D arrays
    u_3d, v_3d, w_3d = read_wind_output('wind_output.asc')
    
    return u_3d, v_3d, w_3d

# Use in coupled simulation
fire = WildfireSolver("inputs.i")
while fire.time < final_time:
    u_3d, v_3d, w_3d = run_wind_solver_and_load(fire.time, fire)
    fire.update_wind_3d(u_3d, v_3d, w_3d, nz, zmin, zmax)
    fire.step()
```

#### Option 3: WRF Wind Fields

Extract wind from WRF output files:

```python
from netCDF4 import Dataset
import wrf
from wildfire_solver import WildfireSolver

def extract_wrf_wind(wrfout_file, time_idx):
    """Extract 3D wind from WRF output"""
    ncfile = Dataset(wrfout_file)
    u = wrf.getvar(ncfile, 'ua', timeidx=time_idx)  # m/s
    v = wrf.getvar(ncfile, 'va', timeidx=time_idx)
    w = wrf.getvar(ncfile, 'wa', timeidx=time_idx)
    
    # Interpolate to fire grid (implementation needed)
    u_interp, v_interp, w_interp = interpolate_to_fire_grid(u, v, w)
    
    return u_interp, v_interp, w_interp

# Use in simulation
fire = WildfireSolver("inputs.i")
for time_idx in range(num_wrf_times):
    u_3d, v_3d, w_3d = extract_wrf_wind('wrfout_d01', time_idx)
    fire.update_wind_3d(u_3d, v_3d, w_3d, nz, zmin, zmax)
    fire.step()
```

### Wind Data Requirements

For any wind solver, the data must satisfy:

| Requirement | Description |
|-------------|-------------|
| **Array shape** | `(nz, ny, nx)` - vertical levels × latitude × longitude |
| **Array order** | Fortran (column-major) order |
| **Units** | Velocities in m/s |
| **Coordinate system** | Same as fire solver (typically UTM) |
| **Domain overlap** | Wind domain must cover fire domain |
| **Vertical extent** | `zmin` to `zmax` in meters above ground |

Example:
```python
# Create wind arrays
nz, ny, nx = 10, 50, 50  # 10 vertical levels, 50x50 horizontal
u_3d = np.zeros((nz, ny, nx), dtype=np.float64, order='F')
v_3d = np.zeros((nz, ny, nx), dtype=np.float64, order='F')
w_3d = np.zeros((nz, ny, nx), dtype=np.float64, order='F')

# Fill with wind data (your wind solver output)
# ...

# Pass to fire solver
zmin, zmax = 0.0, 100.0  # vertical domain
fire.update_wind_3d(u_3d, v_3d, w_3d, nz, zmin, zmax)
```

## Current Limitations

1. **Single-level fire solver**: The fire solver operates on a 2D horizontal domain. The 3D wind field is column-averaged before being used.

2. **Domain matching**: Wind and fire domains must overlap. The fire solver will interpolate wind data to its grid.

3. **Time synchronization**: The user is responsible for ensuring wind and fire solvers are synchronized in time.

4. **One-way coupling only**: Currently, fire does not feed back to the wind solver (no heat release effects on wind). Two-way coupling requires additional implementation.

5. **MPI support**: MPI-parallel simulations are not yet fully tested with Python bindings.

6. **GPU support**: GPU builds work, but data transfer between GPU and Python is not optimized.

## Future Enhancements

Planned features for the Python API:

- [ ] Two-way coupling (fire heat release → wind buoyancy)
- [ ] Sub-cycling support (different timesteps for wind/fire)
- [ ] MPI-parallel Python interface
- [ ] GPU-aware Python bindings (zero-copy on GPU)
- [ ] Simplified wind solver discovery/loading
- [ ] Pre-built interfaces for common wind solvers (WRF, WindNinja, etc.)

## See Also

- [PYTHON_API_IMPLEMENTATION.md](../../../PYTHON_API_IMPLEMENTATION.md) - Detailed API documentation
- [massconsistent_amr](https://github.com/hgopalan/massconsistent_amr) - Mass-consistent wind solver
- [Full Documentation](https://hgopalan.github.io/wildfire_levelset/) - Complete wildfire_levelset docs
