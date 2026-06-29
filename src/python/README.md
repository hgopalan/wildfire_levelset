# pyWildfire Python Bindings

Python bindings for wildfire_levelset that enable:
1. **Direct wind data loading** from numpy arrays without file I/O
2. **Fire solver control** for time-stepping simulations from Python
3. **Coupled wind-fire simulations** with external wind solvers

## New in v0.2.0: Fire Solver Control

The Python bindings now include complete fire solver control, allowing you to:
- Initialize and run wildfire simulations entirely from Python
- Couple with external wind solvers (e.g., massconsistent_amr)
- Extract fire state (ROS, intensity, flame length) at each timestep
- Update wind fields dynamically during simulation

## Features

- Load 3D wind velocity fields from numpy arrays
- Initialize and control the fire solver from Python
- Automatic computation of 2D column-averaged wind for fire spotting
- Time-stepping with automatic CFL condition
- State extraction (phi, ROS, intensity, flame length, wind)
- Wind field updates for coupled simulations
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

### Fire Solver Control (New in v0.2.0)

#### Quick Start: Simple Fire Simulation

```python
from wildfire_solver import WildfireSolver

# Initialize solver from inputs file
fire = WildfireSolver("inputs.i")

# Run for 10 timesteps
for i in range(10):
    fire.step()
    state = fire.get_state()
    print(f"Step {i+1}: t={state['time']:.2f} s, "
          f"burned area={np.sum(state['phi'] <= 0) * fire.dx * fire.dy:.1f} m²")

# Clean up
fire.finalize()
```

#### Coupled Wind-Fire Simulation

```python
from wildfire_solver import WildfireSolver
import numpy as np

# Initialize fire solver
fire = WildfireSolver("inputs.i")

# Time loop with external wind solver
for n in range(num_steps):
    # 1. Solve wind field (e.g., massconsistent_amr)
    # u_3d, v_3d, w_3d = wind_solver.solve(fire.time)
    
    # For demonstration, use synthetic wind
    u_3d = np.full((nz, fire.ny, fire.nx), 5.0)  # 5 m/s easterly
    v_3d = np.zeros((nz, fire.ny, fire.nx))
    w_3d = np.zeros((nz, fire.ny, fire.nx))
    
    # 2. Pass wind to fire solver
    fire.update_wind_3d(u_3d, v_3d, w_3d, nz, zmin=0.0, zmax=100.0)
    
    # 3. Advance fire simulation
    fire.step()
    
    # 4. Optionally extract state for analysis
    if n % 10 == 0:
        state = fire.get_state()
        fire.write_plotfile(f"plt{n:05d}")

fire.finalize()
```

#### Using the Low-Level API

```python
import pyWildfire

# Initialize solver
result = pyWildfire.initialize("inputs.i")
if result['success']:
    print(f"Initialized {result['nx']}×{result['ny']} grid")

# Time-stepping loop
for step in range(10):
    # Advance one timestep
    step_result = pyWildfire.advance()
    
    # Extract state
    state = pyWildfire.get_state()
    phi = state['phi']              # Level set (< 0 = burned)
    ros = state['ros']              # Rate of spread (m/s)
    intensity = state['intensity']  # Fire line intensity (kW/m)
    flame_length = state['flame_length']  # Flame length (m)
    
    # Update wind if needed
    u_new = np.full((state['ny'], state['nx']), 10.0)
    v_new = np.zeros((state['ny'], state['nx']))
    pyWildfire.update_wind(u_new, v_new)

# Write final plotfile
pyWildfire.write_plotfile("plt_final")

# Clean up
pyWildfire.finalize()
```

### Wind Data Loading (Original Functionality)

#### Basic Example: Loading from NumPy Arrays

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

### Fire Solver Control Functions

#### `pyWildfire.initialize(inputs_file)`

Initialize the wildfire solver from an inputs file.

**Parameters:**
- `inputs_file` (str): Path to the inputs file

**Returns:**
- dict with keys:
  - `success` (bool): Whether initialization succeeded
  - `nx, ny` (int): Grid dimensions
  - `xmin, xmax, ymin, ymax` (float): Domain bounds (meters)
  - `dx, dy` (float): Cell sizes (meters)

#### `pyWildfire.advance()`

Advance the fire simulation by one timestep.

**Returns:**
- dict with keys:
  - `success` (bool): Whether timestep succeeded
  - `dt` (float): Timestep used (seconds)
  - `time` (float): Current simulation time (seconds)
  - `step` (int): Current timestep number

#### `pyWildfire.get_state()`

Extract the current state of the fire simulation.

**Returns:**
- dict with keys:
  - `time` (float): Current simulation time (seconds)
  - `step` (int): Current timestep number
  - `nx, ny` (int): Grid dimensions
  - `xmin, xmax, ymin, ymax` (float): Domain bounds (meters)
  - `dx, dy` (float): Cell sizes (meters)
  - `phi` (ndarray): Level set field (ny, nx), < 0 burned, > 0 unburned
  - `ros` (ndarray): Rate of spread (m/s), shape (ny, nx)
  - `intensity` (ndarray): Fire line intensity (kW/m), shape (ny, nx)
  - `flame_length` (ndarray): Flame length (m), shape (ny, nx)
  - `u_wind, v_wind` (ndarray): Wind components (m/s), shape (ny, nx)
  - `arrival_time` (ndarray): Fire arrival time (s), shape (ny, nx)

#### `pyWildfire.update_wind(u_wind, v_wind)`

Update the wind field from external 2D arrays.

**Parameters:**
- `u_wind` (ndarray): U-component of wind (m/s), shape (ny, nx)
- `v_wind` (ndarray): V-component of wind (m/s), shape (ny, nx)

**Returns:**
- bool: True if wind was updated successfully

#### `pyWildfire.update_wind_3d(nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax, u_array, v_array, w_array)`

Update the wind field from external 3D arrays.

**Parameters:**
- `nx, ny, nz` (int): Grid dimensions
- `xmin, xmax, ymin, ymax, zmin, zmax` (float): Domain bounds (meters)
- `u_array, v_array, w_array` (ndarray): 3D wind velocity components (m/s), flattened in Fortran order

**Returns:**
- bool: True if wind was updated successfully

#### `pyWildfire.write_plotfile(plotfile_name)`

Write the current state to an AMReX plotfile.

**Parameters:**
- `plotfile_name` (str): Name/path of the plotfile to write

**Returns:**
- bool: True if plotfile was written successfully

#### `pyWildfire.finalize()`

Clean up and finalize the fire solver.

#### `pyWildfire.is_initialized()`

Check if the fire solver is currently initialized.

**Returns:**
- bool: True if solver is initialized

### Wind Data Loading Functions

#### `pyWildfire.load_wind_from_arrays()`

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

## Examples

### Complete Examples Provided

1. **`test_fire_solver_api.py`** - Comprehensive test suite for the fire solver API
   ```bash
   PYTHONPATH=build/python python3 src/python/test_fire_solver_api.py
   ```

2. **`coupled_wind_fire_example.py`** - Demonstration of coupled wind-fire simulation
   ```bash
   PYTHONPATH=build/python python3 src/python/coupled_wind_fire_example.py inputs.i
   ```

3. **`wildfire_solver.py`** - High-level Python class (can be run standalone)
   ```bash
   PYTHONPATH=build/python python3 src/python/wildfire_solver.py inputs.i
   ```

### Running Examples

```bash
# Set PYTHONPATH to find the module
export PYTHONPATH=$PWD/build/python:$PYTHONPATH

# Run comprehensive test suite
python3 src/python/test_fire_solver_api.py

# Run coupled wind-fire example
python3 src/python/coupled_wind_fire_example.py \
    regtest/surface_spread/farsite_ellipse/inputs.i

# Run simple fire simulation
python3 src/python/wildfire_solver.py \
    regtest/surface_spread/farsite_ellipse/inputs.i
```

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

## Wind-Fire Coupling

The Python API enables flexible coupling with external wind solvers, supporting both
one-way and two-way coupling modes. This section describes the coupling capabilities
and how to use them.

### Coupling Modules

Two modules are available for wind-fire coupling:

1. **`levelset_coupling.py`** - High-level module for dedicated wind-fire coupling
   - Recommended for most users
   - Mirrors the design from massconsistent_amr
   - Handles all coupling workflows automatically
   - See `example_two_way_coupling.py` for usage

2. **`coupled_solver.py`** - Lower-level coupling framework
   - Provides building blocks for custom coupling workflows
   - More flexible for advanced use cases

### Coupling Modes

#### 1. One-Way Coupling (wind → fire)

Wind field is computed independently; fire responds to wind but does NOT affect wind.

**Use cases:**
- Pre-computed wind fields from weather models
- Coupling with atmospheric simulators
- Scenarios where fire-induced wind effects are negligible
- Faster simulations

**Workflow:**
```python
from wildfire_solver import WildfireSolver
import numpy as np

fire = WildfireSolver("fire_inputs.i")

for step in range(num_steps):
    # 1. Get wind field (external source, e.g., massconsistent_amr)
    u_3d, v_3d, w_3d = wind_solver.get_velocity()  # shape: (nz, ny, nx)
    
    # 2. Update fire with wind
    fire.update_wind_3d(u_3d, v_3d, w_3d, nz, zmin, zmax)
    
    # 3. Advance fire
    fire.step()

fire.finalize()
```

**Example:**
```bash
PYTHONPATH=build/python python3 src/python/example_one_way_coupling.py \
    regtest/surface_spread/farsite_ellipse/inputs.i
```

#### 2. Two-Way Coupling (wind ↔ fire)

Wind field is computed WITH fire heating effects; fire heating is extracted
and fed back to wind solver for fire-induced wind changes.

**Use cases:**
- Fire-atmosphere interaction research
- Understanding fire-induced updrafts and flow deflection
- Scenarios where fire creates significant atmospheric changes
- Detailed fire behavior studies

**Workflow:**
```python
from levelset_coupling import CoupledWindFireSimulation

# Initialize with two-way coupling
coupled = CoupledWindFireSimulation(
    wind_inputs="wind_inputs.i",
    fire_inputs="fire_inputs.i",
    coupling_mode='two_way'
)

# Run simulation
results = coupled.run(
    final_time=3600.0,      # 1 hour
    plot_interval=600.0,    # Write plots every 10 minutes
    callback=lambda step, result: print(f"Step {step}: Heat={result.get('fire_state', {}).get('phi', 0).sum()}")
)

coupled.finalize()
```

**Detailed workflow:**
```
Loop:
  1. Solve wind field with fire heat source from previous step
  2. Extract 3D wind velocity
  3. Update fire with wind field
  4. Advance fire simulation
  5. Extract heat release from fire
  6. Add heat to wind solver (for next iteration)
  Repeat
```

**Example:**
```bash
# Requires massconsistent_amr to be available
PYTHONPATH=build/python:$MASSCONSISTENT_PYTHONPATH python3 \
    src/python/example_two_way_coupling.py wind_inputs.i fire_inputs.i
```

### Heat Flux Extraction for Two-Way Coupling

The fire solver provides methods to extract surface heat flux:

```python
from wildfire_solver import WildfireSolver

fire = WildfireSolver("fire_inputs.i")

# ... simulate for a while ...

# Extract surface heat flux
heat_data = fire.get_surface_fluxes()
heat_flux = heat_data['heat_flux']  # shape: (ny, nx), units: kW/m²

# Or compute from state
state = fire.get_state()
heat_release = fire.compute_heat_release(state)
surface_flux = heat_release['surface_flux']  # kW/m²

# Pass to wind solver
grid_info = {
    'xmin': fire.xmin,
    'xmax': fire.xmax,
    'ymin': fire.ymin,
    'ymax': fire.ymax,
    'dx': fire.dx,
    'dy': fire.dy
}
wind_solver.add_heat_source(heat_flux, grid_info)
```

### Integration with massconsistent_amr

The primary use case is coupling with massconsistent_amr (3D mass-consistent wind solver).

#### Current Workflow (File-Based)
```
massconsistent_amr → plotfile (disk) → plt_wind_reader.H → wildfire_levelset
```

#### New Workflow (Memory-Based via Python)
```
massconsistent_amr ↔ pyWindSolver ↔ levelset_coupling ↔ pyWildfire ↔ wildfire_levelset
```

#### Example: Integration with massconsistent_amr

```python
from levelset_coupling import CoupledWindFireSimulation

# Use levelset_coupling module which automatically handles both solvers
coupled = CoupledWindFireSimulation(
    wind_inputs="amr_inputs.i",
    fire_inputs="fire_inputs.i",
    coupling_mode='two_way'  # Enable two-way coupling
)

# Run coupled simulation
results = coupled.run(final_time=3600.0)

# Access results
print(f"Final fire time: {results['final_time']:.1f} s")
print(f"Final fire state: {results['fire_state']}")
```

### Benefits of Python-Based Coupling

- **Eliminates disk I/O overhead** - No intermediate plotfiles needed
- **Faster data transfer** - Direct memory access via numpy arrays
- **Enables coupled simulations** - Single Python script controls both solvers
- **Flexible coupling strategies** - One-way, two-way, sub-cycling, ensemble
- **Easier workflow integration** - Python-based analysis, visualization, ML
- **Interoperability** - Works with massconsistent_amr and other solvers

### Coupling Examples

Three complete examples are provided:

1. **`example_one_way_coupling.py`** - One-way wind→fire coupling
   - Uses synthetic wind or external wind solver
   - Simpler, faster
   - Run: `python3 example_one_way_coupling.py fire_inputs.i`

2. **`example_two_way_coupling.py`** - Two-way wind↔fire coupling
   - Requires wind solver support for `add_heat_source()`
   - More realistic fire-atmosphere interaction
   - Run: `python3 example_two_way_coupling.py wind_inputs.i fire_inputs.i`

3. **`coupled_wind_fire_example.py`** - Advanced coupling patterns
   - Demonstrates sub-cycling, callbacks, diagnostics
   - Run: `python3 coupled_wind_fire_example.py fire_inputs.i wind_inputs.i`

## Roadmap

### Completed (v0.2.0)
- ✅ Fire solver initialization from Python
- ✅ Time-stepping control
- ✅ State extraction (phi, ROS, intensity, etc.)
- ✅ Wind field updates (2D and 3D)
- ✅ Plotfile writing
- ✅ High-level WildfireSolver class
- ✅ Coupled simulation example
- ✅ Heat flux extraction (for two-way coupling)
- ✅ levelset_coupling module for dedicated wind-fire coupling
- ✅ One-way coupling (wind → fire)
- ✅ Two-way coupling framework (wind ↔ fire with heat feedback)
- ✅ Comprehensive coupling examples

### Planned
- ⬜ Direct pyAMReX MultiFab support (zero-copy)
- ⬜ Advanced diagnostics (burn probability, isochrones)
- ⬜ Parallel execution (MPI support)
- ⬜ GPU acceleration via Python
- ⬜ Sub-cycling (different timesteps for wind and fire)

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
