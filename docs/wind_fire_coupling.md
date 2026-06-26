# Wind-Fire Coupling Interface

## Overview

This directory contains the Python interface for coupled wind-fire simulations, enabling two-way interaction between the wildfire_levelset fire solver and the massconsistent_amr wind solver.

## Files

### Core Modules
- **`wildfire_solver.py`** - Main fire solver wrapper with heat release and flux computation
- **`coupled_solver.py`** - High-level coupled simulation manager with one-way/two-way modes
- **`coupled_wind_fire_example.py`** - Example showing 3D wind coupling with heat release

### Examples & Demos
- **`coupling_modes_example.py`** - Comparison of one-way vs two-way coupling modes
- **`example_integration.py`** - Wind field integration example
- **`test_fire_solver_api.py`** - Unit tests for fire solver API
- **`test_pywildfire.py`** - Unit tests for pyWildfire bindings

## Coupling Modes

### One-Way Coupling (Wind → Fire)

Wind field is computed independently and affects fire spread, but fire does NOT affect wind.

**Use cases:**
- Pre-computed wind fields
- Offline coupling with weather models
- Simpler, faster simulations

**Workflow:**
```python
from wildfire_solver import WildfireSolver

fire = WildfireSolver("inputs.i")

while fire.time < final_time:
    # Solve wind independently (e.g., from massconsistent_amr)
    u_3d, v_3d, w_3d = wind_solver.get_velocity()
    
    # Update fire with wind
    fire.update_wind_3d(u_3d, v_3d, w_3d, nz, zmin, zmax)
    
    # Advance fire
    fire.step()
    
    # Get heat release (for information, not fed back to wind)
    state = fire.get_state()
    heat_data = fire.compute_heat_release(state)

fire.finalize()
```

### Two-Way Coupling (Wind ↔ Fire)

Wind field is computed with fire heating effects. Fire heating is extracted and fed back to the wind solver.

**Use cases:**
- Detailed fire-atmosphere interaction
- Fire-induced wind changes
- Realistic fire propagation

**Workflow:**
```python
from coupled_solver import CoupledWindFireSolver

# Create coupled system in two-way mode
coupled = CoupledWindFireSolver(
    fire_inputs="fire_inputs.i",
    wind_inputs="wind_inputs.i",
    coupling_mode='two_way'
)

# Run simulation
coupled.run(
    final_time=3600.0,
    wind_update_interval=5,  # Update wind every 5 fire steps
    callback=progress_callback
)

coupled.finalize()
```

The framework will:
1. Solve wind field (with heating from previous step)
2. Extract 3D velocity
3. Update fire solver with wind
4. Advance fire solver
5. Extract heat release and surface fluxes
6. Feed heat back to wind solver via `wind.add_heat_source()`
7. Repeat

**Note:** Two-way coupling requires the wind solver to implement `add_heat_source()` method (currently in development for massconsistent_amr).

## New Methods in WildfireSolver

### Heat Release Computation

```python
fire = WildfireSolver("inputs.i")
fire.step()

# Get comprehensive heat release data
state = fire.get_state()
heat_data = fire.compute_heat_release(state)

print(f"Total heat release: {heat_data['total_heat_release']} kW")
print(f"Fire perimeter: {heat_data['perimeter']} m")
print(f"Mean intensity: {heat_data['mean_intensity']} kW/m")
print(f"Max surface flux: {heat_data['surface_flux'].max()} kW/m²")
```

**Returns dict with:**
- `total_heat_release` (kW) - Total convective heat release
- `heat_release_rate` (kW/m) - Heat release per unit perimeter
- `mean_intensity` (kW/m) - Mean Byram intensity
- `max_intensity` (kW/m) - Maximum fire intensity
- `max_ros` (m/s) - Maximum rate of spread
- `perimeter` (m) - Estimated fire perimeter
- `burned_area` (m²) - Total burned area
- `surface_flux` (array, kW/m²) - Surface heat flux field

### Surface Flux Extraction

```python
fire.step()

# Get surface fluxes for wind solver
flux_data = fire.get_surface_fluxes()

heat_flux = flux_data['heat_flux']  # (ny, nx) array in kW/m²
grid_info = flux_data['grid_info']  # Grid metadata

# Pass to wind solver (when add_heat_source is available)
# wind.add_heat_source(heat_flux, grid_info)
```

### Area and Perimeter Calculations

```python
fire.step()

# Get burned area
burned_area = fire.get_burned_area()  # m²

# Get fire perimeter
perimeter = fire.get_fire_perimeter()  # m
```

## CoupledWindFireSolver Class

High-level manager for coupled simulations.

### Initialization

```python
from coupled_solver import CoupledWindFireSolver

# One-way coupling
coupled = CoupledWindFireSolver(
    fire_inputs="fire.i",
    wind_inputs="wind.i",
    coupling_mode='one_way'  # default
)

# Two-way coupling (awaiting wind.add_heat_source())
coupled = CoupledWindFireSolver(
    fire_inputs="fire.i",
    wind_inputs="wind.i",
    coupling_mode='two_way'
)
```

### Running Simulations

```python
# Run for 1 hour with wind updates every 5 fire steps
coupled.run(
    final_time=3600.0,
    wind_update_interval=5,
    plot_interval=600.0,  # Plotfiles every 10 minutes
    callback=my_callback_function
)

# Or run for fixed number of steps
coupled.run(
    num_steps=100,
    wind_update_interval=10
)
```

### Callbacks

```python
def my_callback(step, result):
    """Called after each coupled step"""
    fire_state = result['fire_state']
    heat_data = result['heat_data']
    wind_info = result['wind_info']
    
    print(f"Step {step}: "
          f"time={fire_state['time']:.1f}s, "
          f"heat={heat_data['total_heat_release']:.0f}kW")

coupled.run(final_time=3600.0, callback=my_callback)
```

## Example Usage

### Simple One-Way Coupling

```python
from wildfire_solver import WildfireSolver
import numpy as np

# Initialize
fire = WildfireSolver("regtest/surface_spread/farsite_ellipse/inputs.i")

# 3D wind parameters (from wind solver or model)
nz, zmin, zmax = 8, 0.0, 100.0

# Run coupled loop
for step in range(50):
    # Get wind from external source
    u_3d, v_3d, w_3d = get_wind_from_model(fire.time)
    
    # Update fire with wind
    fire.update_wind_3d(u_3d, v_3d, w_3d, nz, zmin, zmax)
    
    # Advance fire
    result = fire.step()
    state = fire.get_state()
    
    # Get heat data (for analysis or logging)
    heat = fire.compute_heat_release(state)
    burned_area = heat['burned_area']
    
    print(f"Step {step}: area={burned_area:.0f}m²")

fire.finalize()
```

### Full Two-Way Coupling (when wind supports add_heat_source)

```python
from coupled_solver import CoupledWindFireSolver

# Create two-way coupled system
coupled = CoupledWindFireSolver(
    fire_inputs="fire.i",
    wind_inputs="wind.i",
    coupling_mode='two_way'
)

def progress_callback(step, result):
    state = result['fire_state']
    heat = result['heat_data']
    wind = result['wind_info']
    
    print(f"Step {step}: "
          f"t={state['time']:.1f}s, "
          f"heat={heat['total_heat_release']:.0f}kW, "
          f"wind_iters={wind.get('iters', 0)}")

# Run simulation
coupled.run(
    final_time=3600.0,
    wind_update_interval=10,
    plot_interval=300.0,
    callback=progress_callback
)

# Output
coupled.write_plotfiles("fire_plt", "wind_plt")
coupled.finalize()
```

## Running Examples

### One-Way Coupling Example
```bash
PYTHONPATH=build/python python3 src/python/coupled_wind_fire_example.py \
    regtest/surface_spread/farsite_ellipse/inputs.i
```

### Compare Coupling Modes
```bash
PYTHONPATH=build/python python3 src/python/coupling_modes_example.py \
    regtest/surface_spread/farsite_ellipse/inputs.i
```

## Integration with massconsistent_amr

### Current Status (One-Way Ready)

```python
from wind_solver import WindSolver
from wildfire_solver import WildfireSolver

# Initialize both solvers
wind = WindSolver("wind_inputs.i")
fire = WildfireSolver("fire_inputs.i")

# Verify domain compatibility
assert wind.nx == fire.nx and wind.ny == fire.ny

# One-way coupling loop
for step in range(num_steps):
    # Solve wind
    wind.solve()
    
    # Extract velocity
    vel = wind.get_velocity()
    u_3d, v_3d, w_3d = vel['u'], vel['v'], vel['w']
    
    # Update fire with wind
    fire.update_wind_3d(u_3d, v_3d, w_3d, wind.nz, wind.zmin, wind.zmax)
    
    # Advance fire
    fire.step()
    
    # Get heat data for analysis
    state = fire.get_state()
    heat_data = fire.compute_heat_release(state)

wind.finalize()
fire.finalize()
```

### Future (Two-Way Ready)

When `massconsistent_amr` implements `add_heat_source()`:

```python
from coupled_solver import CoupledWindFireSolver

# Two-way coupling becomes automatic
coupled = CoupledWindFireSolver(
    fire_inputs="fire.i",
    wind_inputs="wind.i",
    coupling_mode='two_way'
)

coupled.run(final_time=3600.0)
coupled.finalize()
```

## Building with Python Bindings

### wildfire_levelset
```bash
cd wildfire_levelset
cmake -S . -B build -DLEVELSET_BUILD_PYTHON_BINDINGS=ON
cmake --build build -j
export PYTHONPATH=build/python:$PYTHONPATH
```

### massconsistent_amr
```bash
cd massconsistent_amr
cmake -S . -B build -DMASSCONSISTENT_BUILD_PYTHON_BINDINGS=ON
cmake --build build -j
export PYTHONPATH=build/python:$PYTHONPATH
```

## Technical Details

### Heat Release Model

Surface flux is computed from fire intensity and flame length:
```
surface_flux[i,j] = intensity[i,j] / flame_length[i,j]
```

Total heat release is approximated as:
```
total_heat_release = mean_intensity × perimeter
```

This is suitable for coupling with wind solvers that account for convective heating effects.

### Array Ordering

- Fire solver: NumPy arrays in shape (ny, nx) for 2D, (nz, ny, nx) for 3D
- Wind solver: 3D arrays expected in shape (nz, ny, nx)
- Fortran-order flattening used for interfacing with C++ backends

### Grid Compatibility

For coupling to work correctly:
- Same nx, ny for both solvers (required for first implementation)
- Same domain bounds (xmin, xmax, ymin, ymax)
- Can have different vertical grid (nz, zmin, zmax)

## Performance Notes

- One-way coupling: Minimal overhead, wind and fire can use different timesteps
- Two-way coupling: More expensive, wind and fire steps must be synchronized
- Wind update interval: Can be tuned to balance accuracy and speed
  - Every step: Most accurate but slowest
  - Every N steps: Faster, suitable if wind changes slowly

## Testing

Run tests:
```bash
pytest src/python/test_fire_solver_api.py -v
pytest src/python/test_pywildfire.py -v
```

Run regression tests:
```bash
cd build && ctest -L regtest --output-on-failure
```

## References

- [Problem Statement](../../README.md#two-way-wind-fire-coupling)
- [Fire Solver API](fire_solver_api.H)
- [Wind Solver Integration](https://github.com/hgopalan/massconsistent_amr)
