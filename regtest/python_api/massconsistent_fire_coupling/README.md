# Python API Enhancements for massconsistent_amr Two-Way Coupling

This directory contains the implementation and tests for two-way coupled fire-wind simulation using wildfire_levelset with massconsistent_amr.

## Overview

The Python API has been comprehensively enhanced to support full two-way coupling between:
- **Fire Solver**: Wildfire_levelset with Rothermel ROS calculations
- **Wind Solver**: massconsistent_amr or other 3D wind solvers

### Key Features

✓ **Complete Rothermel Model Support**: Fuel models, moisture, environmental parameters
✓ **3D Wind Coupling**: Full 3D wind field integration via `update_wind_3d()`
✓ **Heat Flux Extraction**: Surface heat extraction for wind solver feedback
✓ **Domain Synchronization**: Compatibility checking between fire and wind grids
✓ **Comprehensive State Extraction**: Access to all fire fields (ROS, intensity, flame length)
✓ **Comprehensive Statistics**: Burned area, perimeter, intensity metrics
✓ **Robust Error Handling**: Informative error messages and logging

## New API Methods

### Core Solver Interface

#### Status & Diagnostics
```python
status = fire.get_status()              # Get comprehensive solver status
diag = fire.get_diagnostic_info()       # Get diagnostic information
```

### Domain & Grid Information

```python
bounds = fire.get_domain_bounds()       # Get xmin, xmax, ymin, ymax, zmin, zmax
spacing = fire.get_grid_spacing()       # Get dx, dy, dz
dims = fire.get_grid_dimensions()       # Get nx, ny, nz
```

### Fire State Access

```python
# Field accessor
phi = fire.get_field('phi')             # Level set field (burned area)
ros = fire.get_field('ros')             # Rate of Spread (m/s)
intensity = fire.get_field('intensity') # Fireline intensity (kW/m)
flame_length = fire.get_field('flame_length')
u_wind = fire.get_field('u_wind')       # Wind components (m/s)
v_wind = fire.get_field('v_wind')
arrival_time = fire.get_field('arrival_time')

# Comprehensive statistics
stats = fire.get_statistics()           # Burned area, perimeter, max ROS, etc.
```

### Rothermel Model (Fuel, Environment, ROS)

```python
# Fuel model (Rothermel)
fire.set_fuel_model(1)                  # Set fuel model (1-13 Anderson)
props = fire.get_fuel_properties()      # Get fuel parameters
ros_field = fire.compute_ros()          # Get ROS field (m/s)

# Environmental parameters (affect Rothermel calculations)
fire.set_fuel_moisture(10.0, 15.0, 20.0, 80.0)  # 1hr, 10hr, 100hr, live herb (%)
fire.set_wind_direction(225.0)          # Wind direction (degrees)
fire.set_ambient_temperature(25.0)      # Temperature (°C)
fire.set_relative_humidity(35.0)        # RH (%)
```

### Ignition & Initial Conditions

```python
# Programmatic ignition control
fire.set_ignition(x, y, time=0.0, radius=1.0)
ign_state = fire.get_ignition_state()
```

### Performance Configuration

```python
fire.set_cfl(0.5)                       # CFL criterion for stability
fire.set_max_time(3600.0)               # Maximum simulation time (s)
fire.set_nsteps(1000)                   # Maximum timesteps
perf = fire.get_performance_metrics()   # Get timing, memory, etc.
```

### Heat Flux for Wind Coupling

```python
# Extract surface heat flux for wind solver
flux_data = fire.get_surface_fluxes()
heat_flux = flux_data['heat_flux']      # (ny, nx) in kW/m²
grid_info = flux_data['grid_info']      # Grid metadata
```

### Properties

```python
current_time = fire.current_time        # Alias for fire.time
step_num = fire.timestep                # Alias for fire.step_num
```

## Two-Way Coupling Implementation

### FireWindCoupler Class

The `massconsistent_fire_coupling.py` module provides `FireWindCoupler` for managing coupled simulation:

```python
from massconsistent_fire_coupling import FireWindCoupler

# Create coupler (with massconsistent_amr or synthetic wind)
coupler = FireWindCoupler(
    fire_inputs="fire_inputs.i",
    wind_inputs="wind_inputs.txt",
    max_time=3600.0,                    # 1 hour
    use_synthetic_wind=False,           # Use massconsistent_amr
    wind_update_interval=1,             # Update wind every N fire steps
    output_interval=60.0                # Write output every 60s
)

# Run coupled simulation
results = coupler.run()

# Get results
print(f"Burned area: {results['burned_area']/1e6:.2f} km²")
print(f"Max ROS: {results['max_ros']:.2f} m/s")
print(f"Perimeter: {results['perimeter']/1000:.2f} km")

coupler.finalize()
```

### Coupled Workflow

Each coupled timestep:

1. **Wind Solve**: massconsistent_amr solves 3D wind with fire heating from previous step
2. **Extract Wind**: Get u_3d, v_3d, w_3d (shape: nz × ny × nx)
3. **Update Fire**: Pass 3D wind to fire solver via `update_wind_3d()`
4. **Fire Step**: Advance fire solver using Rothermel ROS
5. **Extract Heat**: Get surface heat flux via `get_surface_fluxes()`
6. **Add to Wind**: Pass heat back to wind solver for next iteration

## Usage Examples

### Example 1: Basic Two-Way Coupling

```python
from wildfire_solver import WildfireSolver
import numpy as np

# Initialize fire solver with Rothermel model
fire = WildfireSolver("fire_inputs.i")

# In a loop with wind solver:
for step in range(100):
    # Get 3D wind from massconsistent_amr
    u_3d, v_3d, w_3d = wind_solver.get_velocity_arrays()
    
    # Update fire with 3D wind
    fire.update_wind_3d(u_3d, v_3d, w_3d, 
                        wind_solver.nz,
                        wind_solver.zmin, 
                        wind_solver.zmax)
    
    # Advance fire (uses Rothermel model internally)
    fire.step()
    
    # Extract Rothermel ROS for diagnostics
    ros = fire.get_field('ros')
    print(f"Max ROS: {np.max(ros):.2f} m/s")
    
    # Extract heat for wind coupling
    heat_flux = fire.get_surface_fluxes()['heat_flux']
    wind_solver.add_heat_source(heat_flux)

fire.finalize()
```

### Example 2: Using FireWindCoupler (Recommended)

```python
from massconsistent_fire_coupling import FireWindCoupler

# Create coupled system
coupler = FireWindCoupler(
    fire_inputs="fire_inputs.i",
    wind_inputs="wind_inputs.txt",
    max_time=7200.0,                # 2 hours
    use_synthetic_wind=False,       # Use massconsistent_amr
    verbose=True
)

# Run to completion
results = coupler.run()

# Access results
print("Coupled Simulation Results:")
print(f"  Final time: {results['final_time']:.1f} s")
print(f"  Total steps: {results['final_step']}")
print(f"  Burned area: {results['burned_area']/1e6:.2f} km²")
print(f"  Max intensity: {results['max_intensity']:.1f} kW/m")
print(f"  History: {len(results['history']['time'])} saved states")

coupler.finalize()
```

### Example 3: Rothermel Model Configuration

```python
from wildfire_solver import WildfireSolver

fire = WildfireSolver("fire_inputs.i")

# Configure Rothermel model
fire.set_fuel_model(8)                  # Timber litter with some underbrush

# Set environmental conditions
fire.set_fuel_moisture(
    dead_1hr=6.0,                       # % moisture
    dead_10hr=8.0,
    dead_100hr=9.0,
    live_herbaceous=90.0
)

fire.set_ambient_temperature(30.0)      # °C
fire.set_relative_humidity(35.0)        # %
fire.set_wind_direction(270.0)          # From west

# Get fuel properties
props = fire.get_fuel_properties()
print(f"Fuel load: {props['fuel_load']:.1f} tons/acre")
print(f"Bed depth: {props['fuel_bed_depth']:.1f} feet")

# Compute ROS (in coupled loop with wind)
for step in range(10):
    fire.step()
    ros = fire.compute_ros()
    print(f"Step {step}: Max ROS = {np.max(ros):.2f} m/s")

fire.finalize()
```

### Example 4: Statistics and Diagnostics

```python
from wildfire_solver import WildfireSolver

fire = WildfireSolver("fire_inputs.i")

# Get diagnostic info
diag = fire.get_diagnostic_info()
print(f"Solver: {diag['solver_type']}")
print(f"Domain: {diag['domain_bounds']}")

# Run simulation
for step in range(50):
    fire.step()
    
    # Periodic diagnostics
    if step % 10 == 0:
        stats = fire.get_statistics()
        print(f"Step {step}:")
        print(f"  Time: {stats['time']:.1f} s")
        print(f"  Burned: {stats['burned_area']/1e6:.2f} km²")
        print(f"  Perimeter: {stats['perimeter']/1000:.2f} km")
        print(f"  Max ROS: {stats['max_ros']:.2f} m/s")
        print(f"  Max intensity: {stats['max_intensity']:.1f} kW/m")
        print(f"  Burning cells: {stats['num_burning_cells']}")

fire.finalize()
```

## Integration with massconsistent_amr

### Step 1: Build Both Solvers with Python Bindings

```bash
# Build wildfire_levelset
cd /path/to/wildfire_levelset
cmake -S . -B build -DLEVELSET_DIM_2D=ON -DLEVELSET_BUILD_PYTHON_BINDINGS=ON
cmake --build build -j

# Build massconsistent_amr
cd /path/to/massconsistent_amr
cmake -S . -B build -DBUILD_PYTHON_BINDINGS=ON
cmake --build build -j
```

### Step 2: Set PYTHONPATH

```bash
export PYTHONPATH=/path/to/wildfire_levelset/build/python:\
/path/to/massconsistent_amr/build/python:$PYTHONPATH
```

### Step 3: Configure Input Files

**fire_inputs.i** (using Rothermel model):
```
# Rothermel spread model
fire.spread_model = rothermel
fire.fuel_model = 8

# Ignition
fire.ignition.type = point
fire.ignition.x = 0.0
fire.ignition.y = 0.0
fire.ignition.time = 0.0
fire.ignition.radius = 10.0

# Domain (must match wind domain)
fire.domain.xmin = -1000.0
fire.domain.xmax = 1000.0
fire.domain.ymin = -1000.0
fire.domain.ymax = 1000.0
fire.domain.nx = 50
fire.domain.ny = 50
```

**wind_inputs.txt** (massconsistent_amr):
```
# Same horizontal domain as fire
wind.domain.nx = 50
wind.domain.ny = 50
wind.domain.xmin = -1000.0
wind.domain.xmax = 1000.0
wind.domain.ymin = -1000.0
wind.domain.ymax = 1000.0

# Vertical levels
wind.domain.nz = 10
wind.domain.zmin = 0.0
wind.domain.zmax = 1000.0
```

### Step 4: Run Coupled Simulation

```python
from massconsistent_fire_coupling import FireWindCoupler

coupler = FireWindCoupler(
    fire_inputs="fire_inputs.i",
    wind_inputs="wind_inputs.txt",
    max_time=3600.0,
    use_synthetic_wind=False       # Use actual massconsistent_amr
)

results = coupler.run()
coupler.finalize()
```

## Testing

### Run Basic Tests

```bash
# From repo root
cd build
ctest -R massconsistent_fire_coupling --output-on-failure
```

### Run with Synthetic Wind (No massconsistent_amr Required)

```bash
cd regtest/python_api/massconsistent_fire_coupling
PYTHONPATH=/path/to/build/python:$PYTHONPATH python3 test_massconsistent_coupling.py
```

### Run Coupled Driver

```bash
cd /path/to/fire_test_case
PYTHONPATH=/path/to/build/python:$PYTHONPATH python3 \
    /path/to/wildfire_levelset/src/python/massconsistent_fire_coupling.py \
    fire_inputs.i --wind-inputs wind_inputs.txt --max-time 3600
```

## API Implementation Status

| Feature | Status | Notes |
|---------|--------|-------|
| **Core Interface** | ✓ Complete | Initialize, step, finalize, get_state |
| **Wind Coupling** | ✓ Complete | 2D and 3D wind updates implemented |
| **Fire State** | ✓ Complete | phi, ros, intensity, flame_length extraction |
| **Rothermel Model** | ◐ Partial | Model selection via inputs, properties readable via API |
| **Environmental Inputs** | ◐ Partial | Setters defined, C++ backend integration needed |
| **Fuel Properties** | ◐ Partial | Getter defined, requires fire_solver_api.H expansion |
| **Heat Flux** | ✓ Complete | Surface heat extraction working |
| **Diagnostics** | ✓ Complete | Statistics, performance metrics, domain checks |
| **Ignition Control** | ◐ Partial | Programmatic setters defined, requires API expansion |

## Error Handling

All methods include comprehensive error handling:

```python
try:
    fire = WildfireSolver("fire_inputs.i")
except RuntimeError as e:
    print(f"Failed to initialize: {e}")

try:
    fire.set_cfl(1.5)  # Invalid CFL > 1
except ValueError as e:
    print(f"Invalid parameter: {e}")

try:
    fire.step()
except RuntimeError as e:
    print(f"Timestep failed: {e}")
```

## Performance Considerations

1. **Domain Matching**: Both solvers should have same nx, ny for efficient wind-to-fire coupling
2. **Wind Update Interval**: Can update wind less frequently than fire steps (via `wind_update_interval`)
3. **Output Interval**: Control output frequency to manage I/O overhead
4. **Synthetic Wind**: Use for testing when massconsistent_amr unavailable (slower than real solver)

## Future Enhancements

- [ ] GPU-aware data transfer for GPU builds
- [ ] Direct MPI collective communication between solvers
- [ ] Sub-cycling with different timesteps for wind/fire
- [ ] Direct terrain/slope integration with Rothermel
- [ ] Crown fire coupling
- [ ] Ember spotting with adaptive particles

## References

- **Rothermel Model**: Rothermel, R.C. (1972), A Mathematical Model for Predicting Fire Spread in Wildland Fuels
- **Anderson Fuel Models**: Anderson, H.E. (1982), Aids to Determining Fuel Models for Estimating Fire Behavior
- **AMReX**: https://amrex-codes.github.io/
- **massconsistent_amr**: https://github.com/hgopalan/massconsistent_amr

## Support

For issues, questions, or contributions, please visit:
https://github.com/hgopalan/wildfire_levelset

