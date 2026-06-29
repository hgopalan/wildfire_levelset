# Wind-Fire Coupling Guide

This guide explains how to use the two-way coupling capability in wildfire_levelset to couple with external wind solvers like massconsistent_amr.

## Overview

Two-way coupling enables fire simulations where:
1. **Wind affects fire**: The fire spread responds to wind dynamics
2. **Fire affects wind**: The fire heating creates updrafts and wind changes

This is critical for realistic fire-atmosphere interaction studies.

## Quick Start

### Option 1: Using levelset_coupling Module (Recommended)

The easiest way to set up wind-fire coupling is using the dedicated `levelset_coupling` module:

```python
from levelset_coupling import CoupledWindFireSimulation

# Initialize with two-way coupling
coupled = CoupledWindFireSimulation(
    wind_inputs="wind_inputs.i",
    fire_inputs="fire_inputs.i",
    coupling_mode='two_way'  # Use 'one_way' for simpler coupling
)

# Run simulation
results = coupled.run(
    final_time=3600.0,        # Run for 1 hour
    plot_interval=600.0       # Write plots every 10 minutes
)

# Cleanup
coupled.finalize()

print(f"Final time: {results['final_time']:.1f} s")
print(f"Final state: {results['fire_state']}")
```

### Option 2: Manual Control (More Flexibility)

For finer control over the coupling process:

```python
from wildfire_solver import WildfireSolver
from wind_solver import WindSolver  # from massconsistent_amr

# Initialize solvers
fire = WildfireSolver("fire_inputs.i")
wind = WindSolver("wind_inputs.i")

# Time loop
final_time = 3600.0
while fire.time < final_time:
    # Step 1: Solve wind (includes heat source from previous iteration)
    wind.solve()
    
    # Step 2: Extract wind velocity
    vel = wind.get_velocity()
    u_3d = vel['u']  # shape (nz, ny, nx)
    v_3d = vel['v']
    w_3d = vel['w']
    
    # Step 3: Update fire with wind
    fire.update_wind_3d(u_3d, v_3d, w_3d, wind.nz, wind.zmin, wind.zmax)
    
    # Step 4: Advance fire
    fire.step()
    
    # Step 5: Extract heat release (for next wind solve)
    heat_data = fire.get_surface_fluxes()
    
    # Step 6: Add heat to wind solver (for next iteration)
    grid_info = {
        'xmin': fire.xmin,
        'xmax': fire.xmax,
        'ymin': fire.ymin,
        'ymax': fire.ymax,
        'dx': fire.dx,
        'dy': fire.dy
    }
    wind.add_heat_source(heat_data['heat_flux'], grid_info)

# Cleanup
fire.finalize()
wind.finalize()
```

## Coupling Modes

### One-Way Coupling (wind → fire)

Wind affects fire, but fire does NOT affect wind.

**When to use:**
- Pre-computed wind fields
- Coupling with weather models
- Fast simulations
- Fire effects on atmosphere are negligible

**Example:**
```python
coupled = CoupledWindFireSimulation(
    wind_inputs="wind_inputs.i",
    fire_inputs="fire_inputs.i",
    coupling_mode='one_way'
)
results = coupled.run(final_time=3600.0)
```

### Two-Way Coupling (wind ↔ fire)

Fire affects wind through heat release; wind affects fire through dynamics.

**When to use:**
- Fire-atmosphere interaction studies
- Fire-induced wind effects important
- Research simulations
- Accurate fire behavior with atmosphere feedback

**Example:**
```python
coupled = CoupledWindFireSimulation(
    wind_inputs="wind_inputs.i",
    fire_inputs="fire_inputs.i",
    coupling_mode='two_way'  # Enable heat feedback
)
results = coupled.run(final_time=3600.0)
```

## Integration with massconsistent_amr

### Prerequisites

1. **massconsistent_amr** built with Python bindings:
   ```bash
   cmake -S . -B build -DMASSCONSISTENT_BUILD_PYTHON_BINDINGS=ON
   cmake --build build -j
   ```

2. **wildfire_levelset** built with Python bindings:
   ```bash
   cmake -S . -B build -DLEVELSET_BUILD_PYTHON_BINDINGS=ON
   cmake --build build -j
   ```

3. **Python path** includes both:
   ```bash
   export PYTHONPATH=$MASSCONSISTENT_BUILD/python:$WILDFIRE_BUILD/python:$PYTHONPATH
   ```

### Setting Up Coupled Simulation

1. **Prepare input files**:
   - `wind_inputs.i` - massconsistent_amr configuration
   - `fire_inputs.i` - wildfire_levelset configuration
   - Both should have compatible domains

2. **Check domain compatibility**:
   ```python
   from levelset_coupling import CoupledWindFireSimulation
   
   coupled = CoupledWindFireSimulation(
       wind_inputs="wind_inputs.i",
       fire_inputs="fire_inputs.i",
       coupling_mode='two_way'
   )
   
   # Prints compatibility info automatically
   # Logs warnings if domains don't match
   ```

3. **Run simulation**:
   ```python
   results = coupled.run(final_time=3600.0)
   coupled.finalize()
   ```

## Heat Flux Interface

### Extracting Heat Flux from Fire

```python
from wildfire_solver import WildfireSolver

fire = WildfireSolver("fire_inputs.i")

# Simulate for a while...
for step in range(100):
    fire.step()

# Extract surface heat flux
heat_data = fire.get_surface_fluxes()
heat_flux = heat_data['heat_flux']  # shape (ny, nx), units: kW/m²

# Grid information
grid_info = heat_data['grid_info']
print(f"Grid: {grid_info['nx']} × {grid_info['ny']}")
print(f"Cell size: {grid_info['dx']:.2f} × {grid_info['dy']:.2f} m")
```

### Adding Heat to Wind Solver

```python
from wind_solver import WindSolver

wind = WindSolver("wind_inputs.i")

# After extracting heat from fire
grid_info = {
    'xmin': fire.xmin,
    'xmax': fire.xmax,
    'ymin': fire.ymin,
    'ymax': fire.ymax,
    'dx': fire.dx,
    'dy': fire.dy,
    'scaling_factor': 1.0,      # Optional: scale heat effect
    'temporal_decay': 1.0        # Optional: decay factor (0-1)
}

wind.add_heat_source(heat_flux, grid_info)

# Heat will be used in the next wind.solve() call
wind.solve()
```

## Example Scripts

### Complete One-Way Example

See `example_one_way_coupling.py`:
```bash
PYTHONPATH=build/python python3 src/python/example_one_way_coupling.py \
    regtest/surface_spread/farsite_ellipse/inputs.i
```

### Complete Two-Way Example

See `example_two_way_coupling.py`:
```bash
PYTHONPATH=build/python:$MASSCONSISTENT_PYTHONPATH \
    python3 src/python/example_two_way_coupling.py \
    wind_inputs.i fire_inputs.i
```

### Using Callbacks for Custom Processing

```python
def my_callback(step, result):
    """Called after each coupled step."""
    fire_state = result['fire_state']
    print(f"Step {step}: Fire time = {fire_state['time']:.1f} s")
    print(f"  Burned area = {np.sum(fire_state['phi'] <= 0)} cells")

coupled = CoupledWindFireSimulation(
    wind_inputs="wind_inputs.i",
    fire_inputs="fire_inputs.i",
    coupling_mode='two_way'
)

results = coupled.run(
    final_time=3600.0,
    callback=my_callback
)
```

## Domain Compatibility

For successful coupling, fire and wind domains should be compatible:

1. **Horizontal domain** should match (within tolerance):
   ```
   fire.xmin ≈ wind.xmin
   fire.xmax ≈ wind.xmax
   fire.ymin ≈ wind.ymin
   fire.ymax ≈ wind.ymax
   ```

2. **Grid spacing** should match:
   ```
   fire.dx ≈ wind.dx
   fire.dy ≈ wind.dy
   ```

If domains don't match exactly, the coupling will still work but may produce warnings.

## Troubleshooting

### Import Errors

If you get `ImportError: No module named 'wind_solver'`:
1. Make sure massconsistent_amr is built with Python bindings
2. Set PYTHONPATH correctly
3. Check that `pyWindSolver.*.so` exists

### Array Size Mismatches

If you get size mismatch errors:
1. Verify grids are compatible
2. Check array shapes before passing
3. Ensure Fortran-order flattening if needed

### Coupling Not Working

If two-way coupling isn't working:
1. Check that wind solver has `add_heat_source()` method
2. Verify domains are compatible
3. Check that heat_flux is not empty/NaN

### Performance Issues

For faster coupling:
1. Use one-way coupling if possible
2. Reduce wind update frequency: `wind_update_interval=N`
3. Use smaller domains
4. Enable GPU acceleration if available

## Advanced Topics

### Sub-Cycling

Run fire multiple steps for each wind solve:

```python
coupled = CoupledWindFireSimulation(...)

results = coupled.run(
    final_time=3600.0,
    wind_update_interval=10  # Update wind every 10 fire steps
)
```

### Ensemble Simulations

Run multiple scenarios with parameter variations:

```python
for wind_speed in [5, 10, 15]:
    coupled = CoupledWindFireSimulation(...)
    # Modify parameters as needed
    results = coupled.run(final_time=3600.0)
    # Analyze results
```

### Adaptive Refinement

Refine grids near fire front:

```python
# This is a future capability
# Currently requires manual implementation
```

## Performance Benchmarks

Typical performance (on modern CPU):

- **One-way coupling**: ~100-200 fire steps/minute
- **Two-way coupling**: ~10-50 wind-fire iterations/minute
  (depends on wind solver complexity)

## References

- **massconsistent_amr**: https://github.com/hgopalan/massconsistent_amr
- **wildfire_levelset**: https://github.com/hgopalan/wildfire_levelset
- **AMReX**: https://amrex-codes.github.io/amrex/

## License

Same as wildfire_levelset (see LICENSE file)
