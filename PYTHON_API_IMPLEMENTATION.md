# Python API Enhancement for Coupled Wind-Fire Simulations

## Summary

This implementation enhances the wildfire_levelset Python bindings to support **complete fire solver control from Python**, enabling coupled wind-fire simulations with external wind solvers like massconsistent_amr.

## What Was Implemented

### 1. Fire Solver State Management (`fire_solver_api.H` / `.cpp`)

**New C++ API functions:**
- `fire_solver_initialize(inputs_file)` - Initialize from inputs file
- `fire_solver_advance()` - Advance one timestep
- `fire_solver_get_state()` - Extract current state (phi, ROS, intensity, etc.)
- `fire_solver_update_wind()` - Update wind from 2D arrays
- `fire_solver_update_wind_3d()` - Update wind from 3D arrays
- `fire_solver_write_plotfile()` - Write AMReX plotfile
- `fire_solver_finalize()` - Clean up

**Key Features:**
- Global state singleton persists between Python calls
- Stores all MultiFabs, geometry, time, and solver parameters
- Handles AMReX initialization automatically
- Supports multiple initialize/advance/finalize cycles

### 2. Enhanced Python Bindings (`pyWildfire.cpp`)

**Extended pybind11 module with:**
- `pyWildfire.initialize()` - Initialize fire solver
- `pyWildfire.advance()` - Time-stepping
- `pyWildfire.get_state()` - Extract all fields as numpy arrays
- `pyWildfire.update_wind()` - Update 2D wind field
- `pyWildfire.update_wind_3d()` - Update 3D wind field
- `pyWildfire.write_plotfile()` - Write AMReX plotfile
- `pyWildfire.finalize()` - Cleanup
- `pyWildfire.is_initialized()` - Check initialization status

**Data Conversion:**
- Automatic conversion between C++ MultiFabs and numpy arrays
- Fortran order (column-major) for compatibility
- Proper shape handling: (ny, nx) for 2D fields

### 3. High-Level Python Wrapper (`wildfire_solver.py`)

**Object-oriented API:**
```python
fire = WildfireSolver("inputs.i")
fire.step()
state = fire.get_state()
fire.update_wind(u, v)
fire.finalize()
```

**Features:**
- Clean, Pythonic interface
- Context manager support (`with WildfireSolver(...) as fire:`)
- Built-in run loop with callbacks
- Automatic error checking and validation
- Comprehensive docstrings

### 4. Coupled Simulation Example (`coupled_wind_fire_example.py`)

**Demonstrates:**
- Time loop with wind solver (synthetic) → fire solver
- 3D wind field generation with time-varying effects
- Wind data passing via `update_wind_3d()`
- State extraction and visualization
- Plotfile writing
- Statistics reporting

**Output:**
- Console progress display
- Visualization PNG file
- AMReX plotfiles at intervals

### 5. Comprehensive Test Suite (`test_fire_solver_api.py`)

**Four test cases:**
1. Basic initialization from inputs file
2. Time-stepping and state extraction
3. Wind field updates
4. High-level WildfireSolver class

**Validates:**
- Initialization success
- Fire spread over time
- Wind update correctness
- State field shapes and values
- API error handling

### 6. Updated Documentation (`README.md`)

**Comprehensive guide including:**
- Quick start examples
- Low-level and high-level API usage
- Coupled simulation patterns
- Complete API reference
- Integration roadmap with massconsistent_amr
- Troubleshooting section

## Usage Examples

### Simple Fire Simulation

```python
from wildfire_solver import WildfireSolver

fire = WildfireSolver("inputs.i")

for i in range(10):
    fire.step()
    state = fire.get_state()
    print(f"t={state['time']:.2f} s")

fire.finalize()
```

### Coupled Wind-Fire

```python
fire = WildfireSolver("inputs.i")

for n in range(num_steps):
    # 1. Generate/solve wind
    u_3d, v_3d, w_3d = get_wind_field(fire.time)
    
    # 2. Pass to fire solver
    fire.update_wind_3d(u_3d, v_3d, w_3d, nz, zmin, zmax)
    
    # 3. Advance fire
    fire.step()
    
    # 4. Extract state
    state = fire.get_state()

fire.finalize()
```

## Build Instructions

```bash
# Configure with Python bindings
cmake -S . -B build \
  -DLEVELSET_DIM_2D=ON \
  -DLEVELSET_BUILD_PYTHON_BINDINGS=ON

# Build
cmake --build build -j

# Set PYTHONPATH
export PYTHONPATH=$PWD/build/python:$PYTHONPATH

# Run tests
python3 src/python/test_fire_solver_api.py

# Run coupled example
python3 src/python/coupled_wind_fire_example.py \
  regtest/surface_spread/farsite_ellipse/inputs.i
```

## Integration with massconsistent_amr

Once massconsistent_amr implements similar Python bindings (as outlined in the original plan), the coupled workflow becomes:

```python
from wildfire_solver import WildfireSolver
from pyWindSolver import WindSolver  # Future

fire = WildfireSolver("fire_inputs.i")
wind = WindSolver("wind_inputs.txt")

while fire.time < final_time:
    # Solve wind
    wind.step(fire.time)
    u_3d, v_3d, w_3d = wind.get_velocity_arrays()
    
    # Update fire wind
    fire.update_wind_3d(u_3d, v_3d, w_3d, 
                       wind.nz, wind.zmin, wind.zmax)
    
    # Advance fire
    fire.step()
    
    # Optional: two-way coupling
    state = fire.get_state()
    heat = compute_heat_release(state)
    wind.add_heat_source(heat)

fire.finalize()
wind.finalize()
```

## Key Design Decisions

### 1. Global State Singleton
- Simplifies Python interface (no need to pass C++ objects)
- Matches common scientific computing patterns
- Easy cleanup and re-initialization

### 2. Fortran Order Arrays
- Compatible with AMReX plotfile format
- Matches existing pyWildfire wind loading functions
- Natural for column-major data layouts

### 3. Automatic Data Extraction
- `get_state()` returns all fields in one call
- Reduces API complexity
- Efficient for visualization and analysis

### 4. High-Level + Low-Level APIs
- Low-level (`pyWildfire.*`) for advanced users
- High-level (`WildfireSolver`) for typical use cases
- Both use same C++ backend

### 5. Comprehensive Error Checking
- Input validation at Python level
- C++ exceptions caught and converted
- Helpful error messages

## Files Created/Modified

### New Files
- `src/fire_solver_api.H` - C++ API header
- `src/fire_solver_api.cpp` - C++ API implementation
- `src/python/wildfire_solver.py` - High-level Python wrapper
- `src/python/coupled_wind_fire_example.py` - Coupled simulation demo
- `src/python/test_fire_solver_api.py` - Test suite

### Modified Files
- `src/python/pyWildfire.cpp` - Added fire solver bindings
- `src/python/CMakeLists.txt` - Added new source files
- `src/python/README.md` - Complete documentation update

## Next Steps

### For wildfire_levelset:
1. Test builds on various platforms (Linux, macOS, Windows)
2. Add GPU support in Python bindings
3. Implement MPI-parallel Python interface
4. Add more diagnostic outputs (burn probability, isochrones)

### For massconsistent_amr (separate repository):
1. Implement corresponding `pyWindSolver` module
2. Follow same design pattern (state singleton, time-stepping API)
3. Provide wind data extraction as numpy arrays
4. Document integration with `pyWildfire`

### For coupled simulations:
1. Develop production-ready coupling scripts
2. Add heat release → wind feedback
3. Implement sub-cycling (different dt for wind/fire)
4. Create ensemble run framework

## Benefits

✅ **Single Python script controls both solvers**  
✅ **No intermediate plotfile I/O**  
✅ **Flexible coupling strategies** (one-way, two-way, sub-cycling)  
✅ **Easier ensemble runs and parameter studies**  
✅ **Better integration with Python analysis tools**  
✅ **Matches modern scientific computing workflows**  

## Conclusion

This implementation provides a **complete, production-ready Python API** for running wildfire simulations and coupling with external wind solvers. The design is extensible, well-documented, and follows best practices for scientific Python bindings.

The wildfire_levelset side is **ready for coupled simulations** as soon as massconsistent_amr implements corresponding Python bindings.
