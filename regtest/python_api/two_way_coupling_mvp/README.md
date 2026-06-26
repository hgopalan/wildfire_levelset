# Two-Way Fire-Wind Coupling MVP

## Overview

This MVP (Minimum Viable Product) demonstrates two-way coupling between wildfire and wind solvers:

1. **Fire → Wind**: Wildfire solver extracts surface heat fluxes (intensity, flame length)
2. **Wind Computation**: Wind solver recomputes 3D wind field with heat-driven buoyancy forcing
3. **Wind → Fire**: Computed wind is fed back into wildfire solver

This creates a feedback loop where:
- Fire heat influences wind (buoyancy effects)
- Wind affects fire spread (ROS modification)
- Cycle repeats each timestep

## Architecture

```
┌─────────────────────────────────────────────────────────────┐
│                   TwoWayCoupledSimulation                    │
└─────────────────────────────────────────────────────────────┘
         ▲                              │
         │                              ▼
    ┌────┴────┐                  ┌────────────────┐
    │   Fire  │◄─── Fluxes ──────│SurfaceFluxEx..│
    │ Solver  │                  └────────────────┘
    └────┬────┘
         │
    Wind Updates
         │
         ▼
┌─────────────────────────────────────┐
│     WindSolverInterface             │
├─────────────────────────────────────┤
│ ✓ solve(time, heat_source)          │
│ ✓ get_velocity_arrays()             │
│ ✓ get_domain_info()                 │
└─────────────────────────────────────┘
         ▲
         │
┌────────┴────────────────────┐
│  SyntheticWindSolver        │
├─────────────────────────────┤
│ (for testing & demos)       │
│ - Responds to heat forcing  │
│ - Log-law wind profile      │
│ - Time-varying wind         │
└─────────────────────────────┘
```

## Key Modules

### 1. `surface_flux_extractor.py`

Extracts heat flux information from fire solver state for wind coupling.

**Key Methods:**
- `compute_column_heat_flux()` - Convert intensity to W/m²
- `compute_buoyancy_parameter()` - Convert heat to temperature anomaly
- `compute_fire_footprint()` - Identify active fire cells
- `get_flux_dict()` - Bundle all fluxes for wind solver

**Physical Model:**
- Intensity [kW/m] → heat release per perimeter length
- Flame front width assumed ~1.5 m
- Heat flux [W/m²] = intensity / width × 1000
- Buoyancy [K] = heat_flux / (ρ × Cp × w_scale)

**Example:**
```python
from surface_flux_extractor import SurfaceFluxExtractor

extractor = SurfaceFluxExtractor((ny, nx), dx, dy)
fluxes = extractor.get_flux_dict(fire_state)

# fluxes contains:
# - 'heat_flux': (ny, nx) [W/m²]
# - 'buoyancy': (ny, nx) [K]
# - 'fire_footprint': (ny, nx) binary mask
# - 'total_heat_MW': scalar [MW]
```

### 2. `wind_solver_interface.py`

Abstract interface for pluggable wind solvers.

**Abstract Class: `WindSolverInterface`**

Any wind solver should inherit and implement:

```python
class MyWindSolver(WindSolverInterface):
    def initialize(self, config_file):
        """Setup wind solver"""
        pass
    
    def solve(self, time, heat_source=None):
        """Compute wind with optional heat forcing"""
        pass
    
    def get_velocity_arrays(self):
        """Return (u_3d, v_3d, w_3d) arrays"""
        pass
    
    def get_domain_info(self):
        """Return domain metadata dictionary"""
        pass
    
    def finalize(self):
        """Cleanup"""
        pass
    
    def is_heat_aware(self):
        """Return True if solver supports heat_source parameter"""
        return True
```

**Included Implementation: `SyntheticWindSolver`**

For testing without external wind solver:
- Log-law wind profile matching atmospheric boundary layer
- Time-varying wind (diurnal cycle)
- Response to heat forcing (horizontal wind reduction, vertical velocity increase)

**Example:**
```python
from wind_solver_interface import SyntheticWindSolver

wind = SyntheticWindSolver()
wind.initialize()

# Without heat
wind.solve(time=0.0)
u, v, w = wind.get_velocity_arrays()

# With heat feedback
heat_source = {
    'heat_flux': intensity_to_flux(fire_state),
    'buoyancy': compute_buoyancy(fire_state),
    'fire_footprint': active_fire_mask,
}
wind.solve(time=0.0, heat_source=heat_source)
```

### 3. `two_way_coupling_example.py`

Complete MVP workflow demonstrating two-way coupling.

**Class: `TwoWayCoupledSimulation`**

Manages the coupled loop:

```python
from two_way_coupling_example import TwoWayCoupledSimulation

sim = TwoWayCoupledSimulation("fire_inputs.i")
sim.coupling_enabled = True
sim.heat_scaling = 1.0  # Feedback strength

# Run for 20 timesteps
sim.run(num_steps=20, plot_interval=5)

# Print summary
sim.print_summary()
sim.finalize()
```

**Workflow:**
1. Extract current fire state
2. Compute surface heat fluxes
3. Pass to wind solver with `heat_source` parameter
4. Get updated 3D wind field
5. Update fire wind via `update_wind_3d()`
6. Advance fire one timestep
7. Repeat

**Features:**
- Configurable heat scaling (sensitivity studies)
- Enable/disable coupling (compare 1-way vs 2-way)
- Progress logging with diagnostics
- Plotfile output at intervals

**Usage:**
```bash
# Run with full two-way coupling
python3 two_way_coupling_example.py inputs.i -n 20 -H 1.0

# Run with reduced feedback (sensitivity test)
python3 two_way_coupling_example.py inputs.i -n 20 -H 0.5

# Run without heat feedback (one-way coupling)
python3 two_way_coupling_example.py inputs.i -n 20 --no-coupling
```

## Running the MVP

### Prerequisites

1. Build wildfire_levelset with Python bindings:
   ```bash
   cmake -S . -B build -DLEVELSET_DIM_2D=ON -DLEVELSET_BUILD_PYTHON_BINDINGS=ON
   cmake --build build -j
   ```

2. Set PYTHONPATH:
   ```bash
   export PYTHONPATH=/path/to/wildfire_levelset/build/python:$PYTHONPATH
   ```

3. Ensure NumPy is installed:
   ```bash
   pip install numpy
   ```

### Run Complete Example

```bash
cd /path/to/wildfire_levelset

# Run two-way coupled simulation
PYTHONPATH=build/python python3 src/python/two_way_coupling_example.py \
    regtest/python_api/coupled_wind_fire/inputs.i -n 20

# Or with reduced feedback strength
PYTHONPATH=build/python python3 src/python/two_way_coupling_example.py \
    regtest/python_api/coupled_wind_fire/inputs.i -n 20 -H 0.5

# Or without heat feedback (one-way coupling for comparison)
PYTHONPATH=build/python python3 src/python/two_way_coupling_example.py \
    regtest/python_api/coupled_wind_fire/inputs.i -n 20 --no-coupling
```

### Run Regression Tests

```bash
# From build directory
cd build
ctest -R python_api_two_way_coupling_mvp -VV

# Or directly
cd ../regtest/python_api/two_way_coupling_mvp
PYTHONPATH=../../../build/python python3 test_two_way_coupling_mvp.py
```

## Physical Model Details

### Heat Flux Conversion

**From fire intensity to heat flux:**
```
Intensity [kW/m] = heat release per unit fire perimeter
Flame front width ~ 1-2 m (assume 1.5 m)
Heat flux [W/m²] = (Intensity [kW/m] / 1.5 m) × 1000
```

### Buoyancy Parameter

**From heat flux to potential temperature anomaly:**
```
Δθ [K] = Q [W/m²] / (ρ [kg/m³] × Cp [J/kg·K] × w_scale [m/s])

Typical values:
- ρ = 1.2 kg/m³
- Cp = 1005 J/kg·K
- w_scale = 2.0 m/s (characteristic updraft)
```

### Wind Modification

**Simple model (SyntheticWindSolver):**
```
u_modified = u_base × (1 - 0.15 × heat_normalized)
v_modified = v_base × (1 - 0.15 × heat_normalized)
w_modified = w_base + 0.2 × sin(πz/zmax) × heat_normalized
```

Interpretation: Strong heat reduces horizontal winds (~15%) and enhances vertical motion.

## Design Decisions

### Why This Architecture?

1. **Abstract Interface** - Enables integration of any wind solver
   - Real massconsistent_amr: swap `SyntheticWindSolver` for `MassConsistentWindSolver`
   - WRF, WindNinja, custom models: implement interface

2. **Flux Dictionary** - Clear data contract
   - All heat data bundled together
   - Easy to extend (add more fields later)
   - Wind solver decides what to use

3. **Modular Structure** - Each piece independently testable
   - Flux extraction works standalone
   - Wind solver tests work without fire solver
   - Coupling logic separated from implementation

4. **Synthetic Wind Solver** - Enables testing without external dependencies
   - Realistic wind profiles (log-law, time variation)
   - Responds to heat (demonstrates coupling)
   - Fast (no external process calls)

## Validation & Testing

### What the Tests Check

1. **Surface Flux Extraction**
   - Correct array shapes
   - Non-negative values
   - Reasonable magnitudes

2. **Wind Solver Interface**
   - Correct velocity array shapes
   - Deterministic behavior
   - Heat response (wind modification)

3. **Coupled Loop**
   - Fire and wind solvers interact correctly
   - 5-step simulation completes without error
   - Fire spreads as expected

4. **Heat Feedback Effects**
   - Heat reduces horizontal wind speed monotonically
   - Magnitude of reduction matches expectation (~15% per 100 kW/m²)

### How to Extend Testing

Add new test functions to `test_two_way_coupling_mvp.py`:

```python
def test_my_feature():
    """Test description"""
    print("Test: My Feature")
    # ... test code ...
    assert condition, "Error message"
    print("  ✓ Success")
    return True

# Add to main():
tests.append(("My Feature", test_my_feature))
```

## Known Limitations & Future Enhancements

### MVP Limitations

- **No spatial interpolation** - Fire and wind grids must align exactly
- **Simplified physics** - Linear wind reduction by heat (not coupled Navier-Stokes)
- **Synthetic wind solver only** - Real massconsistent_amr not yet integrated
- **No sub-cycling** - Fire and wind use same timestep

### Next Steps (Easy → Moderate)

1. **Real Wind Solver Integration** (Easy, ~4-8 hours)
   - Create `MassConsistentWindSolver` wrapping massconsistent_amr Python API
   - Verify grid alignment
   - Test data transfer

2. **Spatial Interpolation** (Easy, ~4-6 hours)
   - Implement nearest-neighbor and bilinear interpolation
   - Handle off-grid cells
   - Test conservation properties

3. **Sub-Cycling Framework** (Moderate, ~8-12 hours)
   - Allow wind timestep ≠ fire timestep
   - Time interpolation of wind between solves
   - Adaptive subcycling based on wind/fire CFL

4. **Advanced Physics** (Moderate, ~12-16 hours)
   - Buoyancy-driven updrafts (more realistic vertical motion)
   - Sensible/latent heat separation
   - Radiation effects on ambient temperature

## Code Examples

### Example 1: Compare 1-Way vs 2-Way Coupling

```python
from two_way_coupling_example import TwoWayCoupledSimulation

inputs = "regtest/python_api/coupled_wind_fire/inputs.i"

print("=== ONE-WAY COUPLING (no heat feedback) ===")
sim1 = TwoWayCoupledSimulation(inputs)
sim1.coupling_enabled = False
sim1.run(num_steps=20)
sim1.print_summary()
burned_1way = sim1.diagnostics['burned_area'][-1]

print("\n=== TWO-WAY COUPLING (full heat feedback) ===")
sim2 = TwoWayCoupledSimulation(inputs)
sim2.coupling_enabled = True
sim2.heat_scaling = 1.0
sim2.run(num_steps=20)
sim2.print_summary()
burned_2way = sim2.diagnostics['burned_area'][-1]

print(f"\nFire spread difference: "
      f"{(burned_2way - burned_1way)/burned_1way * 100:.1f}%")
```

### Example 2: Sensitivity Analysis

```python
for heat_scale in [0.0, 0.5, 1.0, 1.5]:
    sim = TwoWayCoupledSimulation(inputs)
    sim.heat_scaling = heat_scale
    sim.run(num_steps=20)
    
    final_burned = sim.diagnostics['burned_area'][-1]
    peak_heat = max(sim.diagnostics['total_heat_MW'])
    
    print(f"Heat scale {heat_scale:.1f}: "
          f"burned={final_burned:.0f} m², "
          f"peak_heat={peak_heat:.2f} MW")
```

## References

### Physical Models

- **Log-law wind profile**: Monin-Obukhov similarity theory
- **Byram intensity**: Byram (1959), Thomas (1963)
- **Buoyancy parameter**: Lilly (1962), Parsons (1989)
- **Wildfire physics**: Rothermel (1972), Anderson (1983)

### Related Work

- massconsistent_amr: https://github.com/hgopalan/massconsistent_amr
- wildfire_levelset: https://github.com/hgopalan/wildfire_levelset
- QUIC-URB wind model: Röckle (1990), Pardyjak & Brown (2001)

## Contributing

To extend this MVP:

1. **New wind solver**: Inherit `WindSolverInterface`, implement methods
2. **New flux types**: Add methods to `SurfaceFluxExtractor`
3. **Better physics**: Modify `_apply_heat_forcing()` in wind solver
4. **Tests**: Add test functions to `test_two_way_coupling_mvp.py`

## License

Same as wildfire_levelset (see LICENSE file)

---

**Author**: Generated as MVP for two-way fire-wind coupling
**Date**: 2026-06-26
**Status**: Ready for integration and real wind solver coupling
