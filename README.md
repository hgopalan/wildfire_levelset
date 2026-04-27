# wildfire_levelset

An AMReX-based C++ level-set solver for wildfire front propagation modeling with FARSITE elliptical expansion, Rothermel fire spread equations, and terrain effects.

## Documentation

📚 **[Read the full documentation](https://hgopalan.github.io/wildfire_levelset/)**

The documentation includes:
- Overview and mathematical models
- Detailed equations for Rothermel, FARSITE, and terrain effects
- Code structure and API reference
- Build and usage instructions
- Input parameter reference

## Features

- **Level-set advection** for fire front tracking
- **Rothermel fire spread model** with NFFL fuel database (FM1-FM13)
- **FARSITE elliptical expansion** (Richards 1990) with Anderson L/W ratio
- **Terrain effects** including slope and aspect corrections
- **Firebrand spotting** with probability-based stochastic model
- **Crown fire initiation** (Van Wagner 1977)
- **Bulk fuel consumption** modeling
- **2D and 3D** simulation capabilities
- **Embedded Boundary (EB)** support for complex geometries
- **AMReX-based** with adaptive mesh refinement

## Prerequisites

- C++17 compiler
- CMake (3.20+)
- Git

## Quick Start

### Clone

```bash
git clone --recurse-submodules https://github.com/hgopalan/wildfire_levelset.git
cd wildfire_levelset
```

If you already cloned without submodules:

```bash
git submodule update --init --recursive
```

### Build

Default 3D build:
```bash
cmake -S . -B build
cmake --build build -j
```

For 2D:
```bash
cmake -S . -B build -DLEVELSET_DIM_2D=ON
cmake --build build -j
```

With Embedded Boundary support:
```bash
cmake -S . -B build -DLEVELSET_ENABLE_EB=ON
cmake --build build -j
```

### Run

```bash
./build/levelset
```

Override parameters from command line:
```bash
./build/levelset nsteps=200 plot_int=20 cfl=0.4 u_x=0.3 u_y=0.1 u_z=0.0
```

### Example Configurations

**Basic model (no terrain effects):**
```bash
./build/levelset u_x=0.5 u_y=0.0
```

**Rothermel with terrain:**
```bash
./build/levelset u_x=0.5 u_y=0.0 use_terrain_effects=true terrain_slope=15.0 terrain_aspect=90.0
```

**FARSITE with Anderson L/W ratio:**
```bash
./build/levelset u_x=0.5 u_y=0.0 use_terrain_effects=true use_farsite_model=true terrain_slope=15.0 terrain_aspect=90.0
```

**With firebrand spotting:**
```bash
./build/levelset skip_levelset=1 farsite.enable=1 spotting.enable=1 spotting.P_base=0.03 spotting.d_mean=0.15 u_x=0.4
```

## Key Runtime Parameters

- **Grid/Domain**: `n_cell=64`, `prob_lo_x/y/z=0.0`, `prob_hi_x/y/z=1.0`
- **Time stepping**: `nsteps=300`, `cfl=0.5`, `plot_int=50`
- **Reinitialization**: `reinit_int=20`, `reinit_iters=20`
- **Velocity**: `u_x=0.25`, `u_y=0.0`, `u_z=0.0`, `velocity_file="wind.csv"`
- **Time-dependent wind**: `use_time_dependent_wind=0`, `wind_time_spacing=60.0`
- **Terrain**: `use_terrain_effects=false`, `terrain_slope=0.0`, `terrain_aspect=0.0`
- **FARSITE**: `farsite.enable=0`, `farsite.length_to_width_ratio=3.0`
- **Fuel model**: `rothermel.fuel_model=FM4` (FM1-FM13 available)
- **Spotting**: `spotting.enable=0`, `spotting.P_base=0.02`
- **Crown fire**: `crown.enable=0`, `crown.CBH=4.0`, `crown.CBD=0.15`

See the [online documentation](https://hgopalan.github.io/wildfire_levelset/) for complete parameter reference.

## Testing

Run regression tests:
```bash
cd build
ctest -L regtest --output-on-failure
```

Or use the custom target:
```bash
make regtest
```

Available regression tests:
- `basic_levelset` - Basic level-set advection
- `farsite_ellipse` - FARSITE elliptical expansion
- `rothermel_fuel` - Rothermel with fuel models
- `anderson_lw` - Anderson dynamic L/W ratio
- `reinitialization` - Level-set reinitialization
- `ellipse_sdf` - Elliptical SDF initial conditions
- `eb_implicit` - Embedded boundary capabilities
- `firebrand_spotting` - Firebrand spotting model
- `crown_initiation` - Crown fire initiation
- `bulk_fuel_consumption` - Fuel consumption modeling
- `3d_sphere` - Full 3D simulation (3D builds only)
- `terrain_wind` - External terrain and wind (2D builds only)

## Output

Plotfiles are written as `plt####` directories containing:
- `phi` - Level-set function (signed distance or indicator)
- `velx/y/z` - Velocity field components
- `farsite_dx/dy/dz` - FARSITE spread displacements
- `R` - Rothermel rate of spread
- `spot_prob/count/dist/active` - Spotting fields (if enabled)
- `fuel_consumption` - Fuel consumption fraction (if enabled)
- `crown_fraction` - Crown fire fraction (if enabled)

View with ParaView or other AMReX-compatible visualization tools.

## Limitations and Future Work

### Current Limitations

1. **Simplified Fuel Modeling**
   - Limited to NFFL fuel models (FM1-FM13)
   - No comprehensive fuel and moisture databases
   - Constant fuel properties across domain (no spatial variation except through fuel model selection)

2. **Wind Field Constraints**
   - Support for constant, prescribed, and time-dependent wind fields (2D only)
   - Spatially-varying wind supported via CSV input files
   - Time-dependent wind with automatic interpolation between time snapshots
   - One-way coupling only (fire does not affect wind)

3. **Terrain Representation**
   - 2D terrain data only (no 3D terrain in 3D simulations)
   - Terrain effects use simplified slope/aspect model
   - No automated terrain preprocessing or DEM file readers

4. **Physical Sub-models**
   - Simplified crown fire model (Van Wagner empirical criteria only)
   - Stochastic firebrand spotting without detailed physics-based transport
   - No radiation or convection heat transfer models
   - No explicit fuel moisture dynamics

5. **Validation and Calibration**
   - Limited validation against real wildfire data
   - Model parameters not extensively calibrated for diverse fuel types and conditions
   - No automated parameter optimization or data assimilation

### Scope for Future Work

1. **Enhanced Coupling**
   - Two-way coupling with atmospheric models (ERF) including surface heat fluxes
   - Dynamic wind field updates based on fire-induced convection
   - Integration with weather prediction systems

2. **Advanced Fuel Modeling**
   - Comprehensive fuel and moisture databases
   - Spatially-varying fuel properties
   - Dynamic fuel moisture evolution
   - Multi-layer fuel structure (surface, ladder, canopy)

3. **Improved Physical Fidelity**
   - Physics-based firebrand transport (ballistic trajectories, buoyancy)
   - Radiation heat transfer to unburned fuels
   - Detailed crown fire propagation models
   - Fuel consumption rate modeling

4. **3D Capabilities**
   - Full 3D terrain representation
   - 3D atmospheric boundary layer effects
   - Vertical fire spread and plume dynamics

5. **Data Integration**
   - Automated readers for standard terrain formats (GeoTIFF, DEM)
   - Weather data ingestion (GRIB, NetCDF)
   - Real-time data assimilation for operational forecasting

6. **Computational Enhancements**
   - GPU acceleration for large-scale simulations
   - Improved parallel scaling
   - Adaptive mesh refinement optimization

7. **Validation and Verification**
   - Systematic comparison with experimental burns
   - Validation against historical wildfire data
   - Uncertainty quantification
   - Sensitivity analysis frameworks

## References

- **Richards, G.D. (1990)**. "An elliptical growth model of forest fire fronts and its numerical solution." International Journal of Numerical Methods in Engineering.
- **Anderson, H.E. (1983)**. "Predicting wind-driven wild land fire size and shape." Research Paper INT-305, USDA Forest Service.
- **Rothermel, R.C. (1972)**. "A mathematical model for predicting fire spread in wildland fuels." Research Paper INT-115, USDA Forest Service.
- **Van Wagner, C.E. (1977)**. "Conditions for the start and spread of crown fire." Canadian Journal of Forest Research.

## License

See [LICENSE](LICENSE) file.

## Contact

For questions or issues, please open an issue on the GitHub repository.
