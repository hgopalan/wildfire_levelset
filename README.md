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

## Tools

Python utilities live in the `tools/` directory.

### `terrain_wind_preprocess.py` – Unified terrain, landscape, and wind preprocessing

The primary preprocessing tool. Downloads SRTM elevation data and LANDFIRE
fuel/slope/aspect rasters for a lat/lon bounding box, extracts wind from a
WRF-style netCDF file, and automatically writes a ready-to-run `inputs.i` file
for a FARSITE simulation (no firebrand spotting, no crown fire).

```bash
# SRTM terrain + LANDFIRE landscape from a lat/lon bounding box
python3 tools/terrain_wind_preprocess.py \
    --lat-min 40 --lat-max 40.5 \
    --lon-min -105 --lon-max -104.5

# Add WRF wind (bbox is read automatically from the WRF file)
python3 tools/terrain_wind_preprocess.py \
    --wrf-file wrfout_d01 \
    --wind wind.csv

# Extract a range of WRF time steps and compute final_time automatically
python3 tools/terrain_wind_preprocess.py \
    --wrf-file wrfout_d01 \
    --wind wind.csv \
    --time-range 0:5

# Skip SRTM download; extract wind only
python3 tools/terrain_wind_preprocess.py \
    --wrf-file wrfout_d01 \
    --wind wind.csv \
    --no-terrain --no-landscape

# Use local rasters for landscape (no LANDFIRE API call)
python3 tools/terrain_wind_preprocess.py \
    --lat-min 40 --lat-max 40.5 \
    --lon-min -105 --lon-max -104.5 \
    --elev-file elev.tif \
    --slope-file slope.tif \
    --aspect-file aspect.tif \
    --fuel-file fbfm13.tif
```

Key options:
- `--wrf-file FILE` – WRF netCDF; bounding box is read from file automatically
- `--time-range T1:TN` – extract WRF time steps T1…TN (inclusive); the tool
  sets `final_time = (TN − T1) × wind_time_spacing` in the generated `inputs.i`
- `--time-index N` – single WRF time step (default: 0)
- `--interpolate-wind` – interpolate WRF wind to the SRTM terrain grid
- `--subsample N` – keep every N-th point (default: 1)
- `--inputs FILE` – output `inputs.i` filename (default: `inputs.i`)
- `--no-inputs` – skip `inputs.i` generation
- `--no-terrain` / `--no-landscape` / `--no-wind` – skip individual steps
- `--vintage YEAR` – LANDFIRE vintage year (default: 2020)
- `--keep-nonburnable` – retain non-burnable LANDFIRE pixels
- `--cache-dir DIR`, `--timeout N` – LANDFIRE download options
- `--elev-file`, `--slope-file`, `--aspect-file`, `--fuel-file` – local rasters

Requires: `pip install elevation rasterio numpy pyproj requests netCDF4`

The tool writes:
1. `terrain.xyz` – UTM terrain (`X Y Z`) for `rothermel.terrain_file`
2. `landscape.lcp` – FARSITE landscape (`X Y ELEVATION SLOPE ASPECT FUEL_MODEL`)
   for `rothermel.landscape_file`
3. `wind.csv` (or `wind.csv`, `wind_1.csv`, … for multiple time steps) – wind
   field for `velocity_file`
4. `inputs.i` – ready-to-run FARSITE inputs file (no spotting, no crown)

Use the generated files in a simulation:
```
rothermel.terrain_file    = terrain.xyz
rothermel.landscape_file  = landscape.lcp
velocity_file             = wind.csv
```

> **Deprecated tools**: The following scripts have been moved to
> `tools/deprecated/` and are superseded by `terrain_wind_preprocess.py`:
> `dem_to_xyz.py`, `landfire_to_lcp.py`, `srtm_to_xyz_stl.py`,
> `srtm_landfire_to_terrain.py`, `wrf_to_terrain_wind.py`, and
> `utm_convert.py`.

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
- `landfire_farsite` - FARSITE with auto-downloaded LANDFIRE landscape (requires Python3 + `srtm_landfire_to_terrain`)

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
   - No comprehensive fuel moisture databases
   - Constant fuel properties across domain (no spatial variation except through fuel model selection)
   - Spatially-varying fuel models are supported via landscape files but require LANDFIRE data

2. **Wind Field Constraints**
   - Support for constant, prescribed, and time-dependent wind fields (2D only)
   - Spatially-varying wind supported via CSV input files
   - Time-dependent wind with automatic interpolation between time snapshots
   - One-way coupling only (fire does not affect wind)

3. **Terrain Representation**
   - 2D terrain data only (no 3D terrain in 3D simulations)
   - Terrain effects use simplified slope/aspect model
   - Automated DEM preprocessing supported via `srtm_landfire_to_terrain.py` (SRTM download + GeoTIFF conversion)
   - LANDFIRE landscape files (elevation, slope, aspect, fuel model) downloadable via `srtm_landfire_to_terrain.py`

4. **Physical Sub-models**
   - Simplified crown fire model (Van Wagner empirical criteria only)
   - Stochastic firebrand spotting without detailed physics-based transport
   - No radiation or convection heat transfer models
   - No explicit fuel moisture dynamics

5. **Validation and Calibration**
   - Limited validation against real wildfire data
   - Model parameters not extensively calibrated for diverse fuel types and conditions
   - No automated parameter optimization or data assimilation

6. **LANDFIRE Integration**
   - Automated download via LFPS API requires internet connectivity
   - Non-burnable LANDFIRE classes (urban, water, agriculture, snow/ice, barren) are filtered out by default; use `--keep-nonburnable` to retain these pixels (written as fuel code 0)

### Scope for Future Work

1. **Enhanced Coupling**
   - Two-way coupling with atmospheric models (ERF) including surface heat fluxes
   - Dynamic wind field updates based on fire-induced convection
   - Integration with weather prediction systems

2. **Advanced Fuel Modeling**
   - Comprehensive fuel and moisture databases
   - Support for Scott/Burgan 40-class fuel models (FBFM40) from LANDFIRE
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
   - Automated readers for standard terrain formats (GeoTIFF, DEM) — implemented via `srtm_landfire_to_terrain.py`
   - LANDFIRE data download for arbitrary extents — implemented via `srtm_landfire_to_terrain.py`
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
