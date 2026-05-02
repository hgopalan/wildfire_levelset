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
- **Rothermel fire spread model** with Anderson 13 (FM1-FM13) and Scott & Burgan 40 (GR1–GR9, GS1–GS4, SH1–SH9, TU1–TU5, TL1–TL9, SB1–SB4) fuel databases
- **FARSITE elliptical expansion** (Richards 1990) with Anderson L/W ratio
- **Terrain effects** including slope and aspect corrections via constant values, terrain files, or FARSITE landscape files
- **FARSITE landscape files** (`.lcp`) with per-cell elevation, slope, aspect, and fuel model (landscape file takes precedence over terrain file or constant values)
- **Stochastic firebrand spotting** with probability-based model
- **Physics-based firebrand spotting** using Albini (1983) lofting height and 2-D trajectory integration
- **Crown fire initiation** (Van Wagner 1977)
- **Bulk fuel consumption** modeling
- **Multi-point ignition** from a CSV file of coordinates
- **2D and 3D** simulation capabilities
- **Embedded Boundary (EB)** support for complex geometries
- **AMReX-based** with GPU-ready kernels (CUDA/HIP/SYCL via AMReX) and optional MPI parallelism

## Prerequisites

- C++17 compiler (GCC 9+, Clang 10+, MSVC 2019+)
- CMake (3.20+)
- Git
- **For GPU builds**: CUDA Toolkit 11+ (CUDA), ROCm 5+ (HIP), or oneAPI (SYCL)
- **For MPI builds**: An MPI implementation (OpenMPI, MPICH, etc.)

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

Default 3D build (CPU-only):
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

With MPI parallelism:
```bash
cmake -S . -B build -DLEVELSET_ENABLE_MPI=ON
cmake --build build -j
```

With CUDA GPU acceleration:
```bash
cmake -S . -B build -DLEVELSET_GPU_BACKEND=CUDA
cmake --build build -j
```

With AMD GPU (HIP/ROCm):
```bash
cmake -S . -B build -DLEVELSET_GPU_BACKEND=HIP
cmake --build build -j
```

With Intel GPU (SYCL/oneAPI):
```bash
cmake -S . -B build -DLEVELSET_GPU_BACKEND=SYCL
cmake --build build -j
```

Combined (e.g., 2D + MPI + CUDA):
```bash
cmake -S . -B build -DLEVELSET_DIM_2D=ON -DLEVELSET_ENABLE_MPI=ON -DLEVELSET_GPU_BACKEND=CUDA
cmake --build build -j
```

#### CMake Build Options

| Option | Default | Description |
|--------|---------|-------------|
| `LEVELSET_DIM_2D` | `OFF` | Build for 2D instead of 3D |
| `LEVELSET_ENABLE_EB` | `OFF` | Enable Embedded Boundary support |
| `LEVELSET_ENABLE_MPI` | `OFF` | Enable MPI parallelism |
| `LEVELSET_GPU_BACKEND` | `NONE` | GPU backend: `NONE`, `CUDA`, `HIP`, or `SYCL` |
| `LEVELSET_BUILD_DOCS` | `OFF` | Build Sphinx documentation |
| `LEVELSET_USE_VENDORED_AMREX` | `ON` | Use the bundled AMReX submodule |

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
- **Time stepping**: `nsteps=300`, `cfl=0.5`, `plot_int=50`, `final_time=0.0` (overrides `nsteps` when > 0)
- **Reinitialization**: `reinit_int=20`, `reinit_iters=20`
- **Velocity**: `u_x=0.25`, `u_y=0.0`, `u_z=0.0`, `velocity_file="wind.csv"`
- **Time-dependent wind**: `use_time_dependent_wind=0`, `wind_time_spacing=60.0`
- **Terrain**: `use_terrain_effects=false`, `terrain_slope=0.0`, `terrain_aspect=0.0`
- **Terrain/Landscape files**: `rothermel.terrain_file=""`, `rothermel.landscape_file=""` (landscape file takes precedence over terrain file and constant slope/aspect)
- **FARSITE**: `farsite.enable=0`, `farsite.length_to_width_ratio=3.0`, `farsite.use_anderson_LW=0`
- **Bulk fuel consumption**: `farsite.use_bulk_fuel_consumption=0`, `farsite.tau_residence=60.0`
- **Fuel model**: `rothermel.fuel_model=FM4` — Anderson 13 (FM1–FM13) and Scott & Burgan 40 (GR1–GR9, GS1–GS4, SH1–SH9, TU1–TU5, TL1–TL9, SB1–SB4)
- **Stochastic spotting**: `spotting.enable=0`, `spotting.P_base=0.02`, `spotting.d_mean=0.1`, `spotting.distance_model=lognormal`
- **Albini spotting**: `albini_spotting.enable=0`, `albini_spotting.terminal_velocity=1.0`, `albini_spotting.P_base=0.01`, `albini_spotting.I_B_min=10.0`
- **Crown fire**: `crown.enable=0`, `crown.CBH=4.0`, `crown.CBD=0.15`, `crown.FMC=100.0`
- **Multi-point ignition**: `fire_points_file=""` (CSV with `X Y [Z]` columns), `fire_gaussian_sigma=0.0` (≤0 = auto: 3 cells)

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
- `firebrand_spotting` - Stochastic firebrand spotting model
- `albini_spotting` - Albini (1983) physics-based firebrand spotting
- `crown_initiation` - Crown fire initiation
- `bulk_fuel_consumption` - Fuel consumption modeling
- `3d_sphere` - Full 3D simulation (3D builds only)
- `terrain_wind` - External terrain and wind (2D builds only)
- `time_dependent_wind` - Time-varying wind fields
- `terrain_wind_preprocess` - Preprocessing tool integration
- `landfire_farsite` - FARSITE with auto-downloaded LANDFIRE landscape (requires Python3 + `terrain_wind_preprocess.py`)

## Output

Plotfiles are written as `plt####` directories containing:
- `phi` - Level-set function (signed distance or indicator)
- `velx/y/z` - Velocity field components
- `farsite_dx/dy/dz` - FARSITE spread displacements
- `R` - Rothermel rate of spread
- `spot_prob/count/dist/active` - Stochastic spotting fields (if enabled)
- `albini_spot_count/active` - Albini spotting fields (if enabled)
- `fuel_consumption` - Fuel consumption fraction (if enabled)
- `crown_fraction` - Crown fire fraction (if enabled)

View with ParaView or other AMReX-compatible visualization tools.

## Limitations and Known Constraints

1. **Simplified fuel modeling**
   - Anderson 13 (FM1–FM13) and Scott & Burgan 40 fuel models are supported; custom or dynamic fuel models are not.
   - Fuel moisture is a user-specified constant — no dynamic moisture transport.
   - Spatially-varying fuel is supported through FARSITE landscape files (LANDFIRE data via `terrain_wind_preprocess.py`).

2. **Wind field**
   - Spatially-varying wind supported via CSV files; fully 3-D spatially-varying wind fields are not.
   - Time-dependent wind uses linear interpolation between snapshots (2D only).
   - One-way coupling only — fire does not modify the wind field.

3. **Terrain**
   - 2D slope/aspect representation; true 3D terrain geometry is not modelled in 3D simulations.
   - Landscape files provide per-cell slope and aspect for 2D simulations.

4. **Physical sub-models**
   - Crown fire uses Van Wagner (1977) empirical criteria only; no mechanistic crown fire spread model.
   - Firebrand transport uses either a stochastic model or the Albini (1983) 2-D trajectory; no full 3-D plume transport.
   - No radiation or convective heat transfer.
   - No post-frontal smouldering or ember cast accumulation.

5. **Parallel execution**
   - GPU kernels are implemented via AMReX `ParallelFor`/`AMREX_GPU_DEVICE` macros and will execute on GPU when built with `LEVELSET_GPU_BACKEND=CUDA/HIP/SYCL`. Some serial CPU fallback paths (e.g., spotting, landscape I/O) remain.
   - MPI domain decomposition is provided by AMReX when built with `LEVELSET_ENABLE_MPI=ON`.

6. **Validation**
   - Limited validation against real wildfire events; parameters are calibrated to standard fuel models.

## Comparison with FARSITE and WRF-SFIRE

| Capability | wildfire_levelset | FARSITE | WRF-SFIRE |
|---|---|---|---|
| **Fire spread model** | Rothermel ROS + FARSITE ellipse (level-set) | Rothermel ROS + Huygens wavelet | Rothermel ROS (level-set) |
| **Terrain representation** | 2-D slope/aspect (constant or per-cell) | Full 2-D landscape (.lcp) | Full 3-D terrain from WRF grid |
| **Wind coupling** | One-way (prescribed; WRF output converted to CSV via `terrain_wind_preprocess.py`) | One-way (prescribed; gridded wind) | Two-way (fire ↔ atmosphere) |
| **Atmospheric model** | None | None | Full WRF mesoscale atmosphere |
| **Spatial resolution** | User-defined (AMReX grid) | Landscape raster resolution | WRF grid resolution |
| **Fuel models** | Anderson 13 + Scott & Burgan 40 | Anderson 13 + Scott & Burgan 40 | Anderson 13 |
| **Fuel moisture** | User-specified constant | Dead/live moisture per size class | Dead/live moisture (prescribed) |
| **Crown fire** | Van Wagner (1977) initiation | Van Wagner (1977) initiation | Van Wagner (1977) initiation |
| **Firebrand spotting** | Stochastic + Albini (1983) 2-D trajectory | Albini empirical spotting | Not included by default |
| **Bulk fuel burnout** | Residence-time model | Full burnout tracking per cell | Post-frontal fuel consumption |
| **Multi-ignition** | CSV ignition file | Interactive ignition map | WPS/WRF fire restart |
| **GPU acceleration** | Yes (AMReX CUDA/HIP/SYCL kernels) | No | No |
| **MPI parallelism** | Yes (AMReX domain decomposition) | No | Yes (WRF MPI) |
| **Embedded boundaries** | Yes (AMReX EB) | No | No |
| **Operational use** | Research / prototype | Operational (US Forest Service) | Research / operational |
| **Post-processing** | AMReX plotfiles (ParaView) | GIS shapefiles / rasters | WRF standard output (NCO/Python) |

### Key differences from FARSITE

- **Level-set vs. Huygens wavelets**: This solver tracks the fire perimeter as the zero contour of a smooth signed-distance function; FARSITE explicitly propagates individual vertex wavelets. The level-set approach handles complex topology changes (merging fronts, islands) automatically.
- **Embedded Boundary**: AMReX EB allows irregular obstacles (buildings, fuel breaks) on the Cartesian grid without remeshing.
- **GPU and MPI**: The AMReX backend supports CUDA/HIP/SYCL and MPI for large domains; FARSITE is serial and CPU-only.
- **No interactive GUI**: This solver is driven by an input text file and command-line arguments; FARSITE ships with a Windows GUI.

### Key differences from WRF-SFIRE

- **No atmospheric coupling**: This solver uses prescribed wind fields (constant, CSV, or WRF output); WRF-SFIRE fully couples the fire with the WRF atmospheric model, including fire-induced wind, heat flux, and smoke transport.
- **Albini spotting vs. none**: WRF-SFIRE does not include firebrand spotting by default.
- **Simpler setup**: Running this solver requires only CMake and an AMReX build; WRF-SFIRE requires a full WRF stack (NetCDF, MPI, WPS, WRF pre-processing).
- **GPU-native kernels**: WRF-SFIRE is primarily CPU-MPI; this solver's hotspot loops use AMReX GPU kernels.



## References

- **Richards, G.D. (1990)**. "An elliptical growth model of forest fire fronts and its numerical solution." International Journal of Numerical Methods in Engineering.
- **Anderson, H.E. (1983)**. "Predicting wind-driven wild land fire size and shape." Research Paper INT-305, USDA Forest Service.
- **Anderson, H.E. (1982)**. "Aids to determining fuel models for estimating fire behavior." Research Paper INT-122, USDA Forest Service.
- **Rothermel, R.C. (1972)**. "A mathematical model for predicting fire spread in wildland fuels." Research Paper INT-115, USDA Forest Service.
- **Van Wagner, C.E. (1977)**. "Conditions for the start and spread of crown fire." Canadian Journal of Forest Research.
- **Albini, F.A. (1983)**. "Potential spotting distance from wind-driven surface fires." Research Paper INT-309, USDA Forest Service.
- **Scott, J.H. and Burgan, R.E. (2005)**. "Standard Fire Behavior Fuel Models: A Comprehensive Set for Use with Rothermel's Surface Fire Spread Model." Technical Report RMRS-GTR-153, USDA Forest Service.

## License

See [LICENSE](LICENSE) file.

## Contact

For questions or issues, please open an issue on the GitHub repository.
