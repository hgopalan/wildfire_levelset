> ⚠️ **AI-Generated Code Disclaimer**
> A major portion of this codebase was written with the assistance of AI tools (GitHub Copilot / large language models).
> It has not been exhaustively validated against operational wildfire prediction systems.
> **Use with caution** — review all outputs carefully before applying to real-world fire management decisions.

# Wildfire-AMR

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
- **Andrews (2018) wind adjustments** for the Rothermel model: Wind Adjustment Factor (WAF) converting 20-ft open wind to midflame height (`rothermel.use_waf = 1`); maximum effective wind speed (MEWS) limit (`rothermel.use_wind_limit = 1`); both are per-cell when a landscape file is used
- **Balbi (2009) physical fire spread model** – radiation-driven, physics-based ROS; fully replaces Rothermel; auto-generates a Balbi fuel parameter table at startup; accepts all LCP inputs and adds thermal parameters via parmparse (`balbi.*`); fully compatible with all six Viegas wind-terrain options (including Option 2 Viegas-ROS override which uses the Balbi amplitude baseline)
- **Cheney & Gould (1995 / 1998) grassland fire spread model** – empirical model calibrated for open Australian grasslands; piecewise wind-speed formula with moisture and curing corrections; activated by `fire_spread_model = cheney_gould`; configured via `cheney_gould.*`
- **FARSITE elliptical expansion** (Richards 1990) with Anderson L/W ratio; Eulerian level-set implementation of the Huygens wavelet principle
- **Terrain effects** including slope and aspect corrections via constant values, terrain files, or FARSITE landscape files
- **FARSITE landscape files** (`.lcp`) with per-cell elevation, slope, aspect, and fuel model (landscape file takes precedence over terrain file or constant values)
- **Per size-class fuel moisture**: 1-hr, 10-hr, 100-hr dead fuel; live herbaceous and live woody moisture inputs for the multi-class Rothermel (1972) reaction intensity
- **Fire behavior diagnostics**: Byram (1959) fireline intensity [kW/m] and flame length [m] written to every plotfile
- **Weise & Biging (1996) fire whirl model** – optional diagnostic sub-model computing flame tilt, whirl height/radius, angular velocity, and tangential velocity from fireline intensity and wind; enabled via `weise_biging.enable = 1`
- **Viegas (2004) eruptive fire diagnostics** – optional parallel diagnostic that computes the Viegas exponential slope enhancement factor, eruptive-regime flag (critical slope), Viegas ROS, ROS excess vs primary model, and flame-tilt angle for hazard assessment; enabled via `viegas.enable = 1`; works with both Rothermel (R₀ baseline) and Balbi (amplitude A baseline); by default diagnostic-only, but can be coupled to fire propagation via the **wind-terrain feedback model** (see below)
- **Wind-terrain feedback models** – seven selectable options (`wind_terrain.model`) for how terrain-induced or fire-feedback winds modify fire spread: (1) default, no modification; (2) Viegas ROS as actual spread rate (works with both Rothermel and Balbi); (3) Viegas-induced buoyancy wind as velocity perturbation (eruptive cells only); (4) Rothermel (1983) canyon wind amplification; (5) Viegas & Neto (1994) buoyancy-driven upslope wind; (6) Pimont et al. (2009) exponential slope correction; (7) **WindNinja ridge/canyon empirical speed-up** based on wind-slope alignment: ridge acceleration when wind climbs upslope (`f = 1 + k_ridge × tan φ × alignment`), canyon channeling when wind descends (`f = 1 + k_canyon_wn × tan φ × |alignment|`)
- **Heat flux MultiFab** – spatially-varying or uniform fire heat release rate [W/m²] initialised from a value (`heat_flux.value`) or an ASCII file (`heat_flux.file`); drives two WindNinja-style fire-induced wind corrections: (a) upward convective velocity from fire plume buoyancy (`w_up = k_upward × sqrt(g × Q × H / (ρ Cp T_a))`); (b) induced horizontal inflow directed toward the fire perimeter; for the Balbi model the buoyancy velocity `v_b` is additionally augmented by the fire heat flux (`v_b_eff = sqrt(v_b_fuel² + v_b_Q²)`)
- **GIS output**: `tools/plotfile_to_geotiff.py` converts plotfiles to GeoTIFF rasters and GeoJSON fire-perimeter contours
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

**Basic level-set (no fire spread model, constant velocity):**
```bash
./build/levelset u_x=0.5 u_y=0.0
```

**Rothermel + level-set with terrain file:**
```bash
./build/levelset fire_spread_model=rothermel propagation_method=levelset \
    rothermel.terrain_file=terrain.xyz u_x=0.5 u_y=0.0
```

**Rothermel + FARSITE ellipse propagation with Anderson L/W ratio:**
```bash
./build/levelset fire_spread_model=rothermel propagation_method=farsite \
    farsite.use_anderson_LW=1 u_x=0.5 u_y=0.0
```

**Balbi + FARSITE propagation:**
```bash
./build/levelset fire_spread_model=balbi propagation_method=farsite \
    farsite.use_anderson_LW=1 u_x=0.5 u_y=0.0
```

**With firebrand spotting (FARSITE propagation):**
```bash
./build/levelset propagation_method=farsite spotting.enable=1 spotting.P_base=0.03 spotting.d_mean=0.15 u_x=0.4
```

## Key Runtime Parameters

- **Grid/Domain**: `n_cell=64`, `prob_lo_x/y/z=0.0`, `prob_hi_x/y/z=1.0`
- **Time stepping**: `nsteps=300`, `cfl=0.5`, `plot_int=50`, `final_time=0.0` (overrides `nsteps` when > 0)
- **Reinitialization**: `reinit_int=20`, `reinit_iters=20`
- **Velocity**: `u_x=0.25`, `u_y=0.0`, `u_z=0.0`, `velocity_file="wind.csv"`
- **Time-dependent wind**: `use_time_dependent_wind=0`, `wind_time_spacing=60.0`
- **Fire spread model**: `fire_spread_model=rothermel` — choose `rothermel` (default), `balbi`, or `cheney_gould`
- **Propagation method**: `propagation_method=levelset` — choose `levelset` (default) or `farsite`
- **Terrain/Landscape files**: `rothermel.terrain_file=""`, `rothermel.landscape_file=""` (landscape file takes precedence over terrain file and constant slope/aspect)
- **FARSITE parameters**: `farsite.length_to_width_ratio=3.0`, `farsite.use_anderson_LW=0`
- **Bulk fuel consumption**: `farsite.use_bulk_fuel_consumption=0`, `farsite.tau_residence=60.0`
- **Fuel model**: `rothermel.fuel_model=FM4` — Anderson 13 (FM1–FM13) and Scott & Burgan 40 (GR1–GR9, GS1–GS4, SH1–SH9, TU1–TU5, TL1–TL9, SB1–SB4)
- **Andrews (2018) wind adjustments**: `rothermel.use_waf=0` (WAF on/off), `rothermel.use_wind_limit=0` (MEWS cap on/off)
- **Balbi parameters**: `balbi.T_a=300.0`, `balbi.T_f=1000.0`, `balbi.T_i=600.0` (active when `fire_spread_model=balbi`)
- **Cheney–Gould parameters**: `cheney_gould.moisture=10.0` [%], `cheney_gould.curing=1.0` (active when `fire_spread_model=cheney_gould`)
- **Stochastic spotting**: `spotting.enable=0`, `spotting.P_base=0.02`, `spotting.d_mean=0.1`, `spotting.distance_model=lognormal`
- **Albini spotting**: `albini_spotting.enable=0`, `albini_spotting.terminal_velocity=1.0`, `albini_spotting.P_base=0.01`, `albini_spotting.I_B_min=10.0`
- **Crown fire**: `crown.enable=0`, `crown.CBH=4.0`, `crown.CBD=0.15`, `crown.FMC=100.0`
- **Weise & Biging fire whirl**: `weise_biging.enable=0`, `weise_biging.c_r=0.1`, `weise_biging.I_B_min=1.0`
- **Viegas eruptive fire diagnostics**: `viegas.enable=0`, `viegas.a_V=1.83`, `viegas.tan_phi_c=0.4`, `viegas.T_a=300.0`, `viegas.T_f=1000.0`
- **Wind-terrain feedback model**: `wind_terrain.model=none` — choose `none` (default), `viegas_ros`, `viegas_wind`, `canyon_wind`, `viegas_neto`, or `pimont`; `wind_terrain.k_canyon=1.0` (Option 4); `wind_terrain.k_pimont=0.5` (Option 6)
- **Multi-point ignition**: `fire_points_file=""` (CSV with `X Y [Z]` columns), `fire_gaussian_sigma=0.0` (≤0 = auto: 3 cells)

See the [online documentation](https://hgopalan.github.io/wildfire_levelset/) for complete parameter reference.

## Balbi (2009) Physical Fire Spread Model

A radiation-driven, physics-based ROS model that replaces Rothermel when `fire_spread_model = balbi`. The flame tilt angle is derived from wind speed and a buoyancy velocity, and the rate of spread follows from a radiation balance between the tilted flame and unburned fuel ahead. Thermal parameters are set via `balbi.*` parmparse keys; fuel geometry is shared with Rothermel. The solver auto-generates a per-fuel Balbi parameter table at startup.

See the [full documentation](https://hgopalan.github.io/wildfire_levelset/mathematical_models.html#balbi-2009-physical-fire-spread-model) for equations, parameter tables, and examples.

## Andrews (2018) Wind Adjustments for Rothermel

Two optional wind adjustments from Andrews (2018) improve physical accuracy when wind input is from NWP or WRF (20-ft / 10-m height): the **Wind Adjustment Factor** (WAF, `rothermel.use_waf = 1`) converts 20-ft open wind to midflame height; the **Maximum Effective Wind Speed** (MEWS, `rothermel.use_wind_limit = 1`) caps the wind factor at 90% of the reaction intensity to prevent runaway ROS at high wind speeds. Both are per-cell when a landscape file is active.

See the [full documentation](https://hgopalan.github.io/wildfire_levelset/mathematical_models.html#andrews-2018-wind-adjustments-for-rothermel) for equations and examples.

## Cheney & Gould (1995 / 1998) Grassland Fire Spread Model

A purely empirical ROS model calibrated on Australian grassland fires, activated with `fire_spread_model = cheney_gould`. Head-fire ROS is a piecewise-linear function of 10-m wind speed, scaled by a moisture correction factor (`cheney_gould.moisture`) and a curing factor (`cheney_gould.curing`). Terrain slope is intentionally omitted.

See the [full documentation](https://hgopalan.github.io/wildfire_levelset/mathematical_models.html#cheney-gould-1995-1998-grassland-fire-spread-model) for equations, parameters, and examples.

## Weise & Biging (1996) Fire Whirl Model

An optional diagnostic sub-model (`weise_biging.enable = 1`) that computes fire whirl geometry and kinematics (flame tilt, whirl height/radius, angular velocity, tangential velocity) from Byram fireline intensity and wind speed. It runs alongside any primary spread model without modifying fire propagation, and writes six diagnostic fields to every plotfile.

See the [full documentation](https://hgopalan.github.io/wildfire_levelset/mathematical_models.html#weise-biging-1996-fire-whirl-model) for equations, parameters, and examples.

## Viegas (2004) Eruptive Fire Diagnostics

An optional parallel diagnostic (`viegas.enable = 1`) that computes an exponential slope enhancement factor and identifies eruptive-fire cells (where Rothermel's quadratic slope factor under-predicts spread). Works alongside Rothermel or Balbi without modifying propagation by default; coupling to actual spread is controlled via `wind_terrain.model` (see below). Writes five diagnostic fields to every plotfile.

See the [full documentation](https://hgopalan.github.io/wildfire_levelset/mathematical_models.html#viegas-2004-eruptive-fire-diagnostics) for equations, parameters, and examples.

## Wind-Terrain Feedback Models

Seven selectable options (`wind_terrain.model`) control how terrain-induced or fire-feedback winds modify the effective wind seen by the spread model. The default (`none`) preserves existing behaviour. Options 2–6 couple Viegas or terrain-channel physics into propagation (Viegas ROS override, buoyancy upslope wind, Rothermel 1983 canyon wind, Viegas & Neto 1994 wind, Pimont et al. 2009 exponential correction); Option 7 (`windninja_ridge_canyon`) applies WindNinja ridge/canyon empirical speed-up.

See the [full documentation](https://hgopalan.github.io/wildfire_levelset/mathematical_models.html#wind-terrain-feedback-models-options-1-6) for equations, all options, parameters, and examples.

## Future Fire Spread Models (TODO)

The following spread models are candidates for future integration into the `fire_spread_model` interface. Each would be added as a new option (e.g., `fire_spread_model = <name>`) and plugged into the same ROS compute step, making it interchangeable with Rothermel and Balbi.

- [ ] **Cruz & Alexander (2010) crown fire ROS** – Physics-informed empirical model for active crown fire spread combining Rothermel surface ROS with a crown fire transition criterion. Reference: Cruz, M.G. & Alexander, M.E. (2010). *Assessing crown fire potential in coniferous forests of western North America.* Int. J. Wildland Fire 19, 8–21.
- [ ] **Linn (FIRETEC) simplified surrogate** – FIRETEC is a coupled fire–atmosphere CFD model; a reduced-order surrogate or table-lookup could expose its ROS as a `fire_spread_model = firetec` option. Reference: Linn, R.R. et al. (2002). *Studying wildfire behavior using FIRETEC.* Int. J. Wildland Fire 11, 233–246.
- [ ] **Data-driven / ML surrogate** – Neural-network or Gaussian-process ROS emulator trained on historical fire perimeters (e.g., MTBS/GeoMAC/NIFC) or high-fidelity simulation data from FIRETEC/FDS-Fire. Could be loaded as a saved model artifact at runtime.

To add a new spread model, implement a `compute_<name>_R` function (analogous to `compute_rothermel_R` in `src/compute_rothermel_R.H`) and register it with a new branch in `main.cpp` guarded by `if (inputs.fire_spread_model == "<name>")`.

## Tools

Python utilities live in the `tools/` directory.

### `wrf_wind_reader.py` – WRF wind extraction

Reads U and V wind components from a WRF-style netCDF file and writes one or
more CSV wind files suitable for use as `velocity_file` in the solver.

```bash
# Extract wind at t=0 (bounding box derived automatically from WRF file)
python3 tools/wrf_wind_reader.py --wrf-file wrfout_d01 --wind wind.csv

# Extract a range of WRF time steps
python3 tools/wrf_wind_reader.py \
    --wrf-file wrfout_d01 --wind wind.csv --time-range 0:5

# Clip to an explicit lat/lon bounding box
python3 tools/wrf_wind_reader.py \
    --wrf-file wrfout_d01 --wind wind.csv \
    --lat-min 40 --lat-max 40.5 --lon-min -105 --lon-max -104.5

# Interpolate wind onto an existing terrain XYZ grid
python3 tools/wrf_wind_reader.py \
    --wrf-file wrfout_d01 --wind wind.csv \
    --terrain-file terrain.xyz

# Interpolate wind onto a regular grid at 30 m resolution
python3 tools/wrf_wind_reader.py \
    --wrf-file wrfout_d01 --wind wind.csv \
    --grid-resolution 30
```

Key options:
- `--wrf-file FILE` – WRF netCDF file (required)
- `--time-range T1:TN` – extract WRF time steps T1…TN (inclusive)
- `--time-index N` – single WRF time step (default: 0)
- `--lat-min/max`, `--lon-min/max` – optional bounding box override
- `--terrain-file PATH` – existing terrain XYZ file; wind is interpolated
  onto its (x,y) grid
- `--grid-resolution M` – target grid spacing in metres for regular-grid
  interpolation
- `--subsample N` – keep every N-th grid point (default: 1)
- `--inputs FILE` – output `inputs.i` filename (default: `inputs.i`)
- `--no-inputs` – skip `inputs.i` generation

Requires: `pip install netCDF4 numpy pyproj`
Optional: `pip install scipy` (for `--terrain-file` or `--grid-resolution`)

---

### `srtm_terrain_reader.py` – SRTM terrain download

Downloads SRTM1 elevation data for a lat/lon bounding box and writes a UTM
terrain XYZ file suitable for `rothermel.terrain_file`.

```bash
# Download SRTM and write terrain.xyz
python3 tools/srtm_terrain_reader.py \
    --lat-min 40 --lat-max 40.5 \
    --lon-min -105 --lon-max -104.5

# Save intermediate GeoTIFF and subsample
python3 tools/srtm_terrain_reader.py \
    --lat-min 40 --lat-max 40.5 \
    --lon-min -105 --lon-max -104.5 \
    --tif srtm_raw.tif --subsample 2
```

Key options:
- `--lat-min/max`, `--lon-min/max` – latitude/longitude bounds (required)
- `--terrain FILE` – output terrain XYZ file (default: `terrain.xyz`)
- `--tif FILE` – keep intermediate SRTM GeoTIFF at this path
- `--subsample N` – keep every N-th point (default: 1)

Requires: `pip install elevation rasterio numpy pyproj scipy`

---

### `landscape_writer.py` – LANDFIRE landscape file writer

Downloads LANDFIRE elevation, slope, aspect, and fuel model rasters for a
lat/lon bounding box and writes a FARSITE landscape file (`.lcp`) for
`rothermel.landscape_file`.  Three remote sources are supported
(`cog` → AWS S3 COGs, `mspc` → Microsoft Planetary Computer, `lfps` → USGS
REST API).  Local GeoTIFFs can also be used to skip any download.

```bash
# Download LANDFIRE via API
python3 tools/landscape_writer.py \
    --lat-min 40 --lat-max 40.5 \
    --lon-min -105 --lon-max -104.5

# Use local rasters (no internet download)
python3 tools/landscape_writer.py \
    --lat-min 40 --lat-max 40.5 \
    --lon-min -105 --lon-max -104.5 \
    --elev-file elev.tif --slope-file slope.tif \
    --aspect-file aspect.tif --fuel-file fbfm13.tif

# Local fuel raster + SRTM-derived terrain
python3 tools/landscape_writer.py \
    --lat-min 40 --lat-max 40.5 \
    --lon-min -105 --lon-max -104.5 \
    --fuel-file fuel.tif --srtm-slope-aspect

# Scott & Burgan 40 fuel model system
python3 tools/landscape_writer.py \
    --lat-min 40 --lat-max 40.5 \
    --lon-min -105 --lon-max -104.5 \
    --fuel-model-type 40
```

Key options:
- `--lat-min/max`, `--lon-min/max` – latitude/longitude bounds (required)
- `--landscape FILE` – output `.lcp` file (default: `landscape.lcp`)
- `--subsample N` – keep every N-th point (default: 1)
- `--vintage YEAR` – LANDFIRE vintage year (default: 2020)
- `--fuel-model-type {13,40}` – Anderson 13 or Scott & Burgan 40
- `--keep-nonburnable` – retain non-burnable LANDFIRE pixels
- `--cache-dir DIR`, `--timeout N` – LANDFIRE download options
- `--sources LIST` – comma-separated source priority (`cog,mspc,lfps`)
- `--use-lfps` – force the legacy LFPS REST API
- `--no-vintage-fallback` – disable automatic vintage fallback
- `--elev-file`, `--slope-file`, `--aspect-file`, `--fuel-file` – local rasters
- `--srtm-slope-aspect` – derive terrain from SRTM when only `--fuel-file` is given

Requires: `pip install rasterio numpy pyproj requests`

---

Use the generated files in a simulation:
```
rothermel.terrain_file    = terrain.xyz
rothermel.landscape_file  = landscape.lcp
velocity_file             = wind.csv
```

> **Deprecated tools**: The following scripts have been moved to
> `tools/deprecated/` and are superseded by the new split tools above:
> `terrain_wind_preprocess.py`, `dem_to_xyz.py`, `landfire_to_lcp.py`,
> `srtm_to_xyz_stl.py`, `srtm_landfire_to_terrain.py`,
> `wrf_to_terrain_wind.py`, and `utm_convert.py`.

### `plotfile_to_geotiff.py` – GIS export of fire behavior fields

Converts AMReX 2-D plotfiles to GeoTIFF rasters and GeoJSON fire-perimeter
contours for import into QGIS, ArcGIS, or any GIS tool.

**Multi-level AMR plotfiles** (e.g. from external AMReX-based codes with
`finest_level > 0`) are now supported.  All available levels are read and
composited onto the Level_0 base grid — finer-level data replaces coarser
values where AMR patches exist — so the output GeoTIFF always covers the
full domain.

```bash
# Export all fire-behaviour variables from one plotfile (simulation units)
python3 tools/plotfile_to_geotiff.py plt0100 --outdir gis_out

# Export specific variables
python3 tools/plotfile_to_geotiff.py plt0100 \
    -v phi R fireline_intensity flame_length elevation \
    --outdir gis_out

# Georeference with a UTM origin (easting, northing) and EPSG code
python3 tools/plotfile_to_geotiff.py plt0100 \
    --utm-origin 450000 4200000 \
    --epsg 32613 \
    --outdir gis_out

# Convert every plt#### directory in the current directory
python3 tools/plotfile_to_geotiff.py --all --outdir gis_out

# Convert a multi-level AMR plotfile from an external AMReX-based code
python3 tools/plotfile_to_geotiff.py plt_amr0050 --outdir gis_out
```

Key options:
- `-v VAR [VAR …]` – export only specified variables (default: fire-behaviour subset)
- `--all-vars` – export every variable in the plotfile
- `--all` – process all `plt####` directories in the current working directory
- `--outdir DIR` – output directory (default: `gis_out`)
- `--utm-origin EASTING NORTHING` – UTM origin to add to simulation coordinates
- `--epsg CODE` – EPSG CRS code (e.g. `32613` for UTM zone 13N)

Requires: `pip install rasterio numpy` (plus `matplotlib` for GeoJSON perimeters)

The tool writes, for each plotfile and variable:
- `plt####_<variable>.tif` – single-band GeoTIFF (float32, deflate-compressed)
- `plt####_fire_perimeter.geojson` – phi = 0 contour as GeoJSON LineString (requires matplotlib)

**Multi-level compositing**: when `finest_level > 0` is present in the AMReX
Header, each level's FAB data is read from `Level_0/Cell_H`, `Level_1/Cell_H`,
etc.  Finer-level cells are block-averaged to the Level_0 resolution and
overlaid where valid, giving the best available data at each location while
preserving full-domain coverage.

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
- `cheney_gould_grassfire` - Cheney & Gould (1995/1998) grassland fire spread
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
- `terrain_wind_preprocess` - Preprocessing tool integration (uses `wrf_wind_reader.py`)
- `landfire_farsite` - FARSITE with auto-downloaded LANDFIRE landscape (requires Python3 + `landscape_writer.py`)

## Output

Plotfiles are written as `plt####` directories containing:
- `phi` - Level-set function (signed distance or indicator)
- `velx/y/z` - Velocity field components
- `farsite_dx/dy/dz` - FARSITE spread displacements
- `R` - Rothermel rate of spread [m/s]
- `fireline_intensity` - Byram (1959) fireline intensity [kW/m]
- `flame_length` - Byram (1959) flame length [m]
- `weise_flame_height` - Weise & Biging vertical flame height [m] (if enabled)
- `weise_flame_tilt` - Weise & Biging flame tilt angle [rad] (if enabled)
- `weise_whirl_height` - Weise & Biging fire whirl height [m] (if enabled)
- `weise_whirl_radius` - Weise & Biging fire whirl core radius [m] (if enabled)
- `weise_angular_velocity` - Weise & Biging angular velocity [rad/s] (if enabled)
- `weise_max_tang_vel` - Weise & Biging maximum tangential velocity [m/s] (if enabled)
- `spot_prob/count/dist/active` - Stochastic spotting fields (if enabled)
- `albini_spot_count/active` - Albini spotting fields (if enabled)
- `fuel_consumption` - Fuel consumption fraction (if enabled)
- `crown_fraction` - Crown fire fraction (if enabled)
- `elevation` - Terrain elevation [m]
- `slope` - Terrain slope [degrees]
- `aspect` - Terrain aspect [degrees]
- `fuel_model` - Per-cell fuel model code

View with ParaView or other AMReX-compatible visualization tools.  For GIS
import, use `tools/plotfile_to_geotiff.py` to convert plotfiles to GeoTIFF
rasters and GeoJSON fire-perimeter contours (see Tools section below).

## Limitations and Known Constraints

1. **Simplified fuel modeling**
   - Anderson 13 (FM1–FM13) and Scott & Burgan 40 fuel models are supported; custom or dynamic fuel models are not.
   - Fuel moisture is specified per size class (1-hr, 10-hr, 100-hr dead; live herbaceous and live woody) using user-supplied constants; no dynamic moisture transport in time.
   - Spatially-varying fuel is supported through FARSITE landscape files (LANDFIRE data via `landscape_writer.py`).

2. **Wind field**
   - Spatially-varying wind supported via CSV files; fully 3-D spatially-varying wind fields are not.
   - Time-dependent wind uses linear interpolation between snapshots (2D only).
   - One-way coupling by default — fire does not modify the wind field. **Wind-terrain feedback models** (`wind_terrain.model`) allow terrain-induced wind effects to modify the effective wind seen by Rothermel, but atmospheric feedback from the fire itself is not modelled.

3. **Terrain**
   - Full 2-D landscape representation: per-cell elevation, slope, aspect, and fuel model from FARSITE `.lcp` landscape files or terrain XYZ files (consistent with the approach used in FARSITE).  True 3-D terrain geometry (for 3-D simulations) and canopy fuel layers can be added in future.
   - Landscape files provide per-cell slope and aspect for spatially-varying terrain effects.

4. **Physical sub-models**
   - Crown fire uses Van Wagner (1977) empirical criteria only; no mechanistic crown fire spread model.  Canopy fuel parameters (base height, bulk density, foliar moisture) are currently user-specified constants; spatially-varying canopy data can be added in future.
   - Firebrand transport uses either a stochastic model or the Albini (1983) 2-D trajectory; no full 3-D plume transport.
   - No radiation or convective heat transfer.
   - No post-frontal smouldering or ember cast accumulation.

5. **Parallel execution**
   - GPU kernels are implemented via AMReX `ParallelFor`/`AMREX_GPU_DEVICE` macros and will execute on GPU when built with `LEVELSET_GPU_BACKEND=CUDA/HIP/SYCL`. Some serial CPU fallback paths (e.g., spotting, landscape I/O) remain.
   - MPI domain decomposition is provided by AMReX when built with `LEVELSET_ENABLE_MPI=ON`.

6. **Validation**
   - Limited validation against real wildfire events; parameters are calibrated to standard fuel models.

## Comparison with FARSITE and WRF-Fire (WRF-SFIRE)

> **Note on naming**: WRF-Fire and WRF-SFIRE refer to the same coupled fire–atmosphere system built on the WRF mesoscale model. The official package name is WRF-SFIRE; "WRF-Fire" is a common shorthand used in the community.

| Capability | wildfire_levelset | FARSITE | WRF-SFIRE |
|---|---|---|---|
| **Fire spread model** | Rothermel (1972) ROS with Andrews (2018) wind adjustments (WAF, MEWS); Balbi (2009) physical model; Cheney & Gould (1995) grassland empirical model; elliptical directional spread (Richards 1990); Eulerian level-set implementation of Huygens wavelet principle | Rothermel ROS + explicit Huygens wavelet | Rothermel ROS (level-set) |
| **Wind adjustment** | Optional WAF (20-ft → midflame) and MEWS cap via `rothermel.use_waf` / `rothermel.use_wind_limit`; wind-terrain feedback models (`wind_terrain.model`) with 6 options including Rothermel (1983) canyon wind, Viegas & Neto (1994) buoyancy wind, and Pimont et al. (2009) slope correction | WAF applied internally | WAF not explicitly applied; wind from atmospheric model |
| **Terrain representation** | Full 2-D landscape: per-cell elevation, slope, aspect, and fuel model from FARSITE `.lcp` files or terrain XYZ files; 3-D terrain and canopy layer can be added in future | Full 2-D landscape (.lcp) | Full 3-D terrain from WRF grid |
| **Wind coupling** | One-way (prescribed; WRF output converted to CSV via `wrf_wind_reader.py`) | One-way (prescribed; gridded wind) | Two-way (fire ↔ atmosphere) |
| **Atmospheric model** | None | None | Full WRF mesoscale atmosphere |
| **Spatial resolution** | User-defined (AMReX grid) | Landscape raster resolution | WRF grid resolution |
| **Fuel models** | Anderson 13 + Scott & Burgan 40 | Anderson 13 + Scott & Burgan 40 | Anderson 13 |
| **Fuel moisture** | Per size class (1-hr, 10-hr, 100-hr dead; live herbaceous and live woody) — user-specified constants; no dynamic transport | Dead/live moisture per size class | Dead/live moisture (prescribed) |
| **Crown fire** | Van Wagner (1977) initiation | Van Wagner (1977) initiation | Van Wagner (1977) initiation |
| **Firebrand spotting** | Stochastic + Albini (1983) 2-D trajectory | Albini empirical spotting | Not included by default |
| **Bulk fuel burnout** | Residence-time model | Full burnout tracking per cell | Post-frontal fuel consumption |
| **Multi-ignition** | CSV ignition file | Interactive ignition map | WPS/WRF fire restart |
| **Fire behavior output** | `R` [m/s], fireline intensity [kW/m], flame length [m] in plotfiles | Spread rate, intensity in output files | Spread rate in WRF output |
| **GPU acceleration** | Yes (AMReX CUDA/HIP/SYCL kernels) | No | No |
| **MPI parallelism** | Yes (AMReX domain decomposition) | No | Yes (WRF MPI) |
| **Embedded boundaries** | Yes (AMReX EB) | No | No |
| **Operational use** | Research / prototype | Operational (US Forest Service) | Research / operational |
| **Post-processing** | AMReX plotfiles (ParaView) + GeoTIFF/GeoJSON via `plotfile_to_geotiff.py` | GIS shapefiles / rasters | WRF standard output (NCO/Python) |

### Key differences from FARSITE

- **Eulerian level-set vs. explicit Huygens wavelets**: This solver implements the same elliptical directional spread (Richards 1990) as FARSITE, but embeds it in an Eulerian level-set framework.  The fire perimeter is tracked as the zero contour of a smooth signed-distance function rather than a chain of individually advected vertex wavelets.  The Eulerian approach handles complex topology changes (merging fronts, islands) automatically without explicit connectivity management.
- **Andrews (2018) wind adjustments are optional**: FARSITE applies the Wind Adjustment Factor (WAF) and maximum effective wind speed (MEWS) cap internally. This solver exposes both via `rothermel.use_waf` and `rothermel.use_wind_limit` so users can enable them when wind input is at 20-ft (NWP or WRF) height rather than midflame height.
- **Full 2-D landscape terrain**: Both solvers support per-cell elevation, slope, aspect, and fuel model from FARSITE `.lcp` landscape files (LANDFIRE data). Canopy fuel layers can be added to this solver in future without architectural changes.
- **Per size-class fuel moisture**: Dead fuel moisture is specified independently for 1-hr, 10-hr, and 100-hr size classes; live herbaceous and live woody moisture are specified separately.  This mirrors FARSITE's moisture inputs and drives the multi-class Rothermel (1972) reaction intensity calculation.
- **Embedded Boundary**: AMReX EB allows irregular obstacles (buildings, fuel breaks) on the Cartesian grid without remeshing.
- **GPU and MPI**: The AMReX backend supports CUDA/HIP/SYCL and MPI for large domains; FARSITE is serial and CPU-only.
- **GIS output**: `tools/plotfile_to_geotiff.py` converts plotfiles directly to GeoTIFF rasters and GeoJSON fire-perimeter contours; FARSITE produces GIS shapefiles/rasters natively.
- **No interactive GUI**: This solver is driven by an input text file and command-line arguments; FARSITE ships with a Windows GUI.

### Key differences from WRF-Fire (WRF-SFIRE)

- **No atmospheric coupling**: This solver uses prescribed wind fields (constant, CSV, or WRF output); WRF-SFIRE fully couples the fire with the WRF atmospheric model, including fire-induced wind, heat flux, and smoke transport.
- **Wind adjustment control**: When using WRF output wind via `wrf_wind_reader.py`, enable `rothermel.use_waf = 1` to convert the NWP wind to midflame height. WRF-SFIRE handles this internally as part of the coupled framework.
- **Albini spotting vs. none**: WRF-SFIRE does not include firebrand spotting by default.
- **Simpler setup**: Running this solver requires only CMake and an AMReX build; WRF-SFIRE requires a full WRF stack (NetCDF, MPI, WPS, WRF pre-processing).
- **GPU-native kernels**: WRF-SFIRE is primarily CPU-MPI; this solver's hotspot loops use AMReX GPU kernels.

## References

- **Andrews, P.L. (2018)**. *The Rothermel Surface Fire Spread Model and Associated Developments: A Comprehensive Explanation.* Gen. Tech. Rep. RMRS-GTR-371. USDA Forest Service. <https://doi.org/10.2737/RMRS-GTR-371>
- **Byram, G.M. (1959)**. "Combustion of forest fuels." In: Davis, K.P. (ed.) *Forest Fire: Control and Use*. McGraw-Hill, New York. pp. 61–89.
- **Richards, G.D. (1990)**. "An elliptical growth model of forest fire fronts and its numerical solution." International Journal of Numerical Methods in Engineering.
- **Anderson, H.E. (1983)**. "Predicting wind-driven wild land fire size and shape." Research Paper INT-305, USDA Forest Service.
- **Anderson, H.E. (1982)**. "Aids to determining fuel models for estimating fire behavior." Research Paper INT-122, USDA Forest Service.
- **Rothermel, R.C. (1972)**. "A mathematical model for predicting fire spread in wildland fuels." Research Paper INT-115, USDA Forest Service.
- **Rothermel, R.C. (1983)**. *How to Predict the Spread and Intensity of Forest and Range Fires.* Gen. Tech. Rep. INT-143, USDA Forest Service.
- **Van Wagner, C.E. (1977)**. "Conditions for the start and spread of crown fire." Canadian Journal of Forest Research.
- **Albini, F.A. (1983)**. "Potential spotting distance from wind-driven surface fires." Research Paper INT-309, USDA Forest Service.
- **Scott, J.H. and Burgan, R.E. (2005)**. "Standard Fire Behavior Fuel Models: A Comprehensive Set for Use with Rothermel's Surface Fire Spread Model." Technical Report RMRS-GTR-153, USDA Forest Service.
- **Viegas, D.X. & Neto, L.P.S. (1994)**. "Wind tunnel study of the convective air flow of slope fires." Annual Report, Project COMOESTAS.
- **Viegas, D.X. (2004)**. "Slope and wind effects on fire propagation." *International Journal of Wildland Fire*, 13(2), 143–156. <https://doi.org/10.1071/WF03046>
- **Pimont, F., Dupuy, J.-L., Linn, R.R. & Dupont, S. (2009)**. "Validation of FIRETEC wind-flows over a canopy and fuel-break." *International Journal of Wildland Fire*, 18(7), 775–790. <https://doi.org/10.1071/WF07130>

## License

See [LICENSE](LICENSE) file.

## Contact

For questions or issues, please open an issue on the GitHub repository.
