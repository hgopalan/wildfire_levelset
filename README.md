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
- **Alternative fire shape models** for the FARSITE propagation path (`farsite.fire_shape_model`): Catchpole & de Mestre (1986) true double-ellipse, Wilson (1988) single ellipse from rear focus, and Alexander et al. lemniscate (Limaçon) — see [mathematical models documentation](https://hgopalan.github.io/wildfire_levelset/mathematical_models.html#alternative-fire-shape-models)
- **Terrain effects** including slope and aspect corrections via constant values, terrain files, or FARSITE landscape files
- **Binary and ASCII FARSITE landscape files** (`.lcp`) — the solver auto-detects binary vs ASCII format: binary `.lcp` files (FARSITE 4 format) are read natively, including **spatially-varying crown fuel layers** (CBH, CBD, canopy cover) when present (`CrownFuels != 0`); ASCII `X Y ELEV SLOPE ASPECT FUEL` format is also supported; all formats provide per-cell elevation, slope, aspect, and fuel model (landscape file takes precedence over terrain file or constant values)
- **Spatial crown fuel layers** from binary LCP — per-cell crown base height (CBH), crown bulk density (CBD), and canopy cover are populated from the landscape file when `use_spatial_crown = 1` (default); written to every plotfile as `cbh`, `cbd`, and `canopy_cover` fields; the global `crown.*` scalars serve as fallbacks when no crown layers are present
- **Fuel adjustment file (`.adj`)** — per-fuel-model rate-of-spread multipliers loaded from a FARSITE `.adj` file (`fuel_adj_file`); applied once at startup to the Rothermel lookup table; compatible with any landscape fuel type
- **Time-varying fuel moisture (`.fmd`)** — diurnal or scenario-driven moisture schedules loaded from a FARSITE-compatible `.fmd` file (`fmd_file`); the solver interpolates linearly between timestamps and updates `M_d1`, `M_d10`, `M_d100`, `M_lh`, `M_lw` at every time step; rebuilds the Rothermel / Balbi lookup tables on each update
- **Per size-class fuel moisture**: 1-hr, 10-hr, 100-hr dead fuel; live herbaceous and live woody moisture inputs for the multi-class Rothermel (1972) reaction intensity
- **Fire behavior diagnostics**: Byram (1959) fireline intensity [kW/m] and flame length [m] written to every plotfile
- **Scott / Albini (1979) maximum spotting distance table** — built-in lookup table of maximum typical spotting distances [m] by fuel model code (FBFM13 and FBFM40), optionally scaled by ambient wind speed; printed as a diagnostic alongside each Albini spotting step; functions `get_max_spot_dist_m()` and `get_max_spot_dist_scaled_m()` available in `src/scott_spotting_table.H`
- **Fire perimeter CSV / GeoJSON writer** — at every plotfile write, the fire front boundary (cells with φ < 0 adjacent to φ ≥ 0) is extracted and written as `perimeter_NNNN.csv` (comma-separated X, Y; `write_perimeter_csv = 1`) and/or `perimeter_NNNN.geojson` (GeoJSON Polygon, `write_perimeter_geojson = 1`)
- **Burned area / perimeter time series** — `fire_stats.csv` is created at startup and a row is appended at every plotfile step with: `step`, `time_s`, `burned_area_ha`, `perimeter_km`, `active_front_cells`, mean CO₂/CO/PM₂.₅ emissions (`fire_stats_file` parameter)
- **Weise & Biging (1996) fire whirl model** – optional diagnostic sub-model computing flame tilt, whirl height/radius, angular velocity, and tangential velocity from fireline intensity and wind; enabled via `weise_biging.enable = 1`
- **Viegas (2004) eruptive fire diagnostics** – optional parallel diagnostic that computes the Viegas exponential slope enhancement factor, eruptive-regime flag (critical slope), Viegas ROS, ROS excess vs primary model, and flame-tilt angle for hazard assessment; enabled via `viegas.enable = 1`; works with both Rothermel (R₀ baseline) and Balbi (amplitude A baseline); by default diagnostic-only, but can be coupled to fire propagation via the **wind-terrain feedback model** (see below)
- **Wind-terrain feedback models** – seven selectable options (`wind_terrain.model`) for how terrain-induced or fire-feedback winds modify fire spread: (1) default, no modification; (2) Viegas ROS as actual spread rate (works with both Rothermel and Balbi); (3) Viegas-induced buoyancy wind as velocity perturbation (eruptive cells only); (4) Rothermel (1983) canyon wind amplification; (5) Viegas & Neto (1994) buoyancy-driven upslope wind; (6) Pimont et al. (2009) exponential slope correction; (7) **WindNinja ridge/canyon empirical speed-up** based on wind-slope alignment: ridge acceleration when wind climbs upslope (`f = 1 + k_ridge × tan φ × alignment`), canyon channeling when wind descends (`f = 1 + k_canyon_wn × tan φ × |alignment|`)
- **Turbulent wind perturbation** (`turb_wind.model`) – stochastic wind variability for sensitivity studies and ensemble runs: **Ornstein-Uhlenbeck (OU) process** with user-defined temporal decorrelation time (`1/theta`) and perturbation amplitude (`sigma`); optional **Gaussian spatial correlation kernel** (`L_c` [m]) driving per-cell OU states via Gaussian-smoothed white noise; **Random Fourier Feature (RFF) spectral noise** (`spectral_noise`) combining OU temporal evolution with a physically correct Gaussian power spectrum — `N_modes` wavenumbers sampled at init, scalar OU amplitudes per mode evolved on CPU, field reconstructed on GPU as cosine superposition; or **direction random walk** for bounded direction fluctuations with preserved wind speed
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
| `LEVELSET_BUILD_WIND_SOLVER` | `ON` | Build the terrain-following mass-consistent wind solver (3-D only) |
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

## Terrain-Following Mass-Consistent Wind Solver (`wind_solver`)

A stand-alone 3-D terrain-following wind diagnostic executable is built
alongside `levelset` in default 3-D builds.  It is inspired by **QUIC-URB**
and the mass-consistent models of Sherman (1978) and Mathiesen (1987).

### What it does

1. Reads a terrain file (`X Y Z`, metres) — the horizontal domain extents are
   derived automatically from the data.
2. Constructs a log-law wind profile using the height above local terrain
   (interpolated column-by-column via IDW).
3. Solves the anisotropic Poisson equation for the Lagrange multiplier λ
   using **AMReX MLMG** (`MLABecLaplacian`) on a single-level Cartesian
   grid (level 0 only), enforcing ∇·**u** = 0.
4. Writes an AMReX plotfile with the corrected 3-D wind field and diagnostics.

### Quick start

```bash
# Prepare a terrain file (or use the supplied Gaussian hill example)
cp regtest/gaussian_hill_wind_solver/terrain.csv .

# Create inputs.i
cat > inputs.i << 'EOF'
terrain_file  = terrain.csv
U_ref         = 10.0      # 10 m/s westerly at z_ref
V_ref         = 0.0
z_ref         = 10.0      # reference height [m]
z0            = 0.1       # surface roughness [m]
dx            = 30.0
dy            = 30.0
dz            = 10.0      # finer vertical spacing for near-surface resolution
domain_height = 300.0
plot_file     = plt_wind
EOF

./build/wind_solver inputs.i
```

### Key input parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `terrain_file` | `terrain.csv` | Terrain point cloud (X Y Z) |
| `U_ref` | `10.0` | x-wind at reference height [m/s] |
| `V_ref` | `0.0` | y-wind at reference height [m/s] |
| `z_ref` | `10.0` | Reference height above local terrain [m] |
| `z0` | `0.1` | Surface roughness length [m] |
| `dx` | `30.0` | Grid spacing x [m] |
| `dy` | `30.0` | Grid spacing y [m] |
| `dz` | `30.0` | Grid spacing z [m] — reduce for vertical refinement |
| `domain_height` | `300.0` | Vertical domain extent [m] |
| `alpha_h` | `1.0` | Horizontal Lagrange anisotropy factor |
| `alpha_v` | `1.0` | Vertical Lagrange anisotropy factor |
| `plot_file` | `plt_wind` | Output plotfile prefix |

See the [full documentation](https://hgopalan.github.io/wildfire_levelset/wind_solver.html)
for the complete parameter reference, physical model description, and GIS export workflow.


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

A radiation-driven, physics-based ROS model that replaces Rothermel when `fire_spread_model = balbi`. The flame tilt angle is derived from wind speed and a buoyancy velocity, and the rate of spread follows from a radiation balance between the tilted flame and unburned fuel ahead. Thermal parameters are set via `balbi.*` ParmParse keys; fuel geometry is shared with Rothermel. The solver auto-generates a per-fuel Balbi parameter table at startup.

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

## Tools

Python utilities live in the `tools/` directory.
See the 📚 [full documentation](https://hgopalan.github.io/wildfire_levelset/tools.html) for
complete option references and worked examples.

| Script | Purpose |
|--------|---------|
| `wrf_wind_reader.py` | Extract U/V wind from a WRF netCDF file → solver wind CSV |
| `srtm_terrain_reader.py` | Download SRTM1 elevation data → UTM terrain XYZ file |
| `landscape_writer.py` | Download LANDFIRE rasters → FARSITE landscape file (`.lcp`) |
| `farsite_weather_reader.py` | Parse FARSITE `.wtr` RAWS weather files → time-stamped wind CSVs |
| `fuel_moisture_from_weather.py` | Estimate equilibrium dead-fuel moisture from temperature and RH |
| `plotfile_to_geotiff.py` | Convert AMReX plotfiles → GeoTIFF rasters and GeoJSON fire-perimeter contours |
| `farsite_adj_reader.py` | Inspect, generate, and apply FARSITE `.adj` fuel adjustment files |
| `farsite_fmd_reader.py` | Parse, convert, query, and generate FARSITE `.fmd` fuel moisture schedules |

### `farsite_adj_reader.py` — FARSITE fuel adjustment files

Reads per-fuel-model ROS multipliers from a FARSITE `.adj` file (``fuel_adj_file`` parameter), generates template files, and can apply adjustments to a Rothermel parameter CSV for inspection.

```bash
# Inspect an existing .adj file
python3 tools/farsite_adj_reader.py --adj fire.adj

# Generate a template for all FBFM13 models with adj_factor = 1.0
python3 tools/farsite_adj_reader.py --generate --fuel-system 13 \
    --output template_13.adj

# Generate a template for specific models
python3 tools/farsite_adj_reader.py --generate --models 4 9 10 \
    --output chaparral_adj.adj

# Apply adjustments to a Rothermel CSV for calibration inspection
python3 tools/farsite_adj_reader.py --adj fire.adj \
    --apply-csv rothermel_params.csv \
    --output rothermel_adjusted.csv
```

**Fuel adjustment file format** (`fire.adj`):
```
# Comments start with # or !
1
4  1.35    # FM4 chaparral: 35% ROS increase
10 0.85    # FM10: 15% ROS decrease
```

Integrate with the solver via `inputs.i`:
```
fuel_adj_file  = fire.adj
fuel_adj_model = 4        # global model code for single-model runs (no landscape)
```

### `farsite_fmd_reader.py` — FARSITE fuel moisture schedules

Reads, converts, queries, and generates FARSITE-compatible `.fmd` fuel moisture schedule files.  The solver reads these via `fmd_file` and updates moisture at every time step.

```bash
# Inspect a .fmd file
python3 tools/farsite_fmd_reader.py --fmd fire.fmd

# Convert to a flat CSV for analysis
python3 tools/farsite_fmd_reader.py --fmd fire.fmd --output moisture.csv

# Query interpolated moisture at t = 3600 s for FM4
python3 tools/farsite_fmd_reader.py --fmd fire.fmd \
    --query-time 3600 --fuel-model 4

# Generate a 24-hour constant-moisture template for FBFM13 models
python3 tools/farsite_fmd_reader.py --generate \
    --models 4 9 10 \
    --start-month 7 --start-day 15 --hours 24 \
    --M-d1 8 --M-d10 10 --M-d100 15 \
    --M-lh 90 --M-lw 120 \
    --output summer_moisture.fmd
```

**FMD file format** (`fire.fmd`):
```
# MONTH DAY HOUR PRECIP NUM_MODELS
7 15  800 0 2
4  8 10 15 90 120
9  9 11 16 88 115
7 15 1400 0 2
4 12 14 18 95 125
9 13 16 20 92 120
```

Integrate with the solver via `inputs.i`:
```
fmd_file        = summer_moisture.fmd
fmd_start_month = 7
fmd_start_day   = 15
fmd_fuel_model  = 0    # 0 = apply global schedule to all models
```

**Quick examples:**

```bash
# Download terrain and landscape data
python3 tools/srtm_terrain_reader.py \
    --lat-min 40 --lat-max 40.5 --lon-min -105 --lon-max -104.5

python3 tools/landscape_writer.py \
    --lat-min 40 --lat-max 40.5 --lon-min -105 --lon-max -104.5

# Extract wind from a WRF forecast
python3 tools/wrf_wind_reader.py --wrf-file wrfout_d01 --wind wind.csv

# Export simulation results for GIS
python3 tools/plotfile_to_geotiff.py --all --outdir gis_out \
    --utm-origin 450000 4200000 --epsg 32613
```

Use the generated files in a simulation:
```
rothermel.terrain_file    = terrain.xyz
rothermel.landscape_file  = landscape.lcp
velocity_file             = wind.csv
```

> **Deprecated tools**: `tools/deprecated/` contains `terrain_wind_preprocess.py` and related
> legacy scripts superseded by the split tools above.

## Testing

Run regression tests:
```bash
cd build
ctest -L regtest --output-on-failure
```

Run only the timing benchmark:
```bash
ctest -L benchmark --output-on-failure -V
```

Or use the custom target:
```bash
make regtest
```

Available regression tests:

**Core functionality**
- `basic_levelset` - Basic level-set advection
- `farsite_ellipse` - FARSITE elliptical expansion
- `rothermel_fuel` - Rothermel with fuel models
- `anderson_lw` - Anderson dynamic L/W ratio
- `cheney_gould_grassfire` - Cheney & Gould (1995/1998) grassland fire spread
- `catchpole_demestre` - Catchpole & de Mestre (1986) double-ellipse shape
- `wilson_spread` - Wilson (1988) single ellipse from rear focus
- `alexander_lemniscate` - Alexander et al. lemniscate (Limaçon) shape

**Advanced features (uncomment in `regtest/CMakeLists.txt` to enable)**
- `reinitialization` - Level-set reinitialization
- `ellipse_sdf` - Elliptical SDF initial conditions
- `eb_implicit` - Embedded boundary capabilities
- `firebrand_spotting` - Stochastic firebrand spotting model
- `albini_spotting` - Albini (1983) physics-based firebrand spotting
- `crown_initiation` - Crown fire initiation
- `bulk_fuel_consumption` - Fuel consumption modeling

**New feature tests (2025)**
- `fuel_adj_file` — FARSITE `.adj` fuel adjustment file reader (requires Python3)
- `fmd_moisture` — Time-varying `.fmd` fuel moisture schedule reader (requires Python3)
- `fire_perimeter_output` — CSV + GeoJSON perimeter writers and `fire_stats.csv` time series
- `fire_perimeter_output_validate` — Validates perimeter file contents after the solver run
- `burned_area_timeseries` — `fire_stats.csv` burned area / perimeter / emissions time series

**Dimension-specific tests**
- `3d_sphere` - Full 3D simulation (3D builds only)
- `terrain_wind` - External terrain and wind (2D builds only)
- `time_dependent_wind` - Time-varying wind fields (2D builds only)

**External data tests (requires Python3)**
- `landfire_farsite` - FARSITE with auto-downloaded LANDFIRE landscape (requires `landscape_writer.py`)

**Timing benchmark**
- `timing_benchmark` *(label: `benchmark`)* — Multi-resolution wall-clock timing benchmark; runs the solver at several grid sizes for both `levelset` and `farsite` scenarios; checks monotonically increasing runtime and estimates scaling exponent α; writes `timing_results.csv`. Run manually:

```bash
# From the build directory:
python3 ../regtest/timing_benchmark/run_benchmark.py \
    --exe ./levelset --dim 2 --nsteps 20 \
    --resolutions 32 64 128 256

# Dry-run: preview generated inputs without running solver
python3 ../regtest/timing_benchmark/run_benchmark.py \
    --exe ./levelset --dim 2 --dry-run

# 3D benchmark with custom resolutions
python3 ../regtest/timing_benchmark/run_benchmark.py \
    --exe ./levelset --dim 3 --nsteps 15 \
    --resolutions 16 32 48 64
```

The `timing_results.csv` written to the CTest working directory contains:
`scenario, n_cells, total_cells, wall_time_s, nsteps, steps_per_second, cells_per_step_per_s, returncode`

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
- `cbh` - Crown base height [m] (from binary LCP crown layers or global `crown.CBH`)
- `cbd` - Crown bulk density [kg/m³] (from binary LCP or global `crown.CBD`)
- `canopy_cover` - Canopy cover [%] (from binary LCP layer 4)
- `scorch_height` - Van Wagner (1973) scorch height [m] (always present)
- `prob_ignition` - Anderson (1970) probability of sustained ignition [-] (always present)
- `tree_mortality` - Ryan-Reinhardt (1988) style tree mortality fraction [-] (always present)
- `crown_activity` - Crown activity class: 0=surface, 1=passive, 2=active (always present)
- `co2_emissions` / `co_emissions` / `pm25_emissions` - Fuel emissions [kg/m²] (always present)
- `arrival_time` - Simulation time when each cell first burned [s]
- `heat_per_unit_area` - Heat per unit area [kJ/m²]
- `vorticity_z` - Vertical vorticity component [1/s]

**Additional output files** (written alongside each plotfile):
- `phi_negative_NNNN.dat` — x,y coordinates of all burned cells (φ < 0)
- `phi_envelope_NNNN.dat` — convex hull of the burned region
- `perimeter_NNNN.csv` — fire front boundary cells (X,Y columns), enabled via `write_perimeter_csv = 1`
- `perimeter_NNNN.geojson` — fire front as a GeoJSON Polygon, enabled via `write_perimeter_geojson = 1`
- `fire_stats.csv` — time series of burned area [ha], perimeter [km], active front cells, and mean emissions; appended at every plotfile step (configured via `fire_stats_file`)
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
- `scorch_height` - Van Wagner (1973) scorch height [m] (always present)
- `prob_ignition` - Anderson (1970) probability of sustained ignition [-] (always present)
- `tree_mortality` - Ryan-Reinhardt (1988) style tree mortality fraction [-] (always present)
- `crown_activity` - Crown activity class: 0=surface, 1=passive, 2=active (always present)
- `co2_emissions` / `co_emissions` / `pm25_emissions` - Fuel emissions [kg/m²] (always present)

Burned area [ha] and perimeter [km] statistics are printed to stdout at each plot interval.

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

## Comparison with Other Wildfire Simulation Tools

A detailed comparison of `wildfire_levelset` with FARSITE, WRF-SFIRE, FlamMap,
BehavePlus, FSPRO, QUIC-Fire, and FIRETEC — including a capability table and
key differences — is available in the documentation:

📄 **[Comparison with Other Wildfire Simulation Tools](https://hgopalan.github.io/wildfire_levelset/comparison.html)**

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
- **Van Wagner, C.E. (1973)**. "Height of crown scorch in forest fires." *Canadian Journal of Forest Research*, 3(3), 373–378.
- **Anderson, H.E. (1970)**. "Forest fuel ignitibility." *Fire Technology*, 6(4), 312–319.
- **Ryan, K.C. & Reinhardt, E.D. (1988)**. "Predicting postfire mortality of seven western conifers." *Canadian Journal of Forest Research*, 18(10), 1291–1297.
- **Scott, J.H. & Reinhardt, E.D. (2001)**. *Assessing Crown Fire Potential by Linking Models of Surface and Crown Fire Behavior.* Research Paper RMRS-RP-29, USDA Forest Service.
- **Nelson, R.M. Jr. (2000)**. "Prediction of diurnal change in 10-h fuel stick moisture content." *Canadian Journal of Forest Research*, 30(7), 1071–1087.
- **Seiler, W. & Crutzen, P.J. (1980)**. "Estimates of gross and net fluxes of carbon between the biosphere and the atmosphere from biomass burning." *Climatic Change*, 2(3), 207–247.
- **Finney, M.A. (2004)**. *FARSITE: Fire Area Simulator – Model Development and Evaluation.* Research Paper RMRS-RP-4, USDA Forest Service.
- **Finney, M.A. (2006)**. "An overview of FlamMap fire modeling capabilities." In: Andrews, P.L. & Butler, B.W. (comps.) *Fuels Management—How to Measure Success*. USDA Forest Service Proceedings RMRS-P-41. <https://www.firelab.org/project/flammap>
- **Andrews, P.L. (2014)**. "Current status and future needs of the BehavePlus Fire Modeling System." *International Journal of Wildland Fire*, 23(1), 21–33. <https://doi.org/10.1071/WF12167>
- **Finney, M.A., McHugh, C.W., Grenfell, I.C., Riley, K.L. & Short, K.C. (2011)**. "A simulation of probabilistic wildfire risk components for the continental United States." *Stochastic Environmental Research and Risk Assessment*, 25(7), 973–1000. <https://doi.org/10.1007/s00477-011-0462-z>
- **Linn, R.R. (1997)**. *A transport model for prediction of wildfire behavior.* PhD Dissertation / Technical Report LA-13334-T, Los Alamos National Laboratory.
- **Linn, R.R., Reisner, J., Colman, J.J. & Winterkamp, J. (2002)**. "Studying wildfire behavior using FIRETEC." *International Journal of Wildland Fire*, 11(4), 233–246. <https://doi.org/10.1071/WF02007>
- **Linn, R.R., Cunningham, P., Goodrick, S. & Mell, W. (2012)**. "Introduction to special issue on numerical simulation of wildland fire." *International Journal of Wildland Fire*, 21(6), i–ii.
- **Parsons, R.A., Mell, W.E. & McCauley, P. (2011)**. "Linking 3D spatial models of fuels and fire: effects of spatial heterogeneity on fire behavior." *Ecological Modelling*, 222(3), 679–691. <https://doi.org/10.1016/j.ecolmodel.2010.10.023>

## License

See [LICENSE](LICENSE) file.

## Contact

For questions or issues, please open an issue on the GitHub repository.
