> ⚠️ **AI-Generated Code Disclaimer**
> A major portion of this codebase was written with the assistance of AI tools (GitHub Copilot / large language models).
> It has not been exhaustively validated against operational wildfire prediction systems.
> **Use with caution** — review all outputs carefully before applying to real-world fire management decisions.

# Wildfire-AMR

A unified AMReX-based C++ wildfire front propagation framework providing a single interface to operational fire behaviour tools — FARSITE, BehavePlus-style spread models, and physics-based alternatives — with a path toward future two-way coupling with the Energy Research and Forecasting (ERF) atmospheric model.

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
- **Rothermel fire spread model** with Anderson 13 (FM1–FM13) and Scott & Burgan 40 fuel databases
- **Andrews (2018) wind adjustments**: Wind Adjustment Factor (WAF) and Maximum Effective Wind Speed (MEWS) cap
- **Balbi (2009) physical fire spread model** – radiation-driven, physics-based rate of spread
- **Cheney & Gould (1995/1998) grassland fire spread model** – empirical Australian grassland calibration
- **FARSITE elliptical expansion** (Richards 1990) with Anderson L/W ratio
- **Alternative fire shape models**: Catchpole & de Mestre (1986) double-ellipse, Wilson (1988) single ellipse, Alexander et al. lemniscate
- **Minimum Travel Time (MTT) propagation** via Dijkstra fast-marching sweep
- **Barrier polygon / firebreak CSV files** for roads, fuel breaks, and firebreaks
- **Terrain effects** – slope and aspect corrections via constants, terrain files, or FARSITE landscape files
- **Binary and ASCII FARSITE landscape files** (`.lcp`) with spatially-varying crown fuel layers
- **Fuel adjustment file (`.adj`)** – per-fuel-model ROS multipliers
- **Time-varying fuel moisture (`.fmd`)** – diurnal or scenario-driven moisture schedules
- **Per size-class fuel moisture**: 1-hr, 10-hr, 100-hr dead; live herbaceous and live woody
- **Fire behavior diagnostics**: Byram fireline intensity and flame length
- **Fire perimeter CSV / GeoJSON writer** and burned area / perimeter time series (`fire_stats.csv`)
- **Weise & Biging (1996) fire whirl model** – optional diagnostic for flame tilt and whirl kinematics
- **Viegas (2004) eruptive fire diagnostics** – slope enhancement factor and eruptive-regime classification
- **Wind-terrain feedback models** – seven options including Viegas ROS, canyon wind, and WindNinja ridge/canyon speed-up
- **Turbulent wind perturbation** – Ornstein-Uhlenbeck process, spectral (RFF) noise, and direction random walk
- **Heat flux MultiFab** – fire-induced plume buoyancy and horizontal inflow corrections
- **GIS output**: plotfile → GeoTIFF rasters and GeoJSON perimeter contours
- **Ensemble burn probability** – FSPro-style LHS/random driver with optional MPI parallelism
- **Stochastic and physics-based (Albini 1983) firebrand spotting**
- **Crown fire initiation** (Van Wagner 1977) and **bulk fuel consumption**
- **Multi-point ignition** from a CSV file
- **2D and 3D** simulation capabilities; **Embedded Boundary (EB)** support
- **AMReX-based** with GPU-ready kernels (CUDA/HIP/SYCL) and optional MPI parallelism

See the 📚 [full documentation](https://hgopalan.github.io/wildfire_levelset/) for detailed model descriptions, parameter references, and worked examples.

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

See the 📚 [full documentation](https://hgopalan.github.io/wildfire_levelset/) for the complete runtime parameter reference, fire spread model equations, and tool usage guides.

## Tools

Python utilities live in the `tools/` directory for terrain/landscape download, weather file parsing, GIS export, and ensemble analysis.
See the 📚 [full documentation](https://hgopalan.github.io/wildfire_levelset/tools.html) for complete option references and worked examples.

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

Simulation results are written as AMReX plotfiles (`plt####`) containing the level-set function, velocity field, rate of spread, fireline intensity, terrain fields, and optional diagnostics (spotting, crown fire, fire whirl, eruptive fire, emissions). Additional files include burned-cell coordinates, convex hull envelopes, fire perimeter CSV/GeoJSON files, and a `fire_stats.csv` time series of burned area, perimeter, and emissions.

View plotfiles with ParaView or use `tools/plotfile_to_geotiff.py` to export GeoTIFF rasters and GeoJSON contours for GIS. See the 📚 [full documentation](https://hgopalan.github.io/wildfire_levelset/usage.html#output) for a complete field listing.

## Limitations and Known Constraints

1. **Fuel modeling** — Anderson 13 and Scott & Burgan 40 fuel models only; no custom or dynamic fuel models. Fuel moisture is user-supplied per size class; no dynamic moisture transport. Spatially-varying fuel requires a FARSITE landscape file.

2. **Wind field** — Spatially-varying wind via CSV; fully 3-D wind fields are not supported. Time-dependent wind uses linear interpolation (2D only). Fire–atmosphere coupling is one-way by default; wind-terrain feedback options adjust effective wind but atmospheric feedback from the fire is not modelled.

3. **Terrain** — Per-cell slope/aspect from FARSITE `.lcp` or XYZ files. True 3-D terrain geometry for 3-D simulations is not yet supported.

4. **Physical sub-models** — Crown fire uses Van Wagner (1977) empirical criteria; no mechanistic crown spread. Firebrand transport is stochastic or Albini 2-D; no full 3-D plume. No radiation/convective heat transfer or post-frontal smouldering.

5. **Parallel execution** — GPU via AMReX `ParallelFor`; some paths (spotting, landscape I/O) have serial CPU fallbacks. MPI decomposition provided by AMReX.

6. **Validation** — Limited validation against real wildfire events; parameters calibrated to standard fuel models.

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
