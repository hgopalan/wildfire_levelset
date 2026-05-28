> ⚠️ **AI-Generated Code Disclaimer**
> A major portion of this codebase was written with the assistance of AI tools (GitHub Copilot / large language models).
> It has not been exhaustively validated against operational wildfire prediction systems.
> **Use with caution** — review all outputs carefully before applying to real-world fire management decisions.

# Wildfire-AMR

## CI / Build Status

| Configuration | Status |
|---------------|--------|
| Linux & macOS CPU (GCC/Clang, Release + Debug) | [![CMake Build](https://github.com/hgopalan/wildfire_levelset/actions/workflows/cmake_build.yml/badge.svg)](https://github.com/hgopalan/wildfire_levelset/actions/workflows/cmake_build.yml) |
| Windows CPU (MSVC, Release + Debug) | [![CMake Build](https://github.com/hgopalan/wildfire_levelset/actions/workflows/cmake_build.yml/badge.svg)](https://github.com/hgopalan/wildfire_levelset/actions/workflows/cmake_build.yml) |
| Linux GPU — CUDA 12.6 | [![CMake Build](https://github.com/hgopalan/wildfire_levelset/actions/workflows/cmake_build.yml/badge.svg)](https://github.com/hgopalan/wildfire_levelset/actions/workflows/cmake_build.yml) |
| Linux GPU — HIP/ROCm 6.2 | [![CMake Build](https://github.com/hgopalan/wildfire_levelset/actions/workflows/cmake_build.yml/badge.svg)](https://github.com/hgopalan/wildfire_levelset/actions/workflows/cmake_build.yml) |
| Linux GPU — SYCL/oneAPI 2025.x | [![CMake Build](https://github.com/hgopalan/wildfire_levelset/actions/workflows/cmake_build.yml/badge.svg)](https://github.com/hgopalan/wildfire_levelset/actions/workflows/cmake_build.yml) |
| Documentation | [![Docs](https://github.com/hgopalan/wildfire_levelset/actions/workflows/docs.yml/badge.svg)](https://github.com/hgopalan/wildfire_levelset/actions/workflows/docs.yml) |

An AMReX-based C++ wildfire front propagation framework providing a unified interface to FARSITE-compatible fire behaviour models with GPU-ready kernels (CUDA/HIP/SYCL) and optional MPI parallelism.

## 📚 Full Documentation

**[Read the full documentation](https://hgopalan.github.io/wildfire_levelset/)**

The documentation covers mathematical models, input parameters, tools, worked examples, and a comparison with FARSITE, FlamMap, WRF-SFIRE, and other simulators.

### Implementation Details

For detailed information about specific features and implementations:

* **[FARSITE Features Implementation](https://hgopalan.github.io/wildfire_levelset/farsite_features_implementation.html)** - Fire intensity classification, spot fire ignition delay, Pasquill-Gifford atmospheric stability, and FARSITE temporal fire acceleration model
* **[Fire Acceleration Model](https://hgopalan.github.io/wildfire_levelset/implementation_summary.html)** - Full FARSITE temporal acceleration model with wind-onset time-lag capability
* **[Cell Size Correction](https://hgopalan.github.io/wildfire_levelset/cell_size_correction.html)** - Empirical correction factors for grid resolution effects on fire spread rates
* **[Spatial T/RH Interpolation](https://hgopalan.github.io/wildfire_levelset/implementation_spatial_trh.html)** - Multi-station weather interpolation for spatially-varying fuel moisture

## Quick Start

```bash
git clone --recurse-submodules https://github.com/hgopalan/wildfire_levelset.git
cd wildfire_levelset
cmake -S . -B build -DLEVELSET_DIM_2D=ON
cmake --build build -j
./build/levelset regtest/surface_spread/farsite_ellipse/inputs.i
```

### Key Build Options

| Option | Default | Description |
|--------|---------|-------------|
| `LEVELSET_DIM_2D` | `OFF` | 2D build |
| `LEVELSET_ENABLE_EB` | `OFF` | Embedded Boundary support |
| `LEVELSET_ENABLE_MPI` | `OFF` | MPI parallelism |
| `LEVELSET_GPU_BACKEND` | `NONE` | `CUDA`, `HIP`, or `SYCL` |
| `LEVELSET_BUILD_PYTHON_BINDINGS` | `OFF` | Python API (pyWildfire) |
| `LEVELSET_BUILD_DOCS` | `OFF` | Sphinx documentation |

See the [full build guide](https://hgopalan.github.io/wildfire_levelset/building.html) for GPU, MPI, and Windows instructions.

## Python API for Coupled Simulations

The Python API (`pyWildfire`) enables **full programmatic control** of the fire solver from Python, designed for coupled wind-fire simulations with external wind solvers like [massconsistent_amr](https://github.com/hgopalan/massconsistent_amr).

### Quick Start with Python API

Build with Python bindings:
```bash
cmake -S . -B build -DLEVELSET_DIM_2D=ON -DLEVELSET_BUILD_PYTHON_BINDINGS=ON
cmake --build build -j
export PYTHONPATH=$PWD/build/python:$PYTHONPATH
```

Run a fire simulation from Python:
```python
from wildfire_solver import WildfireSolver

# Initialize and run
fire = WildfireSolver("inputs.i")
for i in range(100):
    fire.step()
    state = fire.get_state()
    print(f"t={state['time']:.1f}s, burned={sum(state['phi']<=0)*fire.dx*fire.dy:.0f}m²")
fire.finalize()
```

### Coupled Wind-Fire Simulation

Integrate with external 3D wind solvers:
```python
from wildfire_solver import WildfireSolver
# from pyWindSolver import WindSolver  # Example: massconsistent_amr

fire = WildfireSolver("fire_inputs.i")
# wind = WindSolver("wind_inputs.txt")  # External wind solver

while fire.time < final_time:
    # 1. Solve/update 3D wind field
    # u_3d, v_3d, w_3d = wind.get_velocity_arrays()
    u_3d, v_3d, w_3d = generate_wind_field(...)  # Your wind solver
    
    # 2. Pass wind to fire solver
    fire.update_wind_3d(u_3d, v_3d, w_3d, nz, zmin, zmax)
    
    # 3. Advance fire
    fire.step()

fire.finalize()
```

**Key Features:**
- Complete fire solver state management (initialization, time-stepping, finalization)
- Extract all fire fields as NumPy arrays (`phi`, `ros`, `intensity`, `flame_length`, etc.)
- Update wind from 2D or 3D fields during simulation
- Write AMReX plotfiles from Python
- Zero-copy data transfer between C++ and Python

**Applications:**
- Two-way coupled atmosphere-fire simulations
- Ensemble runs with varying wind scenarios
- Machine learning training data generation
- Custom fire-weather coupling strategies
- Integration with WRF, WRF-Fire, or other atmospheric models

**Documentation:**
- [Python API Guide](https://hgopalan.github.io/wildfire_levelset/python_api.html) - Complete API reference and examples
- [PYTHON_API_IMPLEMENTATION.md](PYTHON_API_IMPLEMENTATION.md) - Implementation details
- [Python API Regression Tests](regtest/python_api/README.md) - How to run tests and integrate wind solvers
- [massconsistent_amr](https://github.com/hgopalan/massconsistent_amr) - Compatible 3D wind solver

## Core Capabilities

| Category | Feature | Key Reference |
|----------|---------|---------------|
| **Fire Spread Models** | Rothermel (1972) surface fire with Anderson 13 & Scott-Burgan 40 fuel databases | `rothermel.*` |
| | FARSITE elliptical expansion (Richards 1990) + Anderson L/W ratio | `farsite.*` |
| | Alternative models: Balbi (2009), Cheney-Gould (1995), Cruz-Alexander-Wakimoto (2005) | `balbi.*`, `cheney.*`, `cruz.*` |
| | Canadian FBP System (O1a/O1b grass, S1/S2/S3 slash) | `fbp.*` |
| | Lautenberger (2013) physics-based spread | `lautenberger.*` |
| **Fire Behavior** | Wind Adjustment Factor: Andrews logarithmic or BehavePlus linear | `rothermel.waf_formula` |
| | Crown fire: Van Wagner (1977) + Cruz et al. (2005) + Rothermel (1991) | `crown.*` |
| | Per-fuel burnout time (Rothermel 1983 residence time) | Auto from landscape |
| | Burn-period daytime window (FARSITE/FSPro concept) | `burn_period.enable` |
| **Spotting & Embers** | Firebrand spotting: Albini (1983) trajectory + torching-tree (1979) | `firebrand.*` |
| | Vorticity-enhanced spotting: Weise & Biging (1996) fire whirl effects | `weise_biging.enhance_spotting` |
| | Ember cascade: Gaussian flux field (Sardoy 2007 approach) | `ember_cascade.*` |
| | GPU-accelerated 3-D wind interpolation (CUDA/HIP/SYCL) | Optional from massconsistent_amr |
| | Post-fire fuel adjustment for re-entry spots | `fuel_depletion.adjust_spotting_reentry` |
| **Terrain & Weather** | Slope/aspect from FARSITE LCP or XYZ; 8-direction horizon scan | `terrain.*` |
| | Time-varying, turbulent (OU/spectral), direction-schedule winds | `wind.*` |
| | Multi-station weather with IDW interpolation | `multi_wtr_file` |
| **Fuel Moisture** | FMD schedule, diurnal (Nelson 2000), precipitation wetting, FMC phenology | `fuel_moisture.*` |
| | Spatial moisture from FARSITE .fms files | `.fms` support |
| | Moisture fields (d1/d10/d100/lh/lw) in every plotfile | Auto output |
| **Ignition & Suppression** | Point CSV, polygon, polyline ignitions | `ignition.*` |
| | Retardant drop zones (zero ROS & spotting) | `retardant.*` |
| | Satellite fire detection (GOES/VIIRS) assimilation | `satellite.*` |
| **Diagnostics & Output** | Fire ecology: scorch height, tree mortality, TI/CI ratios | `ecology.*` |
| | Smoke plume rise (Briggs 1965/1969) | `smoke_plume.enable` |
| | Vorticity field: vertical component for fire whirl identification | Auto in plotfiles |
| | Fire whirl characteristics: Weise & Biging (1996) model | `weise_biging.enable` |
| | Fire emissions: CO₂, CO, PM₂.₅ (WRF-Fire) | Auto output |
| | KML perimeter export (UTM → WGS-84) | `write_perimeter_kml` |
| | FARSITE .fsa/.pst files | `fsa_file`, `pst_file` |
| | Flame length exceedance raster | `fl_exceedance` in `plot_vars` |
| | Post-frontal fuel consumption raster | Auto in plotfiles |
| | Simulation date/time display | `sim_datetime.*` |
| | Fire line intensity classification (Byram I–VI) in HTML report | `fire_report_file` |
| **Technical** | AMReX-based: GPU kernels (CUDA/HIP/SYCL), AMR, MPI | Build options |
| | MTT (Minimum Travel Time) Dijkstra propagation | `mtt.*` |
| | Non-burnable masking (codes 91–99 / NB1–NB9) | Auto from fuel codes |
| | ROS stall threshold (FARSITE-compatible 1×10⁻⁴ m/s) | Built-in |
| | Vectorial slope/wind combination | `rothermel.use_slope_wind_vectors` |
| | Conditional weather ERC/BI/SC trigger | `conditional_weather_trigger` |

## Ensemble / FSim-Style Probabilistic Simulation

`tools/ensemble_burn_probability.py` implements an FSim-style Monte Carlo
driver.  Beyond the original wind-speed / direction / moisture perturbations it
now also supports:

| Feature | Flag | Description |
|---|---|---|
| **Probabilistic ignition locations** | `--ignition-prob-csv` | Weighted random draw of ignition centre from a GIS raster (CSV with `x_m`, `y_m`, `weight`). |
| **Containment probability** | `--containment-prob P` `--suppression-file F` | Bernoulli draw per run: enables/disables suppression schedule with probability P. |
| **Crown fire probability map** | `--crown-fire-prob` `--crown-fire-out F` | Accumulates P(crown fire) per cell from the `crown_fraction` plotfile field. |
| **Burn probability map** | `--out` | P_burn per cell (existing). |
| **Flame-length exceedance** | `--fl-thresholds` | P(FL > M) per cell for each threshold (existing). |
| **Area exceedance CCDF** | `--area-exceedance` | P(burned area ≥ A) curve across runs (existing). |

Quick example combining all new features:
```bash
python3 tools/ensemble_burn_probability.py \
    --exe ./wildfire_levelset \
    --inputs my_scenario/inputs.i \
    --n-runs 200 \
    --ignition-prob-csv ignition_risk.csv \
    --containment-prob 0.35 \
    --suppression-file suppression_lines.csv \
    --crown-fire-prob \
    --area-exceedance \
    --fl-thresholds 1.0 2.0 4.0 \
    --seed 42
```

## Regression Tests

Tests are organised into sub-folders under `regtest/`, all using UTM Zone 11N coordinates
(Southern California reference: 330000 E, 3775000 N):

| Sub-folder | Tests |
|---|---|
| `surface_spread/` | basic_levelset, farsite_ellipse, rothermel_fuel, anderson_lw, catchpole_demestre, wilson_spread, alexander_lemniscate, ellipse_sdf, reinitialization, fbp_o1a_grassfire, fbp_s1_slash, lautenberger_spread |
| `crown_fire/` | crown_initiation, cruz_crown_continental_us, fmc_seasonal, rothermel1991_crown |
| `spotting/` | firebrand_spotting, albini_spotting, **albini_spotting_3d_wind** *(new)*, **ember_cascade_flux** *(new)*, **vorticity_enhanced_spotting** *(new)* |
| `terrain/` | terrain_wind, balbi_viegas_heatflux, windninja_ridge_canyon, **solar_horizon_shading** *(new)*, **terrain_gradient_correction** *(new)* |
| `moisture/` | fmd_moisture, cheney_gould_grassfire, precip_wetting, spatial_moisture_output, **wtr_diurnal** *(new)*, **wtr_rain_wetting** *(new)* |
| `fuel/` | fuel_adj_file |
| `ignition/` | barrier_polygons, **satellite_assimilation** *(new)*, fire_perimeter_output, polygon_ignition, polyline_ignition, **kml_perimeter** *(new)* |
| `wind/` | time_dependent_wind, turb_wind, wind_dir_schedule, **waf_andrews**, **waf_behaviorplus** |
| `diagnostics/` | scott_reinhardt_indices, scott_reinhardt_full_ti_ci, **fl_exceedance** *(new)*, **farsite_fsa_pst** *(new)*, **conditional_weather_bi** *(new)* |
| `misc/` | 3d_sphere, eb_implicit, mtt_propagation, bulk_fuel_consumption, landfire_farsite, **nonburnable_mask**, **smoke_plume_rise** *(new)*, **reentry_spotting** *(new)* |
| `python_api/` | **basic_fire_solver**, **coupled_wind_fire** *(Python bindings)* |

The `albini_spotting_3d_wind` and `ember_cascade_flux` tests require a Python pre-step to generate the synthetic plt wind file:

```bash
cd regtest/spotting/albini_spotting_3d_wind
python3 generate_plt_wind.py   # creates plt_wind_3d/ directory

cd ../ember_cascade_flux
python3 generate_plt_wind.py   # creates plt_wind_3d/ directory
```

```bash
cmake -S . -B build -DLEVELSET_DIM_2D=ON
cmake --build build -j
cd build && ctest -L regtest --output-on-failure
```

## Tools

Python utilities in `tools/` for terrain download, weather parsing, GIS export, ensemble analysis, and post-processing.

| Category | Tool | Description |
|----------|------|-------------|
| **Ensemble** | `ensemble_burn_probability.py` | FSim-style Monte Carlo driver |
| | `plot_burn_probability.py` | Visualise burn probability maps |
| | `values_at_risk.py` | FSPro-style values-at-risk overlay |
| **Fire Analysis** | `fire_size_summary.py` | Fire statistics with **percentile analysis** |
| | `fire_period_analysis.py` | **Day/night burn classification** |
| | `minimum_travel_path.py` | **MTT path extraction** |
| | `isochrone_extractor.py` | Arrival-time isochrones with **visualization** |
| | `farsite_fsa_pst_reader.py` | FARSITE .fsa/.pst file reader |
| **GIS/Export** | `plotfile_to_geotiff.py` | AMReX plotfiles → GeoTIFF/GeoJSON |
| | `perimeter_to_shapefile.py` | Fire perimeter → Esri Shapefile |
| **Input Prep** | `landscape_writer.py` | LANDFIRE → FARSITE .lcp |
| | `srtm_terrain_reader.py` | SRTM elevation → terrain CSV |
| | `farsite_weather_reader.py` | FARSITE .wtr → wind CSV |
| | `historical_wildfires.py` | 29 major US fires (2009–2024) database |
| **Worksheets** | `surface_fire_worksheet.py` | BehavePlus-style surface fire |
| | `crown_fire_worksheet.py` | Van Wagner crown fire |
| | `ignition_probability_table.py` | Anderson P_ignition |
| | `behavior_matrix.py` | Rothermel behavior matrices |

**New FARSITE-parity features:** `fire_size_summary.py` (percentile statistics), `isochrone_extractor.py` (visualization with time labels), `minimum_travel_path.py` (MTT path extraction), `fire_period_analysis.py` (day/night burn classification). See [tools documentation](https://hgopalan.github.io/wildfire_levelset/tools.html) and `tools/NEW_FARSITE_FEATURES.md`.

## References

- Albini, F.A. (1979). Spot fire distance from burning trees — a predictive model. USDA For. Serv. Gen. Tech. Rep. INT-56.
- Albini, F.A. (1983). Potential spotting distance from wind-driven surface fires. USDA For. Serv. Res. Pap. INT-309.
- Andrews, P.L. (2018). The Rothermel surface fire spread model and associated developments. USDA For. Serv. Gen. Tech. Rep. RMRS-GTR-371.
- Briggs, G.A. (1965). A plume rise model compared with observations. JAPCA 15(9):433–438.
- Byram, G.M. (1959). Combustion of forest fuels. In: Davis, K.P. (ed.) Forest Fire: Control and Use. McGraw-Hill.
- Cheney, N.P. & Gould, J.S. (1995). Fire growth in grassland fuels. Int. J. Wildland Fire 5(4):237–247.
- Cruz, M.G., Alexander, M.E. & Wakimoto, R.H. (2005). Development and testing of models for predicting crown fire rate of spread in conifer forest stands. Can. J. For. Res. 35:1626–1639.
- Finney, M.A. (2004). FARSITE: Fire Area Simulator — Model Development and Evaluation. USDA For. Serv. Res. Pap. RMRS-RP-4 Revised.
- Lautenberger, C. (2013). Wildland fire modeling with an Eulerian level set method and automated calibration. Fire Safety J. 62:477–485.
- Nelson, R.M. Jr. (2000). Prediction of diurnal change in 10-h fuel stick moisture content. Can. J. For. Res. 30:1071–1087.
- Richards, G.D. (1990). An elliptical growth model of forest fire fronts and its numerical solution. Int. J. Numer. Meth. Eng. 30:1163–1179.
- Rothermel, R.C. (1972). A mathematical model for predicting fire spread in wildland fuels. USDA For. Serv. Res. Pap. INT-115.
- Rothermel, R.C. (1991). Predicting behavior and size of crown fires in the Northern Rocky Mountains. USDA For. Serv. Res. Pap. INT-438.
- Scott, J.H. & Reinhardt, E.D. (2001). Assessing crown fire potential by linking models of surface and crown fire behavior. USDA For. Serv. Res. Pap. RMRS-RP-29.
- Van Wagner, C.E. (1977). Conditions for the start and spread of crown fire. Can. J. For. Res. 7:23–34.

## Satellite Fire Detection Assimilation

Active-fire detections from GOES-16/17/18 (public NOAA AWS S3), VIIRS (NASA FIRMS REST API),
or a pre-downloaded CSV can be ingested to constrain the initial fire perimeter or to correct
the simulated perimeter during a running simulation.  Enable with `satellite.enable = 1` and
set `satellite.source` to `file`, `viirs`, or `goes`.  See the
[Python Tools documentation](docs/tools.rst) for the full parameter reference and usage examples.

## License

See [LICENSE](LICENSE) file.
