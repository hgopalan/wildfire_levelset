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
| `LEVELSET_BUILD_DOCS` | `OFF` | Sphinx documentation |

See the [full build guide](https://hgopalan.github.io/wildfire_levelset/building.html) for GPU, MPI, and Windows instructions.

## Core Capabilities

- **Rothermel (1972)** surface fire spread — Anderson 13 and Scott-Burgan 40 fuel databases
- **WAF formula selection** — Andrews (2018) logarithmic WAF or BehavePlus linear WAF (`rothermel.waf_formula`); exponential Beer–Lambert canopy sheltering for forest fuels (`rothermel.waf_canopy_alpha`)
- **FARSITE elliptical expansion** (Richards 1990) + Anderson (1983) L/W ratio
- **Alternative spread models**: Balbi (2009), Cheney-Gould (1995), Cruz-Alexander-Wakimoto (2005)
- **Canadian FBP System** — O1a/O1b grass and S1/S2/S3 slash fuel types
- **Lautenberger (2013) physics-based** fire spread model
- **Non-burnable cell masking** — codes 91–99 / NB1–NB9 zero ROS (water, rock, urban, bare ground)
- **ROS stall threshold** — FARSITE-compatible floor (1×10⁻⁴ m/s)
- **Crown fire** — Van Wagner (1977) initiation + Cruz et al. (2005) + Rothermel (1991) multiplier + passive blend
- **Per-fuel burnout time** — Rothermel (1983) residence time from SAV / particle density when landscape file is present
- **Firebrand spotting** — Albini (1983) 2-D trajectory + Albini (1979) torching-tree; optional 3-D wind from [massconsistent_amr](https://github.com/hgopalan/massconsistent_amr) plotfiles; **GPU-accelerated 3-D wind interpolation** (CUDA/HIP/SYCL via `amrex::ParallelFor` on device-side `PltWindData::d_u2d`/`d_v2d`) with CPU fallback intact
- **Retardant suppression** — ROS and spotting probability zeroed inside active drop zones
- **Terrain effects** — per-cell slope/aspect from FARSITE LCP files or XYZ terrain; FARSITE full topographic horizon scan (8-direction ridge shading)
- **Wind models** — time-varying, turbulent (OU/spectral), compact direction schedule
- **Multi-station weather** — IDW spatial interpolation of per-station .wtr files (`multi_wtr_file`)
- **Fuel moisture** — FMD schedule, diurnal (Nelson 2000), precipitation wetting, FMC seasonal phenology; live fuel conditioning ramp
- **Spatial moisture output** — moisture_d1 / d10 / d100 / lh / lw in every plotfile
- **Ignition types** — point CSV, closed polygon, polyline (line fire)
- **Per-cell spatial moisture** from FARSITE .fms scenario files
- **Fire ecology diagnostics** — scorch height, tree mortality, Scott-Reinhardt TI/CI ratios
- **Full TI/CI bisection** — exact Scott & Reinhardt (2001) Torching/Crowning Index
- **Fire emissions** — CO₂, CO, PM₂.₅ (WRF-Fire convention)
- **MTT propagation** — Minimum Travel Time Dijkstra fast-marching
- **Smoke plume-rise model** — Briggs (1965/1969) buoyancy-dominated plume rise per fire-front cell; written to plotfile as `plume_rise_m`; domain maximum printed in logs (`smoke_plume.enable`)
- **KML perimeter export** — Google Earth–compatible `.kml` file alongside CSV/GeoJSON at every plotfile step; full UTM → WGS-84 conversion (`write_perimeter_kml`, `kml_utm_zone`)
- **Simulation date/time display** — calendar date/time (YYYY-MM-DD HH:MM) shown in per-step log output and HTML fire report when `sim_datetime.year/month/day` or `solar_radiation` start date is set
- **Post-fire fuel adjustment for re-entry spots** — firebrand spot fires landing in previously-burned cells have their ignition probability scaled by residual fuel fraction; spots below a configurable threshold are suppressed (`fuel_depletion.adjust_spotting_reentry`)
- **Real-time satellite fire detection assimilation** — ingests GOES-16/17/18 (NOAA, public AWS S3) or VIIRS (NASA FIRMS REST API) active-fire detections; applied as initial conditions and/or mid-simulation re-marking events with configurable confidence threshold, detection radius, and UTM coordinate projection (`satellite.*`)
- **Vectorial slope/wind combination** — FARSITE-style vector sum of φ_w and φ_s before computing ROS, capturing alignment effects; selectable via `rothermel.use_slope_wind_vectors = 1` (separate from the existing additive cross-term option `use_slope_wind_cross`)
- **Burn-period (daytime burning window) gate** — zeroes R_mf outside a configurable local clock window, pausing all spread paths during inactive hours while moisture and diagnostics continue; mirrors FARSITE / FSPro burn-period concept (`burn_period.enable`)
- **Post-frontal fuel consumption raster** — `residual_fuel` (exponential burnout fraction behind the fire front) and `fuel_consumption` (per-class consumption during FARSITE spread) written to every plotfile and included in the default GeoTIFF export set
- **AMReX-based** — GPU kernels, AMR-ready, MPI parallelism

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
| `spotting/` | firebrand_spotting, albini_spotting, **albini_spotting_3d_wind** *(new)* |
| `terrain/` | terrain_wind, balbi_viegas_heatflux, windninja_ridge_canyon, **solar_horizon_shading** *(new)*, **terrain_gradient_correction** *(new)* |
| `moisture/` | fmd_moisture, cheney_gould_grassfire, precip_wetting, spatial_moisture_output, **wtr_diurnal** *(new)*, **wtr_rain_wetting** *(new)* |
| `fuel/` | fuel_adj_file |
| `ignition/` | barrier_polygons, **satellite_assimilation** *(new)*, fire_perimeter_output, polygon_ignition, polyline_ignition, **kml_perimeter** *(new)* |
| `wind/` | time_dependent_wind, turb_wind, wind_dir_schedule, **waf_andrews**, **waf_behaviorplus** |
| `diagnostics/` | scott_reinhardt_indices, scott_reinhardt_full_ti_ci |
| `misc/` | 3d_sphere, eb_implicit, mtt_propagation, bulk_fuel_consumption, landfire_farsite, **nonburnable_mask**, **smoke_plume_rise** *(new)*, **reentry_spotting** *(new)* |

The `albini_spotting_3d_wind` test requires a Python pre-step to generate the synthetic plt wind file:

```bash
cd regtest/spotting/albini_spotting_3d_wind
python3 generate_plt_wind.py   # creates plt_wind_3d/ directory
```

```bash
cmake -S . -B build -DLEVELSET_DIM_2D=ON
cmake --build build -j
cd build && ctest -L regtest --output-on-failure
```

## Tools

Python utilities in `tools/` for terrain download, weather file parsing, GIS export, ensemble analysis, and satellite detection preparation.

Key ensemble tools:
- `ensemble_burn_probability.py` – FSim-style Monte Carlo driver (burn probability, crown fire probability, ignition sampling, containment probability, area exceedance, flame-length exceedance)
- `plot_burn_probability.py` – Visualise burn probability maps (PNG/GeoTIFF)
- `values_at_risk.py` – FSPro-style values-at-risk overlay on burn probability maps
- `fire_size_summary.py` – Fire area / perimeter / emissions statistics vs. time
- `ignition_probability_table.py` – Anderson (1970) / Rothermel (1983) P_ignition worksheet
- `plotfile_to_geotiff.py` – Export solver plotfiles to GeoTIFF; default export now includes `residual_fuel` (post-frontal burnout raster) and `fuel_consumption`

See the [tools documentation](https://hgopalan.github.io/wildfire_levelset/tools.html).

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
