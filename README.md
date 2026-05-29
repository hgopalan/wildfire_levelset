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
| `LEVELSET_BUILD_PYTHON_BINDINGS` | `OFF` | Python API (pyWildfire) |
| `LEVELSET_BUILD_DOCS` | `OFF` | Sphinx documentation |

See the [full build guide](https://hgopalan.github.io/wildfire_levelset/building.html) for GPU, MPI, and Windows instructions.

## Python API for Coupled Simulations

The Python API (`pyWildfire`) enables full programmatic control of the fire solver from Python for coupled wind-fire simulations, ensemble forecasting, and machine learning applications.

**Quick example:**
```python
from wildfire_solver import WildfireSolver

fire = WildfireSolver("inputs.i")
for i in range(100):
    fire.step()
fire.finalize()
```

**Key capabilities:** State extraction as NumPy arrays • 2D/3D wind updates • AMReX plotfile writing • Zero-copy data transfer

**See the [Python API documentation](https://hgopalan.github.io/wildfire_levelset/python_api.html) for complete API reference, coupled simulation examples, wind solver integration patterns, and applications.**

## Core Capabilities

| Category | Feature | Key Reference |
|----------|---------|---------------|
| **Fire Spread Models** | Rothermel (1972) surface fire with Anderson 13 & Scott-Burgan 40 fuel databases | `rothermel.*` |
| | FARSITE elliptical expansion (Richards 1990) + Anderson L/W ratio | `farsite.*` |
| | Alternative models: Balbi (2009), Cheney-Gould (1995), Cruz-Alexander-Wakimoto (2005) | `balbi.*`, `cheney.*`, `cruz.*` |
| | Canadian FBP System (O1a/O1b grass, S1/S2/S3 slash) | `fbp.*` |
| | Lautenberger (2013) physics-based spread | `lautenberger.*` |
| **Fire Behavior** | Wind Adjustment Factor: Andrews logarithmic or BehavePlus linear | `rothermel.waf_formula` |
| | Wind-terrain feedback: 8 models (Viegas, canyon, ridge, WindNinja, FARSITE) | `wind_terrain.model` |
| | Wind-fuel interaction: canopy sheltering via LAI-based attenuation | `wind_fuel.*` |
| | Crown fire: Van Wagner (1977) + Cruz et al. (2005) + Rothermel (1991) | `crown.*` |
| | Per-fuel burnout time (Rothermel 1983 residence time) | Auto from landscape |
| | Intensity-dependent flame residence time (Byram 1959) | `intensity_residence.*` |
| | Burn-period daytime window (FARSITE/FSPro concept) | `burn_period.enable` |
| | Radiation-driven preheating distance ahead of fire front | `preheating.*` |
| | Fuel particle temperature evolution for ignition timing | `fuel_temperature.*` |
| | Critical heat flux for moisture-dependent ignition | `critical_heat_flux.*` |
| | Fuel moisture of extinction gradient (SAV-dependent) | `M_x.*` |
| | Flame intermittency factor for heat transfer efficiency | `intermittency.*` |
| **Spotting & Embers** | Firebrand spotting: Albini (1983) trajectory + torching-tree (1979) | `firebrand.*` |
| | Vorticity-enhanced spotting: Weise & Biging (1996) fire whirl effects | `weise_biging.enhance_spotting` |
| | Ember cascade: Gaussian flux field (Sardoy 2007 approach) | `ember_cascade.*` |
| | GPU-accelerated 3-D wind interpolation (CUDA/HIP/SYCL) | Optional from massconsistent_amr |
| | Post-fire fuel adjustment for re-entry spots | `fuel_depletion.adjust_spotting_reentry` |
| | Plume entrainment momentum feedback on surface wind | `plume_momentum.*` |
| **Terrain & Weather** | Slope/aspect from FARSITE LCP or XYZ; 8-direction horizon scan | `terrain.*` |
| | Time-varying, turbulent (OU/spectral), direction-schedule winds | `wind.*` |
| | Multi-station weather with IDW interpolation | `multi_wtr_file` |
| **Fuel Moisture** | FMD schedule, diurnal (Nelson 2000), precipitation wetting, FMC phenology | `fuel_moisture.*` |
| | Spatial moisture from FARSITE .fms files | `.fms` support |
| | Moisture fields (d1/d10/d100/lh/lw) in every plotfile | Auto output |
| | Spatially varying fuel bed bulk density/loading | `fuel_loading_variation.*` |
| | McArthur temperature/RH-dependent drying (new) | `mcarthur_moisture.*` |
| | Enhanced FMC phenology: sinusoidal & GDD models (new) | `fmc_phenology.*` |
| **Spotting & Embers** | Firebrand spotting: Albini (1983) trajectory + torching-tree (1979) | `firebrand.*` |
| | Vorticity-enhanced spotting: Weise & Biging (1996) fire whirl effects | `weise_biging.enhance_spotting` |
| | Ember cascade: Gaussian flux field (Sardoy 2007 approach) | `ember_cascade.*` |
| | Ember accumulation with probabilistic ignition (new) | `ember_accumulation.*` |
| | GPU-accelerated 3-D wind interpolation (CUDA/HIP/SYCL) | Optional from massconsistent_amr |
| | Post-fire fuel adjustment for re-entry spots | `fuel_depletion.adjust_spotting_reentry` |
| | Plume entrainment momentum feedback on surface wind | `plume_momentum.*` |
| **Terrain & Weather** | Slope/aspect from FARSITE LCP or XYZ; 8-direction horizon scan | `terrain.*` |
| | Time-varying, turbulent (OU/spectral), direction-schedule winds | `wind.*` |
| | Periodic wind gust factor (new) | `wind_gust.*` |
| | Multi-station weather with IDW interpolation | `multi_wtr_file` |
| **Fire Behavior** | Wind Adjustment Factor: Andrews logarithmic or BehavePlus linear | `rothermel.waf_formula` |
| | Wind-terrain feedback: 8 models (Viegas, canyon, ridge, WindNinja, FARSITE) | `wind_terrain.model` |
| | Wind-fuel interaction: canopy sheltering via LAI-based attenuation | `wind_fuel.*` |
| | Crown fire: Van Wagner (1977) + Cruz et al. (2005) + Rothermel (1991) | `crown.*` |
| | Per-fuel burnout time (Rothermel 1983 residence time) | Auto from landscape |
| | Intensity-dependent flame residence time (Byram 1959) | `intensity_residence.*` |
| | Burn-period daytime window (FARSITE/FSPro concept) | `burn_period.enable` |
| | Radiation-driven preheating distance ahead of fire front | `preheating.*` |
| | Slope-dependent flame tilt for radiation (new) | `flame_tilt.*` |
| | Fuel particle temperature evolution for ignition timing | `fuel_temperature.*` |
| | Critical heat flux for moisture-dependent ignition | `critical_heat_flux.*` |
| | Fuel moisture of extinction gradient (SAV-dependent) | `M_x.*` |
| | Flame intermittency factor for heat transfer efficiency | `intermittency.*` |
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
| | Fire line intensity rate of change (dI/dt) for blow-up detection | `intensity_rate_of_change.*` |
| **Technical** | AMReX-based: GPU kernels (CUDA/HIP/SYCL), AMR, MPI | Build options |
| | MTT (Minimum Travel Time) Dijkstra propagation | `mtt.*` |
| | Non-burnable masking (codes 91–99 / NB1–NB9) | Auto from fuel codes |
| | ROS stall threshold (FARSITE-compatible 1×10⁻⁴ m/s) | Built-in |
| | Vectorial slope/wind combination | `rothermel.use_slope_wind_vectors` |
| | Conditional weather ERC/BI/SC trigger | `conditional_weather_trigger` |

## Ensemble / FSim-Style Probabilistic Simulation

The `tools/ensemble_burn_probability.py` tool implements FSim-style Monte Carlo simulation with perturbations of wind, moisture, and ignition parameters. It supports probabilistic ignition locations, containment probability, crown fire probability mapping, and burn probability analysis.

See the [tools documentation](https://hgopalan.github.io/wildfire_levelset/tools.html#ensemble-burn-probability-py-ensemble-burn-probability-driver) for complete feature list and usage examples.

## Regression Tests

Comprehensive regression tests are organized in `regtest/` subdirectories covering surface spread, crown fire, spotting, terrain, moisture, fuel, ignition, wind, diagnostics, and Python API functionality. All tests use UTM Zone 11N coordinates (Southern California reference: 330000 E, 3775000 N).

**Run all tests:**
```bash
cmake -S . -B build -DLEVELSET_DIM_2D=ON
cmake --build build -j
cd build && ctest -L regtest --output-on-failure
```

See the [regression tests documentation](https://hgopalan.github.io/wildfire_levelset/regtests.html) for complete test descriptions and requirements.

## Tools

Python utilities in `tools/` for terrain download, weather parsing, GIS export, ensemble analysis, and post-processing. Categories include ensemble simulation, fire analysis, GIS/export, input preparation, and BehavePlus-style worksheets.

See the [tools documentation](https://hgopalan.github.io/wildfire_levelset/tools.html) for complete tool descriptions and usage examples.

## Documentation and References

Full documentation is available at [https://hgopalan.github.io/wildfire_levelset/](https://hgopalan.github.io/wildfire_levelset/), including:

* [Mathematical Models](https://hgopalan.github.io/wildfire_levelset/mathematical_models.html) - Complete model equations and references
* [Usage Guide](https://hgopalan.github.io/wildfire_levelset/usage.html) - Input parameters and configuration
* [Python API](https://hgopalan.github.io/wildfire_levelset/python_api.html) - Coupled simulation interface
* [Tools](https://hgopalan.github.io/wildfire_levelset/tools.html) - Pre/post-processing utilities
* [Comparison](https://hgopalan.github.io/wildfire_levelset/comparison.html) - Feature comparison with FARSITE, FlamMap, WRF-SFIRE

## Satellite Fire Detection Assimilation

Active-fire detections from GOES-16/17/18, VIIRS, or CSV can be ingested to constrain the initial fire perimeter or correct the simulated perimeter during simulation.

See the [satellite assimilation section in tools documentation](https://hgopalan.github.io/wildfire_levelset/tools.html#satellite-goes-to-csv-py-satellite-fire-detection-assimilation) for parameter reference and usage examples.

## License

See [LICENSE](LICENSE) file.
