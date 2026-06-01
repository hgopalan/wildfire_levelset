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

**Fire Spread Models**: Rothermel (1972) with Anderson 13 & Scott-Burgan 40 fuel databases • FARSITE elliptical expansion • Alternative models (Balbi, Cheney-Gould, Cruz, Canadian FBP, Lautenberger) • Two-fuel model blending (linear, harmonic, maximum, Finney-style ROS)

**Fire Behavior**: Crown fire initiation & spread • Wind-terrain-fuel interactions • Radiation preheating with slope-dependent flame tilt • Byram convective number • Flame tilt angle • Packing ratio diagnostics (β/β_opt) • Flame front depth • McArthur FFDI • Fire acceleration (Anderson temporal & size-based) • Scott & Reinhardt Crown Fire Surface Area (CFSA) for 3-D canopy structure • Backing fire ROS empirical ratios

**Fuel Moisture**: Time-varying schedules (FMD/FMC) • Diurnal cycles • Precipitation wetting • McArthur temperature/RH scaling • Enhanced phenology models (sinusoidal & growing degree day) • Canadian FWI System (FFMC/DMC/DC/BUI/FWI) • Dynamic moisture of extinction (M_x) • Duff moisture with smoldering combustion • Fine fuel moisture time-lag differential equations • Grass curing model

**Fire Danger Indices**: Keetch-Byram Drought Index (KBDI) • Haines Index atmospheric stability • McArthur FFDI • Canadian FWI System • NFDRS Spread Component (SC) • Chandler Burning Index (CBI)

**Spotting**: Albini trajectory model • Ember cascade • Ember accumulation tracking with decay & probabilistic ignition

**Weather & Terrain**: FARSITE LCP terrain support • Turbulent wind models • Periodic gust factor • Multi-station interpolation • Multi-layer canopy wind profile (exponential/logarithmic vertical distribution) • Diurnal weather cycles • Elevation temperature lapse rate

**Technical**: GPU acceleration (CUDA/HIP/SYCL) • MPI parallelism • Python API for coupled simulations

See [full documentation](https://hgopalan.github.io/wildfire_levelset/) for complete feature list, model equations, and parameters.

**Documentation Links:**
- [Mathematical Models](https://hgopalan.github.io/wildfire_levelset/mathematical_models.html)
- [Usage Guide & Parameters](https://hgopalan.github.io/wildfire_levelset/usage.html)
- [Python API](https://hgopalan.github.io/wildfire_levelset/python_api.html)
- [Tools & Utilities](https://hgopalan.github.io/wildfire_levelset/tools.html)
- [Scientific References](https://hgopalan.github.io/wildfire_levelset/references.html)

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
