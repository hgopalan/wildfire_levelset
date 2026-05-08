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
- **FARSITE elliptical expansion** (Richards 1990) + Anderson (1983) L/W ratio
- **Alternative spread models**: Balbi (2009), Cheney-Gould (1995), Cruz-Alexander-Wakimoto (2005)
- **Crown fire** — Van Wagner (1977) initiation + Cruz et al. (2005) active crown ROS
- **Firebrand spotting** — Albini (1983) 2-D trajectory + Albini (1979) torching-tree
- **Terrain effects** — per-cell slope/aspect from FARSITE LCP files or XYZ terrain
- **Wind models** — time-varying, turbulent (OU/spectral), compact direction schedule
- **Fuel moisture** — FMD schedule, diurnal (Nelson 2000), precipitation wetting, FMC seasonal phenology
- **Ignition types** — point CSV, closed polygon, polyline (line fire)
- **Per-cell spatial moisture** from FARSITE .fms scenario files
- **Fire ecology diagnostics** — scorch height, tree mortality, Scott-Reinhardt TI/CI ratios
- **Fire emissions** — CO₂, CO, PM₂.₅ (WRF-Fire convention)
- **MTT propagation** — Minimum Travel Time Dijkstra fast-marching
- **AMReX-based** — GPU kernels, AMR-ready, MPI parallelism

## Regression Tests

Tests are organised into sub-folders under `regtest/`, all using UTM Zone 11N coordinates
(Southern California reference: 330000 E, 3775000 N):

| Sub-folder | Tests |
|---|---|
| `surface_spread/` | basic_levelset, farsite_ellipse, rothermel_fuel, anderson_lw, catchpole_demestre, wilson_spread, alexander_lemniscate, ellipse_sdf, reinitialization |
| `crown_fire/` | crown_initiation, cruz_crown_continental_us, **fmc_seasonal** *(new)* |
| `spotting/` | firebrand_spotting, albini_spotting |
| `terrain/` | terrain_wind, balbi_viegas_heatflux, windninja_ridge_canyon |
| `moisture/` | fmd_moisture, cheney_gould_grassfire, **precip_wetting** *(new)* |
| `fuel/` | fuel_adj_file |
| `ignition/` | barrier_polygons, fire_perimeter_output, **polygon_ignition** *(new)*, **polyline_ignition** *(new)* |
| `wind/` | time_dependent_wind, turb_wind, **wind_dir_schedule** *(new)* |
| `diagnostics/` | **scott_reinhardt_indices** *(new)* |
| `misc/` | 3d_sphere, eb_implicit, mtt_propagation, bulk_fuel_consumption, landfire_farsite |

```bash
cmake -S . -B build -DLEVELSET_DIM_2D=ON
cmake --build build -j
cd build && ctest -L regtest --output-on-failure
```

## Tools

Python utilities in `tools/` for terrain download, weather file parsing, GIS export, and ensemble analysis.
See the [tools documentation](https://hgopalan.github.io/wildfire_levelset/tools.html).

## References

See the [full reference list](https://hgopalan.github.io/wildfire_levelset/references.html) in the documentation.

Key references: Rothermel (1972), Richards (1990), Van Wagner (1977), Cruz et al. (2005),
Albini (1983), Andrews (2018), Nelson (2000), Scott & Reinhardt (2001), Finney (2004).

## License

See [LICENSE](LICENSE) file.
