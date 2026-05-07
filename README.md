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
| `LEVELSET_USE_VENDORED_AMREX` | `ON` | Use the bundled AMReX submodule |

### Run

```bash
./build/levelset
```

Override parameters from command line:
```bash
./build/levelset nsteps=200 plot_int=20 cfl=0.4 u_x=0.3 u_y=0.1 u_z=0.0
```

## Tools

Python utilities live in the `tools/` directory for terrain/landscape download, weather file parsing, GIS export, and ensemble analysis.
See the 📚 [full documentation](https://hgopalan.github.io/wildfire_levelset/tools.html) for complete option references and worked examples.

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
