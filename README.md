# wildfire_levelset

This repository contains a small AMReX-based C++ level-set solver for wildfire-front style advection based on the community fire model. The fuel and moisture databases are not added yet. Work is on-going to add a non-uniform wind field read from files for one way-coupling. A future work will use a two-way coupling with ERF by modifying the surface heat fluxes. 

## Features

### Fire Spread Models

The code now includes advanced fire spread models:

1. **Anderson L/W Ratio (FARSITE Model)**
   - Length-to-Width ratio calculation based on wind speed
   - Based on Anderson (1983): L/W = 0.936 * exp(0.2566 * U) + 0.461 * exp(-0.1548 * U) - 0.397
   - Accounts for elliptical fire shape under wind influence
   - Used when `use_farsite_model=true`

2. **Rothermel Terrain Effects**
   - Slope correction factor based on Rothermel (1972)
   - Accounts for uphill/downhill fire spread acceleration/deceleration
   - φ_s = 5.275 * tan²(slope)
   - Used when `use_terrain_effects=true` and `use_farsite_model=false`

3. **FARSITE Combined Wind and Terrain Model**
   - Vectorial combination of wind and slope effects
   - Accounts for alignment between wind direction and slope aspect
   - Enhanced rate of spread calculations
   - Used when both `use_terrain_effects=true` and `use_farsite_model=true`

## Prerequisites

- C++17 compiler
- CMake (3.20+)
- An installed AMReX build discoverable by CMake (`AMReXConfig.cmake`)

> Note: This project currently assumes a 3D setup (`x/y/z` parameters are used throughout).

## Build

From the repository root:

```bash
cmake -S . -B build -DAMReX_DIR=/path/to/amrex/install/lib/cmake/AMReX
cmake --build build -j
```

If AMReX is already on your `CMAKE_PREFIX_PATH`, the `-DAMReX_DIR=...` flag can be omitted.

## Run

From the repository root:

```bash
./build/levelset
```

You can override runtime parameters directly from the command line (AMReX `ParmParse` style):

```bash
./build/levelset nsteps=200 plot_int=20 cfl=0.4 u_x=0.3 u_y=0.1 u_z=0.0
```

### Fire Spread Model Examples

**Basic model (no terrain effects):**
```bash
./build/levelset u_x=0.5 u_y=0.0
```

**Rothermel model with terrain effects:**
```bash
./build/levelset u_x=0.5 u_y=0.0 use_terrain_effects=true terrain_slope=15.0 terrain_aspect=90.0
```

**FARSITE model with Anderson L/W ratio:**
```bash
./build/levelset u_x=0.5 u_y=0.0 use_terrain_effects=true use_farsite_model=true terrain_slope=15.0 terrain_aspect=90.0
```

## Runtime parameters (defaults)

- Grid / domain:
  - `n_cell=64` (or `n_cell_x`, `n_cell_y`, `n_cell_z`)
  - `max_grid_size=32`
  - `prob_lo_x/y/z=0.0`, `prob_hi_x/y/z=1.0`
- Time stepping:
  - `nsteps=300`
  - `cfl=0.5`
  - `plot_int=50`
- Reinitialization:
  - `reinit_int=20`
  - `reinit_iters=20`
  - `reinit_dtau=0.5`
- Velocity:
  - `u_x=0.25`, `u_y=0.0`, `u_z=0.0`
- Fire spread models:
  - `use_terrain_effects=false` - Enable terrain effects (Rothermel/FARSITE)
  - `use_farsite_model=false` - Use FARSITE model with Anderson L/W ratio
  - `terrain_slope=0.0` - Terrain slope in degrees (0-90)
  - `terrain_aspect=0.0` - Terrain aspect in degrees (0=East, 90=North, 180=West, 270=South)
- Initial source:
  - `source_type=sphere`
  - `sphere_center_x/y/z=0.5`
  - `sphere_radius=0.25`

## Output

- Plotfiles are written in the run directory as `plt####`.
- These can be opened with AMReX/ParaView workflows or other tools that support AMReX plotfile format.
