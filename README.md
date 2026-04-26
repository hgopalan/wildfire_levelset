# wildfire_levelset

This repository contains a small AMReX-based C++ level-set solver for wildfire-front style advection based on the community fire model. The fuel and moisture databases are not added yet. Work is on-going to add a non-uniform wind field read from files for one way-coupling. A future work will use a two-way coupling with ERF by modifying the surface heat fluxes.

The code now includes Richards' (1990) FARSITE (Fire Area Simulator) elliptical expansion model, which computes fire spread using an elliptical pattern with coefficients a, b, and c for directional spread rates based on wind conditions and fuel characteristics. 

## Prerequisites

- C++17 compiler
- CMake (3.20+)
- Git

> Note: This project supports both 2D and 3D configurations. The default is 3D, but you can build for 2D using CMake options (see Build section below).

## Clone

Clone the repository with submodules so the pinned AMReX source is available locally:

```bash
git clone --recurse-submodules https://github.com/hgopalan/wildfire_levelset.git
cd wildfire_levelset
```

If you already cloned the repository without submodules, initialize them from the repository root:

```bash
git submodule update --init --recursive
```

## Build

From the repository root:

```bash
cmake -S . -B build
cmake --build build -j
```

This project prefers the vendored `external/amrex` submodule and builds it as part of the top-level CMake configure. The build is configured for 3D AMReX setup by default.

### Building for 2D

To build for 2D instead of 3D, use the `-DLEVELSET_DIM_2D=ON` option:

```bash
cmake -S . -B build -DLEVELSET_DIM_2D=ON
cmake --build build -j
```

**Note:** Grid tagging (adaptive mesh refinement) is automatically disabled in 2D builds.

### Building with Embedded Boundary Support

To enable AMReX Embedded Boundary (EB) support for complex geometries, use the `-DLEVELSET_ENABLE_EB=ON` option:

```bash
cmake -S . -B build -DLEVELSET_ENABLE_EB=ON
cmake --build build -j
```

Embedded Boundary support allows the solver to handle complex geometries using implicit function representations, enabling features like terrain following coordinates and irregular domain boundaries.

EB support can be combined with 2D builds:

```bash
cmake -S . -B build -DLEVELSET_DIM_2D=ON -DLEVELSET_ENABLE_EB=ON
cmake --build build -j
```

If you intentionally want to use an already installed AMReX instead, configure CMake with `-DLEVELSET_USE_VENDORED_AMREX=OFF` and provide `AMReX_DIR` or `CMAKE_PREFIX_PATH` as needed.

## Run

From the repository root:

```bash
./build/levelset
```

You can override runtime parameters directly from the command line (AMReX `ParmParse` style):

```bash
./build/levelset nsteps=200 plot_int=20 cfl=0.4 u_x=0.3 u_y=0.1 u_z=0.0
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
- Initial source:
  - `source_type=sphere`
  - `sphere_center_x/y/z=0.5`
  - `sphere_radius=0.25`
- Rothermel fire spread model:
  - `rothermel.fuel_model=FM4` (Use standard fuel model from database: FM1-FM13)
    - FM1: Short Grass (1 ft)
    - FM2: Timber (Grass and Understory)
    - FM3: Tall Grass (2.5 ft)
    - FM4: Chaparral (6 ft)
    - FM5: Brush (2 ft)
    - FM6: Dormant Brush, Hardwood Slash
    - FM7: Southern Rough
    - FM8: Closed Timber Litter
    - FM9: Hardwood Litter
    - FM10: Timber (Litter and Understory)
    - FM11: Light Logging Slash
    - FM12: Medium Logging Slash
    - FM13: Heavy Logging Slash
  - Can also use aliases: `1`, `SHORT_GRASS`, `GRASS`, `CHAPARRAL`, `BRUSH`, etc.
  - Default (when no fuel model specified): custom Southern California chaparral
    - w0=0.230, σ=1739.0, δ=2.0, Mx=0.30
  - Individual fuel parameters can be overridden:
    - `rothermel.w0=0.230` (oven-dry fuel load, lb/ft²)
    - `rothermel.sigma=2000.0` (surface-area-to-volume ratio, ft⁻¹)
    - `rothermel.delta=6.0` (fuel bed depth, ft)
    - `rothermel.M_f=0.08` (fuel moisture content, fraction)
    - `rothermel.M_x=0.20` (moisture of extinction, fraction)
    - `rothermel.h_heat=8000.0` (heat content, BTU/lb)
    - `rothermel.S_T=0.0555` (total mineral content, fraction)
    - `rothermel.S_e=0.010` (effective mineral content, fraction)
    - `rothermel.rho_p=32.0` (particle density, lb/ft³)
  - Terrain parameters:
    - `rothermel.slope_x=0.0` (constant terrain slope in x-direction, tan(angle))
    - `rothermel.slope_y=0.0` (constant terrain slope in y-direction, tan(angle))
    - `rothermel.terrain_file=""` (path to terrain file with X Y Z format for spatially-varying slopes)
      - When provided, slopes are computed from terrain data at each grid cell
      - Terrain file format: ASCII file with 3 columns (X Y Z), one data point per line
      - Lines starting with '#' are treated as comments
      - Slopes are computed using inverse distance weighting interpolation
      - **Note:** When terrain_file is specified, it takes precedence over constant slope_x/slope_y values
  - Unit conversion factors:
    - `rothermel.wind_conv=196.85` (converts simulation velocity to ft/min)
    - `rothermel.ros_conv=0.00508` (converts ft/min to simulation units)
- FARSITE ellipse model (Richards 1990):
  - `farsite.enable=1` (1 to enable FARSITE ellipse model, 0 to disable)
  - `farsite.length_to_width_ratio=3.0` (L/W ratio of fire spread ellipse)
  - `farsite.phi_threshold=0.1` (threshold for identifying fire front, cells with |phi| < threshold)
  - `farsite.coeff_a=1.0` (head fire coefficient for elliptical expansion)
  - `farsite.coeff_b=0.5` (flank fire coefficient for elliptical expansion)
  - `farsite.coeff_c=0.2` (backing fire coefficient for elliptical expansion)
- Level set control:
  - `skip_levelset=0` (1 to skip level set advection and use initial phi throughout simulation, 0 for normal operation with evolving level set)


### Example: Creating a terrain file

The terrain file should be an ASCII file with three columns (X, Y, Z) representing the spatial coordinates and elevation:

```
# Terrain data for wildfire simulation
# X Y Z (coordinates and elevation)
0.0 0.0 100.0
0.0 0.5 105.0
0.0 1.0 110.0
0.5 0.0 102.0
0.5 0.5 107.0
0.5 1.0 112.0
1.0 0.0 104.0
1.0 0.5 109.0
1.0 1.0 114.0
```

The terrain slopes (∂z/∂x and ∂z/∂y) are automatically computed at each grid cell using inverse distance weighting interpolation and central differences. The slope factor φ_s in the Rothermel model is then calculated as `φ_s = 5.275 * tan²(slope)` where `tan(slope) = √(slope_x² + slope_y²)`.

### Example: Using FARSITE with fixed initial phi

To run FARSITE ellipse model with a fixed initial phi (skipping level set evolution):

```bash
./build/levelset skip_levelset=1 farsite.enable=1
```

This mode is useful for analyzing FARSITE spread patterns based on initial fire geometry without the complexity of level set advection. When both `skip_levelset=1` and `farsite.enable=1` are set, the phi field is initialized as an indicator function (phi = 1 inside the fire region, phi = 0 outside) instead of a signed distance function. The phi field remains at its initial configuration throughout the simulation, while FARSITE computes directional spread rates at each timestep.

**Note:** In normal level set mode (`skip_levelset=0`), phi is initialized as a signed distance function (negative inside, positive outside). The indicator function initialization (phi = 1 inside, phi = 0 outside) is only used when both FARSITE is enabled and level set is skipped.

## Output

- Plotfiles are written in the run directory as `plt####`.
- These can be opened with AMReX/ParaView workflows or other tools that support AMReX plotfile format.
- Each plotfile contains:
  - `phi`: Level-set function 
    - In normal mode (`skip_levelset=0`): signed distance function (negative inside burned region, positive outside)
    - In FARSITE-only mode (`skip_levelset=1` and `farsite.enable=1`): indicator function (1 inside burned region, 0 outside)
  - `velx`, `vely`, `velz`: Velocity field components
  - `farsite_dx`, `farsite_dy`, `farsite_dz`: FARSITE ellipse spread displacements (when enabled)
