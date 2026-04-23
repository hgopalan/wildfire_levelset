# wildfire_levelset

This repository contains a small AMReX-based C++ level-set solver for wildfire-front style advection based on the community fire model. The fuel and moisture databases are not added yet. Work is on-going to add a non-uniform wind field read from files for one way-coupling. A future work will use a two-way coupling with ERF by modifying the surface heat fluxes. 

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
- Dynamic AMR for negative `phi`:
  - `amr_enable_negative_phi_refine=1`
  - `amr_regrid_int=10`
  - `amr_refine_ratio=2`
  - `amr_max_refinements=1`
  - `amr_tag_phi_threshold=0.0` (cells with `phi < threshold` trigger refinement)
- Velocity:
  - `u_x=0.25`, `u_y=0.0`, `u_z=0.0`
- Initial source:
  - `source_type=sphere`
  - `sphere_center_x/y/z=0.5`
  - `sphere_radius=0.25`

## Output

- Plotfiles are written in the run directory as `plt####`.
- These can be opened with AMReX/ParaView workflows or other tools that support AMReX plotfile format.
