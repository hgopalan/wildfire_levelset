Code Structure
==============

This section describes the organization and architecture of the wildfire level-set solver.

Directory Structure
-------------------

::

    wildfire_levelset/
    ├── CMakeLists.txt              # Main CMake build configuration
    ├── README.md                   # Project README
    ├── LICENSE                     # License file
    ├── docs/                       # Documentation (this directory)
    ├── external/                   # External dependencies
    │   └── amrex/                  # AMReX submodule
    ├── src/                        # Source code
    │   ├── main.cpp                # Main entry point
    │   ├── parse_inputs.H          # Input parameter data structures
    │   ├── parse_inputs.cpp        # Input parameter parsing implementation
    │   ├── advection.H             # Level-set advection (WENO5-Z + RK3)
    │   ├── rothermel_model.H       # Rothermel (1972) fire spread model
    │   ├── compute_rothermel_R.H   # Per-cell Rothermel ROS computation
    │   ├── farsite_ellipse.H       # FARSITE elliptical expansion model
    │   ├── fire_models.H           # Fire model utilities
    │   ├── compute_fire_behavior.H # Byram fireline intensity and flame length
    │   ├── terrain_slope.H         # Terrain slope calculations
    │   ├── landscape_file.H        # FARSITE landscape file (.lcp) reader
    │   ├── spatial_grid.H          # Spatial grid utilities
    │   ├── velocity_field.H        # Wind field management (static + time-dep.)
    │   ├── crown_initiation.H      # Van Wagner crown fire initiation
    │   ├── firebrand_spotting.H    # Probability-based firebrand spotting
    │   ├── albini_spotting.H       # Albini (1983) physics-based spotting
    │   ├── bulk_fuel_consumption.H # Post-frontal fuel consumption
    │   ├── initial_conditions.H    # Initial fire setup (sphere/box/ellipse/EB/CSV)
    │   ├── boundary_conditions.H   # Boundary conditions
    │   ├── numerical_schemes.H     # Numerical discretization utilities
    │   ├── reinitialization.H      # Level-set reinitialization (Sussman)
    │   ├── compute_dt.H            # CFL time step calculation
    │   ├── plot_results.H          # AMReX plotfile output
    │   ├── write_xy_data.H         # XY data and convex hull export
    │   └── fuel_database.H         # Anderson 13 / Scott & Burgan 40 fuel database
    ├── regtest/                    # Regression tests
    │   ├── basic_levelset/         # Basic level-set advection test
    │   ├── 3d_sphere/              # 3D sphere level-set test
    │   ├── anderson_lw/            # Anderson L/W ratio test
    │   ├── farsite_ellipse/        # FARSITE ellipse spread test
    │   ├── ellipse_sdf/            # Elliptical SDF initial condition test
    │   ├── rothermel_fuel/         # Rothermel fuel model test
    │   ├── terrain_wind/           # Terrain slope + wind test
    │   ├── terrain_wind_preprocess/# Terrain preprocessing test
    │   ├── crown_initiation/       # Van Wagner crown fire test
    │   ├── firebrand_spotting/     # Probability-based spotting test
    │   ├── albini_spotting/        # Albini (1983) spotting test
    │   ├── bulk_fuel_consumption/  # Fuel consumption test
    │   ├── reinitialization/       # Level-set reinitialization test
    │   ├── time_dependent_wind/    # Time-dependent wind field test
    │   ├── landfire_farsite/       # LANDFIRE landscape + FARSITE test
    │   └── eb_implicit/            # Embedded boundary initial condition test
    └── tests/                      # Additional unit tests

Core Components
---------------

Main Simulation Loop
^^^^^^^^^^^^^^^^^^^^

**File**: ``src/main.cpp``

The main simulation loop performs the following steps:

1. Parse input parameters from configuration file
2. Set up AMReX geometry and grid hierarchy
3. Initialize level-set function :math:`\phi`, velocity field, terrain slopes, and diagnostic fields
4. Build per-cell Rothermel lookup table from landscape file (if provided)
5. Time stepping loop:

   a. Update time-dependent wind field (if enabled)
   b. Compute Rothermel rate of spread using per-cell fuel data
   c. Compute Byram fireline intensity and flame length diagnostics
   d. Advance level-set function (level-set advection) or apply FARSITE ellipse
   e. Apply probability-based firebrand spotting (if enabled)
   f. Apply Albini physics-based spotting (if enabled)
   g. Apply level-set reinitialization (if enabled and in level-set mode)
   h. Write AMReX plotfiles and XY data files
   i. Advance time :math:`t \leftarrow t + \Delta t`

Input Parsing
^^^^^^^^^^^^^

**Files**: ``src/parse_inputs.H``, ``src/parse_inputs.cpp``

The input parser reads configuration files and stores parameters in structured data types:

* ``RothermelParams``: Fuel properties for Rothermel model, per-class fuel loads and moistures, terrain file paths
* ``InputParameters::FARSITEParams``: FARSITE ellipse parameters and bulk fuel consumption
* ``InputParameters::CrownInitiationParams``: Van Wagner crown fire parameters
* ``InputParameters::SpottingParams``: Probability-based firebrand spotting parameters
* ``InputParameters::AlbiniSpottingParams``: Albini (1983) physics-based spotting parameters

Rothermel Model
^^^^^^^^^^^^^^^

**File**: ``src/rothermel_model.H``

Implements the Rothermel (1972) fire spread equations. Key function:

* ``compute_rothermel_params()``: Computes all Rothermel parameters including:

  - Packing ratio :math:`\beta`
  - Reaction intensity :math:`I_R`
  - No-wind, no-slope ROS :math:`R_0`
  - Wind factor coefficients :math:`C`, :math:`B`, :math:`E`
  - Slope factor :math:`\phi_s`

Supports both single-class (aggregate) and multi-class (per-size-class) fuel
moisture and loading paths. Returns a ``RothermelComputed`` structure with all
computed values.

Per-Cell ROS Computation
^^^^^^^^^^^^^^^^^^^^^^^^^

**File**: ``src/compute_rothermel_R.H``

Computes the Rothermel rate-of-spread field for every grid cell. When a landscape
file is provided, per-cell fuel model data from a pre-built lookup table (indexed
by fuel code) are used, enabling spatially varying fuel properties. When no
landscape file is present, the global ``RothermelParams`` are used uniformly.

FARSITE Ellipse Model
^^^^^^^^^^^^^^^^^^^^^^

**File**: ``src/farsite_ellipse.H``

Implements Richards' (1990) elliptical fire spread model:

* ``compute_farsite_spread()``: Main function that:

  1. Identifies fire front locations (where :math:`\phi \approx 0`)
  2. Computes base ROS using Rothermel model
  3. Applies Van Wagner crown fire modification if enabled
  4. Computes elliptical expansion using Richards' coefficients :math:`a`, :math:`b`, :math:`c`
  5. Updates level-set function with new fire front
  6. Computes bulk fuel consumption fraction if enabled

Fire Behavior Diagnostics
^^^^^^^^^^^^^^^^^^^^^^^^^^

**File**: ``src/compute_fire_behavior.H``

Computes Byram (1959) fire behavior metrics from the Rothermel ROS field:

* ``compute_fire_behavior()``: For each grid cell computes:

  - Fireline intensity :math:`I_B = H \cdot w_a \cdot R` [kW/m]
  - Flame length :math:`L_f = 0.0775 \cdot I_B^{0.46}` [m]

Output fields ``fireline_intensity`` and ``flame_length`` are written to every plotfile.

Level-Set Advection
^^^^^^^^^^^^^^^^^^^

**File**: ``src/advection.H``

Implements the level-set advection equation:

.. math::

   \frac{\partial \phi}{\partial t} + V|\nabla\phi| = 0

* Uses WENO5-Z (5th-order Weighted Essentially Non-Oscillatory) spatial discretization
* Third-order Runge-Kutta (RK3) time integration
* Computes :math:`|\nabla\phi|` using WENO5-Z upwind differences

Numerical Schemes
^^^^^^^^^^^^^^^^^

**File**: ``src/numerical_schemes.H``

Provides numerical discretization utilities:

* WENO5-Z upwind gradient computation
* CFL time step calculation
* Spatial derivative operators

Reinitialization
^^^^^^^^^^^^^^^^

**File**: ``src/reinitialization.H``

Maintains the level-set function as a signed distance function using Sussman's
iterative reinitialization equation:

* Ensures :math:`|\nabla\phi| = 1` away from the interface
* Controlled by ``reinit_int`` (frequency) and ``reinit_iters`` (iterations)

Terrain Slope
^^^^^^^^^^^^^

**File**: ``src/terrain_slope.H``

Computes terrain slopes from elevation data:

* Reads X Y Z terrain data from file
* Interpolates elevation using inverse distance weighting (IDW)
* Computes slopes :math:`\partial z/\partial x` and :math:`\partial z/\partial y` using central differences
* Also computes elevation, slope (degrees), and aspect fields for plotfile output

Landscape File Reader
^^^^^^^^^^^^^^^^^^^^^^

**File**: ``src/landscape_file.H``

Reads FARSITE landscape files (ASCII format) containing X Y ELEVATION SLOPE ASPECT FUEL_MODEL:

* ``read_landscape_file()``: Parses the ASCII landscape file
* ``compute_slopes_from_landscape()``: IDW interpolation to grid for slope/aspect
* ``compute_elevation_from_landscape()``: IDW interpolation for elevation
* ``compute_fuel_model_from_landscape()``: Maps fuel model codes to grid cells
* ``build_fuel_rothermel_table()``: Builds a per-fuel-code Rothermel lookup table
  supporting both Anderson 13 (FBFM13) and Scott & Burgan 40 (FBFM40) fuel systems

Velocity Field
^^^^^^^^^^^^^^

**File**: ``src/velocity_field.H``

Manages wind velocity field:

* Constant (uniform) wind from ``u_x``, ``u_y``, ``u_z`` inputs
* Spatially-varying wind read from CSV file (``velocity_file``)
* Time-dependent wind: loads sequential CSV snapshots and performs temporal
  linear interpolation + spatial IDW interpolation at each time step
  (only available in 2D builds; controlled by ``use_time_dependent_wind``)

Crown Fire Initiation
^^^^^^^^^^^^^^^^^^^^^^

**File**: ``src/crown_initiation.H``

Implements Van Wagner (1977) crown fire model:

* Computes critical surface intensity :math:`I_0` using CBH, FMC parameters
* Checks crown fire initiation criterion :math:`I > I_0`
* Active crown fire criterion using canopy bulk density CBD
* Modifies fire spread rate for active crown fire
* Outputs ``crown_fraction`` diagnostic field

Firebrand Spotting (Probability-Based)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**File**: ``src/firebrand_spotting.H``

Simulates stochastic firebrand spotting:

* Probabilistic spotting model weighted by fire intensity and wind speed
* Landing distance drawn from lognormal or exponential distribution
* Lateral angular dispersion perpendicular to wind direction
* Creates new ignition point ignitions ahead of the main fire

Albini (1983) Physics-Based Spotting
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**File**: ``src/albini_spotting.H``

Simulates firebrand spotting using Albini's thermal-plume lofting formula:

* Byram's fire line intensity computed from Rothermel ROS
* Albini lofting height :math:`H_z = 12.2 I_B^{1/3}` [m]
* 2-D horizontal trajectory integrated forward-Euler using bilinear wind
  field interpolation
* Intensity-weighted launch probability per fire-front cell
* Outputs four diagnostic fields: ``albini_Hz``, ``albini_count``,
  ``albini_dist``, ``albini_active``

Bulk Fuel Consumption
^^^^^^^^^^^^^^^^^^^^^^

**File**: ``src/bulk_fuel_consumption.H``

Computes post-frontal fuel consumption:

* Exponential consumption model with residence time parameter
* Min/max consumption fraction bounds
* Heat release rate based on consumed fuel load
* Outputs ``fuel_consumption`` diagnostic field

Fuel Database
^^^^^^^^^^^^^

**File**: ``src/fuel_database.H``

Standard fuel model property database:

* Anderson 13 fuel models (FBFM13): FM1–FM13
* Scott & Burgan 40 fuel models (FBFM40): FM101–FM256 and others
* Full per-class fuel loads (1-hr, 10-hr, 100-hr dead; live herbaceous; live woody)
* Per-class SAV ratios
* Lookup by name (e.g. ``FM4``, ``GR2``, ``SH7``)

Initial Conditions
^^^^^^^^^^^^^^^^^^

**File**: ``src/initial_conditions.H``

Sets up initial fire configuration:

* Spherical ignition (signed distance function or indicator)
* Box ignition (rectangular region)
* Elliptical ignition
* Embedded boundary implicit function (``eb_type``: sphere, ellipsoid, cylinder, plane)
* CSV fire points: reads ignition coordinates from file and initializes a true
  signed distance function from the union of circular ignition disks

Boundary Conditions
^^^^^^^^^^^^^^^^^^^

**File**: ``src/boundary_conditions.H``

Applies boundary conditions to level-set function:

* Extrapolation (zero-gradient fill into ghost cells)
* Periodic

Output and Visualization
^^^^^^^^^^^^^^^^^^^^^^^^^

**Files**: ``src/plot_results.H``, ``src/write_xy_data.H``

* AMReX plotfile output for visualization with VisIt, ParaView, yt
* All simulation fields written to every plotfile (phi, vel, R, spotting
  diagnostics, Albini diagnostics, terrain fields, fireline intensity, flame length)
* XY data export: ``phi_negative_NNNN.dat`` (cells inside fire) and
  ``phi_envelope_NNNN.dat`` (convex hull of fire perimeter)

Data Structures
---------------

RothermelParams
^^^^^^^^^^^^^^^

Fuel properties for Rothermel model::

    struct RothermelParams {
        Real w0;          // Oven-dry total fuel loading [lb/ft²]
        Real sigma;       // Surface-area-to-volume ratio [ft⁻¹]
        Real delta;       // Fuel bed depth [ft]
        Real M_f;         // Aggregate fuel moisture content [fraction]
        Real M_x;         // Moisture of extinction [fraction]
        Real h_heat;      // Heat content [BTU/lb]
        Real S_T;         // Total mineral content [fraction]
        Real S_e;         // Effective mineral content [fraction]
        Real rho_p;       // Oven-dry particle density [lb/ft³]
        // Per-class moistures
        Real M_d1, M_d10, M_d100;  // Dead fuel moistures [fraction]
        Real M_lh, M_lw;           // Live fuel moistures [fraction]
        // Per-class fuel loads (from database)
        Real w_d1, sigma_d1;       // 1-hr dead
        Real w_d10, w_d100;        // 10-hr, 100-hr dead
        Real w_lh, sigma_lh;       // Live herbaceous
        Real w_lw, sigma_lw;       // Live woody
        // Terrain
        Real slope_x, slope_y;
        std::string terrain_file;
        std::string landscape_file;
        std::string landscape_fuel_type;  // "13" or "40"
        // Unit conversions
        Real wind_conv, ros_conv;
    };

RothermelComputed
^^^^^^^^^^^^^^^^^

Computed Rothermel parameters::

    struct RothermelComputed {
        Real R0;              // No-wind, no-slope ROS [ft/min]
        Real I_R;             // Reaction intensity [BTU/ft²/min]
        Real beta_ratio_E;    // Packing-ratio part of wind factor
        Real C;               // Wind factor coefficient C
        Real B;               // Wind factor coefficient B
        Real phi_s;           // Slope factor
        Real beta;            // Packing ratio
        Real wind_conv;       // Wind unit conversion
        Real ros_conv;        // ROS unit conversion
    };

GPU/CPU Execution
-----------------

The code is designed to run on both CPU and GPU using AMReX's unified memory model:

* GPU kernels are launched using ``ParallelFor``
* Data is transferred automatically between host and device
* CUDA, HIP, and SYCL backends are supported through AMReX

.. note::
   The Albini (1983) spotting model requires serial CPU execution and will abort
   if built with OpenMP or GPU backends enabled. Ensure the build is configured
   without these backends::

       cmake -S . -B build -DAMReX_OMP=OFF -DAMReX_GPU_BACKEND=NONE

Build System
------------

CMake Configuration
^^^^^^^^^^^^^^^^^^^

The project uses CMake (version 3.20+) with the following key options:

* ``LEVELSET_USE_VENDORED_AMREX``: Use the vendored AMReX submodule (default: ON)
* ``LEVELSET_DIM_2D``: Build for 2D instead of 3D (default: OFF)
* ``LEVELSET_ENABLE_EB``: Enable embedded boundary support (default: OFF)

AMReX is configured with:

* No MPI, OpenMP, or GPU backends by default (CPU-only)
* Minimal feature set (no particles, linear solvers, AmrLevel)
* Spatial dimension controlled by ``LEVELSET_DIM_2D``

.. note::
   Time-dependent wind fields are only available in 2D builds
   (``LEVELSET_DIM_2D=ON``).

Testing Framework
-----------------

The project uses CTest for regression testing. Each test case:

1. Runs the solver with specific input parameters
2. Compares output against reference data
3. Reports pass/fail status

Run all tests::

    cd build
    ctest

Run specific test::

    ctest -R basic_levelset
