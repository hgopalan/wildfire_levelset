Code Structure
==============

This section describes the organization and architecture of the wildfire level-set solver.

Directory Structure
-------------------

::

    wildfire_levelset/
    ├── CMakeLists.txt          # Main CMake build configuration
    ├── README.md               # Project README
    ├── LICENSE                 # License file
    ├── docs/                   # Documentation (this directory)
    ├── external/               # External dependencies
    │   └── amrex/              # AMReX submodule
    ├── src/                    # Source code
    │   ├── main.cpp            # Main entry point
    │   ├── parse_inputs.H      # Input parameter parsing
    │   ├── parse_inputs.cpp    # Input parameter implementation
    │   ├── advection.H         # Level-set advection
    │   ├── rothermel_model.H   # Rothermel fire spread model
    │   ├── farsite_ellipse.H   # FARSITE elliptical expansion
    │   ├── fire_models.H       # Fire model utilities
    │   ├── terrain_slope.H     # Terrain slope calculations
    │   ├── velocity_field.H    # Wind field management
    │   ├── crown_initiation.H  # Crown fire initiation
    │   ├── firebrand_spotting.H # Firebrand spotting model
    │   ├── bulk_fuel_consumption.H # Fuel consumption
    │   ├── initial_conditions.H # Initial fire setup
    │   ├── boundary_conditions.H # Boundary conditions
    │   ├── numerical_schemes.H # Numerical discretization
    │   ├── reinitialization.H  # Level-set reinitialization
    │   ├── compute_dt.H        # Time step calculation
    │   ├── plot_results.H      # Output/plotting
    │   ├── write_xy_data.H     # XY data export
    │   ├── compute_rothermel_R.H # Rothermel ROS
    │   └── fuel_database.H     # Fuel property database
    ├── regtest/                # Regression tests
    │   ├── basic_levelset/     # Basic level-set test
    │   ├── anderson_lw/        # Anderson L/W ratio test
    │   ├── farsite_ellipse/    # FARSITE ellipse test
    │   ├── rothermel_fuel/     # Rothermel fuel test
    │   ├── terrain_wind/       # Terrain/wind test
    │   ├── crown_initiation/   # Crown fire test
    │   ├── firebrand_spotting/ # Spotting test
    │   ├── bulk_fuel_consumption/ # Fuel consumption test
    │   └── ...                 # Other tests
    └── tests/                  # Additional test files

Core Components
---------------

Main Simulation Loop
^^^^^^^^^^^^^^^^^^^^

**File**: ``src/main.cpp``

The main simulation loop performs the following steps:

1. Parse input parameters from configuration file
2. Set up AMReX geometry and grid hierarchy
3. Initialize level-set function :math:`\phi` and velocity field
4. Time stepping loop:
   
   a. Compute time step :math:`\Delta t` based on CFL condition
   b. Compute fire spread rate using Rothermel model
   c. Advance level-set function or apply FARSITE ellipse
   d. Apply crown fire initiation if enabled
   e. Apply firebrand spotting if enabled
   f. Compute bulk fuel consumption if enabled
   g. Write output files
   h. Advance time :math:`t \leftarrow t + \Delta t`

Input Parsing
^^^^^^^^^^^^^

**Files**: ``src/parse_inputs.H``, ``src/parse_inputs.cpp``

The input parser reads configuration files and stores parameters in structured data types:

* ``RothermelParams``: Fuel properties for Rothermel model
* ``InputParameters::FARSITEParams``: FARSITE ellipse parameters
* ``InputParameters::CrownInitiationParams``: Crown fire parameters
* ``InputParameters::SpottingParams``: Firebrand spotting parameters
* ``InputParameters::FuelConsumptionParams``: Bulk fuel consumption parameters

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

Returns a ``RothermelComputed`` structure with all computed values.

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

Level-Set Advection
^^^^^^^^^^^^^^^^^^^

**File**: ``src/advection.H``

Implements the level-set advection equation:

.. math::

   \frac{\partial \phi}{\partial t} + V|\nabla\phi| = 0

* Uses upwind finite differences for spatial derivatives
* Explicit time integration
* Computes :math:`|\nabla\phi|` using central or upwind differences

Numerical Schemes
^^^^^^^^^^^^^^^^^

**File**: ``src/numerical_schemes.H``

Provides numerical discretization utilities:

* Upwind gradient computation
* CFL time step calculation
* Spatial derivative operators

Reinitialization
^^^^^^^^^^^^^^^^

**File**: ``src/reinitialization.H``

Maintains the level-set function as a signed distance function:

* Iterative reinitialization equation
* Ensures :math:`|\nabla\phi| = 1` away from the interface

Terrain Slope
^^^^^^^^^^^^^

**File**: ``src/terrain_slope.H``

Computes terrain slopes from elevation data:

* Reads X Y Z terrain data from file
* Interpolates elevation using inverse distance weighting (IDW)
* Computes slopes :math:`\partial z/\partial x` and :math:`\partial z/\partial y` using central differences

Velocity Field
^^^^^^^^^^^^^^

**File**: ``src/velocity_field.H``

Manages wind velocity field:

* Uniform or spatially-varying wind
* Time-dependent wind (future capability)
* Wind direction and magnitude

Crown Fire Initiation
^^^^^^^^^^^^^^^^^^^^^^

**File**: ``src/crown_initiation.H``

Implements Van Wagner (1977) crown fire model:

* Computes critical surface intensity :math:`I_0`
* Checks crown fire initiation criterion :math:`I > I_0`
* Modifies fire spread rate for active crown fire

Firebrand Spotting
^^^^^^^^^^^^^^^^^^

**File**: ``src/firebrand_spotting.H``

Simulates firebrand spotting:

* Probabilistic spotting model
* Distance-dependent spot fire probability
* Creates new ignition points ahead of main fire

Bulk Fuel Consumption
^^^^^^^^^^^^^^^^^^^^^^

**File**: ``src/bulk_fuel_consumption.H``

Computes post-frontal fuel consumption:

* Residence time calculation
* Fuel consumption fraction
* Heat release rate

Fuel Database
^^^^^^^^^^^^^

**File**: ``src/fuel_database.H``

Standard fuel model properties:

* Anderson 13 fuel models
* Custom fuel properties
* Fuel moisture scenarios

Initial Conditions
^^^^^^^^^^^^^^^^^^

**File**: ``src/initial_conditions.H``

Sets up initial fire configuration:

* Circular ignition
* Elliptical ignition
* Custom ignition patterns from file

Boundary Conditions
^^^^^^^^^^^^^^^^^^^

**File**: ``src/boundary_conditions.H``

Applies boundary conditions to level-set function:

* Neumann (zero gradient)
* Dirichlet (fixed value)
* Periodic

Output and Visualization
^^^^^^^^^^^^^^^^^^^^^^^^^

**Files**: ``src/plot_results.H``, ``src/write_xy_data.H``

* AMReX plotfile output for visualization with VisIt, ParaView, yt
* XY data export for time series and 2D plotting

Data Structures
---------------

RothermelParams
^^^^^^^^^^^^^^^

Fuel properties for Rothermel model::

    struct RothermelParams {
        Real w0;        // Oven-dry fuel loading [lb/ft²]
        Real sigma;     // Surface-area-to-volume ratio [ft⁻¹]
        Real delta;     // Fuel bed depth [ft]
        Real M_f;       // Fuel moisture content [fraction]
        Real M_x;       // Moisture of extinction [fraction]
        Real h_heat;    // Heat content [BTU/lb]
        Real S_T;       // Total mineral content [fraction]
        Real S_e;       // Effective mineral content [fraction]
        Real rho_p;     // Oven-dry particle density [lb/ft³]
        Real slope_x;   // Terrain slope in x-direction
        Real slope_y;   // Terrain slope in y-direction
        Real wind_conv; // Wind unit conversion factor
        Real ros_conv;  // ROS unit conversion factor
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
