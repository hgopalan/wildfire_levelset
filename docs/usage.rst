Usage Guide
===========

This section describes how to use the wildfire level-set solver to run simulations.

Input Files
-----------

The solver reads configuration from an input file (typically named ``inputs.i``). The input file uses AMReX's ParmParse format with key-value pairs.

Basic Input File Structure
^^^^^^^^^^^^^^^^^^^^^^^^^^^

A minimal input file looks like::

    # Domain and grid
    n_cell_x = 128
    n_cell_y = 128
    max_grid_size = 32
    prob_lo_x = 0.0
    prob_lo_y = 0.0
    prob_hi_x = 100.0
    prob_hi_y = 100.0

    # Time stepping
    nsteps = 100
    cfl = 0.5

    # Initial condition
    source_type = sphere
    sphere_center_x = 50.0
    sphere_center_y = 50.0
    sphere_center_z = 0.5
    sphere_radius = 5.0

    # Fire model
    fire_spread_model = rothermel
    propagation_method = farsite

    # Rothermel fuel properties
    rothermel.fuel_model = FM4
    rothermel.M_f = 0.08

Running Simulations
-------------------

Basic Execution
^^^^^^^^^^^^^^^

Run the solver with an input file::

    ./levelset inputs.i

The solver will:

1. Read the input file
2. Set up the grid and initial conditions
3. Run the time-stepping loop
4. Write output files to the ``plt/`` directory

Output Files
^^^^^^^^^^^^

The solver produces AMReX plotfiles in the ``plt/`` directory:

* ``plt0000/``: Initial condition
* ``plt0001/``: After first ``plot_int`` steps
* ``plt0002/``: After second ``plot_int`` steps
* etc.

Each plotfile directory contains:

* ``Header``: Metadata about the grid and variables
* ``Level_0/``: Cell-centered data files
* ``Cell_H``: Cell header information

Each plotfile contains the following fields:

+------------------------------+-----------------------------------+
| Field name                   | Description                       |
+==============================+===================================+
| ``phi``                      | Level-set function (fire front    |
|                              | where phi = 0)                    |
+------------------------------+-----------------------------------+
| ``velx``, ``vely``           | Wind velocity components [m/s]    |
+------------------------------+-----------------------------------+
| ``farsite_dx``, ``farsite_dy`` | FARSITE spread displacements    |
+------------------------------+-----------------------------------+
| ``R``                        | Rothermel rate of spread [m/s]    |
+------------------------------+-----------------------------------+
| ``spot_prob``, ``spot_count``| Spotting probability and count    |
+------------------------------+-----------------------------------+
| ``spot_dist``, ``spot_active`` | Spotting distance and flag      |
+------------------------------+-----------------------------------+
| ``fuel_consumption``         | Bulk fuel consumption fraction    |
+------------------------------+-----------------------------------+
| ``crown_fraction``           | Crown fire fraction               |
+------------------------------+-----------------------------------+
| ``albini_Hz``                | Albini lofting height [m]         |
+------------------------------+-----------------------------------+
| ``albini_count``             | Firebrands launched per cell      |
+------------------------------+-----------------------------------+
| ``albini_dist``              | Maximum firebrand landing dist [m]|
+------------------------------+-----------------------------------+
| ``albini_active``            | Spot ignition flag                |
+------------------------------+-----------------------------------+
| ``elevation``                | Terrain elevation [m]             |
+------------------------------+-----------------------------------+
| ``slope``                    | Terrain slope [degrees]           |
+------------------------------+-----------------------------------+
| ``aspect``                   | Terrain aspect [degrees]          |
+------------------------------+-----------------------------------+
| ``fuel_model``               | Fuel model number (from landscape)|
+------------------------------+-----------------------------------+
| ``fireline_intensity``       | Byram fireline intensity [kW/m]   |
+------------------------------+-----------------------------------+
| ``flame_length``             | Byram flame length [m]            |
+------------------------------+-----------------------------------+

In addition, two plain-text files are written at each output step:

* ``phi_negative_NNNN.dat``: X, Y coordinates of all cells where ``phi < 0`` (inside the fire)
* ``phi_envelope_NNNN.dat``: X, Y coordinates of the convex hull of the fire perimeter

Input Parameters Reference
---------------------------

Domain and Grid
^^^^^^^^^^^^^^^

**n_cell** (default: 64)
  Number of cells in all directions (cube). Overridden by per-direction values.

  Example: ``n_cell = 64``

**n_cell_x**, **n_cell_y**, **n_cell_z** (default: n_cell)
  Number of cells in each direction independently.

  Example: ``n_cell_x = 128``, ``n_cell_y = 128``

**max_grid_size** (default: 32)
  Maximum size of each grid box for domain decomposition.

  Example: ``max_grid_size = 64``

**prob_lo_x**, **prob_lo_y**, **prob_lo_z** (default: 0.0)
  Lower corner of physical domain in each direction.

  Example: ``prob_lo_x = 330000.0``, ``prob_lo_y = 3775000.0``

**prob_hi_x**, **prob_hi_y**, **prob_hi_z** (default: 1.0)
  Upper corner of physical domain in each direction.

  Example: ``prob_hi_x = 331000.0``, ``prob_hi_y = 3776000.0``

Time Stepping
^^^^^^^^^^^^^

**nsteps** (default: 300)
  Maximum number of time steps. Used when ``final_time`` is not set.

  Example: ``nsteps = 500``

**final_time** (default: -1.0)
  Simulation stop time in seconds. When positive, overrides ``nsteps``.

  Example: ``final_time = 3600.0``

**cfl** (default: 0.5)
  CFL number for time step calculation:

  .. math::

     \Delta t = \text{CFL} \times \frac{\Delta x}{V_{max}}

  Example: ``cfl = 0.7``

**plot_int** (default: 50)
  Plot interval (number of steps between output files).

  Example: ``plot_int = 10``

**reinit_int** (default: 20)
  Reinitialization interval (number of steps between level-set reinitialization).
  Set to -1 to disable.

  Example: ``reinit_int = 20``

Initial Conditions
^^^^^^^^^^^^^^^^^^

**source_type** (default: ``sphere``)
  Type of initial fire shape. Options: ``sphere``, ``box``, ``ellipse``, ``eb``.

  Example: ``source_type = box``

Sphere ignition
~~~~~~~~~~~~~~~

**center_x**, **center_y**, **center_z** (default: 0.5, 0.5, 0.5)
  Center coordinates of the initial spherical fire.

  Example: ``center_x = 330050.0``, ``center_y = 3775050.0``

**sphere_radius**
  Radius of the initial spherical fire [physical units].

  Example: ``sphere_radius = 10.0``

Box ignition
~~~~~~~~~~~~

**box_xmin**, **box_xmax**, **box_ymin**, **box_ymax**, **box_zmin**, **box_zmax**
  Bounds of the initial rectangular fire region.

  Example: ``box_xmin = 10.0``, ``box_xmax = 15.0``, ``box_ymin = 30.0``, ``box_ymax = 70.0``

Ellipse ignition
~~~~~~~~~~~~~~~~

**ellipse_center_x**, **ellipse_center_y**, **ellipse_center_z**
  Center coordinates of the initial elliptical fire.

  Example: ``ellipse_center_x = 0.4``, ``ellipse_center_y = 0.5``

**ellipse_radius_x**, **ellipse_radius_y**, **ellipse_radius_z**
  Semi-axes of the initial elliptical fire.

  Example: ``ellipse_radius_x = 0.25``, ``ellipse_radius_y = 0.15``

Embedded Boundary (EB) ignition
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use ``source_type = eb`` to initialize the fire from an implicit-function shape.
This is useful when the EB feature is enabled in the build.

**eb_type**
  Implicit function type. Options: ``sphere``, ``ellipsoid``, ``cylinder``, ``plane``.

  Example: ``eb_type = ellipsoid``

**eb_param1** – **eb_param6**
  Parameters for the chosen EB implicit function (center, radii, etc.).

  For an ellipsoid: ``eb_param1`` = center_x, ``eb_param2`` = center_y,
  ``eb_param3`` = center_z, ``eb_param4`` = radius_x, ``eb_param5`` = radius_y,
  ``eb_param6`` = radius_z.

  Example::

      eb_type = ellipsoid
      eb_param1 = 0.3
      eb_param2 = 0.5
      eb_param3 = 0.5
      eb_param4 = 0.20
      eb_param5 = 0.15
      eb_param6 = 0.10

CSV Fire Points ignition
~~~~~~~~~~~~~~~~~~~~~~~~

**fire_points_file** (default: "")
  Path to a CSV file containing ignition point coordinates (X Y [Z]), one point
  per row. Takes precedence over ``source_type`` when non-empty. The level-set
  field is initialized as a signed distance function from the union of small disks
  centred on each point.

  Example: ``fire_points_file = ignition_points.csv``

**fire_gaussian_sigma** (default: auto)
  Radius of each ignition disk in physical units. When ≤ 0, the radius is set
  automatically to 3 times the minimum cell spacing.

  Example: ``fire_gaussian_sigma = 15.0``

Fire Model Selection
^^^^^^^^^^^^^^^^^^^^

**fire_spread_model** (default: ``rothermel``)
  Select the rate-of-spread model. Currently supported values:

  - ``rothermel`` – Rothermel (1972) empirical fire spread model (default)
  - ``balbi`` – Balbi (2009) radiation-driven physics-based model
  - ``cheney_gould`` – Cheney & Gould (1995/1998) empirical grassland fire spread model

  Example: ``fire_spread_model = cheney_gould``

**propagation_method** (default: ``levelset``)
  Select how the fire perimeter is propagated. Currently supported values:

  - ``levelset`` – WENO5-Z level-set advection (default); Rothermel/Balbi ROS drives
    the signed-distance evolution equation.
  - ``farsite`` – FARSITE elliptical Huygens-wavelet expansion (Richards 1990);
    uses the FARSITE spread parameters (``farsite.*``) and pre-computed ROS.

  Example: ``propagation_method = farsite``

Rothermel Fuel Properties
^^^^^^^^^^^^^^^^^^^^^^^^^^

**rothermel.fuel_model** (default: custom chaparral)
  Fuel model name from the database (e.g. ``FM1``–``FM13`` for Anderson 13
  models, or ``FM101``–``FM256`` for Scott & Burgan 40 models).

  Example: ``rothermel.fuel_model = FM4``

**rothermel.M_f** (default: 0.08)
  Aggregate fuel moisture content (fraction, 0.0–1.0). Used as the default for
  all dead size classes when per-class moistures are not specified.

  Example: ``rothermel.M_f = 0.08``

**rothermel.M_x** (default: 0.30)
  Moisture of extinction (fraction).

  Example: ``rothermel.M_x = 0.25``

**rothermel.w0** (default: 0.230)
  Oven-dry total fuel loading [lb/ft²]. Overrides the database value.

  Example: ``rothermel.w0 = 0.3``

**rothermel.sigma** (default: 1739.0)
  Surface-area-to-volume ratio [ft⁻¹]. Overrides the database value.

  Example: ``rothermel.sigma = 1500.0``

**rothermel.delta** (default: 2.0)
  Fuel bed depth [ft].

  Example: ``rothermel.delta = 1.5``

**rothermel.h_heat** (default: 8000.0)
  Heat content [BTU/lb].

**rothermel.S_T** (default: 0.0555)
  Total mineral content (fraction).

**rothermel.S_e** (default: 0.010)
  Effective mineral content (fraction).

**rothermel.rho_p** (default: 32.0)
  Oven-dry particle density [lb/ft³].

**rothermel.slope_x**, **rothermel.slope_y** (default: 0.0)
  Constant terrain slope (tan of slope angle) in x and y directions,
  when no terrain file is used.

  Example: ``rothermel.slope_x = 0.1``

Per-Class Fuel Moisture (Multi-Class Rothermel Path)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When a fuel model from the database is selected (e.g.
``rothermel.fuel_model = FM4``), per-class fuel loads and SAV ratios are
populated automatically and the full multi-class Rothermel (1972) formulation
is used for reaction intensity and moisture damping.

The following per-class moisture inputs allow independent moisture contents
for each fuel size class. Dead-class defaults match ``rothermel.M_f``; live
defaults are physical equilibrium values.

**rothermel.M_d1** (default: equal to ``rothermel.M_f``)
  1-hr dead fuel moisture (fraction).  Fine sticks and fine herbaceous litter.

  Example: ``rothermel.M_d1 = 0.06``

**rothermel.M_d10** (default: equal to ``rothermel.M_f``)
  10-hr dead fuel moisture (fraction).  Small woody material (0.25–1 in).

  Example: ``rothermel.M_d10 = 0.08``

**rothermel.M_d100** (default: equal to ``rothermel.M_f``)
  100-hr dead fuel moisture (fraction).  Medium woody material (1–3 in).

  Example: ``rothermel.M_d100 = 0.10``

**rothermel.M_lh** (default: 0.90)
  Live herbaceous fuel moisture (fraction).  Green grass and forbs.

  Example: ``rothermel.M_lh = 0.60``

**rothermel.M_lw** (default: 1.20)
  Live woody fuel moisture (fraction).  Shrubs and small trees.

  Example: ``rothermel.M_lw = 1.00``

Per-Class Fuel Loads (Advanced Override)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

These parameters override the per-class loads provided by the fuel model
database entry.  They are not normally needed; use ``rothermel.fuel_model`` to
select a standard model instead.

* ``rothermel.w_d1``    — 1-hr dead fuel load [lb/ft²]
* ``rothermel.sigma_d1``— 1-hr dead SAV [ft⁻¹]
* ``rothermel.w_d10``   — 10-hr dead fuel load [lb/ft²]  (SAV = 109 ft⁻¹)
* ``rothermel.w_d100``  — 100-hr dead fuel load [lb/ft²] (SAV = 30 ft⁻¹)
* ``rothermel.w_lh``    — live herbaceous fuel load [lb/ft²]
* ``rothermel.sigma_lh``— live herbaceous SAV [ft⁻¹]
* ``rothermel.w_lw``    — live woody fuel load [lb/ft²]
* ``rothermel.sigma_lw``— live woody SAV [ft⁻¹]

When all per-class loads are zero the solver falls back to the single-class
aggregate path (``rothermel.w0`` / ``rothermel.sigma`` / ``rothermel.M_f``)
preserving backward compatibility.

Cheney & Gould (1995 / 1998) Grassland Spread Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

These parameters are used when ``fire_spread_model = cheney_gould``. The model
is an empirical piecewise-linear formula calibrated for open Australian
grasslands.

**cheney_gould.moisture** (default: 10.0)
  Dead fine fuel moisture content [%]. Used in the exponential moisture
  correction factor :math:`f_{MC} = \exp(-0.108 \times MC)`.

  Example: ``cheney_gould.moisture = 8.0``

**cheney_gould.curing** (default: 1.0)
  Degree of curing of the grass [0-1; 1 = fully cured dry standing grass,
  0 = completely green]. Values outside [0, 1] are clamped at runtime.

  Example: ``cheney_gould.curing = 0.90``

Wind Parameters
^^^^^^^^^^^^^^^

**u_x**, **u_y**, **u_z** (default: 0.25, 0.0, 0.0)
  Constant wind velocity components [m/s].

  Example: ``u_x = 2.0``, ``u_y = 0.5``

**velocity_file** (default: "")
  Path to a spatially-varying wind field CSV file (X Y U V columns).
  When set, overrides ``u_x``/``u_y``/``u_z``.

  Example: ``velocity_file = wind_data.csv``

**use_time_dependent_wind** (default: 0)
  Enable time-dependent wind fields (1=yes, 0=no). When enabled, the solver
  loads a sequence of wind field files and interpolates between them in time.
  Only available in 2D builds.

  Example: ``use_time_dependent_wind = 1``

**wind_time_spacing** (default: 60.0)
  Time spacing in seconds between consecutive wind field files. Used only when
  ``use_time_dependent_wind = 1``.

  Example: ``wind_time_spacing = 60.0``

Time-Dependent Wind Fields
~~~~~~~~~~~~~~~~~~~~~~~~~~~

When ``use_time_dependent_wind = 1``, the solver expects a series of wind field files
following this naming convention:

* Base file: ``<velocity_file>`` (time index 0)
* First update: ``<base>_1.csv`` (time index 1)
* Second update: ``<base>_2.csv`` (time index 2)
* And so on...

The solver will:

1. Determine which two wind field files bracket the current simulation time
2. Load both files (if not already loaded)
3. Perform spatial interpolation using inverse distance weighting for each field
4. Perform temporal linear interpolation between the two fields based on the current time

For example, with ``wind_time_spacing = 60.0``:

* At t=0s: Uses index-0 and index-1 files
* At t=30s: Interpolates 50% between indices 0 and 1
* At t=60s: Uses index-1 and index-2 files
* At t=120s: Uses index-2 and index-3 files

Each wind field file should be in CSV format with X, Y, U, V columns::

    # X Y U V
    0 0 1.0 0.0
    1000 0 1.0 0.0
    2000 0 1.0 0.0
    0 1000 1.2 0.1
    1000 1000 1.2 0.1
    2000 1000 1.2 0.1

Terrain Parameters
^^^^^^^^^^^^^^^^^^

**rothermel.terrain_file** (default: "")
  Path to a terrain elevation CSV file (X Y Z format). Computes per-cell slopes
  using inverse distance weighting interpolation and central differences.

  Example: ``rothermel.terrain_file = terrain.xyz``

**rothermel.landscape_file** (default: "")
  Path to a FARSITE landscape file (ASCII, X Y ELEVATION SLOPE ASPECT FUEL_MODEL).
  Provides elevation, slope, aspect, and fuel model data. Takes precedence over
  ``rothermel.terrain_file`` for slope and elevation; fuel model always comes from
  the landscape file when present.

  Example: ``rothermel.landscape_file = socal_chaparral_landscape.lcp``

**rothermel.landscape_fuel_type** (default: ``"13"``)
  Fuel model system used in the landscape file.  Must be ``"13"`` (Anderson 13 /
  FBFM13) or ``"40"`` (Scott & Burgan 40 / FBFM40).

  Example: ``rothermel.landscape_fuel_type = 40``

FARSITE Landscape File Format
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The landscape file is an ASCII text file with the following format::

    # Comments start with #
    X Y ELEVATION SLOPE ASPECT FUEL_MODEL

Where:

* **X, Y**: Coordinates in meters
* **ELEVATION**: Elevation above sea level in meters
* **SLOPE**: Slope angle in degrees (0-90)
* **ASPECT**: Slope aspect in degrees (0-360, where 0=North, 90=East, 180=South, 270=West)
* **FUEL_MODEL**: NFFL fuel model number (optional, defaults to 0)

Example landscape file::

    # Southern California chaparral terrain
    # X Y ELEVATION SLOPE ASPECT FUEL_MODEL
    0.0 0.0 100.0 15.0 225.0 4
    10.0 0.0 101.5 16.0 220.0 4
    20.0 0.0 103.5 18.0 215.0 4

The code uses inverse distance weighting (IDW) interpolation to map landscape data
points to the simulation grid. When a landscape file is specified, slope and elevation
from any terrain file are ignored.

FARSITE Parameters
^^^^^^^^^^^^^^^^^^

**farsite.use_anderson_LW** (default: 0)
  Use Anderson (1983) length-to-width ratio derived from wind speed (1=yes, 0=no).
  When disabled, uses ``farsite.length_to_width_ratio`` directly.

  Example: ``farsite.use_anderson_LW = 1``

**farsite.length_to_width_ratio** (default: 3.0)
  Fixed ellipse length-to-width ratio. Used when ``farsite.use_anderson_LW = 0``.

  Example: ``farsite.length_to_width_ratio = 3.0``

**farsite.coeff_a** (default: 1.0)
  Richards' head fire coefficient (maximum spread direction).

  Example: ``farsite.coeff_a = 1.0``

**farsite.coeff_b** (default: 0.5)
  Richards' flank fire coefficient (perpendicular to wind).

  Example: ``farsite.coeff_b = 0.5``

**farsite.coeff_c** (default: 0.2)
  Richards' backing fire coefficient (minimum spread, upwind direction).

  Example: ``farsite.coeff_c = 0.2``

**farsite.phi_threshold** (default: 0.1)
  Level-set value threshold for identifying fire front cells.

  Example: ``farsite.phi_threshold = 0.1``

Crown Fire Parameters
^^^^^^^^^^^^^^^^^^^^^

**crown.enable** (default: 0)
  Enable Van Wagner crown fire initiation model (1=yes, 0=no).

  Example: ``crown.enable = 1``

**crown.CBH** (default: 4.0)
  Canopy base height [m]. Must be > 0.

  Example: ``crown.CBH = 4.0``

**crown.CBD** (default: 0.15)
  Canopy bulk density [kg/m³]. Must be > 0.

  Example: ``crown.CBD = 0.15``

**crown.FMC** (default: 100.0)
  Foliar moisture content [%]. Must be between 50 and 300.

  Example: ``crown.FMC = 100.0``

**crown.crown_fraction_weight** (default: 1.0)
  Crown fire weighting factor (0.0–2.0). Scales the crown fire contribution to
  the spread rate.

  Example: ``crown.crown_fraction_weight = 1.0``

**crown.use_metric_units** (default: 1)
  Unit system for crown fire calculations. 1 = metric (m, kW/m), 0 = imperial.

  Example: ``crown.use_metric_units = 1``

Bulk Fuel Consumption Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The bulk fuel consumption model tracks the fraction of fuel consumed behind the
fire front. It is configured under the ``farsite`` namespace.

**farsite.use_bulk_fuel_consumption** (default: 0)
  Enable bulk fuel consumption calculation (1=yes, 0=no).

  Example: ``farsite.use_bulk_fuel_consumption = 1``

**farsite.tau_residence** (default: 60.0)
  Residence time for fuel consumption [seconds]. Must be > 0.

  Example: ``farsite.tau_residence = 60.0``

**farsite.f_consumed_max** (default: 0.9)
  Maximum fuel consumption fraction (0.0–1.0). Applies to slow-moving,
  high-intensity fires.

  Example: ``farsite.f_consumed_max = 0.9``

**farsite.f_consumed_min** (default: 0.5)
  Minimum fuel consumption fraction (0.0–1.0). Applies to fast-moving fires.
  Must be ≤ ``farsite.f_consumed_max``.

  Example: ``farsite.f_consumed_min = 0.5``

Probability-Based Firebrand Spotting Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**spotting.enable** (default: 0)
  Enable the stochastic firebrand spotting model (1=yes, 0=no).

  Example: ``spotting.enable = 1``

**spotting.P_base** (default: 0.02)
  Base probability of generating a firebrand per fire-front cell (0.0–1.0).

  Example: ``spotting.P_base = 0.05``

**spotting.k_wind** (default: 0.3)
  Wind speed coefficient for spotting probability scaling.

  Example: ``spotting.k_wind = 0.3``

**spotting.I_critical** (default: 1000.0)
  Critical fire intensity threshold for spotting [kW/m or BTU/ft/s].

  Example: ``spotting.I_critical = 800.0``

**spotting.d_mean** (default: 0.1)
  Mean spotting distance parameter [physical units].

  Example: ``spotting.d_mean = 15.0``

**spotting.d_sigma** (default: 0.5)
  Standard deviation for lognormal spotting distance distribution.

  Example: ``spotting.d_sigma = 40.0``

**spotting.d_lambda** (default: 10.0)
  Decay rate for exponential spotting distance distribution [1/unit].

  Example: ``spotting.d_lambda = 0.08``

**spotting.distance_model** (default: ``"lognormal"``)
  Distribution model for spotting distance. Options: ``"lognormal"`` or
  ``"exponential"``.

  Example: ``spotting.distance_model = lognormal``

**spotting.lateral_spread_angle** (default: 15.0)
  Angular spread perpendicular to wind direction [degrees].

  Example: ``spotting.lateral_spread_angle = 20.0``

**spotting.spot_radius** (default: 0.02)
  Radius of new spot-fire ignition zones [physical units].

  Example: ``spotting.spot_radius = 3.0``

**spotting.random_seed** (default: 0)
  Seed for the random number generator. 0 uses the system clock (non-reproducible).

  Example: ``spotting.random_seed = 12345``

**spotting.check_interval** (default: 5)
  Run the spotting model every N time steps.

  Example: ``spotting.check_interval = 3``

Albini (1983) Physics-Based Spotting Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The Albini spotting model uses a physics-based lofting height from Byram's fire
line intensity and integrates a 2-D horizontal firebrand trajectory through the
wind field.

**albini_spotting.enable** (default: 0)
  Enable the Albini (1983) firebrand spotting model (1=yes, 0=no).

  Example: ``albini_spotting.enable = 1``

**albini_spotting.terminal_velocity** (default: 1.0)
  Firebrand terminal descent velocity [m/s]. Must be > 0.

  Example: ``albini_spotting.terminal_velocity = 5.0``

**albini_spotting.P_base** (default: 0.01)
  Maximum probability of launching a firebrand per fire-front cell (0.0–1.0).
  The actual probability is weighted by fire line intensity.

  Example: ``albini_spotting.P_base = 0.05``

**albini_spotting.I_B_min** (default: 10.0)
  Minimum Byram fire line intensity to allow firebrand generation [kW/m].

  Example: ``albini_spotting.I_B_min = 500.0``

**albini_spotting.spot_radius** (default: 5.0)
  Radius of new spot-fire ignition zone [m].

  Example: ``albini_spotting.spot_radius = 15.0``

**albini_spotting.random_seed** (default: 0)
  Seed for the random number generator. 0 uses the system clock.

  Example: ``albini_spotting.random_seed = 42``

**albini_spotting.check_interval** (default: 5)
  Run the Albini spotting model every N time steps.

  Example: ``albini_spotting.check_interval = 5``

**albini_spotting.n_traj_steps** (default: 100)
  Number of forward-Euler sub-steps for the 2-D firebrand trajectory integration.

  Example: ``albini_spotting.n_traj_steps = 200``

Output Parameters
^^^^^^^^^^^^^^^^^

**plot_int** (default: 50)
  Number of time steps between output files.

  Example: ``plot_int = 10``

Wind-Terrain Feedback Model Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**wind_terrain.model** (default: ``"none"``)
  Selects the wind-terrain feedback model. Seven options are available:

  - ``none`` – Option 1: default behaviour; no wind modification
  - ``viegas_ros`` – Option 2: override spread rate with
    :math:`\max(R_{\text{primary}}, R_V)` in eruptive cells; works with
    both Rothermel and Balbi (uses Balbi amplitude baseline when
    ``fire_spread_model = balbi``); auto-enables ``viegas.enable = 1``
  - ``viegas_wind`` – Option 3: add Viegas buoyancy-induced upslope wind
    only in eruptive cells (:math:`\tan\varphi > \tan\varphi_c`)
  - ``canyon_wind`` – Option 4: Rothermel (1983) canyon wind amplification
    :math:`U_{\text{eff}} = U (1 + k_{\text{canyon}} \tan\varphi)`
  - ``viegas_neto`` – Option 5: Viegas & Neto (1994) buoyancy-driven
    upslope wind at all cells
  - ``pimont`` – Option 6: Pimont et al. (2009) exponential slope correction
    :math:`U_{\text{eff}} = U \exp(k_{\text{pimont}} \tan\varphi)`
  - ``windninja_ridge_canyon`` – Option 7: WindNinja empirical ridge/canyon
    speed-up based on wind-slope alignment (see below)

  Example: ``wind_terrain.model = windninja_ridge_canyon``

**wind_terrain.k_canyon** (default: 1.0)
  Terrain channeling coefficient for Option 4 (``canyon_wind``). Must be > 0.

  Example: ``wind_terrain.k_canyon = 1.5``

**wind_terrain.k_pimont** (default: 0.5)
  Exponential slope correction coefficient for Option 6 (``pimont``). Must be > 0.

  Example: ``wind_terrain.k_pimont = 0.5``

**wind_terrain.k_ridge** (default: 1.0)
  Ridge speed-up coefficient for Option 7 (``windninja_ridge_canyon``).
  Scales the speed-up when wind climbs a slope:
  :math:`f = 1 + k_{\text{ridge}} \tan\varphi \cdot \text{alignment}`. Must be > 0.

  Example: ``wind_terrain.k_ridge = 1.5``

**wind_terrain.k_canyon_wn** (default: 0.5)
  Canyon channeling coefficient for Option 7 (``windninja_ridge_canyon``).
  Scales the amplification when wind flows downslope:
  :math:`f = 1 + k_{\text{canyon\_wn}} \tan\varphi \cdot |\text{alignment}|`. Must be > 0.

  Example: ``wind_terrain.k_canyon_wn = 0.8``

Heat Flux Parameters
^^^^^^^^^^^^^^^^^^^^^

The heat flux MultiFab provides a spatially-varying (or uniform) fire heat
release rate :math:`Q` [W/m²] that drives WindNinja-style fire-induced wind
corrections and augments the Balbi buoyancy velocity.

**heat_flux.value** (default: 0.0)
  Uniform heat flux [W/m²] applied to the entire domain. Set to a positive
  value to activate. Requires at least one of ``enable_upward`` or
  ``enable_induced`` to also be set to 1.

  Example: ``heat_flux.value = 5000.0``

**heat_flux.file** (default: "")
  Path to an ASCII file with columns ``X Y Q`` for spatially-varying heat flux.
  Uses the same Wendland C² IDW interpolation as terrain and wind files.
  Overrides ``heat_flux.value`` when non-empty. Only available in 2D builds.

  Example: ``heat_flux.file = fire_heat_flux.csv``

**heat_flux.rho_air** (default: 1.2)
  Air density [kg/m³].

  Example: ``heat_flux.rho_air = 1.1``

**heat_flux.Cp_air** (default: 1005.0)
  Specific heat of air [J/(kg·K)].

  Example: ``heat_flux.Cp_air = 1005.0``

**heat_flux.T_a** (default: 300.0)
  Ambient temperature [K] for the buoyancy velocity calculation.

  Example: ``heat_flux.T_a = 310.0``

**heat_flux.ref_height** (default: 10.0)
  Reference height [m] for the buoyancy velocity calculation.
  Physically represents the plume height over which thermal energy is
  distributed.

  Example: ``heat_flux.ref_height = 20.0``

**heat_flux.k_upward** (default: 1.0)
  Upward velocity coefficient. Scales the magnitude of the fire-plume
  vertical velocity:
  :math:`w_{\uparrow} = k_{\uparrow} \sqrt{g Q H_{\text{ref}} / (\rho C_p T_a)}`.

  Example: ``heat_flux.k_upward = 0.8``

**heat_flux.k_induced** (default: 0.5)
  Induced horizontal inflow coefficient. The inflow toward the fire
  perimeter is :math:`U_{\text{ind}} = k_{\text{ind}} \, w_{\uparrow}`.

  Example: ``heat_flux.k_induced = 0.3``

**heat_flux.enable_upward** (default: 0)
  Set to 1 to enable the upward convective velocity term (adds to the
  wind field; for 3D adds to z-component, for 2D projects as inflow
  opposite to ambient wind direction).

  Example: ``heat_flux.enable_upward = 1``

**heat_flux.enable_induced** (default: 0)
  Set to 1 to enable the induced horizontal inflow term directed toward
  the fire perimeter (requires phi MultiFab to be available).

  Example: ``heat_flux.enable_induced = 1``

Example Input Files
-------------------

Basic Level-Set Simulation
^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    # Basic circular fire with level-set advection
    n_cell_x = 64
    n_cell_y = 64
    max_grid_size = 32
    prob_lo_x = 0.0
    prob_lo_y = 0.0
    prob_hi_x = 100.0
    prob_hi_y = 100.0

    final_time = 600.0
    cfl = 0.5
    plot_int = 20
    reinit_int = 20

    source_type = sphere
    center_x = 50.0
    center_y = 50.0
    center_z = 0.5
    sphere_radius = 5.0

    rothermel.fuel_model = FM1
    rothermel.M_f = 0.05
    u_x = 0.25

    propagation_method = levelset

FARSITE Ellipse Simulation
^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    # FARSITE elliptical expansion with Anderson L/W
    n_cell_x = 100
    n_cell_y = 100
    max_grid_size = 50
    prob_lo_x = 330000.0
    prob_lo_y = 3775000.0
    prob_hi_x = 330200.0
    prob_hi_y = 3775200.0

    final_time = 3600.0
    cfl = 0.5
    plot_int = 10

    source_type = sphere
    center_x = 330100.0
    center_y = 3775100.0
    center_z = 0.5
    sphere_radius = 10.0

    fire_spread_model = rothermel
    propagation_method = farsite
    farsite.use_anderson_LW = 1
    farsite.phi_threshold = 0.1

    rothermel.fuel_model = FM4
    rothermel.M_f = 0.08
    u_x = 2.0
    u_y = 0.5

Terrain and Crown Fire Simulation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    # Complex simulation with terrain, wind, and crown fire
    n_cell_x = 100
    n_cell_y = 100
    max_grid_size = 50
    prob_lo_x = 330000.0
    prob_lo_y = 3775000.0
    prob_hi_x = 331000.0
    prob_hi_y = 3776000.0

    final_time = 7200.0
    cfl = 0.5
    plot_int = 25

    source_type = sphere
    center_x = 330500.0
    center_y = 3775500.0
    center_z = 0.5
    sphere_radius = 20.0

    fire_spread_model = rothermel
    propagation_method = farsite
    farsite.use_anderson_LW = 1

    rothermel.fuel_model = FM10
    rothermel.M_f = 0.06
    rothermel.terrain_file = elevation.xyz
    u_x = 0.5
    u_y = 0.2

    crown.enable = 1
    crown.CBH = 4.0
    crown.CBD = 0.2
    crown.FMC = 100.0
    crown.crown_fraction_weight = 1.0

    spotting.enable = 1
    spotting.P_base = 0.05
    spotting.d_mean = 50.0
    spotting.distance_model = lognormal

FARSITE with Landscape File
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    # FARSITE simulation with landscape file (Southern California terrain)
    n_cell_x = 50
    n_cell_y = 50
    max_grid_size = 32
    prob_lo_x = 0.0
    prob_lo_y = 0.0
    prob_hi_x = 100.0
    prob_hi_y = 100.0

    final_time = 3000.0
    cfl = 0.5
    plot_int = 10

    # Line fire ignition (simulating Santa Ana wind-driven fire)
    source_type = box
    box_xmin = 10.0
    box_xmax = 15.0
    box_ymin = 30.0
    box_ymax = 70.0

    # Wind from west (Santa Ana conditions)
    u_x = 5.0
    u_y = 0.0

    # FARSITE model with Anderson L/W ratio
    fire_spread_model = rothermel
    propagation_method = farsite
    farsite.use_anderson_LW = 1
    farsite.phi_threshold = 0.1

    # Use landscape file for terrain (slope, aspect, elevation) and fuel model
    rothermel.fuel_model = FM4
    rothermel.M_f = 0.08
    # Optional: differentiate dead/live moisture by size class
    rothermel.M_d1   = 0.06   # 1-hr dead moisture
    rothermel.M_d10  = 0.08   # 10-hr dead moisture
    rothermel.M_d100 = 0.10   # 100-hr dead moisture
    rothermel.M_lh   = 0.60   # live herbaceous moisture
    rothermel.M_lw   = 1.00   # live woody moisture
    rothermel.landscape_file = socal_chaparral_landscape.lcp

Albini Spotting Simulation
^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    # FARSITE with Albini (1983) physics-based firebrand spotting
    n_cell_x = 100
    n_cell_y = 100
    max_grid_size = 50
    prob_lo_x = 330000.0
    prob_lo_y = 3775000.0
    prob_hi_x = 331000.0
    prob_hi_y = 3776000.0

    final_time = 1000.0
    cfl = 0.5
    plot_int = 25
    reinit_int = -1

    source_type = box
    box_xmin = 330050.0
    box_xmax = 330100.0
    box_ymin = 3775300.0
    box_ymax = 3775700.0

    u_x = 2.0
    u_y = 0.5

    rothermel.fuel_model = FM4
    rothermel.M_f = 0.08

    fire_spread_model = rothermel
    propagation_method = farsite
    farsite.use_anderson_LW = 1
    farsite.phi_threshold = 0.1

    albini_spotting.enable = 1
    albini_spotting.terminal_velocity = 5.0
    albini_spotting.P_base = 0.05
    albini_spotting.I_B_min = 500.0
    albini_spotting.spot_radius = 15.0
    albini_spotting.random_seed = 42
    albini_spotting.check_interval = 5
    albini_spotting.n_traj_steps = 200

Balbi + Viegas + Heat Flux Simulation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    # Balbi physical fire spread with Viegas slope diagnostics and heat flux
    # wind corrections on a Gaussian hill terrain.
    n_cell_x = 100
    n_cell_y = 100
    prob_lo_x = 330000.0
    prob_lo_y = 3775000.0
    prob_hi_x = 331000.0
    prob_hi_y = 3776000.0

    final_time = 600.0
    cfl = 0.5
    plot_int = 20
    reinit_int = 20

    source_type = box
    box_xmin = 330100.0
    box_xmax = 330120.0
    box_ymin = 3775200.0
    box_ymax = 3775800.0

    velocity_file = gaussian_hill_wind.csv

    # Balbi (2009) physical fire spread model
    fire_spread_model = balbi
    propagation_method = levelset

    rothermel.fuel_model = FM4
    rothermel.M_f = 0.08
    rothermel.terrain_file = gaussian_hill_terrain.csv

    balbi.T_a   = 300.0
    balbi.T_f   = 1000.0
    balbi.T_i   = 600.0

    # Viegas Option 2: Viegas-Balbi ROS overrides Balbi in eruptive cells
    # viegas.enable is auto-set to 1 by wind_terrain.model = viegas_ros
    wind_terrain.model = viegas_ros
    viegas.a_V       = 1.83
    viegas.tan_phi_c = 0.4

    # Heat flux wind corrections
    heat_flux.value         = 5000.0    # W/m^2
    heat_flux.ref_height    = 10.0      # m
    heat_flux.k_upward      = 1.0
    heat_flux.enable_upward = 1         # adds upward velocity to wind field

WindNinja Ridge/Canyon Simulation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    # Rothermel with WindNinja empirical ridge/canyon wind speed-up
    n_cell_x = 100
    n_cell_y = 100
    prob_lo_x = 330000.0
    prob_lo_y = 3775000.0
    prob_hi_x = 331000.0
    prob_hi_y = 3776000.0

    final_time = 600.0
    cfl = 0.5
    plot_int = 20

    velocity_file = gaussian_hill_wind.csv
    rothermel.terrain_file = gaussian_hill_terrain.csv

    source_type = box
    box_xmin = 330100.0
    box_xmax = 330120.0
    box_ymin = 3775200.0
    box_ymax = 3775800.0

    rothermel.fuel_model = FM4
    rothermel.M_f = 0.08

    fire_spread_model = rothermel
    propagation_method = levelset

    # Option 7: ridge speed-up on windward face, canyon channeling on lee face
    wind_terrain.model       = windninja_ridge_canyon
    wind_terrain.k_ridge     = 1.5    # ridge amplification
    wind_terrain.k_canyon_wn = 0.8   # canyon channeling

Visualization
-------------

VisIt
^^^^^

Open plotfiles in VisIt::

    visit -o plt00000/Header

ParaView
^^^^^^^^

Open plotfiles in ParaView using the AMReX reader.

yt
^^

Use yt to analyze plotfiles in Python::

    import yt
    ds = yt.load("plt00000/Header")
    slc = yt.SlicePlot(ds, 'z', 'phi')
    slc.save()

Python Scripts
^^^^^^^^^^^^^^

Read XY data files with NumPy/Pandas::

    import numpy as np
    data = np.loadtxt('xy_data.txt')
    x = data[:, 0]
    y = data[:, 1]
    phi = data[:, 2]

Troubleshooting
---------------

Simulation Crashes or Diverges
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Reduce the CFL number (try ``cfl = 0.3``)
* Use a finer grid (increase ``n_cell``)
* Check that fuel properties are physical (positive values, moisture < 1)

Output Files Not Created
^^^^^^^^^^^^^^^^^^^^^^^^^

* Check that ``plot_int`` or ``plot_dt`` is set to a positive value
* Ensure the simulation runs for at least one output interval
* Check disk space and write permissions

Performance Issues
^^^^^^^^^^^^^^^^^^

* Increase ``max_grid_size`` for better cache performance
* Reduce grid resolution if simulation is too slow
* Use Release build type for optimization
