Usage Guide
===========

This section describes how to use the wildfire level-set solver to run simulations.

Input Files
-----------

The solver reads configuration from an input file (typically named ``inputs``). The input file uses AMReX's ParmParse format with key-value pairs.

Basic Input File Structure
^^^^^^^^^^^^^^^^^^^^^^^^^^^

A minimal input file looks like::

    # Domain and grid
    n_cell = 128 128 1
    max_grid_size = 32
    prob_lo = 0.0 0.0 0.0
    prob_hi = 1.0 1.0 1.0
    
    # Time stepping
    max_step = 100
    stop_time = 10.0
    cfl = 0.5
    
    # Initial condition
    init_type = sphere
    init_center = 0.5 0.5 0.5
    init_radius = 0.1
    
    # Fire model parameters
    use_farsite_model = 0
    use_terrain_effects = 0
    
    # Rothermel fuel properties
    fuel_model = 1
    fuel_moisture = 0.05

Running Simulations
-------------------

Basic Execution
^^^^^^^^^^^^^^^

Run the solver with an input file::

    ./levelset inputs

The solver will:

1. Read the input file
2. Set up the grid and initial conditions
3. Run the time-stepping loop
4. Write output files to the ``plt/`` directory

Output Files
^^^^^^^^^^^^

The solver produces AMReX plotfiles in the ``plt/`` directory:

* ``plt00000/``: Initial condition
* ``plt00001/``: First output time
* ``plt00002/``: Second output time
* etc.

Each plotfile directory contains:

* ``Header``: Metadata about the grid and variables
* ``Level_0/``: Cell-centered data files
* ``Cell_H``: Cell header information

Input Parameters Reference
---------------------------

Domain and Grid
^^^^^^^^^^^^^^^

**n_cell** (required)
  Number of cells in each direction. For 2D, set z-dimension to 1.
  
  Example: ``n_cell = 128 128 1``

**max_grid_size** (default: 32)
  Maximum size of each grid box for domain decomposition.
  
  Example: ``max_grid_size = 64``

**prob_lo** (required)
  Lower corner of physical domain.
  
  Example: ``prob_lo = 0.0 0.0 0.0``

**prob_hi** (required)
  Upper corner of physical domain.
  
  Example: ``prob_hi = 100.0 100.0 10.0``

**is_periodic** (default: 0 0 0)
  Periodic boundaries in each direction (1=periodic, 0=not periodic).
  
  Example: ``is_periodic = 0 0 1``

Time Stepping
^^^^^^^^^^^^^

**max_step** (default: 1000)
  Maximum number of time steps.
  
  Example: ``max_step = 500``

**stop_time** (default: 1.0)
  Simulation stop time.
  
  Example: ``stop_time = 3600.0``

**cfl** (default: 0.5)
  CFL number for time step calculation:
  
  .. math::
  
     \Delta t = \text{CFL} \times \frac{\Delta x}{V_{max}}
  
  Example: ``cfl = 0.7``

**dt_fixed** (optional)
  Fixed time step (overrides CFL-based calculation).
  
  Example: ``dt_fixed = 0.01``

Initial Conditions
^^^^^^^^^^^^^^^^^^

**init_type** (required)
  Type of initial fire shape: ``sphere``, ``ellipse``, ``plane``, or ``file``.
  
  Example: ``init_type = sphere``

**init_center** (required for sphere/ellipse)
  Center of initial fire.
  
  Example: ``init_center = 50.0 50.0 0.5``

**init_radius** (required for sphere)
  Radius of initial circular/spherical fire.
  
  Example: ``init_radius = 5.0``

**init_semiaxes** (required for ellipse)
  Semi-axes of elliptical fire (a, b, c).
  
  Example: ``init_semiaxes = 10.0 5.0 1.0``

Fire Model Selection
^^^^^^^^^^^^^^^^^^^^

**use_farsite_model** (default: 0)
  Enable FARSITE elliptical expansion model (1) or use level-set advection (0).
  
  Example: ``use_farsite_model = 1``

**use_terrain_effects** (default: 0)
  Enable terrain slope effects on fire spread (1=yes, 0=no).
  
  Example: ``use_terrain_effects = 1``

Rothermel Fuel Properties
^^^^^^^^^^^^^^^^^^^^^^^^^^

**fuel_model** (default: 1)
  Anderson 13 fuel model number (1-13), or 0 for custom fuel.
  
  Example: ``fuel_model = 3``

**fuel_moisture** (default: 0.05)
  Fuel moisture content (fraction, 0.0-1.0).
  
  Example: ``fuel_moisture = 0.08``

**moisture_extinction** (default: 0.30)
  Moisture of extinction (fraction).
  
  Example: ``moisture_extinction = 0.25``

For custom fuel (``fuel_model = 0``), specify all fuel properties:

* ``fuel_load``: Oven-dry fuel loading [lb/ft²]
* ``fuel_sav``: Surface-area-to-volume ratio [ft⁻¹]
* ``fuel_depth``: Fuel bed depth [ft]
* ``fuel_heat``: Heat content [BTU/lb]
* ``fuel_density``: Particle density [lb/ft³]
* ``fuel_st``: Total mineral content [fraction]
* ``fuel_se``: Effective mineral content [fraction]

Wind Parameters
^^^^^^^^^^^^^^^

**wind_speed** (default: 0.0)
  Midflame wind speed [mph or m/s depending on units].
  
  Example: ``wind_speed = 10.0``

**wind_direction** (default: 0.0)
  Wind direction [degrees, 0=east, 90=north].
  
  Example: ``wind_direction = 270.0``

**wind_from_file** (default: 0)
  Read spatially-varying wind from file (1=yes, 0=no).
  
  Example: ``wind_from_file = 1``

**wind_file** (required if wind_from_file=1)
  Path to wind data file.
  
  Example: ``wind_file = "wind_data.txt"``

**use_time_dependent_wind** (default: 0)
  Enable time-dependent wind fields (1=yes, 0=no). When enabled, the solver will 
  load a sequence of wind field files and interpolate between them in time.
  
  Example: ``use_time_dependent_wind = 1``

**wind_time_spacing** (default: 60.0)
  Time spacing in seconds between consecutive wind field files. This parameter is 
  only used when ``use_time_dependent_wind = 1``.
  
  Example: ``wind_time_spacing = 60.0``

Time-Dependent Wind Fields
~~~~~~~~~~~~~~~~~~~~~~~~~~~

When ``use_time_dependent_wind = 1``, the solver expects a series of wind field files 
following this naming convention:

* Base file: ``turbulent_field_2d.csv`` (time index 0)
* First update: ``turbulent_field_2d_1.csv`` (time index 1)
* Second update: ``turbulent_field_2d_2.csv`` (time index 2)
* And so on...

The solver will:

1. Determine which two wind field files bracket the current simulation time
2. Load both files (if not already loaded)
3. Perform spatial interpolation using inverse distance weighting for each field
4. Perform temporal linear interpolation between the two fields based on the current time

For example, with ``wind_time_spacing = 60.0``:

* At t=0s: Uses ``turbulent_field_2d.csv`` and ``turbulent_field_2d_1.csv``
* At t=30s: Interpolates 50% between indices 0 and 1
* At t=60s: Uses ``turbulent_field_2d_1.csv`` and ``turbulent_field_2d_2.csv``
* At t=120s: Uses ``turbulent_field_2d_2.csv`` and ``turbulent_field_2d_3.csv``

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

**terrain_from_file** (default: 0)
  Read terrain elevation from file (1=yes, 0=no).
  
  Example: ``terrain_from_file = 1``

**terrain_file** (required if terrain_from_file=1)
  Path to terrain data file (X Y Z format).
  
  Example: ``terrain_file = "terrain.xyz"``

FARSITE Parameters
^^^^^^^^^^^^^^^^^^

**farsite.enable** (default: 0)
  Enable FARSITE elliptical model (same as use_farsite_model).
  
  Example: ``farsite.enable = 1``

**farsite.use_anderson_LW** (default: 0)
  Use Anderson (1983) length-to-width ratio (1=yes, 0=no).
  
  Example: ``farsite.use_anderson_LW = 1``

**farsite.coeff_a** (default: 1.0)
  Richards' head fire coefficient.
  
  Example: ``farsite.coeff_a = 1.5``

**farsite.coeff_b** (default: 1.0)
  Richards' flank fire coefficient.
  
  Example: ``farsite.coeff_b = 1.0``

**farsite.coeff_c** (default: 0.5)
  Richards' backing fire coefficient.
  
  Example: ``farsite.coeff_c = 0.3``

**farsite.phi_threshold** (default: 0.1)
  Threshold for identifying fire front (level-set value).
  
  Example: ``farsite.phi_threshold = 0.05``

Crown Fire Parameters
^^^^^^^^^^^^^^^^^^^^^

**crown.enable** (default: 0)
  Enable crown fire initiation (1=yes, 0=no).
  
  Example: ``crown.enable = 1``

**crown.canopy_base_height** (default: 2.0)
  Canopy base height [m].
  
  Example: ``crown.canopy_base_height = 3.0``

**crown.canopy_bulk_density** (default: 0.1)
  Canopy bulk density [kg/m³].
  
  Example: ``crown.canopy_bulk_density = 0.15``

**crown.canopy_moisture** (default: 0.10)
  Canopy moisture content [fraction].
  
  Example: ``crown.canopy_moisture = 0.08``

**crown.crown_consumption_depth** (default: 1.0)
  Crown fuel consumption depth [m].
  
  Example: ``crown.crown_consumption_depth = 1.5``

Firebrand Spotting Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**spotting.enable** (default: 0)
  Enable firebrand spotting (1=yes, 0=no).
  
  Example: ``spotting.enable = 1``

**spotting.probability** (default: 0.01)
  Base probability of firebrand generation per cell.
  
  Example: ``spotting.probability = 0.05``

**spotting.max_distance** (default: 100.0)
  Maximum spotting distance [m].
  
  Example: ``spotting.max_distance = 500.0``

**spotting.wind_speed_threshold** (default: 5.0)
  Minimum wind speed for spotting [m/s].
  
  Example: ``spotting.wind_speed_threshold = 10.0``

Bulk Fuel Consumption Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**fuel_consumption.enable** (default: 0)
  Enable bulk fuel consumption calculation (1=yes, 0=no).
  
  Example: ``fuel_consumption.enable = 1``

**fuel_consumption.residence_time** (default: 60.0)
  Residence time for fuel consumption [seconds].
  
  Example: ``fuel_consumption.residence_time = 120.0``

Output Parameters
^^^^^^^^^^^^^^^^^

**plot_int** (default: -1)
  Plot interval (number of steps between outputs). -1 disables.
  
  Example: ``plot_int = 10``

**plot_dt** (default: -1.0)
  Plot time interval (simulation time between outputs). -1 disables.
  
  Example: ``plot_dt = 60.0``

**plot_file** (default: "plt")
  Prefix for plot file directories.
  
  Example: ``plot_file = "output_plt"``

**write_xy_data** (default: 0)
  Write XY data files for 2D plotting (1=yes, 0=no).
  
  Example: ``write_xy_data = 1``

**xy_data_file** (default: "xy_data.txt")
  Filename for XY data output.
  
  Example: ``xy_data_file = "fire_perimeter.txt"``

Example Input Files
-------------------

Basic Level-Set Simulation
^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    # Basic circular fire with level-set advection
    n_cell = 128 128 1
    max_grid_size = 32
    prob_lo = 0.0 0.0 0.0
    prob_hi = 100.0 100.0 1.0
    
    max_step = 200
    stop_time = 600.0
    cfl = 0.5
    plot_int = 10
    
    init_type = sphere
    init_center = 50.0 50.0 0.5
    init_radius = 5.0
    
    fuel_model = 1
    fuel_moisture = 0.05
    wind_speed = 5.0
    wind_direction = 90.0

FARSITE Ellipse Simulation
^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    # FARSITE elliptical expansion with Anderson L/W
    n_cell = 256 256 1
    max_grid_size = 64
    prob_lo = 0.0 0.0 0.0
    prob_hi = 200.0 200.0 1.0
    
    max_step = 100
    stop_time = 3600.0
    cfl = 0.7
    plot_dt = 300.0
    
    init_type = sphere
    init_center = 100.0 100.0 0.5
    init_radius = 10.0
    
    use_farsite_model = 1
    farsite.enable = 1
    farsite.use_anderson_LW = 1
    farsite.phi_threshold = 0.1
    
    fuel_model = 3
    fuel_moisture = 0.08
    wind_speed = 15.0
    wind_direction = 270.0

Terrain and Crown Fire Simulation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    # Complex simulation with terrain, wind, and crown fire
    n_cell = 512 512 1
    max_grid_size = 64
    prob_lo = 0.0 0.0 0.0
    prob_hi = 1000.0 1000.0 1.0
    
    max_step = 500
    stop_time = 7200.0
    cfl = 0.5
    plot_dt = 600.0
    write_xy_data = 1
    
    init_type = sphere
    init_center = 500.0 500.0 0.5
    init_radius = 20.0
    
    use_terrain_effects = 1
    terrain_from_file = 1
    terrain_file = "elevation.xyz"
    
    wind_speed = 20.0
    wind_direction = 180.0
    
    fuel_model = 10
    fuel_moisture = 0.06
    
    crown.enable = 1
    crown.canopy_base_height = 4.0
    crown.canopy_bulk_density = 0.2
    crown.canopy_moisture = 0.10
    
    spotting.enable = 1
    spotting.probability = 0.05
    spotting.max_distance = 300.0

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
