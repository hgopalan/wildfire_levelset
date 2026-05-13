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

+----------------------------------------+---------------------------------------------------+
| Field name                             | Description                                       |
+========================================+===================================================+
| ``phi``                                | Level-set function (fire front where phi = 0)     |
+----------------------------------------+---------------------------------------------------+
| ``velx``, ``vely`` [, ``velz``]        | Wind velocity components [m/s]                    |
+----------------------------------------+---------------------------------------------------+
| ``farsite_dx``, ``farsite_dy``         | FARSITE spread displacements                      |
| [, ``farsite_dz``]                     |                                                   |
+----------------------------------------+---------------------------------------------------+
| ``R``                                  | Rate of spread [m/s]                              |
+----------------------------------------+---------------------------------------------------+
| ``spot_prob``, ``spot_count``          | Firebrand spotting probability and count          |
+----------------------------------------+---------------------------------------------------+
| ``spot_dist``, ``spot_active``         | Spotting distance [m] and ignition flag           |
+----------------------------------------+---------------------------------------------------+
| ``fuel_consumption``                   | Bulk fuel consumption fraction                    |
+----------------------------------------+---------------------------------------------------+
| ``crown_fraction``                     | Crown fire fraction                               |
+----------------------------------------+---------------------------------------------------+
| ``albini_Hz``                          | Albini lofting height [m]                         |
+----------------------------------------+---------------------------------------------------+
| ``albini_count``                       | Firebrands launched per cell                      |
+----------------------------------------+---------------------------------------------------+
| ``albini_dist``                        | Maximum firebrand landing distance [m]            |
+----------------------------------------+---------------------------------------------------+
| ``albini_active``                      | Albini spot ignition flag                         |
+----------------------------------------+---------------------------------------------------+
| ``elevation``                          | Terrain elevation [m]                             |
+----------------------------------------+---------------------------------------------------+
| ``slope``                              | Terrain slope [degrees]                           |
+----------------------------------------+---------------------------------------------------+
| ``aspect``                             | Terrain aspect [degrees]                          |
+----------------------------------------+---------------------------------------------------+
| ``fuel_model``                         | Fuel model number (from landscape)                |
+----------------------------------------+---------------------------------------------------+
| ``fireline_intensity``                 | Byram fireline intensity [kW/m]                   |
+----------------------------------------+---------------------------------------------------+
| ``flame_length``                       | Byram flame length [m]                            |
+----------------------------------------+---------------------------------------------------+
| ``weise_flame_height``                 | Weise & Biging flame height [m]                   |
+----------------------------------------+---------------------------------------------------+
| ``weise_flame_tilt``                   | Flame tilt angle [rad]                            |
+----------------------------------------+---------------------------------------------------+
| ``weise_whirl_height``                 | Fire whirl height [m]                             |
+----------------------------------------+---------------------------------------------------+
| ``weise_whirl_radius``                 | Fire whirl radius [m]                             |
+----------------------------------------+---------------------------------------------------+
| ``weise_angular_velocity``             | Fire whirl angular velocity [rad/s]               |
+----------------------------------------+---------------------------------------------------+
| ``weise_max_tang_vel``                 | Fire whirl max tangential velocity [m/s]          |
+----------------------------------------+---------------------------------------------------+
| ``viegas_ROS``                         | Viegas eruptive rate of spread [m/s]              |
+----------------------------------------+---------------------------------------------------+
| ``viegas_eruptive_flag``               | Viegas eruptive fire flag (0/1)                   |
+----------------------------------------+---------------------------------------------------+
| ``viegas_ROS_excess``                  | Viegas ROS excess above threshold                 |
+----------------------------------------+---------------------------------------------------+
| ``viegas_flame_tilt``                  | Viegas flame tilt angle [rad]                     |
+----------------------------------------+---------------------------------------------------+
| ``viegas_slope_factor``                | Viegas slope amplification factor                 |
+----------------------------------------+---------------------------------------------------+
| ``scorch_height``                      | Van Wagner scorch height [m]                      |
+----------------------------------------+---------------------------------------------------+
| ``prob_ignition``                      | Anderson probability of ignition [0–1]            |
+----------------------------------------+---------------------------------------------------+
| ``tree_mortality``                     | Ryan–Reinhardt tree mortality [0–1]               |
+----------------------------------------+---------------------------------------------------+
| ``crown_activity``                     | Crown activity class (0=surface, 1=passive,       |
|                                        | 2=active)                                         |
+----------------------------------------+---------------------------------------------------+
| ``torching_ratio``                     | Scott–Reinhardt torching ratio                    |
+----------------------------------------+---------------------------------------------------+
| ``crowning_ratio``                     | Scott–Reinhardt crowning ratio                    |
+----------------------------------------+---------------------------------------------------+
| ``energy_release_component``           | NFDRS energy release component                    |
+----------------------------------------+---------------------------------------------------+
| ``nfdrs_spread_component``             | NFDRS spread component                            |
+----------------------------------------+---------------------------------------------------+
| ``nfdrs_burning_index``                | NFDRS burning index                               |
+----------------------------------------+---------------------------------------------------+
| ``co2_emissions``                      | CO₂ emission rate [kg/m²/s]                       |
+----------------------------------------+---------------------------------------------------+
| ``co_emissions``                       | CO emission rate [kg/m²/s]                        |
+----------------------------------------+---------------------------------------------------+
| ``pm25_emissions``                     | PM₂.₅ emission rate [kg/m²/s]                     |
+----------------------------------------+---------------------------------------------------+
| ``arrival_time``                       | Fire arrival time [s]                             |
+----------------------------------------+---------------------------------------------------+
| ``heat_per_unit_area``                 | Heat per unit area [kJ/m²]                        |
+----------------------------------------+---------------------------------------------------+
| ``vorticity_z``                        | Vertical vorticity [1/s]                          |
+----------------------------------------+---------------------------------------------------+
| ``cbh``                                | Canopy base height [m]                            |
+----------------------------------------+---------------------------------------------------+
| ``cbd``                                | Canopy bulk density [kg/m³]                       |
+----------------------------------------+---------------------------------------------------+
| ``canopy_cover``                       | Canopy cover fraction [0–1]                       |
+----------------------------------------+---------------------------------------------------+
| ``canopy_height``                      | Canopy height [m]                                 |
+----------------------------------------+---------------------------------------------------+
| ``burnout_time``                       | Cell burnout time [s]                             |
+----------------------------------------+---------------------------------------------------+
| ``reaction_intensity``                 | Rothermel reaction intensity [BTU/ft²/min]        |
+----------------------------------------+---------------------------------------------------+
| ``residual_fuel``                      | Residual (unconsumed) fuel fraction               |
+----------------------------------------+---------------------------------------------------+
| ``shade_fraction``                     | Terrain/canopy shade fraction [0–1]               |
+----------------------------------------+---------------------------------------------------+
| ``torching_index_kmh``                 | Scott–Reinhardt torching index [km/h]             |
+----------------------------------------+---------------------------------------------------+
| ``crowning_index_kmh``                 | Scott–Reinhardt crowning index [km/h]             |
+----------------------------------------+---------------------------------------------------+
| ``moisture_d1``                        | 1-hr dead fuel moisture [fraction]                |
+----------------------------------------+---------------------------------------------------+
| ``moisture_d10``                       | 10-hr dead fuel moisture [fraction]               |
+----------------------------------------+---------------------------------------------------+
| ``moisture_d100``                      | 100-hr dead fuel moisture [fraction]              |
+----------------------------------------+---------------------------------------------------+
| ``moisture_lh``                        | Live herbaceous moisture [fraction]               |
+----------------------------------------+---------------------------------------------------+
| ``moisture_lw``                        | Live woody moisture [fraction]                    |
+----------------------------------------+---------------------------------------------------+

.. note::

   The plotfile can become large when all variables are enabled.  Use
   ``plot_vars`` to write only the fields you need (see :ref:`plot-vars`).

In addition, two plain-text files are written at each output step:

* ``phi_negative_NNNN.dat``: X, Y coordinates of all cells where ``phi < 0`` (inside the fire)
* ``phi_envelope_NNNN.dat``: X, Y coordinates of the convex hull of the fire perimeter

Input Parameters Reference
---------------------------

All Parameters Quick Reference
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following table lists every input parameter with its default value and
purpose.  Detailed descriptions and examples follow in the subsections below.

.. list-table::
   :header-rows: 1
   :widths: 35 18 47

   * - Parameter
     - Default
     - Purpose
   * - **Domain and Grid**
     -
     -
   * - ``n_cell``
     - 64
     - Number of cells in all directions (cube); overridden by per-direction values
   * - ``n_cell_x``, ``n_cell_y``, ``n_cell_z``
     - n_cell
     - Number of cells in each Cartesian direction independently
   * - ``max_grid_size``
     - 32
     - Maximum grid-box size for domain decomposition
   * - ``prob_lo_x``, ``prob_lo_y``, ``prob_lo_z``
     - 0.0
     - Lower corner of physical domain [m]
   * - ``prob_hi_x``, ``prob_hi_y``, ``prob_hi_z``
     - 1.0
     - Upper corner of physical domain [m]
   * - **Time Stepping**
     -
     -
   * - ``nsteps``
     - 300
     - Maximum number of time steps (used when ``final_time`` is not set)
   * - ``final_time``
     - -1.0
     - Simulation stop time [s]; overrides ``nsteps`` when positive
   * - ``cfl``
     - 0.5
     - CFL number for adaptive time-step calculation
   * - ``plot_int``
     - 50
     - Number of time steps between output (plot) files
   * - ``plot_vars``
     - "" (all)
     - Space-separated list of variable names to include in plotfiles; when empty all fields are written (see :ref:`plot-vars` for the complete list)
   * - ``reinit_int``
     - 20
     - Steps between level-set reinitialization (-1 = disabled)
   * - **Initial Conditions**
     -
     -
   * - ``source_type``
     - sphere
     - Initial fire shape: ``sphere``, ``box``, ``ellipse``, ``eb``, ``polygon``, ``polyline``
   * - ``center_x``, ``center_y``, ``center_z``
     - 0.5, 0.5, 0.5
     - Center of spherical ignition [m]
   * - ``sphere_radius``
     - (required)
     - Radius of spherical ignition [m]
   * - ``box_xmin`` … ``box_zmax``
     - (required)
     - Bounds of box ignition region [m]
   * - ``ellipse_center_x`` … ``ellipse_radius_z``
     - (required)
     - Center and semi-axes of elliptical ignition [m]
   * - ``fire_points_file``
     - ""
     - CSV of ignition point coordinates (X Y [Z]), one per row
   * - ``fire_gaussian_sigma``
     - auto
     - Ignition disk radius [m] (≤ 0 = 3 × min cell spacing)
   * - ``fire_polygon_file``
     - ""
     - CSV vertex file for ``polygon`` or ``polyline`` ignition
   * - ``polyline_width``
     - 10.0
     - Half-width of line-ignition zone [m]
   * - ``fire_polygon_z_level``
     - 0.5
     - Z-centre of ignition layer in 3-D polygon/polyline runs [m]
   * - **Fire Model Selection**
     -
     -
   * - ``fire_spread_model``
     - rothermel
     - ROS model: ``rothermel``, ``balbi``, ``cheney_gould``, ``cruz_crown``, ``fbp_o1a``, ``fbp_o1b``, ``fbp_s1``, ``fbp_s2``, ``fbp_s3``, ``lautenberger``
   * - ``propagation_method``
     - levelset
     - Propagation method: ``levelset``, ``farsite``, ``mtt``
   * - **Barriers**
     -
     -
   * - ``barrier_files``
     - ""
     - Space-separated list of barrier polyline/polygon CSV files
   * - **Rothermel Fuel Properties**
     -
     -
   * - ``rothermel.fuel_model``
     - custom
     - Fuel model name from database (e.g. ``FM4``, ``FM101``)
   * - ``rothermel.M_f``
     - 0.08
     - Aggregate dead fuel moisture [fraction]
   * - ``rothermel.M_x``
     - 0.30
     - Moisture of extinction [fraction]
   * - ``rothermel.w0``
     - 0.230
     - Oven-dry total fuel loading [lb/ft²]
   * - ``rothermel.sigma``
     - 1739.0
     - Surface-area-to-volume ratio [ft⁻¹]
   * - ``rothermel.delta``
     - 2.0
     - Fuel bed depth [ft]
   * - ``rothermel.h_heat``
     - 8000.0
     - Heat content [BTU/lb]
   * - ``rothermel.S_T``
     - 0.0555
     - Total mineral content [fraction]
   * - ``rothermel.S_e``
     - 0.010
     - Effective mineral content [fraction]
   * - ``rothermel.rho_p``
     - 32.0
     - Oven-dry particle density [lb/ft³]
   * - ``rothermel.slope_x``, ``rothermel.slope_y``
     - 0.0
     - Constant terrain slope (tan θ) in x and y directions
   * - ``rothermel.M_d1``
     - =M_f
     - 1-hr dead fuel moisture [fraction]
   * - ``rothermel.M_d10``
     - =M_f
     - 10-hr dead fuel moisture [fraction]
   * - ``rothermel.M_d100``
     - =M_f
     - 100-hr dead fuel moisture [fraction]
   * - ``rothermel.M_lh``
     - 0.90
     - Live herbaceous fuel moisture [fraction]
   * - ``rothermel.M_lw``
     - 1.20
     - Live woody fuel moisture [fraction]
   * - ``rothermel.w_d1``
     - 0.0
     - 1-hr dead fuel load [lb/ft²] (override)
   * - ``rothermel.sigma_d1``
     - 0.0
     - 1-hr dead SAV [ft⁻¹] (override)
   * - ``rothermel.w_d10``
     - 0.0
     - 10-hr dead fuel load [lb/ft²] (override)
   * - ``rothermel.w_d100``
     - 0.0
     - 100-hr dead fuel load [lb/ft²] (override)
   * - ``rothermel.w_lh``
     - 0.0
     - Live herbaceous fuel load [lb/ft²] (override)
   * - ``rothermel.sigma_lh``
     - 0.0
     - Live herbaceous SAV [ft⁻¹] (override)
   * - ``rothermel.w_lw``
     - 0.0
     - Live woody fuel load [lb/ft²] (override)
   * - ``rothermel.sigma_lw``
     - 0.0
     - Live woody SAV [ft⁻¹] (override)
   * - ``rothermel.use_waf``
     - 0
     - Enable Wind Adjustment Factor (1=yes): converts 20-ft→midflame wind
   * - ``rothermel.waf_formula``
     - ``"andrews"``
     - WAF formula: ``"andrews"`` (logarithmic Albini & Baughman 1979) or
       ``"behaviorplus"`` (linear 0.36 + 0.004×h_in for open/shrub fuels;
       exponential canopy attenuation for forest fuels)
   * - ``rothermel.waf_canopy_alpha``
     - 1.5
     - Canopy attenuation coefficient α_c for the BehavePlus exponential canopy
       WAF: WAF_canopy = WAF_open × exp(−α_c × f_c).  Larger values give
       stronger wind sheltering.  Only used when ``waf_formula = "behaviorplus"``
       and canopy cover data are available.
   * - **Cheney & Gould Grassland Parameters**
     -
     -
   * - ``cheney_gould.moisture``
     - 10.0
     - Dead fine fuel moisture [%]
   * - ``cheney_gould.curing``
     - 1.0
     - Degree of grass curing [0–1; 1 = fully cured]
   * - **Canadian FBP System Parameters**
     -
     -
   * - ``fbp.fuel_type``
     - o1a
     - FBP fuel type (overridden by fire_spread_model): ``o1a``, ``o1b``, ``s1``, ``s2``, ``s3``
   * - ``fbp.moisture``
     - 10.0
     - Dead fine fuel moisture [%]
   * - ``fbp.curing``
     - 80.0
     - Degree of grass curing [%] (O1a/O1b only)
   * - **Lautenberger (2013) Parameters**
     -
     -
   * - ``lautenberger.A_L``
     - 1.05e-5
     - Pre-exponential ROS coefficient [m²/s]
   * - ``lautenberger.B_L``
     - 2.5
     - Wind speed exponent (dimensionless)
   * - ``lautenberger.C_L``
     - 0.40
     - Moisture sensitivity coefficient [(m/s)⁻¹]
   * - ``lautenberger.D_L``
     - 0.50
     - Slope sensitivity coefficient (dimensionless)
   * - **Wind**
     -
     -
   * - ``u_x``, ``u_y``, ``u_z``
     - 0.25, 0.0, 0.0
     - Constant wind velocity components [m/s]
   * - ``velocity_file``
     - ""
     - Spatially-varying wind field CSV (X Y U V)
   * - ``use_time_dependent_wind``
     - 0
     - Enable time-dependent sequence of wind field files (1=yes)
   * - ``wind_time_spacing``
     - 60.0
     - Time spacing between consecutive wind field files [s]
   * - ``wind_dir_schedule_file``
     - ""
     - Three-column CSV wind schedule (time_s, speed_ms, dir_deg)
   * - ``multi_wtr_file``
     - ""
     - Station list CSV for multi-station IDW weather interpolation
   * - ``multi_wtr_idw_power``
     - 2.0
     - IDW exponent for multi-station weather interpolation
   * - **Terrain**
     -
     -
   * - ``rothermel.terrain_file``
     - ""
     - Terrain elevation CSV (X Y Z)
   * - ``rothermel.landscape_file``
     - ""
     - FARSITE landscape file (ASCII: X Y ELEV SLOPE ASPECT FUEL)
   * - ``rothermel.landscape_fuel_type``
     - "13"
     - Fuel model system in landscape file: ``"13"`` (FBFM13) or ``"40"`` (FBFM40)
   * - **FARSITE Parameters**
     -
     -
   * - ``farsite.use_anderson_LW``
     - 0
     - Use Anderson (1983) L/W ratio from wind speed (1=yes)
   * - ``farsite.length_to_width_ratio``
     - 3.0
     - Fixed ellipse length-to-width ratio
   * - ``farsite.coeff_a``
     - 1.0
     - Richards' head-fire (maximum spread) coefficient
   * - ``farsite.coeff_b``
     - 0.5
     - Richards' flank-fire coefficient
   * - ``farsite.coeff_c``
     - 0.2
     - Richards' backing-fire (minimum spread) coefficient
   * - ``farsite.phi_threshold``
     - 0.1
     - Level-set threshold for fire-front cell detection
   * - ``farsite.scale_ellipse_with_crown``
     - 0
     - Scale ellipse L/W during active crown fire (1=yes)
   * - ``farsite.crown_lw_scale``
     - 1.5
     - L/W multiplier applied during active crown fire
   * - **Crown Fire Parameters**
     -
     -
   * - ``crown.enable``
     - 0
     - Enable Van Wagner (1977) crown fire model (1=yes)
   * - ``crown.CBH``
     - 4.0
     - Canopy base height [m]
   * - ``crown.CBD``
     - 0.15
     - Canopy bulk density [kg/m³]
   * - ``crown.FMC``
     - 100.0
     - Foliar moisture content [%]
   * - ``crown.crown_fraction_weight``
     - 1.0
     - Crown fire weighting factor [0–2]
   * - ``crown.use_metric_units``
     - 1
     - Unit system: 1 = metric (m, kW/m), 0 = imperial
   * - ``crown.use_cruz_crown``
     - 0
     - Use Cruz (2005) empirical crown ROS formula instead of Van Wagner proxy (1=yes)
   * - ``cruz_crown.CBD``
     - 0.10
     - Canopy bulk density for Cruz (2005) formula [kg/m³]
   * - ``cruz_crown.MC10``
     - 10.0
     - 10-hr dead fuel moisture for Cruz (2005) formula [%]
   * - **Bulk Fuel Consumption**
     -
     -
   * - ``farsite.use_bulk_fuel_consumption``
     - 0
     - Enable exponential bulk fuel consumption model (1=yes)
   * - ``farsite.tau_residence``
     - 60.0
     - Fuel residence / burnout time [s]
   * - ``farsite.f_consumed_max``
     - 0.9
     - Maximum fuel consumption fraction [0–1]
   * - ``farsite.f_consumed_min``
     - 0.5
     - Minimum fuel consumption fraction [0–1]
   * - **Probability-Based Spotting**
     -
     -
   * - ``spotting.enable``
     - 0
     - Enable stochastic firebrand spotting model (1=yes)
   * - ``spotting.P_base``
     - 0.02
     - Base spotting probability per fire-front cell [0–1]
   * - ``spotting.k_wind``
     - 0.3
     - Wind speed scaling coefficient for spotting probability
   * - ``spotting.I_critical``
     - 1000.0
     - Critical fire intensity threshold for spotting [kW/m]
   * - ``spotting.d_mean``
     - 0.1
     - Mean spotting distance parameter
   * - ``spotting.d_sigma``
     - 0.5
     - Lognormal standard deviation for spotting distance
   * - ``spotting.d_lambda``
     - 10.0
     - Exponential decay rate for spotting distance
   * - ``spotting.distance_model``
     - "lognormal"
     - Spotting distance distribution: ``lognormal`` or ``exponential``
   * - ``spotting.lateral_spread_angle``
     - 15.0
     - Angular spread perpendicular to wind direction [°]
   * - ``spotting.spot_radius``
     - 0.02
     - Spot-fire ignition zone radius [m]
   * - ``spotting.random_seed``
     - 0
     - RNG seed (0 = system clock)
   * - ``spotting.check_interval``
     - 5
     - Run spotting model every N time steps
   * - **Albini (1983) Physics-Based Spotting**
     -
     -
   * - ``albini_spotting.enable``
     - 0
     - Enable Albini physics-based spotting model (1=yes)
   * - ``albini_spotting.terminal_velocity``
     - 1.0
     - Firebrand terminal descent velocity [m/s]
   * - ``albini_spotting.P_base``
     - 0.01
     - Maximum launch probability per front cell [0–1]
   * - ``albini_spotting.I_B_min``
     - 10.0
     - Minimum Byram intensity for firebrand launch [kW/m]
   * - ``albini_spotting.spot_radius``
     - 5.0
     - Spot-fire ignition zone radius [m]
   * - ``albini_spotting.random_seed``
     - 0
     - RNG seed (0 = system clock)
   * - ``albini_spotting.check_interval``
     - 5
     - Run Albini spotting every N time steps
   * - ``albini_spotting.n_traj_steps``
     - 100
     - Forward-Euler sub-steps for 2-D firebrand trajectory
   * - ``albini_spotting.use_3d_wind``
     - 0
     - Use 3-D wind from massconsistent_amr plt file for trajectory (1=yes)
   * - ``albini_spotting.plt_wind_file``
     - ""
     - Path to massconsistent_amr plt directory (required when use_3d_wind=1)
   * - **Torching-Tree Spotting**
     -
     -
   * - ``torching_spotting.enable``
     - 0
     - Enable torching-tree spotting model (1=yes)
   * - ``torching_spotting.k_torch``
     - 4.24
     - Spotting distance coefficient: d = k × U × √H_eff [m]
   * - ``torching_spotting.I_B_min``
     - 100.0
     - Minimum Byram intensity for torching-tree spotting [kW/m]
   * - ``torching_spotting.spot_radius``
     - 5.0
     - Spot-fire ignition zone radius [m]
   * - ``torching_spotting.P_base``
     - 0.05
     - Probability of spotting per torching cell per check [0–1]
   * - ``torching_spotting.check_interval``
     - 5
     - Run torching-tree spotting every N time steps
   * - ``torching_spotting.min_crown_activity``
     - 1
     - Min crown activity to trigger spotting (1=passive+active, 2=active-only)
   * - ``torching_spotting.random_seed``
     - 0
     - RNG seed (0 = system clock)
   * - **Turbulent Wind Perturbation**
     -
     -
   * - ``turb_wind.model``
     - "none"
     - Turbulent wind model: ``none``, ``ou_process``, ``spectral_noise``, ``direction_walk``
   * - ``turb_wind.theta``
     - 0.1
     - Ornstein-Uhlenbeck reversion rate [s⁻¹]
   * - ``turb_wind.sigma``
     - 0.5
     - OU stationary standard deviation [m/s]
   * - ``turb_wind.L_c``
     - 0.0
     - Spatial correlation length [m] (0 = domain-uniform perturbation)
   * - ``turb_wind.N_modes``
     - 32
     - Number of random Fourier modes (``spectral_noise`` model)
   * - ``turb_wind.sigma_theta``
     - 0.1
     - Direction walk angular standard deviation [rad/step]
   * - ``turb_wind.theta_max``
     - 0.5236
     - Max cumulative direction deviation [rad] (≈ ±30°)
   * - ``turb_wind.random_seed``
     - 0
     - RNG seed (0 = system clock)
   * - **Wind-Terrain Feedback**
     -
     -
   * - ``wind_terrain.model``
     - "none"
     - Terrain wind model: ``none``, ``viegas_ros``, ``viegas_wind``, ``canyon_wind``, ``viegas_neto``, ``pimont``, ``windninja_ridge_canyon``
   * - ``wind_terrain.k_canyon``
     - 1.0
     - Canyon wind amplification coefficient (``canyon_wind`` model)
   * - ``wind_terrain.k_pimont``
     - 0.5
     - Exponential slope correction coefficient (``pimont`` model)
   * - ``wind_terrain.k_ridge``
     - 1.0
     - Ridge speed-up coefficient (``windninja_ridge_canyon`` model)
   * - ``wind_terrain.k_canyon_wn``
     - 0.5
     - Canyon channeling coefficient (``windninja_ridge_canyon`` model)
   * - **Heat Flux**
     -
     -
   * - ``heat_flux.value``
     - 0.0
     - Uniform heat flux [W/m²]
   * - ``heat_flux.file``
     - ""
     - Spatially-varying heat flux CSV (X Y Q); overrides ``heat_flux.value``
   * - ``heat_flux.rho_air``
     - 1.2
     - Air density [kg/m³]
   * - ``heat_flux.Cp_air``
     - 1005.0
     - Specific heat of air [J/(kg·K)]
   * - ``heat_flux.T_a``
     - 300.0
     - Ambient temperature for buoyancy calculation [K]
   * - ``heat_flux.ref_height``
     - 10.0
     - Reference plume height for buoyancy calculation [m]
   * - ``heat_flux.k_upward``
     - 1.0
     - Scaling coefficient for fire-plume upward velocity
   * - ``heat_flux.k_induced``
     - 0.5
     - Scaling coefficient for induced horizontal inflow
   * - ``heat_flux.enable_upward``
     - 0
     - Enable upward convective velocity term (1=yes)
   * - ``heat_flux.enable_induced``
     - 0
     - Enable induced horizontal inflow term (1=yes)
   * - **Fire Ecology Coupling**
     -
     -
   * - ``fire_ecology.couple_to_ros``
     - 0
     - Scale ROS by per-cell ignition probability in unburned cells (1=yes)
   * - ``fire_ecology.p_ignition_floor``
     - 0.05
     - Minimum ignition probability used as ROS scale floor
   * - ``fire_ecology.T_a_C``
     - 25.0
     - Ambient temperature for fire ecology calculations [°C]
   * - ``fire_ecology.solar_heating_F``
     - 25.0
     - Solar temperature increment [°F]
   * - ``fire_ecology.tree_height``
     - 10.0
     - Mean stand height for tree mortality sub-model [m]
   * - **Persistent Fuel Depletion**
     -
     -
   * - ``fuel_depletion.enable``
     - 0
     - Track per-cell residual fuel fraction (1=yes)
   * - ``fuel_depletion.tau_burnout``
     - 3600.0
     - Fine-fuel burnout time constant [s]
   * - ``fuel_depletion.couple_to_ros``
     - 0
     - Scale ROS by residual fuel in re-entry cells (1=yes)
   * - **Fire Acceleration**
     -
     -
   * - ``acceleration.enable``
     - 0
     - Enable small-fire acceleration model (1=yes)
   * - ``acceleration.L_acc``
     - 50.0
     - Acceleration length scale [m]
   * - **FMC Seasonal Schedule**
     -
     -
   * - ``fmc_schedule.enable``
     - 0
     - Enable FMC seasonal schedule (1=yes)
   * - ``fmc_schedule.file``
     - ""
     - Two-column CSV (day-of-year, fmc_pct)
   * - ``fmc_schedule.use_farsite_curve``
     - 0
     - Use built-in parametric phenology curve (1=yes)
   * - ``fmc_schedule.start_doy``
     - 1
     - Day-of-year at simulation t = 0
   * - ``fmc_schedule.spring_start``
     - 90
     - DOY when spring green-up begins
   * - ``fmc_schedule.summer_peak``
     - 150
     - DOY when FMC reaches its maximum
   * - ``fmc_schedule.fall_start``
     - 240
     - DOY when autumn curing begins
   * - ``fmc_schedule.fall_end``
     - 300
     - DOY when curing ends
   * - ``fmc_schedule.fmc_min``
     - 85.0
     - Dormant / cured FMC [%]
   * - ``fmc_schedule.fmc_max``
     - 140.0
     - Peak green FMC [%]
   * - **Precipitation Wetting**
     -
     -
   * - ``diurnal_moisture.precip_rain_rate_mm_hr``
     - 0.0
     - Constant rain rate [mm/hr]
   * - ``diurnal_moisture.precip_schedule_file``
     - ""
     - CSV of (time_s, rain_mm_hr) for time-varying rain
   * - ``diurnal_moisture.precip_threshold_mm_hr``
     - 0.25
     - Minimum rain rate to trigger dead-fuel wetting [mm/hr]
   * - ``diurnal_moisture.M_sat``
     - 1.20
     - Saturation moisture content [fraction]
   * - **FMS Moisture Scenario File**
     -
     -
   * - ``fms_file``
     - ""
     - Path to FARSITE ``.fms`` per-fuel-model moisture scenario file

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
  - ``cruz_crown`` – Cruz, Alexander & Wakimoto (2005) algebraic crown fire spread model
  - ``fbp_o1a`` – Canadian Forest Fire Behaviour Prediction (FBP) System: O1a matted grass
  - ``fbp_o1b`` – Canadian FBP System: O1b standing grass
  - ``fbp_s1`` – Canadian FBP System: S1 Jack or Lodgepole Pine slash
  - ``fbp_s2`` – Canadian FBP System: S2 White Spruce-Balsam slash
  - ``fbp_s3`` – Canadian FBP System: S3 Coastal Cedar-Hemlock-Douglas-Fir slash
  - ``lautenberger`` – Lautenberger (2013) semi-empirical firebrand-driven model

  Example: ``fire_spread_model = cheney_gould``

**propagation_method** (default: ``levelset``)
  Select how the fire perimeter is propagated. Currently supported values:

  - ``levelset`` – WENO5-Z level-set advection (default); Rothermel/Balbi ROS drives
    the signed-distance evolution equation.
  - ``farsite`` – FARSITE elliptical Huygens-wavelet expansion (Richards 1990);
    uses the FARSITE spread parameters (``farsite.*``) and pre-computed ROS.
  - ``mtt`` – Minimum Travel Time (Finney 2002); arrival times are pre-computed
    once via a Dijkstra fast-marching sweep using the selected ROS model, then
    ``phi = arrival_time − current_time`` is applied each step.  No
    reinitialization is needed.  Compatible with all fire spread models.

  Example: ``propagation_method = mtt``

Barrier Polygons / Firebreaks
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**barrier_files** (default: empty)
  Space-separated list of CSV files, each defining a barrier polyline or polygon
  as a sequence of ``X Y [Z]`` vertex coordinates (one per line; lines starting
  with ``#`` are ignored).  At every time step the nearest grid cells to each
  barrier vertex are identified and any that are currently burning (``phi < 0``)
  are extinguished by setting ``phi = +0.5 * min(dx,dy)``.

  This models firebreaks, fuel breaks, roads, or other physical barriers without
  using the AMReX Embedded Boundary (EB) framework.

  Example::

      barrier_files = road_break.csv ridge_break.csv

  Each CSV file format::

      # X Y coordinates of barrier vertices (UTM metres)
      330100.0  3775300.0
      330200.0  3775300.0
      330300.0  3775300.0

**burn_period.enable** (default: 0)
  When set to 1, fire spread (rate-of-spread field) is suppressed to zero
  outside a user-specified daily active-burn time window.  This reproduces the
  operational *burn period* concept used in FARSITE and FSPro, where nighttime
  recovery of fuel moisture effectively halts spread.  Moisture evolution and
  all diagnostics continue normally during inactive hours.

  Example::

      burn_period.enable     = 1
      burn_period.start_hour = 10.0   # 10:00 AM – spread begins
      burn_period.end_hour   = 20.0   # 8:00 PM  – spread ends
      burn_period.sim_start_hour = 8.0

**burn_period.start_hour** (default: 10.0)
  Local clock hour (decimal, 0–24) when the daily burn period begins.

**burn_period.end_hour** (default: 20.0)
  Local clock hour (decimal, 0–24) when the daily burn period ends.
  Overnight windows (``end_hour < start_hour``) are handled correctly.

**burn_period.sim_start_hour** (default: inherited from ``solar_radiation.sim_start_hour``)
  Local clock hour at simulation time = 0.  Defaults to
  ``solar_radiation.sim_start_hour`` when solar shading is enabled; otherwise
  defaults to 0.0 (midnight start).

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

Canadian FBP System Spread Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

These parameters are used when ``fire_spread_model`` is one of ``fbp_o1a``,
``fbp_o1b``, ``fbp_s1``, ``fbp_s2``, or ``fbp_s3``.  The Canadian Forest
Fire Behaviour Prediction (FBP) System uses empirical equations derived from
experimental and prescribed-fire data for specific fuel types.

**fbp.fuel_type** (default: "o1a")
  FBP fuel type string.  Automatically overridden by the ``fire_spread_model``
  setting (e.g., ``fire_spread_model = fbp_s2`` sets ``fbp.fuel_type = s2``).
  Valid values: ``o1a``, ``o1b``, ``s1``, ``s2``, ``s3``.

  Example: ``fbp.fuel_type = o1b``

**fbp.moisture** (default: 10.0)
  Dead fine fuel moisture content [%].  Used in the exponential moisture
  dampening factor for grass and slash fuel types.

  Example: ``fbp.moisture = 8.0``

**fbp.curing** (default: 80.0)
  Degree of grass curing [%] (0–100).  Only used for grass fuel types O1a/O1b.
  Higher curing increases the effective ROS.

  Example: ``fbp.curing = 90.0``

Lautenberger (2013) Semi-Empirical Spread Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

These parameters are used when ``fire_spread_model = lautenberger``.  The
model is a semi-empirical rate-of-spread formula that accounts for wind speed,
fuel moisture, and slope through four calibration coefficients.

**lautenberger.A_L** (default: 1.05e-5)
  Pre-exponential rate-of-spread coefficient [m²/s].

  Example: ``lautenberger.A_L = 1.05e-5``

**lautenberger.B_L** (default: 2.5)
  Wind speed exponent (dimensionless).

  Example: ``lautenberger.B_L = 2.5``

**lautenberger.C_L** (default: 0.40)
  Moisture sensitivity coefficient [(m/s)⁻¹].

  Example: ``lautenberger.C_L = 0.40``

**lautenberger.D_L** (default: 0.50)
  Slope sensitivity coefficient (dimensionless).

  Example: ``lautenberger.D_L = 0.50``

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

**farsite.scale_ellipse_with_crown** (default: 0)
  Scale the FARSITE ellipse length-to-width ratio when active crown fire is
  detected (1=yes, 0=no). Off by default for backward compatibility.

  Example: ``farsite.scale_ellipse_with_crown = 1``

**farsite.crown_lw_scale** (default: 1.5)
  Multiplier applied to the ellipse L/W ratio during active crown fire.
  Only used when ``farsite.scale_ellipse_with_crown = 1``.

  Example: ``farsite.crown_lw_scale = 1.5``

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

**crown.use_cruz_crown** (default: 0)
  Use the Cruz, Alexander & Wakimoto (2005) empirical crown ROS formula instead
  of the Van Wagner 3/CBD proxy (1=yes, 0=no).  When enabled, the crown ROS
  is ``R = 11.02 × U10^0.90 × CBD^0.19 × exp(-0.17 × MC10)``.

  Example: ``crown.use_cruz_crown = 1``

**cruz_crown.CBD** (default: 0.10)
  Canopy bulk density used specifically in the Cruz (2005) crown ROS formula
  [kg/m³].  May differ from ``crown.CBD`` used by the Van Wagner threshold.

  Example: ``cruz_crown.CBD = 0.12``

**cruz_crown.MC10** (default: 10.0)
  10-hr dead fuel moisture content used in the Cruz (2005) crown ROS formula
  [%].

  Example: ``cruz_crown.MC10 = 8.0``

Crown Fire State (Passive vs. Active)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Crown fire state is automatically distinguished when ``crown.enable = 1``:

* **Passive crown** (``I_B ≥ I_o`` but ``R < R'_SA``): fire torches
  intermittently but does not propagate through the canopy. Surface ROS is
  used; no crown-ROS enhancement is applied.
* **Active crown** (``I_B ≥ I_o`` AND ``R ≥ R'_SA``): fire propagates
  continuously through the canopy. Crown ROS (Van Wagner or Cruz) is applied.

The critical active-crown ROS is :math:`R'_{SA} = 3.0 / \text{CBD}` [m/min].
No additional inputs are required.

Fire Ecology Coupling Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

These parameters couple the fire ecology model (ignition probability,
tree mortality) to the spread-rate field.

**fire_ecology.couple_to_ros** (default: 0)
  Scale the Rothermel/Balbi ROS by the per-cell ignition probability in
  unburned cells (1=yes, 0=no).

  Example: ``fire_ecology.couple_to_ros = 1``

**fire_ecology.p_ignition_floor** (default: 0.05)
  Minimum ignition probability used as the ROS scale floor.  Prevents ROS
  from being driven to zero in low-probability cells.

  Example: ``fire_ecology.p_ignition_floor = 0.10``

**fire_ecology.T_a_C** (default: 25.0)
  Ambient temperature for fire ecology calculations [°C].

  Example: ``fire_ecology.T_a_C = 30.0``

**fire_ecology.solar_heating_F** (default: 25.0)
  Solar temperature increment added to ambient temperature [°F].  Used in
  fuel moisture and ignition probability computations.

  Example: ``fire_ecology.solar_heating_F = 30.0``

**fire_ecology.tree_height** (default: 10.0)
  Mean stand height [m] used in the tree mortality sub-model.

  Example: ``fire_ecology.tree_height = 15.0``

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

**albini_spotting.use_3d_wind** (default: 0)
  When set to 1, read horizontal wind from a 3-D AMReX plotfile produced by
  `massconsistent_amr <https://github.com/hgopalan/massconsistent_amr>`_ and
  use it for firebrand trajectory integration in place of the 2-D ``vel`` MultiFab.
  The plotfile must contain variables named ``"u"`` and ``"v"`` (and optionally
  ``"w"``).  The 3-D wind is projected to 2-D by column-averaging ``u`` and ``v``
  over all vertical levels at each horizontal grid point.  This allows full
  terrain-following mass-consistent wind fields to drive firebrand trajectories
  while the fire-spread solver remains 2-D.

  Requires ``albini_spotting.plt_wind_file`` to be set.

  Example: ``albini_spotting.use_3d_wind = 1``

**albini_spotting.plt_wind_file** (default: "")
  Path to the AMReX plotfile directory produced by massconsistent_amr.
  The directory must contain a ``Header`` file and a ``Level_0/`` sub-directory
  with ``Cell_H`` and binary data files.  Required when ``use_3d_wind = 1``.

  Example: ``albini_spotting.plt_wind_file = /path/to/plt_wind``

Albini (1979) Torching-Tree Spotting Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Long-range spotting originating from torching trees when crown fire is
present.  Triggered in cells where ``crown_activity ≥ min_crown_activity``.
The spotting distance is ``d = k_torch × U × sqrt(H_eff)`` [m].

**torching_spotting.enable** (default: 0)
  Enable the torching-tree spotting model (1=yes, 0=no).

  Example: ``torching_spotting.enable = 1``

**torching_spotting.k_torch** (default: 4.24)
  Spotting distance coefficient in the formula
  :math:`d = k_{\text{torch}} \times U \times \sqrt{H_{\text{eff}}}` [m].

  Example: ``torching_spotting.k_torch = 4.24``

**torching_spotting.I_B_min** (default: 100.0)
  Minimum Byram fire line intensity required for torching-tree spotting [kW/m].

  Example: ``torching_spotting.I_B_min = 200.0``

**torching_spotting.spot_radius** (default: 5.0)
  Radius of new spot-fire ignition zone [m].

  Example: ``torching_spotting.spot_radius = 10.0``

**torching_spotting.P_base** (default: 0.05)
  Probability of spotting per torching cell per check interval [0–1].

  Example: ``torching_spotting.P_base = 0.10``

**torching_spotting.check_interval** (default: 5)
  Run the torching-tree spotting model every N time steps.

  Example: ``torching_spotting.check_interval = 3``

**torching_spotting.min_crown_activity** (default: 1)
  Minimum crown fire activity level required to trigger spotting.
  1 = passive or active crown fire; 2 = active crown fire only.

  Example: ``torching_spotting.min_crown_activity = 2``

**torching_spotting.random_seed** (default: 0)
  Seed for the random number generator. 0 uses the system clock.

  Example: ``torching_spotting.random_seed = 42``

Persistent Fuel Depletion Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Tracks per-cell residual fuel fraction that decays exponentially after fire
passage: :math:`f_r = \exp(-(t - t_{\text{arrive}}) / \tau_{\text{burnout}})`.
Written to plotfiles as ``residual_fuel``.

**fuel_depletion.enable** (default: 0)
  Track per-cell residual fuel fraction (1=yes, 0=no).

  Example: ``fuel_depletion.enable = 1``

**fuel_depletion.tau_burnout** (default: 3600.0)
  Fine-fuel burnout time constant [s].

  Example: ``fuel_depletion.tau_burnout = 1800.0``

**fuel_depletion.couple_to_ros** (default: 0)
  Scale ROS by residual fuel fraction in cells the fire re-enters (1=yes, 0=no).

  Example: ``fuel_depletion.couple_to_ros = 1``

Fire Acceleration Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Small fires spread slower than their quasi-steady-state ROS.  This model
applies a scaling factor :math:`\alpha = 1 - \exp(-r_{\text{fire}} / L_{\text{acc}})`,
where :math:`r_{\text{fire}}` is the effective fire radius from burned area.
When :math:`r_{\text{fire}} \gg L_{\text{acc}}`, :math:`\alpha \to 1` and large
fires are unaffected.

**acceleration.enable** (default: 0)
  Enable the small-fire acceleration model (1=yes, 0=no). Off by default for
  backward compatibility.

  Example: ``acceleration.enable = 1``

**acceleration.L_acc** (default: 50.0)
  Acceleration length scale [m]. Controls how quickly the fire approaches its
  quasi-steady-state ROS as it grows.

  Example: ``acceleration.L_acc = 30.0``

Output Parameters
^^^^^^^^^^^^^^^^^

**plot_int** (default: 50)
  Number of time steps between output files.

  Example: ``plot_int = 10``

Turbulent Wind Perturbation Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

These parameters control the stochastic wind perturbation model that adds
temporally-correlated and (optionally) spatially-correlated fluctuations to
the ambient wind field at every timestep.  The unperturbed base wind is
always preserved internally; the perturbation is computed on top of it.

**turb_wind.model** (default: ``"none"``)
  Selects the turbulent wind model:

  - ``none`` – no perturbation (default)
  - ``ou_process`` – Ornstein-Uhlenbeck temporally-correlated noise.
    When ``turb_wind.L_c = 0`` all cells receive the same domain-uniform
    perturbation.  When ``turb_wind.L_c > 0`` per-cell spatially correlated
    perturbations are generated via a Gaussian spatial kernel of length scale
    ``L_c`` [m].
  - ``spectral_noise`` – Random Fourier Feature (RFF) spectral noise with
    OU temporal evolution.  ``N_modes`` wavenumber pairs are drawn at
    initialisation from the 2-D isotropic Gaussian power spectrum (length
    scale ``L_c``); scalar OU amplitude coefficients per mode evolve each
    step on the CPU; the perturbation field is reconstructed on the GPU as a
    cosine superposition.  Produces physically correct energy distribution
    across wavenumbers.  Requires ``L_c > 0``.
  - ``direction_walk`` – bounded cumulative random walk of wind direction.
    Wind speed is preserved exactly; only direction fluctuates.

  Example: ``turb_wind.model = spectral_noise``

**turb_wind.theta** (default: 0.1)
  Ornstein-Uhlenbeck reversion rate :math:`\theta` [s⁻¹].  The temporal
  decorrelation time (e-folding time of gust autocorrelation) is
  :math:`\tau = 1/\theta`.  Only used by ``ou_process``.  Must be > 0.

  Example: ``turb_wind.theta = 0.05``  (→ 20 s gust decorrelation)

**turb_wind.sigma** (default: 0.5)
  Stationary standard deviation of the OU perturbation [m/s].  Each cell's
  long-run perturbation has standard deviation exactly ``sigma``.  Only used
  by ``ou_process``.  Must be > 0.

  Example: ``turb_wind.sigma = 1.0``

**turb_wind.L_c** (default: 0.0)
  Spatial correlation length [m].  When > 0, activates the Gaussian kernel
  smoothing of the OU noise field so that cells within distance ``L_c``
  receive correlated perturbations.  Setting ``L_c = 0`` gives domain-uniform
  perturbations (all cells receive the same gust).  Only used by
  ``ou_process``.  Must be ≥ 0.

  The kernel standard deviation in cells is ``sigma_k = L_c / dx``.  For
  ``sigma_k ≫ 1`` the perturbation field is nearly uniform over the domain;
  for ``sigma_k ≈ 1`` each cell is nearly independent.

  Example: ``turb_wind.L_c = 100.0``  (→ sigma_k = 10 cells for dx = 10 m)

**turb_wind.N_modes** (default: 32)
  Number of random Fourier modes for the ``spectral_noise`` model.  More
  modes produce a smoother, more isotropic perturbation field at the cost of
  proportionally more GPU computation per cell per timestep.  Typical values
  are 16–128; 32 is a good default for most wildfire grids.  Only used by
  ``spectral_noise``.  Must be ≥ 1.

  Example: ``turb_wind.N_modes = 64``

**turb_wind.sigma_theta** (default: 0.1)
  Angular standard deviation [rad/step] for the ``direction_walk`` model.
  Controls the magnitude of each step's direction increment.  Must be > 0.

  Example: ``turb_wind.sigma_theta = 0.15``

**turb_wind.theta_max** (default: 0.5236)
  Maximum cumulative directional deviation [rad] for the ``direction_walk``
  model.  The accumulated angle is clamped to :math:`\pm\theta_{\max}`.
  The default is :math:`\pi/6 \approx 30°`.  Must be > 0.

  Example: ``turb_wind.theta_max = 0.785``  (→ ±45°)

**turb_wind.random_seed** (default: 0)
  Seed for the pseudo-random number generator.  ``0`` uses the system clock
  (non-reproducible between runs).  Any positive integer gives a
  reproducible sequence.

  Example: ``turb_wind.random_seed = 42``

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

.. _plot-vars:

Controlling Plotfile Output Variables
--------------------------------------

By default every simulation field is written to each AMReX plotfile.  For
large domains or long runs this can produce very large files.  The
``plot_vars`` input key lets you restrict the output to only the variables
you need.

Usage
^^^^^

Add a space-separated list of variable names to your ``inputs.i`` file::

    # Write only essential fire-spread diagnostics
    plot_vars = phi R arrival_time fireline_intensity flame_length

If ``plot_vars`` is **not** specified (the default), all fields are written.
If a name in the list does not match any known variable a warning is printed
and that entry is ignored.  If the resulting list is empty after filtering,
the full plotfile is written as a fallback.

Complete Variable Name Reference
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following table lists every variable name accepted by ``plot_vars``.
Use the exact string shown in the **Variable name** column.

+------------------------------------+---------------------------------------------------+
| Variable name                      | Description                                       |
+====================================+===================================================+
| ``phi``                            | Level-set function (fire front where phi = 0)     |
+------------------------------------+---------------------------------------------------+
| ``velx``, ``vely`` [, ``velz``]    | Wind velocity components [m/s]                    |
+------------------------------------+---------------------------------------------------+
| ``farsite_dx``, ``farsite_dy``     | FARSITE spread displacements                      |
| [, ``farsite_dz``]                 |                                                   |
+------------------------------------+---------------------------------------------------+
| ``R``                              | Rate of spread [m/s]                              |
+------------------------------------+---------------------------------------------------+
| ``spot_prob``                      | Firebrand spotting probability                    |
+------------------------------------+---------------------------------------------------+
| ``spot_count``                     | Firebrand count per cell                          |
+------------------------------------+---------------------------------------------------+
| ``spot_dist``                      | Firebrand landing distance [m]                    |
+------------------------------------+---------------------------------------------------+
| ``spot_active``                    | Spotting ignition flag                            |
+------------------------------------+---------------------------------------------------+
| ``fuel_consumption``               | Bulk fuel consumption fraction                    |
+------------------------------------+---------------------------------------------------+
| ``crown_fraction``                 | Crown fire fraction                               |
+------------------------------------+---------------------------------------------------+
| ``albini_Hz``                      | Albini lofting height [m]                         |
+------------------------------------+---------------------------------------------------+
| ``albini_count``                   | Firebrands launched per cell                      |
+------------------------------------+---------------------------------------------------+
| ``albini_dist``                    | Max firebrand landing distance [m]                |
+------------------------------------+---------------------------------------------------+
| ``albini_active``                  | Albini spot ignition flag                         |
+------------------------------------+---------------------------------------------------+
| ``elevation``                      | Terrain elevation [m]                             |
+------------------------------------+---------------------------------------------------+
| ``slope``                          | Terrain slope [degrees]                           |
+------------------------------------+---------------------------------------------------+
| ``aspect``                         | Terrain aspect [degrees]                          |
+------------------------------------+---------------------------------------------------+
| ``fuel_model``                     | Fuel model number (from landscape)                |
+------------------------------------+---------------------------------------------------+
| ``fireline_intensity``             | Byram fireline intensity [kW/m]                   |
+------------------------------------+---------------------------------------------------+
| ``flame_length``                   | Byram flame length [m]                            |
+------------------------------------+---------------------------------------------------+
| ``weise_flame_height``             | Weise & Biging flame height [m]                   |
+------------------------------------+---------------------------------------------------+
| ``weise_flame_tilt``               | Flame tilt angle [rad]                            |
+------------------------------------+---------------------------------------------------+
| ``weise_whirl_height``             | Fire whirl height [m]                             |
+------------------------------------+---------------------------------------------------+
| ``weise_whirl_radius``             | Fire whirl radius [m]                             |
+------------------------------------+---------------------------------------------------+
| ``weise_angular_velocity``         | Fire whirl angular velocity [rad/s]               |
+------------------------------------+---------------------------------------------------+
| ``weise_max_tang_vel``             | Fire whirl max tangential velocity [m/s]          |
+------------------------------------+---------------------------------------------------+
| ``viegas_ROS``                     | Viegas eruptive rate of spread [m/s]              |
+------------------------------------+---------------------------------------------------+
| ``viegas_eruptive_flag``           | Viegas eruptive fire flag (0/1)                   |
+------------------------------------+---------------------------------------------------+
| ``viegas_ROS_excess``              | Viegas ROS excess above threshold                 |
+------------------------------------+---------------------------------------------------+
| ``viegas_flame_tilt``              | Viegas flame tilt angle [rad]                     |
+------------------------------------+---------------------------------------------------+
| ``viegas_slope_factor``            | Viegas slope amplification factor                 |
+------------------------------------+---------------------------------------------------+
| ``scorch_height``                  | Van Wagner scorch height [m]                      |
+------------------------------------+---------------------------------------------------+
| ``prob_ignition``                  | Anderson probability of ignition [0–1]            |
+------------------------------------+---------------------------------------------------+
| ``tree_mortality``                 | Ryan–Reinhardt tree mortality [0–1]               |
+------------------------------------+---------------------------------------------------+
| ``crown_activity``                 | Crown activity (0=surface, 1=passive, 2=active)   |
+------------------------------------+---------------------------------------------------+
| ``torching_ratio``                 | Scott–Reinhardt torching ratio                    |
+------------------------------------+---------------------------------------------------+
| ``crowning_ratio``                 | Scott–Reinhardt crowning ratio                    |
+------------------------------------+---------------------------------------------------+
| ``energy_release_component``       | NFDRS energy release component                    |
+------------------------------------+---------------------------------------------------+
| ``nfdrs_spread_component``         | NFDRS spread component                            |
+------------------------------------+---------------------------------------------------+
| ``nfdrs_burning_index``            | NFDRS burning index                               |
+------------------------------------+---------------------------------------------------+
| ``co2_emissions``                  | CO₂ emission rate [kg/m²/s]                       |
+------------------------------------+---------------------------------------------------+
| ``co_emissions``                   | CO emission rate [kg/m²/s]                        |
+------------------------------------+---------------------------------------------------+
| ``pm25_emissions``                 | PM₂.₅ emission rate [kg/m²/s]                     |
+------------------------------------+---------------------------------------------------+
| ``arrival_time``                   | Fire arrival time [s]                             |
+------------------------------------+---------------------------------------------------+
| ``heat_per_unit_area``             | Heat per unit area [kJ/m²]                        |
+------------------------------------+---------------------------------------------------+
| ``vorticity_z``                    | Vertical vorticity [1/s]                          |
+------------------------------------+---------------------------------------------------+
| ``cbh``                            | Canopy base height [m]                            |
+------------------------------------+---------------------------------------------------+
| ``cbd``                            | Canopy bulk density [kg/m³]                       |
+------------------------------------+---------------------------------------------------+
| ``canopy_cover``                   | Canopy cover fraction [0–1]                       |
+------------------------------------+---------------------------------------------------+
| ``canopy_height``                  | Canopy height [m]                                 |
+------------------------------------+---------------------------------------------------+
| ``burnout_time``                   | Cell burnout time [s]                             |
+------------------------------------+---------------------------------------------------+
| ``reaction_intensity``             | Rothermel reaction intensity [BTU/ft²/min]        |
+------------------------------------+---------------------------------------------------+
| ``residual_fuel``                  | Residual (unconsumed) fuel fraction               |
+------------------------------------+---------------------------------------------------+
| ``shade_fraction``                 | Terrain/canopy shade fraction [0–1]               |
+------------------------------------+---------------------------------------------------+
| ``torching_index_kmh``             | Scott–Reinhardt torching index [km/h]             |
+------------------------------------+---------------------------------------------------+
| ``crowning_index_kmh``             | Scott–Reinhardt crowning index [km/h]             |
+------------------------------------+---------------------------------------------------+
| ``moisture_d1``                    | 1-hr dead fuel moisture [fraction]                |
+------------------------------------+---------------------------------------------------+
| ``moisture_d10``                   | 10-hr dead fuel moisture [fraction]               |
+------------------------------------+---------------------------------------------------+
| ``moisture_d100``                  | 100-hr dead fuel moisture [fraction]              |
+------------------------------------+---------------------------------------------------+
| ``moisture_lh``                    | Live herbaceous moisture [fraction]               |
+------------------------------------+---------------------------------------------------+
| ``moisture_lw``                    | Live woody moisture [fraction]                    |
+------------------------------------+---------------------------------------------------+

Example: minimal fire-spread output
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To reduce disk usage to the bare essentials for a fire-spread study::

    plot_vars = phi R arrival_time fireline_intensity flame_length

Example: ecology and emissions only
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To write only post-fire ecology and air-quality diagnostics::

    plot_vars = phi scorch_height prob_ignition tree_mortality \
                co2_emissions co_emissions pm25_emissions

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
