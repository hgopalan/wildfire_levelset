Regression Tests
================

The ``regtest/`` directory contains a comprehensive suite of regression tests that exercise
individual features and capabilities of the wildfire level-set solver.  Each test is
self-contained — it ships with an ``inputs.i`` configuration file and any required data files.

Running the Tests
-----------------

**Run a single test** from the repository root::

    cd regtest/basic_levelset
    ../../build/levelset inputs.i

**Run all regression tests with CTest** (after building)::

    cd build
    ctest -L regtest --output-on-failure

    # Or use the custom target
    make regtest

    # Run a specific test by name
    ctest -R rothermel_fuel --output-on-failure

Build requirements vary by test; see the per-test sections below.

Test Descriptions
-----------------

basic_levelset
^^^^^^^^^^^^^^

**Purpose**: Verifies fundamental level-set advection with a constant velocity field.

The spherical initial condition is advected without distortion for 100 steps on a 64³ grid.
Successful completion confirms that:

* The WENO5-Z + RK3 advection kernel compiles and runs.
* Reinitialization keeps :math:`|\nabla\phi| \approx 1`.
* Plotfile output is written correctly.

**Build**: Default 3D build.

farsite_ellipse
^^^^^^^^^^^^^^^

**Purpose**: Verifies the FARSITE elliptical expansion model (Richards 1990).

A box (line) ignition expands under a constant wind field using fixed Richards coefficients
(:math:`a = 1.0`, :math:`b = 0.4`, :math:`c = 0.2`) and a fixed length-to-width ratio of 3.0.
Confirms that the fire perimeter forms an ellipse with the expected elongation.

**Build**: Default 3D build.

rothermel_fuel
^^^^^^^^^^^^^^

**Purpose**: Tests fire spread with NFFL fuel model FM1 (Short Grass).

The Rothermel (1972) model is evaluated with a standard fuel database entry, 6 % fuel moisture,
and a moderate constant wind.  Confirms that the Rothermel ROS field is computed for every cell
and that different fuel models produce different spread rates.

**Build**: Default 3D build.  Change ``rothermel.fuel_model`` to test FM1–FM13 or any
Scott & Burgan 40 code.

cheney_gould_grassfire
^^^^^^^^^^^^^^^^^^^^^^

**Purpose**: Tests the Cheney & Gould (1995/1998) empirical grassland fire spread model.

A spherical ignition on a flat 1 km × 1 km domain is propagated at 10 m/s wind (36 km/h)
with 8 % dead fuel moisture and full curing (:math:`CF = 1`).  The expected head-fire ROS is
≈ 2.69 m/s from the piecewise formula:

.. math::

   R = (-0.020 + 0.640 \times 36) \times \exp(-0.108 \times 8) \times 1.0 \approx 2.69 \ \text{m/s}

Confirms that the Cheney–Gould model produces an asymmetric elliptical perimeter.

**Build**: Default 3D build.

terrain_wind
^^^^^^^^^^^^

**Purpose**: Tests external terrain (Gaussian hill) and a spatially-varying wind field (2D only).

A 1000 m × 1000 m domain contains a Gaussian hill with σ = 150 m and peak elevation 100 m.
The wind CSV has a base speed of 5 m/s with up to 50 % speed-up over the crest.  Confirms:

* Terrain slopes are correctly computed from an XYZ elevation file.
* IDW interpolation of wind from the CSV to the AMReX grid works.
* The fire accelerates upslope and over the crest.

**Build**: 2D build (``-DLEVELSET_DIM_2D=ON``).

anderson_lw
^^^^^^^^^^^

**Purpose**: Tests dynamic L/W ratio calculation based on wind speed (Anderson 1983).

The Anderson formula is evaluated at ≈ 10 mph wind (expected L/W ≈ 2.5) and the resulting
ellipse coefficients :math:`a`, :math:`b`, :math:`c` are applied to the FARSITE propagation.

.. math::

   L/W = 0.936 \exp(0.2566 U) + 0.461 \exp(-0.1548 U) - 0.397

**Build**: Default 3D build.

reinitialization
^^^^^^^^^^^^^^^^

**Purpose**: Tests aggressive level-set reinitialization to maintain the signed-distance property.

Reinitialization is run every 5 steps for 30 iterations.  Confirms that :math:`|\nabla\phi|`
remains close to 1.0 after prolonged advection.

**Build**: Default 3D build.

ellipse_sdf
^^^^^^^^^^^

**Purpose**: Tests the elliptical SDF initial condition.

An ellipsoidal initial fire region (semi-axes :math:`r_x = 0.25`, :math:`r_y = 0.15`,
:math:`r_z = 0.10`) is advected at constant velocity.  Confirms that the elliptical SDF is
constructed correctly and that shape is preserved under advection.

**Build**: Default 3D build.

eb_implicit
^^^^^^^^^^^

**Purpose**: Tests embedded boundary (EB) implicit function initial conditions.

An EB-defined ellipsoid is used as the initial fire region.  Confirms that the AMReX EB
infrastructure can initialise :math:`\phi` from implicit-function geometries (sphere,
ellipsoid, cylinder, plane).

**Build**: EB build recommended (``-DLEVELSET_ENABLE_EB=ON``).

firebrand_spotting
^^^^^^^^^^^^^^^^^^

**Purpose**: Tests the probability-based stochastic firebrand spotting model.

Firebrands are generated probabilistically from the fire front, with landing distances sampled
from a log-normal distribution.  Confirms that spot ignitions appear ahead of the fire perimeter
and that diagnostic fields (``spot_prob``, ``spot_count``, ``spot_dist``, ``spot_active``) are
written to every plotfile.

**Build**: Default 3D build.

albini_spotting
^^^^^^^^^^^^^^^

**Purpose**: Tests the Albini (1983) physics-based firebrand spotting model.

The Albini model lofts firebrands using:

.. math::

   H_z \;[\text{m}] = 12.2 \; I_B^{1/3}

and then integrates a 2-D horizontal trajectory through the wind field.  Confirms that the
lofting height, landing distance, and diagnostic fields (``albini_Hz``, ``albini_count``,
``albini_dist``, ``albini_active``) are computed and written correctly.

**Build**: CPU or GPU builds (GPU-safe: trajectory runs on the host; device
data is synchronized and copied to host-pinned memory before computation, then
written back to device).

albini_spotting_3d_wind
^^^^^^^^^^^^^^^^^^^^^^^

**Purpose**: Tests the Albini spotting model driven by 3-D wind data from a
`massconsistent_amr <https://github.com/hgopalan/massconsistent_amr>`_ plotfile.

A synthetic AMReX plotfile is generated by ``generate_plt_wind.py`` (an 8×8×4
grid with uniform u = 5 m/s, v = 0.5 m/s wind).  The reader in
``src/plt_wind_reader.H`` loads the plotfile as flat 1-D GPU arrays of
``(x, y, z, u, v, w)``, computes a height-averaged 2-D wind, and bilinearly
interpolates onto the fire-model grid for trajectory integration.

Requires Python 3 to generate the synthetic plt file (handled automatically
by CTest as a fixture setup step).

**Build**: CPU or GPU builds (GPU-safe).

ember_cascade_flux
^^^^^^^^^^^^^^^^^^

**Purpose**: Tests the flux-based ember cascade model (``ember_cascade_flux.H``)
driven by a synthetic 3-D wind plotfile.

Physical scenario: a chaparral ignition (FM4, 7 % moisture) on a 600 m × 600 m
domain burns under a 6 m/s westerly wind.  ``generate_plt_wind.py`` creates a
synthetic plt wind directory (8×8×4 grid, u = 6 m/s, v = 1 m/s) that the
model reads via ``plt_wind_reader.H``.

The test confirms that:

* The Albini (1983) plume height is computed per fire-front cell.
* The Gaussian landing-flux field is non-zero downwind of the fire.
* The plotfile variables ``ember_cascade_flux`` and ``ember_cascade_ignition``
  are written to every output step.
* The model produces spot-fire ignitions (``ember_cascade_ignition`` > 0)
  after the fire front reaches sufficient intensity.

Requires Python 3 for the setup step (same pattern as ``albini_spotting_3d_wind``).

**Build**: CPU or GPU builds (GPU-safe).

crown_initiation
^^^^^^^^^^^^^^^^

**Purpose**: Tests Van Wagner (1977) crown fire initiation.

The model checks whether the Byram fireline intensity exceeds the critical threshold:

.. math::

   I_o \;[\text{kW/m}] = 0.010 \; H_{\text{CBH}} \, (460 + 25.9 \, F_{\text{MC}})

Cells where crown fire initiates are classified as passive (1) or active (2) depending on
whether the ROS exceeds the critical active-crown threshold.  Confirms that ``crown_fraction``
and ``crown_activity`` fields are written to plotfiles.

**Build**: Default 3D build.

bulk_fuel_consumption
^^^^^^^^^^^^^^^^^^^^^

**Purpose**: Tests post-frontal bulk fuel consumption modelling.

The exponential consumption model:

.. math::

   f_c = f_{c,\min} + (f_{c,\max} - f_{c,\min})\left(1 - \exp\!\left(-\frac{t}{\tau_{\text{res}}}\right)\right)

is applied to cells inside the fire perimeter.  Confirms that the ``fuel_consumption`` field
increases monotonically with time behind the fire front.

**Build**: Default 3D build.

3d_sphere
^^^^^^^^^

**Purpose**: Tests full 3D fire spread simulation.

A sphere initial condition is propagated using the Rothermel model in a 48³ cell 3D domain
with a 3D wind vector and an upslope terrain.  Confirms that the 3D FARSITE elliptical
expansion and terrain effects work in three dimensions.

**Build**: Default 3D build.  Visualise results with ParaView or VisIt.

terrain_wind_preprocess
^^^^^^^^^^^^^^^^^^^^^^^

**Purpose**: Integration test for the ``wrf_wind_reader.py`` preprocessing tool.

Generates a terrain XYZ file and a wind CSV using the Python tool, then runs the solver with
those files.  Confirms that the WRF wind extraction and IDW interpolation pipeline produces
valid inputs for the C++ solver.

**Build**: 2D build.  Requires Python 3 with ``netCDF4``, ``numpy``, ``pyproj``.

time_dependent_wind
^^^^^^^^^^^^^^^^^^^

**Purpose**: Tests time-varying wind fields (2D only).

Sequential wind CSV snapshots are loaded and linearly interpolated at each time step.  Confirms
that ``use_time_dependent_wind = 1`` correctly advances through the wind time series.

**Build**: 2D build (``-DLEVELSET_DIM_2D=ON``).

balbi_viegas_heatflux
^^^^^^^^^^^^^^^^^^^^^

**Purpose**: Tests the Balbi (2009) model together with Viegas eruptive-fire diagnostics
and the heat-flux-driven wind correction.

Confirms that ``fire_spread_model = balbi``, ``viegas.enable = 1``, and
``heat_flux.enable_upward = 1`` all interoperate correctly.

**Build**: Default 3D build.

windninja_ridge_canyon
^^^^^^^^^^^^^^^^^^^^^^

**Purpose**: Tests the WindNinja ridge/canyon empirical speed-up (Option 7).

A Gaussian hill terrain is combined with ``wind_terrain.model = windninja_ridge_canyon``.
Confirms that wind is accelerated on the upslope face and channelled in canyon geometry.

**Build**: 2D build.

terrain_gradient_correction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Purpose**: Verifies that the level-set gradient :math:`|\nabla\phi|` is computed on
the terrain surface when a terrain file is present.

When terrain slopes are non-zero, ``godunov_norm_grad_phi`` replaces the flat grid
spacings with terrain arc-length spacings:

.. math::

   \Delta x_\text{eff} = \Delta x\,\sqrt{1 + s_x^2}, \qquad
   \Delta y_\text{eff} = \Delta y\,\sqrt{1 + s_y^2}

where :math:`s_x = \partial z/\partial x` and :math:`s_y = \partial z/\partial y`.

The test uses a steep Gaussian hill (:math:`H = 200` m, :math:`\sigma = 100` m) where
the slope magnitude at the inflection ring reaches :math:`\approx 1.21` (50°), giving an
effective spacing factor of :math:`\sqrt{1 + 1.21^2} \approx 1.57`.  A level-set run
(``propagation_method = levelset``) exercises the corrected WENO3 Godunov gradient across
the full hill profile.  Confirms that the solver completes 100 steps without NaN/Inf and
writes 4 plotfiles.

**Build**: 2D build (``-DLEVELSET_DIM_2D=ON``).

landfire_farsite
^^^^^^^^^^^^^^^^

**Purpose**: Integration test for FARSITE with a real (or synthetic) LANDFIRE landscape.

The setup script ``regtest/landfire_farsite/create_landscape.py`` attempts to download
LANDFIRE data; if offline, a synthetic Southern California chaparral landscape is generated as
a fallback.  FARSITE elliptical spread with Anderson dynamic L/W ratio is then run on the
landscape.

**Build**: Default 3D build.  Requires Python 3 with
``pip install requests rasterio numpy pyproj elevation``.

gaussian_hill_wind_solver (Deprecated)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. deprecated::
   The mass-consistent wind solver has been moved to ``src/deprecated/`` and is no longer
   built.  This regression test is retained for historical reference only.

**Purpose**: Validated the terrain-following mass-consistent wind solver (``wind_solver``
executable) over a Gaussian hill terrain.

A 10 × 10 × 10 cell Cartesian grid (dx = dy = dz = 30 m) covers a 300 m × 300 m ×
300 m domain.  The terrain is a Gaussian hill with peak elevation 50 m centred at
(150, 150) m (σ = 60 m).  A 10 m/s westerly reference wind at z_ref = 10 m with
roughness length z₀ = 0.1 m is applied.

Successful completion confirmed that:

* The log-law profile initialisation runs without error.
* The AMReX MLMG (``MLABecLaplacian``) solver converges to the requested tolerance.
* The maximum divergence after correction is substantially smaller than before.
* An AMReX plotfile ``plt_wind_gaussian`` is written.

**Build**: This test is no longer active (wind solver deprecated).

Build Configurations Summary
-----------------------------

.. list-table::
   :header-rows: 1
   :widths: 45 10 45

   * - Compatible tests
     - Build
     - CMake flags
   * - All except ``3d_sphere``, ``terrain_wind``, ``time_dependent_wind``
     - 3D (default)
     - ``cmake -S . -B build``
   * - ``terrain_wind``, ``terrain_wind_preprocess``, ``time_dependent_wind``,
       ``windninja_ridge_canyon``
     - 2D
     - ``-DLEVELSET_DIM_2D=ON``
   * - ``eb_implicit``
     - 3D + EB
     - ``-DLEVELSET_ENABLE_EB=ON``
   * - ``albini_spotting``, ``albini_spotting_3d_wind``
     - CPU or GPU
     - GPU-safe: synchronizes device, runs host-side trajectory

Adding a New Regression Test
------------------------------

1. Create ``regtest/<test_name>/`` directory.
2. Add ``inputs.i`` with the test parameters.
3. Add a ``README.md`` documenting expected behavior.
4. Register in ``CMakeLists.txt`` so CTest discovers it.
5. Add a description to the main ``regtest/README.md`` and to this page.

New Feature Regression Tests
-----------------------------

The following tests were added with the 2025 feature update. All use
**UTM Zone 11N, Southern California** coordinates
(reference: 330000 E, 3775000 N — Malibu / Santa Monica Mountains area).

fmc_seasonal (crown_fire/)
^^^^^^^^^^^^^^^^^^^^^^^^^^

**Purpose**: Tests the FMC seasonal phenological schedule.

The foliar moisture content (FMC) used by the Van Wagner (1977) crown fire initiation
model is updated each timestep from a built-in parametric curve representing Southern
California chaparral phenology (green-up in spring, peak in summer, curing in autumn).

Key parameters::

    fmc_schedule.enable            = 1
    fmc_schedule.use_farsite_curve = 1
    fmc_schedule.start_doy         = 200   # mid-July (summer peak)

Confirms that ``crown.FMC`` changes each timestep and that crown fire behaviour
responds to the seasonal FMC change.

precip_wetting (moisture/)
^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Purpose**: Tests precipitation-driven dead fuel moisture wetting.

Dead fuel moisture evolves under a constant light rain rate (2 mm/hr).  The
exponential wetting model drives 1-hr fuel moisture towards saturation (120%)
with a 1-hr time constant, and 10-hr fuel with a 2-hr time constant.

Key parameters::

    diurnal_moisture.enable              = 1
    diurnal_moisture.precip_rain_rate_mm_hr = 2.0
    diurnal_moisture.M_sat               = 1.20

Confirms that ``rothermel.M_d1`` increases over time and fire spread rate
decreases as fuel becomes wetter.

polygon_ignition (ignition/)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Purpose**: Tests closed-polygon fire ignition rasterisation.

An irregular 8-vertex polygon (~20 m × 18 m) is rasterised into the level-set
field using the winding-number algorithm.  The interior is set to ``phi < 0``
(burned) and the exterior to the Euclidean signed distance.

Key parameters::

    source_type        = polygon
    fire_polygon_file  = ignition_polygon.csv

Confirms that the rasterised initial perimeter closely follows the input polygon
vertices and that the signed-distance function is correct inside and outside.

polyline_ignition (ignition/)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Purpose**: Tests polyline (line fire) ignition rasterisation.

A 6-vertex east-west polyline is rasterised with a 6 m half-width, creating an
initial line-fire ignition zone.  Wind drives the fire northward from the ignition
line.

Key parameters::

    source_type       = polyline
    fire_polygon_file = ignition_line.csv
    polyline_width    = 6.0

wind_dir_schedule (wind/)
^^^^^^^^^^^^^^^^^^^^^^^^^^

**Purpose**: Tests the compact wind direction / speed schedule.

Wind starts at 270° (westerly) at 3 m/s and backs to 180° (southerly) at 5 m/s
over 2 hours.  The schedule is read from a three-column CSV
(``time_s, speed_ms, dir_deg``).

Key parameters::

    wind_dir_schedule_file = wind_schedule.csv

Confirms that the wind vector changes direction and magnitude at each timestep,
and that the fire footprint reflects the wind rotation.

waf_andrews (wind/)
^^^^^^^^^^^^^^^^^^^^

**Purpose**: Tests the Andrews (2018) logarithmic Wind Adjustment Factor (WAF)
together with the Maximum Effective Wind Speed (MEWS) cap.

The WAF converts the 20-ft reference wind to midflame height using
Albini & Baughman (1979):

.. math::

   \text{WAF} = \frac{1.83}{\ln\!\left(\dfrac{20 + 0.36\,h}{0.13\,h}\right)}

For FM4 chaparral (:math:`h = 6` ft): WAF ≈ 0.36; the 5 m/s ambient wind is
reduced to ≈ 1.8 m/s at midflame height.

Key parameters::

    rothermel.use_waf        = 1
    rothermel.waf_formula    = andrews
    rothermel.use_wind_limit = 1

Confirms that the solver runs with ``waf_formula = "andrews"`` and that fire
spread is slower than without WAF, consistent with midflame wind sheltering.

**Build**: 2D build (``-DLEVELSET_DIM_2D=ON``).

waf_behaviorplus (wind/)
^^^^^^^^^^^^^^^^^^^^^^^^^

**Purpose**: Tests the BehavePlus linear Wind Adjustment Factor for open and
shrub fuel models.

The BehavePlus linear formula:

.. math::

   \text{WAF} = 0.36 + 0.004\,h_\text{in}   \quad [h_\text{in} = 12 h_\text{ft}]

gives WAF = 0.648 for FM4 (:math:`h = 72` in).  Compare with the Andrews
logarithmic formula (WAF ≈ 0.36 for the same fuel): the BehavePlus formula
produces higher midflame wind and therefore faster ROS for deep fuel beds.

Key parameters::

    rothermel.use_waf          = 1
    rothermel.waf_formula      = behaviorplus
    rothermel.waf_canopy_alpha = 1.5
    rothermel.use_wind_limit   = 1

Confirms that the solver runs with ``waf_formula = "behaviorplus"``, that the
WAF value is computed from the linear formula, and that the fire spread rate
is higher than the ``waf_andrews`` test for the same ambient wind.

**Build**: 2D build (``-DLEVELSET_DIM_2D=ON``).

scott_reinhardt_indices (diagnostics/)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Purpose**: Tests the Scott & Reinhardt (2001) TI/CI diagnostic outputs.

A strong 8 m/s wind and low canopy base height (2.5 m) create active crown fire
conditions.  The plotfile ``torching_ratio`` and ``crowning_ratio`` fields
(ecology components 4 and 5) should be ≥ 1.0 in the active fire zone.

.. math::

    \text{torching\_ratio} = I_B / I_o \qquad
    \text{crowning\_ratio} = R / R^\prime_{SA}

These dimensionless ratios are zero for surface fire, cross 1.0 when torching
begins, and exceed 1.0 again for active crown fire, giving an easily interpreted
diagnostic layer without inverting the Rothermel wind function.

Sub-folder Organisation
-----------------------

Starting from the 2025 release all regression tests reside in named sub-folders
grouped by physical category:

+------------------+---------------------------------------------------------+
| Sub-folder       | Tests                                                   |
+==================+=========================================================+
| ``surface_spread`` | basic_levelset, farsite_ellipse, rothermel_fuel,      |
|                  | anderson_lw, catchpole_demestre, wilson_spread,         |
|                  | alexander_lemniscate, ellipse_sdf, reinitialization,    |
|                  | fbp_o1a_grassfire, fbp_o1b_grassfire,                   |
|                  | fbp_s1_slash, fbp_s2_slash, fbp_s3_slash,               |
|                  | lautenberger_spread                                     |
+------------------+---------------------------------------------------------+
| ``crown_fire``   | crown_initiation, cruz_crown_continental_us, fmc_seasonal|
+------------------+---------------------------------------------------------+
| ``spotting``     | firebrand_spotting, albini_spotting,                    |
|                  | albini_spotting_3d_wind,                                |
|                  | **ember_cascade_flux** *(new)*                          |
+------------------+---------------------------------------------------------+
| ``terrain``      | terrain_wind, balbi_viegas_heatflux,                    |
|                  | windninja_ridge_canyon,                                 |
|                  | terrain_gradient_correction                             |
+------------------+---------------------------------------------------------+
| ``moisture``     | fmd_moisture, cheney_gould_grassfire, precip_wetting,   |
|                  | spatial_moisture_output                                 |
+------------------+---------------------------------------------------------+
| ``fuel``         | fuel_adj_file                                           |
+------------------+---------------------------------------------------------+
| ``ignition``     | barrier_polygons, fire_perimeter_output,                |
|                  | polygon_ignition, polyline_ignition                     |
+------------------+---------------------------------------------------------+
| ``wind``         | time_dependent_wind, turb_wind, wind_dir_schedule,      |
|                  | waf_andrews, waf_behaviorplus                           |
+------------------+---------------------------------------------------------+
| ``diagnostics``  | scott_reinhardt_indices                                 |
+------------------+---------------------------------------------------------+
| ``misc``         | 3d_sphere, eb_implicit, mtt_propagation,                |
|                  | bulk_fuel_consumption, timing_benchmark,                |
|                  | landfire_farsite, nonburnable_mask                      |
+------------------+---------------------------------------------------------+
| ``python_api``   | **basic_fire_solver**, **coupled_wind_fire**           |
|                  | (Python bindings, requires pybind11)                   |
+------------------+---------------------------------------------------------+

All tests use UTM Zone 11N, Southern California reference coordinates
(330 000 E, 3 775 000 N).

Python API Tests
-----------------

The ``python_api/`` sub-folder contains regression tests for the Python bindings
that enable programmatic control of the fire solver from Python.

basic_fire_solver (python_api/)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Purpose**: Tests the core Python API functionality for running wildfire simulations
from Python.

The test exercises:

* Initialization from an inputs file using the ``WildfireSolver`` class
* Time-stepping with ``fire.step()``
* State extraction with ``fire.get_state()``
* Fire spread verification (burned area growth)
* Plotfile writing from Python

Key operations:

.. code-block:: python

   from wildfire_solver import WildfireSolver
   
   fire = WildfireSolver("inputs.i")
   for i in range(10):
       fire.step()
       state = fire.get_state()
   fire.finalize()

**Build**: 2D build with Python bindings::

   cmake -S . -B build -DLEVELSET_DIM_2D=ON -DLEVELSET_BUILD_PYTHON_BINDINGS=ON

**Requirements**: Python 3.6+, NumPy, pybind11

**Run**::

   cd build
   ctest -R python_api_basic_fire_solver --output-on-failure

Or directly::

   cd regtest/python_api/basic_fire_solver
   PYTHONPATH=/path/to/build/python:$PYTHONPATH python3 test_fire_solver.py

See the `Python API regression tests README <../regtest/python_api/README.md>`_
for detailed instructions.

coupled_wind_fire (python_api/)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Purpose**: Demonstrates integration with external 3D wind solvers via the Python API.

This test simulates a coupled wind-fire workflow where:

1. A 3D wind field is generated (synthetic data simulating massconsistent_amr output)
2. Wind is passed to the fire solver via ``fire.update_wind_3d()``
3. The fire solver advances one timestep
4. Steps 1-3 repeat in a coupled time loop

Key operations:

.. code-block:: python

   from wildfire_solver import WildfireSolver
   
   fire = WildfireSolver("inputs.i")
   
   while fire.time < final_time:
       # Generate/update 3D wind field (your wind solver here)
       u_3d, v_3d, w_3d, nz, zmin, zmax = generate_wind_field(...)
       
       # Update fire wind
       fire.update_wind_3d(u_3d, v_3d, w_3d, nz, zmin, zmax)
       
       # Advance fire
       fire.step()
   
   fire.finalize()

**Integrating with massconsistent_amr**:

To replace the synthetic wind generation with real wind solver output:

1. Build massconsistent_amr with Python bindings::

      cd /path/to/massconsistent_amr
      cmake -S . -B build -DBUILD_PYTHON_BINDINGS=ON
      cmake --build build -j

2. Set PYTHONPATH to include both modules::

      export PYTHONPATH=/path/to/wildfire_levelset/build/python:\
      /path/to/massconsistent_amr/build/python:$PYTHONPATH

3. Replace the synthetic wind call with actual solver calls (see test comments)

**Wind Solver Requirements**:

Any wind solver can be used if it provides:

* 3D velocity arrays ``u_3d, v_3d, w_3d`` with shape ``(nz, ny, nx)``
* Fortran (column-major) array order
* Velocities in m/s
* Vertical domain specification (``nz, zmin, zmax``)
* Domain overlap with fire solver domain

**Build**: Same as basic_fire_solver

**Run**::

   cd build
   ctest -R python_api_coupled_wind_fire --output-on-failure

**See also**:

* `Python API documentation <python_api.html>`_ - Complete API reference
* `massconsistent_amr <https://github.com/hgopalan/massconsistent_amr>`_ - Compatible wind solver
* ``regtest/python_api/README.md`` - Integration guide for other wind solvers

Replacing the Wind Solver
~~~~~~~~~~~~~~~~~~~~~~~~~~

The coupled_wind_fire test is designed to work with any wind solver. Common patterns:

**Option 1: AMReX-based solver with Python bindings** (like massconsistent_amr):

.. code-block:: python

   from pyWindSolver import WindSolver
   
   wind = WindSolver("wind_inputs.txt")
   wind.solve(fire.time)
   u_3d, v_3d, w_3d = wind.get_velocity_arrays()

**Option 2: External executable** (WindNinja, etc.):

.. code-block:: python

   import subprocess
   
   def run_external_wind_solver(time):
       subprocess.run(['WindNinja_cli', '--time', str(time), ...])
       return load_wind_output('wind_output.asc')

**Option 3: WRF or NetCDF-based models**:

.. code-block:: python

   from netCDF4 import Dataset
   import wrf
   
   u = wrf.getvar(ncfile, 'ua', timeidx=time_idx)
   # Interpolate to fire grid...

**Option 4: Custom Python wind solver**:

.. code-block:: python

   def solve_log_law_wind(fire, time):
       # Implement your wind model
       return u_3d, v_3d, w_3d, nz, zmin, zmax

See the `Python API documentation <python_api.html>`_ for complete integration examples.
