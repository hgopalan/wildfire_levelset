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

**Build**: CPU-only build (Albini spotting requires serial execution).

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

landfire_farsite
^^^^^^^^^^^^^^^^

**Purpose**: Integration test for FARSITE with a real (or synthetic) LANDFIRE landscape.

The setup script ``regtest/landfire_farsite/create_landscape.py`` attempts to download
LANDFIRE data; if offline, a synthetic Southern California chaparral landscape is generated as
a fallback.  FARSITE elliptical spread with Anderson dynamic L/W ratio is then run on the
landscape.

**Build**: Default 3D build.  Requires Python 3 with
``pip install requests rasterio numpy pyproj elevation``.

gaussian_hill_wind_solver
^^^^^^^^^^^^^^^^^^^^^^^^^

**Purpose**: Validates the terrain-following mass-consistent wind solver (``wind_solver``
executable) over a Gaussian hill terrain.

A 10 × 10 × 10 cell Cartesian grid (dx = dy = dz = 30 m) covers a 300 m × 300 m ×
300 m domain.  The terrain is a Gaussian hill with peak elevation 50 m centred at
(150, 150) m (σ = 60 m).  A 10 m/s westerly reference wind at z_ref = 10 m with
roughness length z₀ = 0.1 m is applied.

Successful completion confirms that:

* The log-law profile initialisation runs without error.
* The AMReX MLMG (``MLABecLaplacian``) solver converges to the requested tolerance.
* The maximum divergence after correction is substantially smaller than before.
* An AMReX plotfile ``plt_wind_gaussian`` is written.

**Build**: Default 3D build with ``LEVELSET_BUILD_WIND_SOLVER=ON`` (the default).

**Label**: ``regtest;wind_solver`` — can be run selectively::

    ctest -L wind_solver --output-on-failure

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
   * - ``albini_spotting``
     - CPU-only
     - no GPU/OpenMP flags
   * - ``gaussian_hill_wind_solver``
     - 3D (default)
     - ``-DLEVELSET_BUILD_WIND_SOLVER=ON`` (default)

Adding a New Regression Test
------------------------------

1. Create ``regtest/<test_name>/`` directory.
2. Add ``inputs.i`` with the test parameters.
3. Add a ``README.md`` documenting expected behavior.
4. Register in ``CMakeLists.txt`` so CTest discovers it.
5. Add a description to the main ``regtest/README.md`` and to this page.
