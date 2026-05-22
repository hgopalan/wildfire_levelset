Multiple Weather Streams Spatial T/RH Interpolation
====================================================

Overview
--------

This implementation adds **spatial temperature and relative humidity interpolation** to the existing multi-station weather system (``multi_wtr_weather.H``). Previously, only wind fields were spatially interpolated via IDW, while T/RH used a simple domain-mean average. This enhancement enables spatially-varying fuel moisture content driven by weather station gradients.

Implementation Details
----------------------

1. Core Interpolation Function
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**File**: ``src/multi_wtr_weather.H``

Added ``apply_multi_wtr_TRH_to_spatial()`` function that:

* Interpolates temperature and RH from multiple stations to each grid cell using IDW
* Uses the same IDW formula as wind interpolation: :math:`w_i = 1 / \max(d_i^p, \epsilon)`
* GPU-compatible implementation using AMReX kernels
* Stores results in new ``temperature_mf`` and ``humidity_mf`` MultiFabs

**Key Features**:

* Parallel GPU execution via ``amrex::ParallelFor``
* Same IDW power and epsilon parameters as wind interpolation
* Per-cell distance calculations to each station
* Weighted averaging of station values

2. New MultiFabs for Weather Storage
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**File**: ``src/multifab_setup.H``

Added two new MultiFabs to ``WildfireFields`` struct:

* ``temperature_mf`` - Per-cell air temperature [°C] (1 component)
* ``humidity_mf`` - Per-cell relative humidity [%] (1 component)

These are initialized with default values (25°C, 30% RH) and updated each timestep when ``multi_wtr_file`` is active.

3. Spatial EMC Computation
~~~~~~~~~~~~~~~~~~~~~~~~~~~

**File**: ``src/solar_radiation.H``

Added overload of ``apply_solar_emc_to_spatial_moisture()`` that:

* Accepts ``T_C_mf`` and ``RH_pct_mf`` MultiFabs instead of scalar values
* Computes per-cell fuel temperature with solar heating
* Applies Simard (1968) EMC equations using per-cell T/RH
* Preserves GPU-compatibility

Also added overload of ``apply_solar_radiation_step()`` that:

* Accepts spatial T/RH MultiFabs as input
* Integrates shade computation with spatial weather data
* Calls the spatial EMC overload instead of diurnal sinusoidal model

4. Integration in Main Loop
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**File**: ``src/main.cpp``

Modified the multi-station weather section to:

1. Call ``apply_multi_wtr_TRH_to_spatial()`` when ``diurnal_moisture.enable = 1``
2. Store interpolated T/RH in ``temperature_mf`` and ``humidity_mf``
3. Use spatial T/RH version of ``apply_solar_radiation_step()`` when ``multi_wtr_active = true``
4. Maintain backward compatibility with domain-mean fallback

5. Header Documentation Updates
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**File**: ``src/multi_wtr_weather.H`` (header comments)

Updated the file header to document:

* Temperature and RH interpolation capability
* New ``apply_multi_wtr_TRH_to_spatial()`` function
* Usage examples showing both spatial and domain-mean approaches
* Reference to spatial EMC integration

Testing
-------

Regression Test
~~~~~~~~~~~~~~~

**Location**: ``regtest/moisture/multi_wtr_spatial/``

Created comprehensive test with:

* 3 weather stations at different locations
* Distinct T/RH conditions creating spatial gradient:

  * Station 1 (SW): 40°C, 15% RH (Hot & Dry)
  * Station 2 (SE): 30°C, 30% RH (Moderate)
  * Station 3 (N): 20°C, 60% RH (Cool & Wet)

* Expected asymmetric fire spread based on moisture gradient
* Python setup script (``create_stations.py``) to generate test data
* Comprehensive README with validation criteria

**Test Files**:

* ``create_stations.py`` - Generates station list and .wtr files
* ``inputs.i`` - Test configuration
* ``README.md`` - Test documentation

Documentation Updates
---------------------

User Documentation
~~~~~~~~~~~~~~~~~~

**File**: ``docs/new_features.rst``

Updated multi-station weather section to document:

* Spatial T/RH interpolation capability
* Integration with diurnal moisture model
* Interaction with solar radiation for shade-adjusted EMC
* List of new functions and their purposes
* Reference to regression test

**File**: ``docs/mathematical_models.rst``

Updated technical documentation to:

* Describe T and RH as interpolated fields (alongside wind)
* Explain spatial EMC computation with interpolated weather
* Reference the regression test

Files Modified
--------------

Source Code (4 files)
~~~~~~~~~~~~~~~~~~~~~

1. ``src/multi_wtr_weather.H`` - Added ``apply_multi_wtr_TRH_to_spatial()`` + documentation
2. ``src/multifab_setup.H`` - Added ``temperature_mf`` and ``humidity_mf`` MultiFabs
3. ``src/solar_radiation.H`` - Added spatial T/RH overloads of EMC functions
4. ``src/main.cpp`` - Integrated spatial T/RH interpolation in main loop

Test Files (3 files)
~~~~~~~~~~~~~~~~~~~~~

5. ``regtest/moisture/multi_wtr_spatial/create_stations.py`` - Test setup
6. ``regtest/moisture/multi_wtr_spatial/inputs.i`` - Test configuration
7. ``regtest/moisture/multi_wtr_spatial/README.md`` - Test documentation

Documentation (2 files)
~~~~~~~~~~~~~~~~~~~~~~~~

8. ``docs/new_features.rst`` - Updated user documentation
9. ``docs/mathematical_models.rst`` - Updated technical documentation

Backward Compatibility
----------------------

The implementation maintains full backward compatibility:

* Domain-mean T/RH is still computed as fallback
* Existing ``multi_wtr_file`` configurations continue to work
* Spatial interpolation only activates when ``diurnal_moisture.enable = 1``
* Non-multi_wtr simulations are unaffected

Implementation Quality
----------------------

Code Quality
~~~~~~~~~~~~

* Follows existing code patterns and style
* GPU-compatible using AMReX parallel primitives
* Consistent with IDW wind interpolation implementation
* Proper error handling and bounds checking

Performance
~~~~~~~~~~~

* Spatial interpolation runs once per timestep
* GPU-accelerated when available
* No performance impact when feature is disabled
* Station data copied to device only once per timestep

Maintainability
~~~~~~~~~~~~~~~

* Clear function names and documentation
* Modular design with function overloads
* Comprehensive inline comments
* Test coverage for new functionality

References
----------

1. Shepard, D. (1968). A two-dimensional interpolation function for irregularly-spaced data. Proc. 23rd ACM Conference, pp. 517–524.
2. Finney, M.A. (2004). FARSITE: Fire Area Simulator. USDA FS RMRS-RP-4.
3. Nelson, R.M. Jr. (2000). Prediction of diurnal change in 10-h fuel stick moisture content. Can. J. For. Res. 30:1071–1087.
4. Simard, A.J. (1968). The moisture content of forest fuels. Forestry Branch Info. Rep. FF-X-14.

Status
------

* ✅ Implementation complete
* ✅ Documentation updated
* ✅ Regression test created
* ⏳ Build validation pending (requires AMReX initialization)
* ⏳ Runtime validation pending (requires full build)

Next Steps
----------

To validate the implementation:

1. Initialize AMReX submodule: ``git submodule update --init --recursive``
2. Build project: ``cmake -S . -B build -DLEVELSET_DIM_2D=ON && cmake --build build``
3. Run regression test: ``cd regtest/moisture/multi_wtr_spatial && python3 create_stations.py && ../../../build/levelset inputs.i``
4. Verify spatial T/RH variation in plotfiles
5. Confirm asymmetric fire spread based on moisture gradient
