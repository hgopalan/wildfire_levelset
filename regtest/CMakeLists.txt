# CMakeLists.txt for regression tests
# This file enables building and running regression tests for wildfire_levelset

cmake_minimum_required(VERSION 3.20)

# Include CTest for test framework
enable_testing()

# Define the levelset executable path (assumes build in build/ directory)
set(LEVELSET_EXE "${CMAKE_BINARY_DIR}/levelset" CACHE FILEPATH "Path to levelset executable")

# Path to the optional wind_solver executable
set(WIND_SOLVER_EXE "${CMAKE_BINARY_DIR}/wind_solver" CACHE FILEPATH "Path to wind_solver executable")

# Find Python3 interpreter (needed for the landfire_farsite setup step)
find_package(Python3 COMPONENTS Interpreter QUIET)

# Function to add a regression test
function(add_regression_test TEST_NAME TEST_DIR)
    set(TEST_INPUT "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_DIR}/inputs.i")
    set(TEST_WORK_DIR "${CMAKE_CURRENT_BINARY_DIR}/${TEST_DIR}")
    
    # Create working directory for test
    file(MAKE_DIRECTORY ${TEST_WORK_DIR})
    
    # Add the test
    add_test(NAME ${TEST_NAME}
             COMMAND ${LEVELSET_EXE} ${TEST_INPUT}
             WORKING_DIRECTORY ${TEST_WORK_DIR})
    
    # Set test properties
    set_tests_properties(${TEST_NAME} PROPERTIES
                        TIMEOUT 600  # 10 minute timeout
                        LABELS "regtest")
    
    # Copy any data files if they exist
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_DIR}/gaussian_hill_terrain.csv")
        file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_DIR}/gaussian_hill_terrain.csv"
             DESTINATION ${TEST_WORK_DIR})
    endif()
    
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_DIR}/gaussian_hill_wind.csv")
        file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_DIR}/gaussian_hill_wind.csv"
             DESTINATION ${TEST_WORK_DIR})
    endif()

    # Copy time-dependent wind field snapshots if present
    file(GLOB _TIME_WIND_FILES "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_DIR}/time_wind*.csv")
    foreach(_twf ${_TIME_WIND_FILES})
        file(COPY "${_twf}" DESTINATION ${TEST_WORK_DIR})
    endforeach()

    # Copy landscape files if present (farsite_ellipse landscape test)
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_DIR}/socal_chaparral_landscape.lcp")
        file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_DIR}/socal_chaparral_landscape.lcp"
             DESTINATION ${TEST_WORK_DIR})
    endif()
endfunction()

# ---------------------------------------------------------------------------
# Function to add a regression test that needs a Python setup step first
# (e.g. generating a landscape file from LANDFIRE data or a synthetic fallback)
# ---------------------------------------------------------------------------
function(add_python_setup_regression_test TEST_NAME TEST_DIR SETUP_SCRIPT)
    set(TEST_INPUT "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_DIR}/inputs.i")
    set(TEST_WORK_DIR "${CMAKE_CURRENT_BINARY_DIR}/${TEST_DIR}")

    file(MAKE_DIRECTORY ${TEST_WORK_DIR})

    if(Python3_FOUND)
        # Setup step: run the Python script to prepare test data
        set(SETUP_TEST_NAME "${TEST_NAME}_setup")
        add_test(NAME ${SETUP_TEST_NAME}
                 COMMAND ${Python3_EXECUTABLE}
                         ${CMAKE_CURRENT_SOURCE_DIR}/${TEST_DIR}/${SETUP_SCRIPT}
                 WORKING_DIRECTORY ${TEST_WORK_DIR})
        set_tests_properties(${SETUP_TEST_NAME} PROPERTIES
                            TIMEOUT 300
                            LABELS "regtest"
                            FIXTURES_SETUP "${TEST_NAME}_fixture")

        # Main test step: run the solver
        add_test(NAME ${TEST_NAME}
                 COMMAND ${LEVELSET_EXE} ${TEST_INPUT}
                 WORKING_DIRECTORY ${TEST_WORK_DIR})
        set_tests_properties(${TEST_NAME} PROPERTIES
                            TIMEOUT 600
                            LABELS "regtest"
                            FIXTURES_REQUIRED "${TEST_NAME}_fixture")
    else()
        message(STATUS "Python3 not found – skipping ${TEST_NAME} (needs Python for setup)")
    endif()
endfunction()

# Add core regression tests
# Core functionality tests (always run)
add_regression_test(basic_levelset basic_levelset)
add_regression_test(farsite_ellipse farsite_ellipse)
add_regression_test(rothermel_fuel rothermel_fuel)
add_regression_test(anderson_lw anderson_lw)
add_regression_test(cheney_gould_grassfire cheney_gould_grassfire)
# Alternative fire shape models (Richards 1990 alternatives)
add_regression_test(catchpole_demestre catchpole_demestre)
add_regression_test(wilson_spread wilson_spread)
add_regression_test(alexander_lemniscate alexander_lemniscate)

# Advanced feature tests (can be enabled as needed)
# Uncomment to run additional tests:
# add_regression_test(reinitialization reinitialization)
# add_regression_test(ellipse_sdf ellipse_sdf)
# add_regression_test(eb_implicit eb_implicit)
# add_regression_test(firebrand_spotting firebrand_spotting)
# add_regression_test(crown_initiation crown_initiation)
# add_regression_test(bulk_fuel_consumption bulk_fuel_consumption)

# 3D tests (only if building in 3D mode)
if(NOT LEVELSET_DIM_2D)
    add_regression_test(3d_sphere 3d_sphere)
endif()

# 2D-specific tests (terrain_wind and time_dependent_wind require 2D build with data files)
if(LEVELSET_DIM_2D)
    add_regression_test(terrain_wind terrain_wind)
    add_regression_test(time_dependent_wind time_dependent_wind)
endif()

# LANDFIRE / landscape-file FARSITE test (requires Python3 for setup step)
# The setup script downloads LANDFIRE data or generates a synthetic fallback
add_python_setup_regression_test(landfire_farsite landfire_farsite create_landscape.py)

# Create a custom target to run all regression tests
add_custom_target(regtest
                 COMMAND ${CMAKE_CTEST_COMMAND} -L regtest --output-on-failure
                 COMMENT "Running all regression tests"
                 USES_TERMINAL)

# Display information about configured tests
message(STATUS "Configured core regression tests:")
message(STATUS "  - basic_levelset: Basic level-set advection with constant velocity")
message(STATUS "  - farsite_ellipse: FARSITE elliptical fire expansion (Richards 1990)")
message(STATUS "  - rothermel_fuel: Rothermel model with fuel types")
message(STATUS "  - anderson_lw: Anderson dynamic L/W ratio")
message(STATUS "  - cheney_gould_grassfire: Cheney & Gould (1995/1998) grassland fire spread")
message(STATUS "  - catchpole_demestre: Catchpole & de Mestre (1986) double-ellipse shape")
message(STATUS "  - wilson_spread: Wilson (1988) single ellipse from rear focus")
message(STATUS "  - alexander_lemniscate: Alexander et al. lemniscate (Limaçon) shape")
if(NOT LEVELSET_DIM_2D)
    message(STATUS "  - 3d_sphere: Full 3D fire spread simulation")
endif()
if(LEVELSET_DIM_2D)
    message(STATUS "  - terrain_wind: External terrain and wind field (2D only)")
    message(STATUS "  - time_dependent_wind: Time-dependent wind field with linear interpolation (2D only)")
endif()
message(STATUS "")
message(STATUS "Note: Additional tests are available but commented out for faster testing.")
message(STATUS "      Uncomment tests in regtest/CMakeLists.txt to enable them:")
message(STATUS "      reinitialization, ellipse_sdf, eb_implicit, firebrand_spotting,")
message(STATUS "      crown_initiation, bulk_fuel_consumption")
message(STATUS "")
message(STATUS "Run tests with: ctest -L regtest or make regtest")
message(STATUS "")
if(Python3_FOUND)
    message(STATUS "Python3 found (${Python3_EXECUTABLE}): landfire_farsite regtest enabled")
    message(STATUS "  - landfire_farsite: FARSITE with auto-downloaded LANDFIRE landscape")
    message(STATUS "    (downloads real terrain/fuel data if internet is available;")
    message(STATUS "     falls back to synthetic data if offline)")
else()
    message(STATUS "Python3 not found: landfire_farsite regtest disabled")
    message(STATUS "  Install Python3 and run: pip install requests rasterio numpy pyproj")
endif()

# ---------------------------------------------------------------------------
# New feature tests (added 2025)
# ---------------------------------------------------------------------------

# fuel_adj_file: .adj fuel adjustment file reader
if(Python3_FOUND)
    add_python_setup_regression_test(fuel_adj_file fuel_adj_file create_adj.py)
endif()

# fmd_moisture: time-varying fuel moisture schedule (.fmd)
if(Python3_FOUND)
    add_python_setup_regression_test(fmd_moisture fmd_moisture create_fmd.py)
endif()

# fire_perimeter_output: CSV + GeoJSON perimeter writers + burned-area CSV
add_regression_test(fire_perimeter_output fire_perimeter_output)
# Validation step (checks generated files)
if(Python3_FOUND)
    set(_perim_work "${CMAKE_CURRENT_BINARY_DIR}/fire_perimeter_output")
    add_test(NAME fire_perimeter_output_validate
             COMMAND ${Python3_EXECUTABLE}
                     ${CMAKE_CURRENT_SOURCE_DIR}/fire_perimeter_output/validate.py
             WORKING_DIRECTORY ${_perim_work})
    set_tests_properties(fire_perimeter_output_validate PROPERTIES
                         TIMEOUT 60
                         LABELS "regtest"
                         DEPENDS fire_perimeter_output)
endif()

# burned_area_timeseries: fire_stats.csv time series
add_regression_test(burned_area_timeseries burned_area_timeseries)

# ---------------------------------------------------------------------------
# Multi-resolution timing benchmark
# ---------------------------------------------------------------------------
# Runs the solver at several grid sizes for two scenarios (levelset and farsite)
# and checks for monotonically increasing runtimes and plausible scaling.
# ---------------------------------------------------------------------------
if(Python3_FOUND)
    set(TIMING_WORK_DIR "${CMAKE_CURRENT_BINARY_DIR}/timing_benchmark")
    file(MAKE_DIRECTORY ${TIMING_WORK_DIR})

    # Determine spatial dimension flag for the benchmark script
    if(LEVELSET_DIM_2D)
        set(_TIMING_DIM 2)
        set(_TIMING_RES "32;64;128")   # smaller set for faster CI
    else()
        set(_TIMING_DIM 3)
        set(_TIMING_RES "16;32;48")    # 3-D is expensive; keep small
    endif()

    # Convert semicolon-separated CMake list to space-separated args
    string(REPLACE ";" " " _TIMING_RES_ARGS "${_TIMING_RES}")
    separate_arguments(_TIMING_RES_LIST NATIVE_COMMAND "${_TIMING_RES_ARGS}")

    add_test(NAME timing_benchmark
             COMMAND ${Python3_EXECUTABLE}
                     ${CMAKE_CURRENT_SOURCE_DIR}/timing_benchmark/run_benchmark.py
                     --exe       ${LEVELSET_EXE}
                     --dim       ${_TIMING_DIM}
                     --nsteps    20
                     --resolutions ${_TIMING_RES_LIST}
                     --timeout   300
                     --output    timing_results.csv
                     --skip-scaling   # scaling check is advisory; don't fail CI
             WORKING_DIRECTORY ${TIMING_WORK_DIR})
    set_tests_properties(timing_benchmark PROPERTIES
                         TIMEOUT 1200        # 20 min overall
                         LABELS  "regtest;benchmark"
                         # Re-run on failure to distinguish transient noise
                         RUN_SERIAL TRUE)
    message(STATUS "Timing benchmark enabled:")
    message(STATUS "  - timing_benchmark: multi-resolution wall-clock benchmark "
                   "(dim=${_TIMING_DIM}, resolutions=${_TIMING_RES})")
else()
    message(STATUS "Python3 not found: timing_benchmark disabled")
endif()

message(STATUS "")
message(STATUS "New feature regtests (2025):")
message(STATUS "  - fuel_adj_file       : FARSITE .adj fuel adjustment reader")
message(STATUS "  - fmd_moisture        : Time-varying .fmd fuel moisture schedule")
message(STATUS "  - fire_perimeter_output: CSV + GeoJSON perimeter writers + stats CSV")
message(STATUS "  - burned_area_timeseries: fire_stats.csv time series")
message(STATUS "  - timing_benchmark    : Multi-resolution wall-clock timing (benchmark label)")

# ---------------------------------------------------------------------------
# wind_solver regression tests (3-D builds only, requires linear solver)
# ---------------------------------------------------------------------------
# Function to register a wind_solver regression test
function(add_wind_solver_regression_test TEST_NAME TEST_DIR)
    if(LEVELSET_DIM_2D)
        message(STATUS "  Skipping wind_solver test '${TEST_NAME}' (requires 3-D build)")
        return()
    endif()

    set(TEST_INPUT   "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_DIR}/inputs.i")
    set(TEST_WORK_DIR "${CMAKE_CURRENT_BINARY_DIR}/${TEST_DIR}")
    file(MAKE_DIRECTORY ${TEST_WORK_DIR})

    # Copy terrain file if present
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_DIR}/terrain.csv")
        file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_DIR}/terrain.csv"
             DESTINATION ${TEST_WORK_DIR})
    endif()

    add_test(NAME ${TEST_NAME}
             COMMAND ${WIND_SOLVER_EXE} ${TEST_INPUT}
             WORKING_DIRECTORY ${TEST_WORK_DIR})
    set_tests_properties(${TEST_NAME} PROPERTIES
                         TIMEOUT 600
                         LABELS  "regtest;wind_solver")
endfunction()

if(LEVELSET_BUILD_WIND_SOLVER AND NOT LEVELSET_DIM_2D)
    add_wind_solver_regression_test(gaussian_hill_wind_solver gaussian_hill_wind_solver)
    message(STATUS "")
    message(STATUS "wind_solver regtests:")
    message(STATUS "  - gaussian_hill_wind_solver: log-law + mass-consistent Poisson solve over a Gaussian hill")
else()
    message(STATUS "")
    message(STATUS "wind_solver regtests: disabled (LEVELSET_BUILD_WIND_SOLVER=OFF or 2-D build)")
endif()
