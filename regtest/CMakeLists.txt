# CMakeLists.txt for regression tests
# This file enables building and running regression tests for wildfire_levelset

cmake_minimum_required(VERSION 3.20)

# Include CTest for test framework
enable_testing()

# Define the levelset executable path (assumes build in build/ directory)
set(LEVELSET_EXE "${CMAKE_BINARY_DIR}/levelset" CACHE FILEPATH "Path to levelset executable")

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

# 2D-specific tests (terrain_wind requires 2D build with data files)
if(LEVELSET_DIM_2D)
    add_regression_test(terrain_wind terrain_wind)
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
message(STATUS "  - farsite_ellipse: FARSITE elliptical fire expansion")
message(STATUS "  - rothermel_fuel: Rothermel model with fuel types")
message(STATUS "  - anderson_lw: Anderson dynamic L/W ratio")
if(NOT LEVELSET_DIM_2D)
    message(STATUS "  - 3d_sphere: Full 3D fire spread simulation")
endif()
if(LEVELSET_DIM_2D)
    message(STATUS "  - terrain_wind: External terrain and wind field (2D only)")
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
