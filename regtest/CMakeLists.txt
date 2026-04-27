# CMakeLists.txt for regression tests
# This file enables building and running regression tests for wildfire_levelset

cmake_minimum_required(VERSION 3.20)

# Include CTest for test framework
enable_testing()

# Define the levelset executable path (assumes build in build/ directory)
set(LEVELSET_EXE "${CMAKE_BINARY_DIR}/levelset" CACHE FILEPATH "Path to levelset executable")

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
endfunction()

# Add all regression tests
add_regression_test(basic_levelset basic_levelset)
add_regression_test(farsite_ellipse farsite_ellipse)
add_regression_test(rothermel_fuel rothermel_fuel)
add_regression_test(anderson_lw anderson_lw)
add_regression_test(reinitialization reinitialization)
add_regression_test(ellipse_sdf ellipse_sdf)
add_regression_test(eb_implicit eb_implicit)

# 3D tests (only if building in 3D mode)
if(NOT LEVELSET_DIM_2D)
    add_regression_test(3d_sphere 3d_sphere)
endif()

# 2D-specific tests (terrain_wind requires 2D build with data files)
if(LEVELSET_DIM_2D)
    add_regression_test(terrain_wind terrain_wind)
endif()

# Create a custom target to run all regression tests
add_custom_target(regtest
                 COMMAND ${CMAKE_CTEST_COMMAND} -L regtest --output-on-failure
                 COMMENT "Running all regression tests"
                 USES_TERMINAL)

# Display information about configured tests
message(STATUS "Configured regression tests:")
message(STATUS "  - basic_levelset: Basic level-set advection with constant velocity")
message(STATUS "  - farsite_ellipse: FARSITE elliptical fire expansion")
message(STATUS "  - rothermel_fuel: Rothermel model with fuel types")
message(STATUS "  - anderson_lw: Anderson dynamic L/W ratio")
message(STATUS "  - reinitialization: Level-set reinitialization")
message(STATUS "  - ellipse_sdf: Elliptical SDF initial condition")
message(STATUS "  - eb_implicit: EB implicit function initial condition")
if(NOT LEVELSET_DIM_2D)
    message(STATUS "  - 3d_sphere: Full 3D fire spread simulation")
endif()
if(LEVELSET_DIM_2D)
    message(STATUS "  - terrain_wind: External terrain and wind field (2D only)")
endif()
message(STATUS "Run tests with: ctest -L regtest or make regtest")
