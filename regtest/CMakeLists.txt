# CMakeLists.txt for regression tests
# Tests are organised into sub-folders:
#   surface_spread/  – core fire spread model tests
#   crown_fire/      – crown fire initiation and spread
#   spotting/        – firebrand / torching spotting
#   terrain/         – terrain and wind-terrain interactions
#   moisture/        – fuel moisture models (FMD, diurnal, precipitation)
#   fuel/            – fuel models and adjustment files
#   ignition/        – ignition types (polygon, polyline, multi-point, barriers)
#   wind/            – wind models (time-dependent, turbulent, direction schedule)
#   diagnostics/     – fire ecology, emissions, Scott-Reinhardt TI/CI
#   misc/            – 3-D tests, EB, MTT, benchmarks

cmake_minimum_required(VERSION 3.20)
enable_testing()

set(LEVELSET_EXE   "${CMAKE_BINARY_DIR}/levelset"    CACHE FILEPATH "Path to levelset executable")

find_package(Python3 COMPONENTS Interpreter QUIET)

# ---------------------------------------------------------------------------
# Helper: add a solver regression test
# ---------------------------------------------------------------------------
function(add_regression_test TEST_NAME TEST_DIR)
    set(TEST_INPUT   "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_DIR}/inputs.i")
    set(TEST_WORK_DIR "${CMAKE_CURRENT_BINARY_DIR}/${TEST_DIR}")
    file(MAKE_DIRECTORY ${TEST_WORK_DIR})

    # Copy data files (wind CSVs, terrain files, polygon CSVs, etc.)
    file(GLOB _DATA_FILES
         "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_DIR}/*.csv"
         "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_DIR}/*.lcp"
         "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_DIR}/*.adj"
         "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_DIR}/*.fmd"
         "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_DIR}/*.fms")
    foreach(_df ${_DATA_FILES})
        file(COPY "${_df}" DESTINATION ${TEST_WORK_DIR})
    endforeach()

    add_test(NAME ${TEST_NAME}
             COMMAND ${LEVELSET_EXE} ${TEST_INPUT}
             WORKING_DIRECTORY ${TEST_WORK_DIR})
    set_tests_properties(${TEST_NAME} PROPERTIES
                         TIMEOUT 600
                         LABELS  "regtest")
endfunction()

# ---------------------------------------------------------------------------
# Helper: add a test that needs a Python setup step first
# ---------------------------------------------------------------------------
function(add_python_setup_regression_test TEST_NAME TEST_DIR SETUP_SCRIPT)
    set(TEST_INPUT   "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_DIR}/inputs.i")
    set(TEST_WORK_DIR "${CMAKE_CURRENT_BINARY_DIR}/${TEST_DIR}")
    file(MAKE_DIRECTORY ${TEST_WORK_DIR})

    file(GLOB _DATA_FILES
         "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_DIR}/*.csv"
         "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_DIR}/*.lcp"
         "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_DIR}/*.adj"
         "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_DIR}/*.fmd"
         "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_DIR}/*.fms")
    foreach(_df ${_DATA_FILES})
        file(COPY "${_df}" DESTINATION ${TEST_WORK_DIR})
    endforeach()

    if(Python3_FOUND)
        set(SETUP_TEST_NAME "${TEST_NAME}_setup")
        add_test(NAME ${SETUP_TEST_NAME}
                 COMMAND ${Python3_EXECUTABLE}
                         ${CMAKE_CURRENT_SOURCE_DIR}/${TEST_DIR}/${SETUP_SCRIPT}
                 WORKING_DIRECTORY ${TEST_WORK_DIR})
        set_tests_properties(${SETUP_TEST_NAME} PROPERTIES
                             TIMEOUT 300
                             LABELS  "regtest"
                             FIXTURES_SETUP "${TEST_NAME}_fixture")
        add_test(NAME ${TEST_NAME}
                 COMMAND ${LEVELSET_EXE} ${TEST_INPUT}
                 WORKING_DIRECTORY ${TEST_WORK_DIR})
        set_tests_properties(${TEST_NAME} PROPERTIES
                             TIMEOUT 600
                             LABELS  "regtest"
                             FIXTURES_REQUIRED "${TEST_NAME}_fixture")
    else()
        message(STATUS "Python3 not found – skipping ${TEST_NAME} (needs Python for setup)")
    endif()
endfunction()

# ===========================================================================
# surface_spread/ – core level-set and Rothermel spread models
# ===========================================================================
add_regression_test(basic_levelset           surface_spread/basic_levelset)
add_regression_test(farsite_ellipse          surface_spread/farsite_ellipse)
add_regression_test(rothermel_fuel           surface_spread/rothermel_fuel)
add_regression_test(anderson_lw              surface_spread/anderson_lw)
add_regression_test(catchpole_demestre       surface_spread/catchpole_demestre)
add_regression_test(wilson_spread            surface_spread/wilson_spread)
add_regression_test(alexander_lemniscate     surface_spread/alexander_lemniscate)
add_regression_test(ellipse_sdf              surface_spread/ellipse_sdf)
add_regression_test(reinitialization         surface_spread/reinitialization)
add_regression_test(fbp_o1a_grassfire        surface_spread/fbp_o1a_grassfire)
add_regression_test(fbp_o1b_grassfire        surface_spread/fbp_o1b_grassfire)
add_regression_test(fbp_s1_slash             surface_spread/fbp_s1_slash)
add_regression_test(fbp_s2_slash             surface_spread/fbp_s2_slash)
add_regression_test(fbp_s3_slash             surface_spread/fbp_s3_slash)
add_regression_test(lautenberger_spread      surface_spread/lautenberger_spread)

# ===========================================================================
# crown_fire/ – crown fire initiation, active crown, FMC seasonal schedule
# ===========================================================================
add_regression_test(crown_initiation         crown_fire/crown_initiation)
add_regression_test(cruz_crown_continental_us crown_fire/cruz_crown_continental_us)
add_regression_test(fmc_seasonal             crown_fire/fmc_seasonal)
add_regression_test(rothermel1991_crown      crown_fire/rothermel1991_crown)

# ===========================================================================
# spotting/ – firebrand spotting models
# ===========================================================================
add_regression_test(firebrand_spotting       spotting/firebrand_spotting)
add_regression_test(albini_spotting          spotting/albini_spotting)
if(Python3_FOUND)
    add_python_setup_regression_test(albini_spotting_3d_wind
                                     spotting/albini_spotting_3d_wind
                                     generate_plt_wind.py)
endif()

# ===========================================================================
# terrain/ – terrain and wind-terrain interaction
# ===========================================================================
add_regression_test(terrain_wind             terrain/terrain_wind)
add_regression_test(balbi_viegas_heatflux    terrain/balbi_viegas_heatflux)
add_regression_test(windninja_ridge_canyon   terrain/windninja_ridge_canyon)
if(Python3_FOUND)
    add_python_setup_regression_test(solar_horizon_shading
                                     terrain/solar_horizon_shading
                                     create_terrain.py)
endif()

# ===========================================================================
# moisture/ – fuel moisture models (FMD schedule, diurnal, precipitation)
# ===========================================================================
if(Python3_FOUND)
    add_python_setup_regression_test(fmd_moisture moisture/fmd_moisture create_fmd.py)
endif()
add_regression_test(cheney_gould_grassfire   moisture/cheney_gould_grassfire)
add_regression_test(precip_wetting           moisture/precip_wetting)
add_regression_test(spatial_moisture_output  moisture/spatial_moisture_output)

# ===========================================================================
# fuel/ – fuel models and adjustment files
# ===========================================================================
if(Python3_FOUND)
    add_python_setup_regression_test(fuel_adj_file fuel/fuel_adj_file create_adj.py)
endif()

# ===========================================================================
# ignition/ – ignition types (barriers, polygon, polyline, perimeter output)
# ===========================================================================
if(Python3_FOUND)
    add_python_setup_regression_test(barrier_polygons ignition/barrier_polygons create_barriers.py)
endif()
add_regression_test(fire_perimeter_output    ignition/fire_perimeter_output)
add_regression_test(polygon_ignition         ignition/polygon_ignition)
add_regression_test(polyline_ignition        ignition/polyline_ignition)
if(Python3_FOUND)
    set(_perim_work "${CMAKE_CURRENT_BINARY_DIR}/ignition/fire_perimeter_output")
    add_test(NAME fire_perimeter_output_validate
             COMMAND ${Python3_EXECUTABLE}
                     ${CMAKE_CURRENT_SOURCE_DIR}/ignition/fire_perimeter_output/validate.py
             WORKING_DIRECTORY ${_perim_work})
    set_tests_properties(fire_perimeter_output_validate PROPERTIES
                         TIMEOUT 60
                         LABELS  "regtest"
                         DEPENDS fire_perimeter_output)
endif()

# ===========================================================================
# wind/ – wind models (time-dependent, turbulent, direction schedule)
# ===========================================================================
add_regression_test(time_dependent_wind      wind/time_dependent_wind)
add_regression_test(turb_wind                wind/turb_wind)
add_regression_test(wind_dir_schedule        wind/wind_dir_schedule)

# ===========================================================================
# diagnostics/ – fire ecology, emissions, Scott-Reinhardt TI/CI
# ===========================================================================
add_regression_test(scott_reinhardt_indices  diagnostics/scott_reinhardt_indices)
add_regression_test(scott_reinhardt_full_ti_ci diagnostics/scott_reinhardt_full_ti_ci)

# ===========================================================================
# misc/ – EB, MTT, benchmarks, LANDFIRE
# ===========================================================================
add_regression_test(sphere_2d                misc/3d_sphere)
add_regression_test(eb_implicit              misc/eb_implicit)
add_regression_test(mtt_propagation          misc/mtt_propagation)
add_regression_test(bulk_fuel_consumption    misc/bulk_fuel_consumption)

if(Python3_FOUND)
    add_python_setup_regression_test(landfire_farsite misc/landfire_farsite create_landscape.py)
    add_python_setup_regression_test(nonburnable_mask misc/nonburnable_mask create_landscape.py)
endif()

# ---------------------------------------------------------------------------
# Timing benchmark
# ---------------------------------------------------------------------------
if(Python3_FOUND)
    set(TIMING_WORK_DIR "${CMAKE_CURRENT_BINARY_DIR}/misc/timing_benchmark")
    file(MAKE_DIRECTORY ${TIMING_WORK_DIR})
    set(_TIMING_RES "32;64;128")
    string(REPLACE ";" " " _TIMING_RES_ARGS "${_TIMING_RES}")
    separate_arguments(_TIMING_RES_LIST NATIVE_COMMAND "${_TIMING_RES_ARGS}")
    add_test(NAME timing_benchmark
             COMMAND ${Python3_EXECUTABLE}
                     ${CMAKE_CURRENT_SOURCE_DIR}/misc/timing_benchmark/run_benchmark.py
                     --exe ${LEVELSET_EXE} --dim 2
                     --nsteps 20 --resolutions ${_TIMING_RES_LIST}
                     --timeout 300 --output timing_results.csv --skip-scaling
             WORKING_DIRECTORY ${TIMING_WORK_DIR})
    set_tests_properties(timing_benchmark PROPERTIES
                         TIMEOUT 1200
                         LABELS  "regtest;benchmark"
                         RUN_SERIAL TRUE)
endif()

# ---------------------------------------------------------------------------
# Custom target
# ---------------------------------------------------------------------------
add_custom_target(regtest
                 COMMAND ${CMAKE_CTEST_COMMAND} -L regtest --output-on-failure
                 COMMENT "Running all regression tests"
                 USES_TERMINAL)

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
message(STATUS "")
message(STATUS "Regression tests (sub-folder organisation):")
message(STATUS "  surface_spread/ : basic_levelset, farsite_ellipse, rothermel_fuel,")
message(STATUS "                    anderson_lw, catchpole_demestre, wilson_spread,")
message(STATUS "                    alexander_lemniscate, ellipse_sdf, reinitialization,")
message(STATUS "                    fbp_o1a_grassfire, fbp_o1b_grassfire,")
message(STATUS "                    fbp_s1_slash, fbp_s2_slash, fbp_s3_slash, lautenberger_spread")
message(STATUS "  crown_fire/     : crown_initiation, cruz_crown_continental_us,")
message(STATUS "                    fmc_seasonal, rothermel1991_crown")
message(STATUS "  spotting/       : firebrand_spotting, albini_spotting,")
message(STATUS "                    albini_spotting_3d_wind")
message(STATUS "  terrain/        : terrain_wind, balbi_viegas_heatflux,")
message(STATUS "                    windninja_ridge_canyon, solar_horizon_shading")
message(STATUS "  moisture/       : fmd_moisture, cheney_gould_grassfire,")
message(STATUS "                    precip_wetting, spatial_moisture_output")
message(STATUS "  fuel/           : fuel_adj_file")
message(STATUS "  ignition/       : barrier_polygons, fire_perimeter_output,")
message(STATUS "                    polygon_ignition, polyline_ignition")
message(STATUS "  wind/           : time_dependent_wind, turb_wind, wind_dir_schedule")
message(STATUS "  diagnostics/    : scott_reinhardt_indices, scott_reinhardt_full_ti_ci")
message(STATUS "  misc/           : sphere_2d, eb_implicit, mtt_propagation,")
message(STATUS "                    bulk_fuel_consumption, timing_benchmark,")
message(STATUS "                    landfire_farsite, nonburnable_mask")
message(STATUS "")
message(STATUS "All tests use 2-D (AMReX_SPACEDIM=2) builds only.")
message(STATUS "")
message(STATUS "Run tests with: ctest -L regtest  or  make regtest")
