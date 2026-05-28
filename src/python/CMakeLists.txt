# Python bindings for wildfire_levelset using pybind11
# Provides interface to load wind data from Python (numpy arrays or pyAMReX MultiFab)

cmake_minimum_required(VERSION 3.20)

# Find Python3 and required components
find_package(Python3 REQUIRED COMPONENTS Interpreter Development NumPy)

message(STATUS "Python bindings enabled")
message(STATUS "  Python executable: ${Python3_EXECUTABLE}")
message(STATUS "  Python version: ${Python3_VERSION}")
message(STATUS "  Python include dirs: ${Python3_INCLUDE_DIRS}")
message(STATUS "  Python libraries: ${Python3_LIBRARIES}")

# Find or fetch pybind11
find_package(pybind11 CONFIG QUIET)
if(NOT pybind11_FOUND)
  message(STATUS "pybind11 not found, fetching from GitHub...")
  include(FetchContent)
  FetchContent_Declare(
    pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11.git
    GIT_TAG        v2.11.1
  )
  FetchContent_MakeAvailable(pybind11)
else()
  message(STATUS "Found pybind11: ${pybind11_DIR}")
endif()

# Create Python module
pybind11_add_module(pyWildfire MODULE
  pyWildfire.cpp
  ${CMAKE_SOURCE_DIR}/src/fire_solver_api.cpp
  ${CMAKE_SOURCE_DIR}/src/parse_inputs.cpp
)

# Link with AMReX
target_link_libraries(pyWildfire PRIVATE AMReX::amrex)

# Add include directories for headers
target_include_directories(pyWildfire PRIVATE
  ${CMAKE_SOURCE_DIR}/src
)

# On CUDA builds, compile bindings as CUDA code
if(LEVELSET_GPU_BACKEND STREQUAL "CUDA")
  set_source_files_properties(
    pyWildfire.cpp
    PROPERTIES LANGUAGE CUDA
  )
endif()

# Set output directory for the Python module
set_target_properties(pyWildfire PROPERTIES
  LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/python
)

# Installation support (optional)
install(TARGETS pyWildfire
  LIBRARY DESTINATION python
)
