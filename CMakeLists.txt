cmake_minimum_required(VERSION 3.20)

set(AMREX_SUBMODULE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external/amrex")
option(LEVELSET_USE_VENDORED_AMREX "Prefer building the vendored AMReX submodule" ON)
option(LEVELSET_DIM_2D "Build for 2D instead of 3D" OFF)
option(LEVELSET_ENABLE_EB "Enable AMReX Embedded Boundary support" OFF)
option(LEVELSET_BUILD_DOCS "Build Sphinx documentation" OFF)
option(LEVELSET_ENABLE_MPI "Enable MPI parallelism" OFF)
# The terrain-following mass-consistent wind solver requires the AMReX linear
# solver library (MLMG / MLABecLaplacian) and is inherently 3-D.  It is
# therefore only built when LEVELSET_DIM_2D=OFF.
if(WIN32)
  set(_wind_solver_default OFF)
else()
  set(_wind_solver_default ON)
endif()
option(LEVELSET_BUILD_WIND_SOLVER
       "Build the stand-alone terrain-following mass-consistent wind solver"
       ${_wind_solver_default})
# Internal option used by the ExternalProject sub-build for wind_solver when
# LEVELSET_DIM_2D=ON.  End users should not set this directly.
option(LEVELSET_BUILD_LEVELSET "Build the levelset fire-spread solver" ON)
set(LEVELSET_GPU_BACKEND "NONE" CACHE STRING
    "GPU backend for AMReX: NONE (CPU-only), CUDA, HIP, or SYCL")
set_property(CACHE LEVELSET_GPU_BACKEND PROPERTY STRINGS NONE CUDA HIP SYCL)

project(LevelSet LANGUAGES CXX C)

# Require C++17 (AMReX supports C++17+)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Enable additional language(s) required for the selected GPU backend
if(LEVELSET_GPU_BACKEND STREQUAL "CUDA")
  enable_language(CUDA)
  set(CMAKE_CUDA_STANDARD 17)
  set(CMAKE_CUDA_STANDARD_REQUIRED ON)
  message(STATUS "GPU backend: CUDA")
elseif(LEVELSET_GPU_BACKEND STREQUAL "HIP")
  enable_language(HIP)
  message(STATUS "GPU backend: HIP")
elseif(LEVELSET_GPU_BACKEND STREQUAL "SYCL")
  # SYCL uses a standard C++ compiler (e.g., Intel icpx); no extra language needed
  message(STATUS "GPU backend: SYCL (ensure icpx or DPC++ compiler is active)")
elseif(NOT LEVELSET_GPU_BACKEND STREQUAL "NONE")
  message(WARNING "Unknown GPU backend '${LEVELSET_GPU_BACKEND}'. Supported values: NONE CUDA HIP SYCL")
endif()

# Enable testing for regression tests
enable_testing()

if(LEVELSET_ENABLE_MPI)
  find_package(MPI REQUIRED)
  message(STATUS "MPI parallelism enabled (${MPI_CXX_COMPILER})")
endif()

if(LEVELSET_USE_VENDORED_AMREX AND EXISTS "${AMREX_SUBMODULE_DIR}/CMakeLists.txt")
  message(STATUS "Using vendored AMReX from ${AMREX_SUBMODULE_DIR}")

  if(LEVELSET_DIM_2D)
    set(AMReX_SPACEDIM "2" CACHE STRING "AMReX spatial dimensions for wildfire_levelset" FORCE)
    message(STATUS "Building for 2D (AMReX_SPACEDIM=2)")
  else()
    set(AMReX_SPACEDIM "3" CACHE STRING "AMReX spatial dimensions for wildfire_levelset" FORCE)
    message(STATUS "Building for 3D (AMReX_SPACEDIM=3)")
  endif()
  if(LEVELSET_ENABLE_MPI)
    set(AMReX_MPI ON CACHE BOOL "Build AMReX with MPI" FORCE)
  else()
    set(AMReX_MPI OFF CACHE BOOL "Build AMReX without MPI" FORCE)
  endif()
  set(AMReX_OMP OFF CACHE BOOL "Build AMReX without OpenMP" FORCE)
  set(AMReX_GPU_BACKEND "${LEVELSET_GPU_BACKEND}" CACHE STRING "AMReX GPU backend" FORCE)
  set(AMReX_FORTRAN OFF CACHE BOOL "Build AMReX without Fortran" FORCE)
  set(AMReX_AMRLEVEL OFF CACHE BOOL "Build AMReX without AmrLevel support" FORCE)
  # The wind_solver requires AMReX linear solvers (MLMG).
  # In a 3-D build, enable them now.  In a 2-D build, the wind_solver gets
  # its own 3-D AMReX via ExternalProject, so we keep this build 2-D.
  if(LEVELSET_BUILD_WIND_SOLVER AND NOT LEVELSET_DIM_2D)
    set(AMReX_LINEAR_SOLVERS ON  CACHE BOOL "Build AMReX with linear solvers (wind_solver)" FORCE)
  else()
    set(AMReX_LINEAR_SOLVERS OFF CACHE BOOL "Build AMReX without linear solvers" FORCE)
  endif()
  set(AMReX_PARTICLES OFF CACHE BOOL "Build AMReX without particle support" FORCE)
  set(AMReX_INSTALL OFF CACHE BOOL "Disable AMReX install rules when used as a subproject" FORCE)
  if(LEVELSET_ENABLE_EB)
    set(AMReX_EB ON CACHE BOOL "Enable AMReX Embedded Boundary support" FORCE)
    message(STATUS "Embedded Boundary support enabled")
  else()
    set(AMReX_EB OFF CACHE BOOL "Disable AMReX Embedded Boundary support" FORCE)
  endif()

  add_subdirectory("${AMREX_SUBMODULE_DIR}" "${CMAKE_CURRENT_BINARY_DIR}/external/amrex")
else()
  if(LEVELSET_USE_VENDORED_AMREX)
    message(STATUS "Vendored AMReX submodule not found at ${AMREX_SUBMODULE_DIR}; trying find_package(AMReX)")
  endif()

  find_package(AMReX QUIET CONFIG)

  if(AMReX_FOUND)
    message(STATUS "Using externally installed AMReX from ${AMReX_DIR}")
  else()
    message(FATAL_ERROR
      "AMReX was not found.\n"
      "Expected vendored submodule at: ${AMREX_SUBMODULE_DIR}\n"
      "Run `git submodule update --init --recursive` from the repository root, or configure with an installed AMReX via -DAMReX_DIR=/path/to/lib/cmake/AMReX."
    )
  endif()
endif()

if(LEVELSET_BUILD_LEVELSET)
add_executable(levelset 
               src/main.cpp
               src/parse_inputs.cpp                
              )
# On Windows with CUDA, MSVC is the default host compiler and does not
# understand CUDA device-side keywords (__device__, __host__, etc.) that
# appear in AMReX headers.  Tell CMake to route these .cpp files through
# nvcc (LANGUAGE CUDA) so they are compiled correctly on all platforms.
if(LEVELSET_GPU_BACKEND STREQUAL "CUDA")
  set_source_files_properties(
    src/main.cpp
    src/parse_inputs.cpp
    PROPERTIES LANGUAGE CUDA
  )
endif()
target_link_libraries(levelset PRIVATE AMReX::amrex)
if(LEVELSET_ENABLE_MPI)
  target_link_libraries(levelset PRIVATE MPI::MPI_CXX)
endif()

# Recommended warnings
# When source files are compiled as CUDA (language=CUDA), CMake uses nvcc as the
# compiler driver.  MSVC-style flags (e.g. /W4) must be wrapped with
# -Xcompiler= so nvcc forwards them to the host compiler; passing them bare
# causes nvcc to treat them as input files and abort with
# "A single input file is required for a non-link phase".
if (MSVC)
  target_compile_options(levelset PRIVATE
    $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=/W4>
    $<$<NOT:$<COMPILE_LANGUAGE:CUDA>>:/W4>
  )
else()
  target_compile_options(levelset PRIVATE -Wall -Wextra -Wpedantic)
endif()
endif() # LEVELSET_BUILD_LEVELSET

# ---------------------------------------------------------------------------
# Optional: terrain-following mass-consistent wind solver (always 3-D)
#
# The wind_solver is inherently 3-D (it uses the full 3-D AMReX linear
# solver stack).  When the user requests LEVELSET_DIM_2D=ON the main AMReX
# build uses SPACEDIM=2, so wind_solver cannot share that library.  In that
# case we kick off a separate CMake sub-project (via ExternalProject_Add)
# configured with SPACEDIM=3 and LEVELSET_BUILD_LEVELSET=OFF so that only
# wind_solver and its AMReX-3D dependency are compiled.  The resulting
# binary is copied to CMAKE_BINARY_DIR alongside the 2-D levelset solver.
# ---------------------------------------------------------------------------
if(LEVELSET_BUILD_WIND_SOLVER)
  if(NOT LEVELSET_DIM_2D)
    # ---- Normal 3-D build: wind_solver is a direct target ----
    add_executable(wind_solver src/wind_solver.cpp)
    # Same as levelset: route through nvcc on CUDA builds so MSVC does not
    # choke on CUDA device-side keywords in AMReX headers.
    if(LEVELSET_GPU_BACKEND STREQUAL "CUDA")
      set_source_files_properties(
        src/wind_solver.cpp
        PROPERTIES LANGUAGE CUDA
      )
    endif()
    target_link_libraries(wind_solver PRIVATE AMReX::amrex)
    if(LEVELSET_ENABLE_MPI)
      target_link_libraries(wind_solver PRIVATE MPI::MPI_CXX)
    endif()
    if(MSVC)
      target_compile_options(wind_solver PRIVATE
        $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=/W4>
        $<$<NOT:$<COMPILE_LANGUAGE:CUDA>>:/W4>
      )
    else()
      target_compile_options(wind_solver PRIVATE -Wall -Wextra -Wpedantic)
    endif()
    message(STATUS "wind_solver: enabled (3-D, same AMReX build)")
  else()
    # ---- 2-D + 3-D co-build: wind_solver in an ExternalProject sub-build ----
    include(ExternalProject)
    ExternalProject_Add(wind_solver_3d_subbuild
      SOURCE_DIR  "${CMAKE_CURRENT_SOURCE_DIR}"
      BINARY_DIR  "${CMAKE_CURRENT_BINARY_DIR}/wind_solver_3d"
      CMAKE_ARGS
        "-DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}"
        "-DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}"
        "-DLEVELSET_USE_VENDORED_AMREX:BOOL=ON"
        "-DLEVELSET_DIM_2D:BOOL=OFF"          # force 3-D for wind_solver
        "-DLEVELSET_BUILD_LEVELSET:BOOL=OFF"  # skip the 2-D levelset target
        "-DLEVELSET_BUILD_WIND_SOLVER:BOOL=ON"
        "-DLEVELSET_ENABLE_MPI:BOOL=${LEVELSET_ENABLE_MPI}"
        "-DLEVELSET_ENABLE_EB:BOOL=OFF"
        "-DLEVELSET_GPU_BACKEND:STRING=${LEVELSET_GPU_BACKEND}"
      BUILD_COMMAND
        "${CMAKE_COMMAND}" --build "<BINARY_DIR>" --target wind_solver
      INSTALL_COMMAND
        "${CMAKE_COMMAND}" -E copy_if_different
          "<BINARY_DIR>/wind_solver"
          "${CMAKE_CURRENT_BINARY_DIR}/wind_solver"
      BUILD_BYPRODUCTS "${CMAKE_CURRENT_BINARY_DIR}/wind_solver"
    )
    message(STATUS "wind_solver: will be built as a 3-D ExternalProject sub-build")
    message(STATUS "  (2-D levelset + 3-D wind_solver co-build)")
    message(STATUS "  Sub-build directory: ${CMAKE_CURRENT_BINARY_DIR}/wind_solver_3d")
  endif()
else()
  message(STATUS "wind_solver: disabled (LEVELSET_BUILD_WIND_SOLVER=OFF)")
endif()

# Add regression test subdirectory
add_subdirectory(regtest)

# Add documentation build target
if(LEVELSET_BUILD_DOCS)
  find_package(Python3 COMPONENTS Interpreter)
  if(Python3_FOUND)
    find_program(SPHINX_EXECUTABLE
                 NAMES sphinx-build
                 DOC "Path to sphinx-build executable")
    
    if(SPHINX_EXECUTABLE)
      set(SPHINX_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/docs)
      set(SPHINX_BUILD ${CMAKE_CURRENT_BINARY_DIR}/docs)
      
      add_custom_target(docs
        COMMAND ${SPHINX_EXECUTABLE} -b html
                ${SPHINX_SOURCE} ${SPHINX_BUILD}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating documentation with Sphinx"
        VERBATIM)
      
      message(STATUS "Documentation build target 'docs' configured")
      message(STATUS "  Build with: cmake --build build --target docs")
      message(STATUS "  Output: ${SPHINX_BUILD}/index.html")
    else()
      message(WARNING "sphinx-build not found. Install with: pip install -r docs/requirements.txt")
    endif()
  else()
    message(WARNING "Python3 not found. Documentation build disabled.")
  endif()
endif()
