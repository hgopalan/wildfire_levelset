cmake_minimum_required(VERSION 3.20)

set(AMREX_SUBMODULE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external/amrex")
option(LEVELSET_USE_VENDORED_AMREX "Prefer building the vendored AMReX submodule" ON)
option(LEVELSET_DIM_2D "Build for 2D (only supported dimension)" ON)
option(LEVELSET_ENABLE_EB "Enable AMReX Embedded Boundary support" OFF)
option(LEVELSET_BUILD_DOCS "Build Sphinx documentation" OFF)
option(LEVELSET_ENABLE_MPI "Enable MPI parallelism" OFF)
option(LEVELSET_BUILD_LEVELSET "Build the levelset fire-spread solver" ON)
option(LEVELSET_BUILD_FCOMPARE "Build the AMReX fcompare plotfile comparison tool" OFF)
option(LEVELSET_BUILD_PYTHON_BINDINGS "Build Python bindings using pybind11" OFF)
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
    message(FATAL_ERROR
      "3D builds are not supported. Remove -DLEVELSET_DIM_2D=OFF "
      "(the default is ON). Only AMREX_SPACEDIM=2 is supported by this solver.")
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
  # For Python bindings, we need position-independent code
  if(LEVELSET_BUILD_PYTHON_BINDINGS)
    set(AMReX_PIC ON CACHE BOOL "Build AMReX with position-independent code" FORCE)
    set(CMAKE_POSITION_INDEPENDENT_CODE ON)
  endif()
  # Linear solvers are not required by the levelset fire-spread solver.
  set(AMReX_LINEAR_SOLVERS OFF CACHE BOOL "Build AMReX without linear solvers" FORCE)
  set(AMReX_PARTICLES OFF CACHE BOOL "Build AMReX without particle support" FORCE)
  set(AMReX_INSTALL OFF CACHE BOOL "Disable AMReX install rules when used as a subproject" FORCE)
  if(LEVELSET_BUILD_FCOMPARE)
    set(AMReX_PLOTFILE_TOOLS ON CACHE BOOL "Build AMReX plotfile tools (fcompare, fextract, …)" FORCE)
  else()
    set(AMReX_PLOTFILE_TOOLS OFF CACHE BOOL "Build AMReX plotfile tools (fcompare, fextract, …)" FORCE)
  endif()
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
# On CUDA builds with MSVC as host compiler, /W4 must be wrapped with
# -Xcompiler so nvcc forwards it to the host compiler rather than treating
# it as a bare nvcc argument (which causes a "single input file" fatal error).
if (MSVC)
  target_compile_options(levelset PRIVATE
    $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=/W4>
    $<$<COMPILE_LANGUAGE:CXX>:/W4>
  )
else()
  target_compile_options(levelset PRIVATE -Wall -Wextra -Wno-pedantic)
endif()
endif() # LEVELSET_BUILD_LEVELSET

# Add Python bindings subdirectory
if(LEVELSET_BUILD_PYTHON_BINDINGS)
  add_subdirectory(src/python)
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
