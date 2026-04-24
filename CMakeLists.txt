cmake_minimum_required(VERSION 3.20)
project(LevelSet LANGUAGES CXX C)

# Require C++17 (AMReX supports C++17+)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(AMREX_SUBMODULE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external/amrex")
option(LEVELSET_USE_VENDORED_AMREX "Prefer building the vendored AMReX submodule" ON)
option(LEVELSET_DIM_2D "Build for 2D instead of 3D" OFF)

if(LEVELSET_USE_VENDORED_AMREX AND EXISTS "${AMREX_SUBMODULE_DIR}/CMakeLists.txt")
  message(STATUS "Using vendored AMReX from ${AMREX_SUBMODULE_DIR}")

  if(LEVELSET_DIM_2D)
    set(AMReX_SPACEDIM "2" CACHE STRING "AMReX spatial dimensions for wildfire_levelset" FORCE)
    message(STATUS "Building for 2D (AMReX_SPACEDIM=2)")
  else()
    set(AMReX_SPACEDIM "3" CACHE STRING "AMReX spatial dimensions for wildfire_levelset" FORCE)
    message(STATUS "Building for 3D (AMReX_SPACEDIM=3)")
  endif()
  #set(AMReX_MPI OFF CACHE BOOL "Build AMReX without MPI" FORCE)
  set(AMReX_OMP OFF CACHE BOOL "Build AMReX without OpenMP" FORCE)
  set(AMReX_GPU_BACKEND NONE CACHE STRING "Build AMReX without a GPU backend" FORCE)
  set(AMReX_FORTRAN OFF CACHE BOOL "Build AMReX without Fortran" FORCE)
  set(AMReX_AMRLEVEL OFF CACHE BOOL "Build AMReX without AmrLevel support" FORCE)
  set(AMReX_LINEAR_SOLVERS OFF CACHE BOOL "Build AMReX without linear solvers" FORCE)
  set(AMReX_PARTICLES OFF CACHE BOOL "Build AMReX without particle support" FORCE)
  set(AMReX_INSTALL OFF CACHE BOOL "Disable AMReX install rules when used as a subproject" FORCE)

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

add_executable(levelset 
               src/main.cpp
               src/parse_inputs.cpp                
              )
target_link_libraries(levelset PRIVATE AMReX::amrex)

# Recommended warnings
if (MSVC)
  target_compile_options(levelset PRIVATE /W4)
else()
  target_compile_options(levelset PRIVATE -Wall -Wextra -Wpedantic)
endif()
