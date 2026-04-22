cmake_minimum_required(VERSION 3.20)
project(LevelSet LANGUAGES CXX C)

# Require C++17 (AMReX supports C++17+)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find an installed AMReX (built with DIM=2 or DIM=3 as you need)
# You may need to pass -DAMReX_DIR=/path/to/amrex/install/lib/cmake/amrex
#find_package(AMReX CONFIG REQUIRED)

find_package(AMReX REQUIRED)

add_executable(levelset src/main.cpp)
target_link_libraries(levelset PRIVATE AMReX::amrex)

# Recommended warnings
if (MSVC)
  target_compile_options(levelset PRIVATE /W4)
else()
  target_compile_options(levelset PRIVATE -Wall -Wextra -Wpedantic)
endif()
