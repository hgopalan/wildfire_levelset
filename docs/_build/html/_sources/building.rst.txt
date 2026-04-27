Building the Code
=================

This section describes how to build the wildfire level-set solver.

Prerequisites
-------------

Required Software
^^^^^^^^^^^^^^^^^

* **C++17 compiler**: GCC 7+, Clang 5+, or MSVC 2017+
* **CMake**: Version 3.20 or higher
* **Git**: For cloning the repository and submodules

Optional Software
^^^^^^^^^^^^^^^^^

* **Python 3**: For visualization and post-processing scripts
* **VisIt, ParaView, or yt**: For visualizing AMReX plotfiles

Cloning the Repository
-----------------------

Clone with Submodules
^^^^^^^^^^^^^^^^^^^^^

The project uses AMReX as a git submodule. Clone the repository with::

    git clone --recurse-submodules https://github.com/hgopalan/wildfire_levelset.git
    cd wildfire_levelset

If you already cloned without submodules, initialize them::

    git submodule update --init --recursive

Basic Build (3D, CPU-only)
---------------------------

Default Configuration
^^^^^^^^^^^^^^^^^^^^^

The default build configuration is:

* 3D spatial dimensions
* CPU-only (no MPI, OpenMP, or GPU)
* No embedded boundary support

Build with::

    cmake -S . -B build
    cmake --build build -j

The executable will be created at ``build/levelset``.

Build Options
-------------

2D Build
^^^^^^^^

To build for 2D instead of 3D::

    cmake -S . -B build -DLEVELSET_DIM_2D=ON
    cmake --build build -j

.. note::
   Grid tagging (adaptive mesh refinement) is automatically disabled in 2D builds.

Embedded Boundary Support
^^^^^^^^^^^^^^^^^^^^^^^^^^

To enable AMReX Embedded Boundary (EB) support for complex geometries::

    cmake -S . -B build -DLEVELSET_ENABLE_EB=ON
    cmake --build build -j

EB support allows the solver to handle complex geometries using implicit function representations, enabling features like terrain following coordinates and irregular domain boundaries.

Combine 2D and EB::

    cmake -S . -B build -DLEVELSET_DIM_2D=ON -DLEVELSET_ENABLE_EB=ON
    cmake --build build -j

Using External AMReX
^^^^^^^^^^^^^^^^^^^^

To use an externally installed AMReX instead of the vendored submodule::

    cmake -S . -B build -DLEVELSET_USE_VENDORED_AMREX=OFF -DAMReX_DIR=/path/to/lib/cmake/AMReX
    cmake --build build -j

Build Types
^^^^^^^^^^^

Debug build::

    cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
    cmake --build build -j

Release build (optimized)::

    cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
    cmake --build build -j

RelWithDebInfo (optimized with debug symbols)::

    cmake -S . -B build -DCMAKE_BUILD_TYPE=RelWithDebInfo
    cmake --build build -j

Advanced Build Options
-----------------------

Compiler Selection
^^^^^^^^^^^^^^^^^^

Specify a different compiler::

    cmake -S . -B build -DCMAKE_CXX_COMPILER=g++-11 -DCMAKE_C_COMPILER=gcc-11
    cmake --build build -j

or::

    cmake -S . -B build -DCMAKE_CXX_COMPILER=clang++-14 -DCMAKE_C_COMPILER=clang-14
    cmake --build build -j

Parallel Build
^^^^^^^^^^^^^^

Use multiple cores for building::

    cmake --build build -j8  # Use 8 cores

or::

    cmake --build build -j$(nproc)  # Use all available cores

Verbose Build Output
^^^^^^^^^^^^^^^^^^^^

See the actual compiler commands::

    cmake --build build --verbose

or::

    cmake -S . -B build -DCMAKE_VERBOSE_MAKEFILE=ON
    cmake --build build

Clean Build
^^^^^^^^^^^

Remove build directory and rebuild::

    rm -rf build
    cmake -S . -B build
    cmake --build build -j

Building Documentation
----------------------

The Sphinx documentation can be built locally. First, install the required Python packages::

    pip install -r docs/requirements.txt

Build HTML documentation::

    cd docs
    make html

The generated HTML files will be in ``docs/_build/html/``. Open ``docs/_build/html/index.html`` in a web browser.

Build PDF documentation::

    cd docs
    make latexpdf

The generated PDF will be in ``docs/_build/latex/``.

Clean documentation build::

    cd docs
    make clean

Running Tests
-------------

Build Tests
^^^^^^^^^^^

The regression tests are built automatically when you build the main project. They are defined in ``regtest/CMakeLists.txt``.

Run All Tests
^^^^^^^^^^^^^

From the build directory::

    cd build
    ctest

or::

    ctest --output-on-failure  # Show output for failed tests

Run Specific Test
^^^^^^^^^^^^^^^^^

Run a single test by name::

    ctest -R basic_levelset

Run tests matching a pattern::

    ctest -R farsite

Verbose Test Output
^^^^^^^^^^^^^^^^^^^

See all test output::

    ctest -V

or::

    ctest --verbose

List Available Tests
^^^^^^^^^^^^^^^^^^^^

See all available tests::

    ctest -N

Installation
------------

Install to System Directory
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Install the executable and headers::

    cmake --build build --target install

By default, this installs to ``/usr/local``. To change the installation prefix::

    cmake -S . -B build -DCMAKE_INSTALL_PREFIX=/custom/path
    cmake --build build --target install

.. note::
   Currently, only the ``levelset`` executable is installed. Library installation and header installation are not configured.

Troubleshooting
---------------

AMReX Submodule Not Found
^^^^^^^^^^^^^^^^^^^^^^^^^^

If you see an error about missing AMReX submodule::

    CMake Error: AMReX was not found.
    Expected vendored submodule at: external/amrex

Solution: Initialize the submodule::

    git submodule update --init --recursive

Compiler Not Found
^^^^^^^^^^^^^^^^^^

If CMake cannot find your compiler, specify it explicitly::

    cmake -S . -B build -DCMAKE_CXX_COMPILER=/path/to/g++ -DCMAKE_C_COMPILER=/path/to/gcc

CMake Version Too Old
^^^^^^^^^^^^^^^^^^^^^

If your CMake version is too old::

    CMake Error: CMake 3.20 or higher is required.

Solution: Install a newer CMake from https://cmake.org/download/ or your package manager.

C++17 Not Supported
^^^^^^^^^^^^^^^^^^^

If your compiler doesn't support C++17::

    error: unrecognized command line option '-std=c++17'

Solution: Upgrade to a newer compiler (GCC 7+, Clang 5+, MSVC 2017+).

Platform-Specific Notes
-----------------------

Linux
^^^^^

Most Linux distributions have the required tools in their package repositories::

    # Ubuntu/Debian
    sudo apt-get install build-essential cmake git

    # Fedora/RHEL/CentOS
    sudo dnf install gcc-c++ cmake git

macOS
^^^^^

Install Xcode Command Line Tools and CMake::

    xcode-select --install
    brew install cmake

Windows
^^^^^^^

Use Visual Studio 2017 or later with CMake support, or install CMake separately and use MinGW or MSYS2.

Example with MinGW::

    cmake -S . -B build -G "MinGW Makefiles"
    cmake --build build
