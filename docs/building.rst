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

GPU Builds
----------

CUDA (Linux)
^^^^^^^^^^^^

Install the CUDA toolkit (12.6 or later recommended) before configuring::

    cmake -S . -B build \
      -DLEVELSET_GPU_BACKEND=CUDA \
      -DAMReX_CUDA_ARCH="8.0"
    cmake --build build --parallel

.. note::
   The mass-consistent wind solver has been deprecated and is no longer built.
   The source is retained in ``src/deprecated/wind_solver.cpp`` for reference.

**Known issues encountered during development**

* **``identifier undefined in device code``** — ``static constexpr`` variables at file scope in
  headers have internal linkage and are invisible to NVCC's separate device-code compilation
  pass.  Fixed by changing them to ``inline constexpr`` (C++17 inline variables), which have
  external linkage and can be inlined directly into device code.

* **``-Wpedantic`` triggers NVCC diagnostics** — NVCC emits ``"style of line directive is a
  GCC extension"`` and ``"ISO C++ forbids taking address of function '::main'"`` when
  ``-Wpedantic`` is active.  Fixed by replacing ``-Wpedantic`` with ``-Wno-pedantic`` for
  CUDA builds in ``CMakeLists.txt``.

* **GPU-unsafe lambdas** — AMReX ``ParallelFor`` device lambdas must capture by value
  (``[=]``), not by reference (``[&]``), because host stack pointers are inaccessible from GPU
  threads.  Additionally, ``std::max({a, b, c})`` (initializer-list overload) is not available
  in CUDA device code; use ``amrex::max(amrex::max(a, b), c)`` instead.

* **``std::min``/``std::max`` in device lambdas** — The ``<algorithm>`` overloads are not
  portable across CUDA/HIP/SYCL backends; replace them with ``amrex::min``/``amrex::max``
  inside ``ParallelFor`` kernels.

CUDA (Windows)
^^^^^^^^^^^^^^

Windows CUDA builds require Visual Studio 2019 or later (MSVC) with the CUDA toolkit installed.
Use the ``build_windows_cuda`` CI job as a reference.

.. note::
   The mass-consistent wind solver has been deprecated and is no longer built.

**Known issues encountered during development**

* **MSVC-style ``/W4`` flag passed to NVCC** — When ``.cpp`` files are compiled through
  ``nvcc`` (``LANGUAGE CUDA``), bare MSVC flags like ``/W4`` land as separate NVCC arguments.
  NVCC misinterprets ``/W4`` as an input file and aborts with ``"A single input file is
  required for a non-link phase when an outputfile is specified"``.  Fixed by wrapping the flag
  in a CMake generator expression::

      $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=/W4>
      $<$<NOT:$<COMPILE_LANGUAGE:CUDA>>:/W4>

* **``static`` functions with device lambdas rejected on Windows CUDA** — NVCC on Windows
  rejects extended ``__device__`` / ``__host__ __device__`` lambdas inside functions with
  internal linkage (``static void`` in headers).  Fixed by changing all such free functions in
  headers from ``static void`` to ``inline void``; ``inline`` provides external linkage while
  remaining safe to define in headers included across multiple translation units.

* **Project ``.cpp`` files not routed through NVCC** — On Windows, ``src/main.cpp`` and
  ``src/parse_inputs.cpp`` were compiled by MSVC instead of NVCC.  AMReX headers expose CUDA
  device-side keywords that MSVC cannot parse.  Fixed by marking the affected source files
  with ``PROPERTIES LANGUAGE CUDA`` in ``CMakeLists.txt`` when the CUDA backend is enabled.

* **Missing CUDA sub-packages** — The Windows CUDA CI job failed with a ``curand`` error
  because ``curand_dev`` was absent from the ``Jimver/cuda-toolkit`` sub-package list.  The
  required sub-packages (matching the AMReX reference recipe) are::

      ["nvcc", "cudart", "cuda_profiler_api", "cufft_dev", "cusparse_dev", "curand_dev"]

* **AMReX headers not ready when project sources compile** — Ninja's parallel scheduler can
  start compiling ``src/*.cpp`` before AMReX's generated headers (e.g., ``AMReX_Config.H``)
  have been written to the build tree.  This manifests as the same misleading ``"A single
  input file is required"`` NVCC error.  Fixed by adding an explicit ``add_dependencies`` or
  pre-build step that compiles the ``amrex_3d`` target before the project sources.

* **Wrong Ninja target name for pre-build step** — The AMReX pre-build step must target
  ``amrex_3d`` (the name CMake generates for the 3-D AMReX library), not the alias ``amrex``.
  Requesting the alias aborts the CI job before the main build starts.

HIP / AMD ROCm (Linux)
^^^^^^^^^^^^^^^^^^^^^^

Install ROCm 6.2 (Ubuntu 22.04 / Jammy) before configuring.  The recommended package set
mirrors the one used by AMReX, ERF, and AMR-Wind CI::

    sudo apt-get install -y \
      rocm-dev rocrand-dev rocprim-dev hiprand-dev rocsparse-dev rocm-cmake

Configure using the ROCm Clang compilers (sourced from ``/etc/profile.d/rocm.sh``)::

    source /etc/profile.d/rocm.sh
    cmake -S . -B build \
      -G Ninja \
      -DCMAKE_C_COMPILER=$(which clang) \
      -DCMAKE_CXX_COMPILER=$(which clang++) \
      -DLEVELSET_GPU_BACKEND=HIP \
      -DAMReX_AMD_ARCH="gfx90a"
    cmake --build build --parallel

.. note::
   The mass-consistent wind solver has been deprecated and is no longer built.

**Known issues encountered during development**

* **Missing ROCm device library packages** — An earlier CI recipe installed only ``hipcc``,
  ``hip-dev``, and ``rocm-cmake``.  This was insufficient: ``rocm-dev`` is the critical
  addition because it pulls in ``rocm-device-libs`` — the GPU bitcode files that
  ``hipcc``/``clang++`` requires when targeting HIP.  Without it, the build fails with
  unresolved device-code symbols.

* **Duplicate YAML job keys silently discarded** — ``cmake_build.yml`` previously contained
  two ``build_hip:`` keys.  YAML does not allow duplicate mapping keys; GitHub Actions
  silently drops the first definition, so the stronger job (with full ROCm packages) was
  never running.  Fixed by removing the duplicate weaker definition.

* **Wrong compiler selected** — An earlier recipe passed ``/opt/rocm/bin/hipcc`` as the CXX
  compiler.  The correct approach (matching AMReX/ERF/AMR-Wind CI) is to source the ROCm
  environment and then use ``$(which clang++)`` / ``$(which clang)``.

SYCL / Intel oneAPI (Linux)
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Install the Intel DPC++/C++ compiler and MKL before configuring::

    # Add the Intel oneAPI apt repository first, then:
    sudo apt-get install -y \
      intel-oneapi-compiler-dpcpp-cpp \
      intel-oneapi-mkl-devel

Configure using the Intel compilers (sourced from ``/opt/intel/oneapi/setvars.sh``)::

    source /opt/intel/oneapi/setvars.sh
    cmake -S . -B build \
      -G Ninja \
      -DCMAKE_C_COMPILER=$(which icx) \
      -DCMAKE_CXX_COMPILER=$(which icpx) \
      -DLEVELSET_GPU_BACKEND=SYCL
    cmake --build build --parallel

**Known issues encountered during development**

* **Missing ``intel-oneapi-mkl-devel``** — AMReX includes ``<oneapi/mkl/rng.hpp>`` from
  ``AMReX_RandomEngine.H`` when building for SYCL.  If only the compiler package is installed
  (without ``intel-oneapi-mkl-devel``), the build fails with a missing header error.

* **Duplicate YAML job keys** — The same issue as HIP above: ``cmake_build.yml`` previously
  had two ``build_sycl:`` keys; GitHub Actions silently used only the second (lighter)
  definition, which was missing the MKL package.

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
