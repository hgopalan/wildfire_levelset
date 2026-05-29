Performance and Benchmarking
=============================

This section describes the performance characteristics of the wildfire level-set solver
and provides guidance on benchmarking, profiling, and optimization.

.. contents::
   :local:
   :depth: 2

Overview
--------

The wildfire level-set solver is designed for high-performance wildfire simulations
with scalable algorithms and efficient data structures. Performance depends on:

* **Grid resolution**: Computational cost scales as O(N²) in 2D, O(N³) in 3D
* **Fire spread model**: Physics-based models (Balbi, Lautenberger) are more expensive than empirical models (Rothermel, FBP)
* **Propagation method**: FARSITE elliptical expansion is faster than level-set advection for coarse grids
* **Feature complexity**: Spotting, crown fire, radiation preheating add overhead

Timing Benchmarks
-----------------

The solver includes automated timing benchmarks in ``regtest/misc/timing_benchmark/``
to track performance across code versions and detect regressions.

Running Benchmarks
^^^^^^^^^^^^^^^^^^

From the build directory::

   cd build
   
   # Default benchmarks (level-set and FARSITE)
   python3 ../regtest/misc/timing_benchmark/run_benchmark.py \
       --exe ./levelset --dim 2 --nsteps 30
   
   # All fire spread models
   python3 ../regtest/misc/timing_benchmark/run_benchmark.py \
       --exe ./levelset --dim 2 --nsteps 30 \
       --scenarios levelset farsite balbi cruz_crown cheney_gould fbp_o1a lautenberger
   
   # Custom resolutions
   python3 ../regtest/misc/timing_benchmark/run_benchmark.py \
       --exe ./levelset --dim 2 --nsteps 50 \
       --resolutions 64 128 256 512

Benchmark Output
^^^^^^^^^^^^^^^^

The benchmark script generates:

* **Console output**: Real-time progress and summary table
* **timing_results.csv**: Detailed timing data with columns:

  - ``scenario``: Fire spread model name
  - ``n_cells``: Grid cells per dimension
  - ``total_cells``: Total grid cells (n_cells^dim)
  - ``wall_time_s``: Wall-clock execution time
  - ``nsteps``: Number of timesteps executed
  - ``steps_per_second``: Throughput metric
  - ``cells_per_step_per_s``: Overall performance metric

Example output::

   Scenario        N  Cells  Time(s)  Steps/s  Cells/step/s
   ----------------------------------------------------------------
   levelset       32   1,024      1.23    24.4      24,976
   levelset       64   4,096      4.82    6.2      25,393
   levelset      128  16,384     19.45    1.5      24,784
   levelset      256  65,536     78.12    0.4      25,042
   farsite        32   1,024      0.95    31.6      32,375
   farsite        64   4,096      3.54    8.5      34,790
   ...

Scaling Analysis
^^^^^^^^^^^^^^^^

The benchmark computes an empirical scaling exponent α where T ∝ N^α:

* **Expected range**: [0.8·dim, 1.6·dim]
* **Ideal scaling**: α = dim (linear in total cells)
* **Sub-linear scaling** (α < dim): Cache efficiency, vectorization
* **Super-linear scaling** (α > dim): Memory bandwidth limitations

Example::

   Estimated scaling exponent α = 2.03  (expected [1.6, 3.2])  ✓

Fire Spread Model Performance
------------------------------

Relative computational costs (normalized to Rothermel level-set):

+-------------------------+------------------+-------------------------------+
| Fire Spread Model       | Relative Cost    | Notes                         |
+=========================+==================+===============================+
| Rothermel (level-set)   | 1.0× (baseline)  | WENO5-Z + RK3 advection       |
+-------------------------+------------------+-------------------------------+
| FARSITE ellipse         | 0.7-0.9×         | Faster for coarse grids       |
+-------------------------+------------------+-------------------------------+
| Balbi physics-based     | 1.2-1.5×         | Additional radiation terms    |
+-------------------------+------------------+-------------------------------+
| Cruz crown fire         | 0.9-1.1×         | Simple algebraic model        |
+-------------------------+------------------+-------------------------------+
| Cheney-Gould grassfire  | 0.8-1.0×         | Piecewise-linear formula      |
+-------------------------+------------------+-------------------------------+
| FBP (Canadian)          | 0.9-1.1×         | Lookup tables, efficient      |
+-------------------------+------------------+-------------------------------+
| Lautenberger semi-phys  | 1.1-1.4×         | Physics-based coefficients    |
+-------------------------+------------------+-------------------------------+

**Note**: Costs vary with grid resolution, timestep size, and feature configuration.

Feature Overhead
----------------

Performance impact of optional features (% overhead relative to baseline):

+-------------------------------+------------------+---------------------------------+
| Feature                       | Overhead         | Notes                           |
+===============================+==================+=================================+
| Albini spotting               | 5-15%            | Depends on spot frequency       |
+-------------------------------+------------------+---------------------------------+
| Ember accumulation            | 2-5%             | Additional field updates        |
+-------------------------------+------------------+---------------------------------+
| Radiation preheating          | 10-20%           | View factor calculations        |
+-------------------------------+------------------+---------------------------------+
| Crown fire (Van Wagner)       | 3-8%             | Initiation criteria checks      |
+-------------------------------+------------------+---------------------------------+
| Bulk fuel consumption         | 5-10%            | Post-frontal burnout tracking   |
+-------------------------------+------------------+---------------------------------+
| FMC phenology                 | 1-3%             | Seasonal moisture updates       |
+-------------------------------+------------------+---------------------------------+
| McArthur moisture scaling     | <2%              | Temperature/RH modulation       |
+-------------------------------+------------------+---------------------------------+
| Periodic wind gusts           | <1%              | Sinusoidal wind updates         |
+-------------------------------+------------------+---------------------------------+

**Recommendation**: Enable only features needed for your simulation to minimize overhead.

Optimization Guidelines
-----------------------

Grid Resolution
^^^^^^^^^^^^^^^

Choose grid resolution based on:

* **Feature scale**: Cell size should resolve fire perimeter features (typically 5-50 m)
* **Computational budget**: Runtime scales as O(N²) in 2D
* **Accuracy requirements**: Finer grids reduce numerical diffusion

**Example**: For 1 km × 1 km domain:

* 32×32 grid (31 m cells): Fast, coarse features only
* 64×64 grid (16 m cells): Good balance for most applications
* 128×128 grid (8 m cells): High accuracy, slower
* 256×256 grid (4 m cells): Very high accuracy, expensive

Timestep Selection
^^^^^^^^^^^^^^^^^^

The CFL condition limits timestep size::

   Δt ≤ CFL × Δx / max(ROS)

**Recommendations**:

* **CFL = 0.5**: Safe default, stable for all scenarios
* **CFL = 0.7-0.9**: Faster, use for well-tested scenarios
* **CFL < 0.5**: Use if stability issues occur

Reinitialization Frequency
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Level-set reinitialization maintains the signed-distance property:

* ``reinit_int = -1``: Never reinitialize (fastest, may drift over time)
* ``reinit_int = 10-20``: Good balance for most scenarios
* ``reinit_int = 5``: Frequent reinitialization (most accurate, slower)

**Cost**: Reinitialization typically adds 20-30% overhead per occurrence.

Plotfile Output
^^^^^^^^^^^^^^^

Plotfile I/O can dominate runtime for large grids and frequent output:

* Disable during benchmarking: ``plot_int = -1``
* Minimize output frequency: ``plot_int = 20`` or higher
* Reduce output fields (future feature)

Example: 256×256 grid, plotfile writing can take 1-5 seconds per output.

Parallelization
---------------

Current Status
^^^^^^^^^^^^^^

The solver is primarily **serial** with some parallelization:

* **AMReX-level parallelism**: MPI domain decomposition (future)
* **Thread-level parallelism**: OpenMP for some kernels (partial)

**Future work**: Full MPI+OpenMP hybrid parallelization for HPC systems.

Memory Usage
------------

Approximate memory requirements (per grid cell, 2D):

* **Base solver**: ~100-200 bytes/cell (level-set φ, ROS, intensity, etc.)
* **With spotting**: +50-100 bytes/cell (ember density, spot tracking)
* **With crown fire**: +50 bytes/cell (crown fuel fields)
* **AMReX overhead**: ~50 bytes/cell (metadata, ghost cells)

**Example**: 256×256 grid (65,536 cells):

* Base: ~13-20 MB
* Full features: ~20-30 MB

**3D scaling**: Memory scales as O(N³), limit 3D grids to <100³ cells without HPC resources.

Profiling
---------

AMReX Built-in Profiling
^^^^^^^^^^^^^^^^^^^^^^^^^

AMReX includes performance profiling tools. Build with::

   cmake -S . -B build -DAMReX_TINY_PROFILE=ON

Run simulation, profiling output appears in console::

   TinyProfiler total time across processes [min...avg...max]
   LevelSet::Advect3D        [0.123, 0.125, 0.128]
   LevelSet::Reinit          [0.045, 0.046, 0.047]
   Rothermel::ComputeROS     [0.089, 0.090, 0.091]
   ...

External Profiling Tools
^^^^^^^^^^^^^^^^^^^^^^^^

**Linux perf**::

   perf record -g ./levelset inputs.i
   perf report

**gprof** (build with ``-pg``)::

   gfortran -pg ...  # or gcc -pg
   ./levelset inputs.i
   gprof levelset gmon.out > profile.txt

**Valgrind callgrind**::

   valgrind --tool=callgrind ./levelset inputs.i
   kcachegrind callgrind.out.*

GPU Acceleration
----------------

**Status**: GPU support via AMReX GPU backends is under development.

**Future capabilities**:

* CUDA/HIP backends for NVIDIA/AMD GPUs
* GPU-accelerated fire spread kernels
* Hybrid CPU/GPU execution

**Estimated speedup**: 5-20× for large grids (>128² cells) on modern GPUs.

**Note**: GPU scaling benchmarks are deferred pending infrastructure development.

Best Practices
--------------

For Performance-Critical Simulations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

1. **Profile first**: Identify actual bottlenecks before optimizing
2. **Right-size grid**: Balance accuracy and performance
3. **Disable unused features**: Turn off spotting, crown fire if not needed
4. **Reduce plotfile frequency**: I/O can dominate runtime
5. **Use coarser propagation**: FARSITE can be faster than level-set for operational forecasts

For Accuracy-Critical Simulations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

1. **Grid convergence study**: Test multiple resolutions
2. **Enable reinitialization**: ``reinit_int = 10`` or less
3. **Conservative CFL**: ``cfl = 0.5``
4. **Full physics**: Enable all relevant features

For Development and Testing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

1. **Small grids**: 32×32 or 64×64 for rapid iteration
2. **Short runs**: ``nsteps = 50`` or ``final_time = 300``
3. **Disable I/O**: ``plot_int = -1``
4. **Single model**: Test one fire spread model at a time

Performance Reporting
---------------------

When reporting performance issues or regressions:

1. **System specs**: CPU model, RAM, compiler version
2. **Build configuration**: Debug/Release, compiler flags
3. **Test case**: Input file or scenario description
4. **Timing data**: ``timing_results.csv`` from benchmark script
5. **Profiling data**: If available

Example issue report::

   System: Intel i9-12900K, 64 GB RAM, GCC 11.3, Ubuntu 22.04
   Build: Release mode, -DCMAKE_BUILD_TYPE=Release
   Test: 256×256 Rothermel level-set, 100 timesteps
   Performance: 0.45 steps/s (expected: 0.6 steps/s based on v1.0)
   Regression: 25% slower than previous version

References
----------

* AMReX documentation: https://amrex-codes.github.io/amrex/docs_html/
* WENO schemes: Jiang & Shu (1996), "Efficient Implementation of Weighted ENO Schemes"
* Level-set methods: Osher & Fedkiw (2003), "Level Set Methods and Dynamic Implicit Surfaces"

See Also
--------

* :doc:`regtests` - Regression test suite including timing benchmarks
* :doc:`building` - Build configuration options for performance
* :doc:`usage` - Runtime parameters affecting performance
