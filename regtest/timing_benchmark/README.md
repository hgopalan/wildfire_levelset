# Timing Benchmark Regression Test
# ==================================
#
# This directory contains a multi-resolution wall-clock timing benchmark for
# the wildfire level-set solver.  It is used to detect performance regressions
# and to verify that execution time scales correctly with grid size.
#
# Scenarios
# ----------
# 1. Level-set advection (propagation_method=levelset)
# 2. FARSITE ellipse propagation  (propagation_method=farsite)
#
# Grid sizes benchmarked
# -----------------------
# 2D builds: 32, 64, 128, 256 cells per side
# 3D builds: 16, 32, 64 cells per side
#
# Each resolution runs for a fixed number of time steps (default: 30) with
# plotfile output disabled so disk I/O does not contaminate timing.
#
# Pass criteria
# -------------
# * The solver exits successfully (return code 0) for all sizes.
# * Wall-clock time is monotonically non-decreasing with grid size
#   (5 % jitter tolerance to forgive OS scheduling noise).
# * The empirical scaling exponent α in T ∝ N^α falls in the range
#   [0.8·dim, 1.6·dim] where dim is 2 or 3.
#
# Output
# ------
# timing_results.csv  –  full table of timing data (written to the CTest
#                         working directory; inspect after a run with:
#                         cat build/regtest/timing_benchmark/timing_results.csv)
#
# Running manually
# ----------------
# From the build directory:
#
#   python3 ../regtest/timing_benchmark/run_benchmark.py \
#       --exe ./levelset --dim 2 --nsteps 20
#
# Adjust --nsteps for faster / more accurate timings.
# Use --dry-run to preview the generated inputs without running the solver.
# Use --skip-scaling to suppress the scaling-law check (e.g. on single-core CI).
