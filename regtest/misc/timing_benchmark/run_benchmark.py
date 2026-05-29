#!/usr/bin/env python3
"""
run_benchmark.py – Multi-resolution timing benchmark for the wildfire level-set
solver.

Runs the solver at a set of increasing grid resolutions, measures wall-clock
execution time for a fixed number of time steps, and checks for monotonically
increasing runtimes and plausible scaling behaviour.

Two scenarios are benchmarked in sequence:

1. **Level-set advection** (``propagation_method=levelset``) – exercises the
   WENO5Z / RK3 advection kernel and Rothermel ROS computation.

2. **FARSITE ellipse propagation** (``propagation_method=farsite``) – exercises
   the FARSITE wavelet expansion with Anderson L/W ratio.

Grid resolutions tested
-----------------------
:2D: 32, 64, 128, 256 cells per side.
:3D: 16, 32, 64 cells per side  (smaller to keep CI times reasonable).

All runs use a flat domain, constant wind, sphere ignition, and Rothermel FM4
(chaparral) fuel.  Plotfiles are *disabled* (``plot_int = -1``) so that disk
I/O does not dominate the timing signal.

Output
------
* ``timing_results.csv`` in the working directory with columns:
  ``scenario, n_cells, wall_time_s, steps_per_second, cells_per_step_per_s``
* Printed summary table and pass/fail verdict.

Failure criteria (non-zero exit code)
--------------------------------------
* The solver exits with a non-zero code for any resolution.
* The wall time is not monotonically non-decreasing with grid size for either
  scenario (small timing jitter is forgiven via a 5 % relative tolerance on
  the finer grid).
* Any individual run takes > ``--timeout`` seconds (default 300 s).

Usage
-----
Direct (from build directory)::

    python3 ../regtest/timing_benchmark/run_benchmark.py \\
        --exe ./levelset --dim 2

Via CTest (set by CMakeLists.txt automatically).

Arguments
---------
--exe PATH        Path to the ``levelset`` executable  **[required]**
--dim {2,3}       Spatial dimension of the build (default: 3)
--nsteps N        Fixed number of time steps per run (default: 30)
--resolutions N…  Custom list of cell counts (overrides built-in defaults)
--timeout T       Per-run wall-time limit in seconds (default: 300)
--output FILE     CSV output file (default: timing_results.csv)
--dry-run         Print generated inputs files and exit without running solver
--skip-scaling    Skip the scaling-law check (useful for debug/CI resource limits)
"""

from __future__ import annotations

import argparse
import csv
import os
import subprocess
import sys
import tempfile
import time
import textwrap
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional


# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------

DEFAULT_RESOLUTIONS_2D = [32, 64, 128, 256]
DEFAULT_RESOLUTIONS_3D = [16, 32, 64]
DEFAULT_NSTEPS         = 30
DEFAULT_TIMEOUT        = 300   # seconds
SCALING_TOLERANCE      = 0.05  # 5 % tolerance for monotonicity check


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class RunResult:
    scenario:       str
    n_cells:        int          # cells per side
    total_cells:    int          # total grid cells
    wall_time_s:    float
    nsteps:         int
    returncode:     int
    stdout:         str  = field(repr=False, default="")
    stderr:         str  = field(repr=False, default="")
    error_msg:      str  = ""

    @property
    def ok(self) -> bool:
        return self.returncode == 0

    @property
    def steps_per_second(self) -> float:
        return self.nsteps / self.wall_time_s if self.wall_time_s > 0 else 0.0

    @property
    def cells_per_step_per_s(self) -> float:
        return (self.total_cells * self.steps_per_second
                if self.wall_time_s > 0 else 0.0)


# ---------------------------------------------------------------------------
# Inputs-file templates
# ---------------------------------------------------------------------------

def _make_inputs_2d(n_cells: int, nsteps: int, scenario: str) -> str:
    """Return an inputs string for a 2-D run."""
    cell_size = 1000.0 / n_cells   # domain is always 1 km × 1 km
    radius    = max(2 * cell_size, 30.0)
    centre    = 500.0
    mg        = min(n_cells, 64)   # max_grid_size
    prop      = ("farsite" if scenario == "farsite" else "levelset")
    
    # Base configuration
    base_config = f"""\
        # Timing benchmark – 2D {scenario}  n={n_cells}
        n_cell_x = {n_cells}
        n_cell_y = {n_cells}
        max_grid_size = {mg}
        prob_lo_x = 0.0
        prob_lo_y = 0.0
        prob_lo_z = 0.0
        prob_hi_x = 1000.0
        prob_hi_y = 1000.0
        prob_hi_z = 1.0

        nsteps   = {nsteps}
        cfl      = 0.5
        plot_int = -1

        u_x = 4.0
        u_y = 1.0
        u_z = 0.0

        source_type    = sphere
        center_x       = {centre}
        center_y       = {centre}
        center_z       = 0.5
        sphere_radius  = {radius:.2f}

        reinit_int = -1

        propagation_method   = {prop}
        
        fire_stats_file = ""
        write_perimeter_csv     = 0
        write_perimeter_geojson = 0
    """
    
    # Model-specific configuration
    if scenario == "balbi":
        model_config = """\
        fire_spread_model = balbi
        balbi.T_a = 300.0
        balbi.T_f = 1000.0
        balbi.T_i = 600.0
        rothermel.fuel_model = FM4
        rothermel.M_f = 0.08
        """
    elif scenario == "cruz_crown":
        model_config = """\
        fire_spread_model = cruz_crown
        cruz_crown.CBD  = 0.15
        cruz_crown.MC10 = 8.0
        """
    elif scenario == "cheney_gould":
        model_config = """\
        fire_spread_model = cheney_gould
        cheney_gould.moisture = 8.0
        cheney_gould.curing   = 1.0
        """
    elif scenario == "fbp_o1a":
        model_config = """\
        fire_spread_model = fbp_o1a
        fbp.fuel_type = o1a
        fbp.curing = 80.0
        fbp.moisture = 8.0
        """
    elif scenario == "lautenberger":
        model_config = """\
        fire_spread_model = lautenberger
        rothermel.fuel_model = FM4
        rothermel.M_f = 0.08
        lautenberger.A_L = 1.05e-5
        lautenberger.B_L = 2.5
        lautenberger.C_L = 0.45
        lautenberger.D_L = 0.50
        """
    else:  # levelset or farsite
        model_config = """\
        rothermel.fuel_model = FM4
        rothermel.M_f        = 0.08
        farsite.use_anderson_LW = 1
        farsite.phi_threshold   = 0.1
        """
    
    return textwrap.dedent(base_config + model_config)


def _make_inputs_3d(n_cells: int, nsteps: int, scenario: str) -> str:
    """Return an inputs string for a 3-D run."""
    cell_size = 1000.0 / n_cells
    radius    = max(2 * cell_size, 30.0)
    centre    = 500.0
    mg        = min(n_cells, 32)
    prop      = ("farsite" if scenario == "farsite" else "levelset")
    
    # Base configuration
    base_config = f"""\
        # Timing benchmark – 3D {scenario}  n={n_cells}
        n_cell_x = {n_cells}
        n_cell_y = {n_cells}
        n_cell_z = {n_cells}
        max_grid_size = {mg}
        prob_lo_x = 0.0
        prob_lo_y = 0.0
        prob_lo_z = 0.0
        prob_hi_x = 1000.0
        prob_hi_y = 1000.0
        prob_hi_z = 1000.0

        nsteps   = {nsteps}
        cfl      = 0.5
        plot_int = -1

        u_x = 4.0
        u_y = 1.0
        u_z = 0.0

        source_type    = sphere
        center_x       = {centre}
        center_y       = {centre}
        center_z       = {centre}
        sphere_radius  = {radius:.2f}

        reinit_int = -1

        propagation_method   = {prop}
        
        fire_stats_file = ""
        write_perimeter_csv     = 0
        write_perimeter_geojson = 0
    """
    
    # Model-specific configuration
    if scenario == "balbi":
        model_config = """\
        fire_spread_model = balbi
        balbi.T_a = 300.0
        balbi.T_f = 1000.0
        balbi.T_i = 600.0
        rothermel.fuel_model = FM4
        rothermel.M_f = 0.08
        """
    elif scenario == "cruz_crown":
        model_config = """\
        fire_spread_model = cruz_crown
        cruz_crown.CBD  = 0.15
        cruz_crown.MC10 = 8.0
        """
    elif scenario == "cheney_gould":
        model_config = """\
        fire_spread_model = cheney_gould
        cheney_gould.moisture = 8.0
        cheney_gould.curing   = 1.0
        """
    elif scenario == "fbp_o1a":
        model_config = """\
        fire_spread_model = fbp_o1a
        fbp.fuel_type = o1a
        fbp.curing = 80.0
        fbp.moisture = 8.0
        """
    elif scenario == "lautenberger":
        model_config = """\
        fire_spread_model = lautenberger
        rothermel.fuel_model = FM4
        rothermel.M_f = 0.08
        lautenberger.A_L = 1.05e-5
        lautenberger.B_L = 2.5
        lautenberger.C_L = 0.45
        lautenberger.D_L = 0.50
        """
    else:  # levelset or farsite
        model_config = """\
        rothermel.fuel_model = FM4
        rothermel.M_f        = 0.08
        farsite.use_anderson_LW = 1
        farsite.phi_threshold   = 0.1
        """
    
    return textwrap.dedent(base_config + model_config)


def make_inputs(n_cells: int, nsteps: int, scenario: str, dim: int) -> str:
    if dim == 2:
        return _make_inputs_2d(n_cells, nsteps, scenario)
    return _make_inputs_3d(n_cells, nsteps, scenario)


# ---------------------------------------------------------------------------
# Execution
# ---------------------------------------------------------------------------

def run_solver(exe: str, inputs_text: str, timeout: float) -> tuple[float, int, str, str]:
    """
    Write *inputs_text* to a temporary file, run the solver, and return
    ``(wall_time_s, returncode, stdout, stderr)``.
    """
    with tempfile.NamedTemporaryFile(mode="w", suffix=".i",
                                    delete=False) as tmp:
        tmp.write(inputs_text)
        tmp_path = tmp.name

    try:
        t0 = time.perf_counter()
        result = subprocess.run(
            [exe, tmp_path],
            capture_output=True,
            text=True,
            timeout=timeout,
        )
        elapsed = time.perf_counter() - t0
        return elapsed, result.returncode, result.stdout, result.stderr
    except subprocess.TimeoutExpired:
        return timeout, -1, "", f"TIMEOUT after {timeout} s"
    finally:
        try:
            os.unlink(tmp_path)
        except OSError:
            pass


# ---------------------------------------------------------------------------
# Benchmark driver
# ---------------------------------------------------------------------------

def run_scenario(
    exe: str,
    scenario: str,
    resolutions: list[int],
    nsteps: int,
    dim: int,
    timeout: float,
    dry_run: bool,
) -> list[RunResult]:
    """Run one scenario at each resolution and return RunResult objects."""
    results: list[RunResult] = []
    print(f"\n{'='*64}")
    print(f"Scenario: {scenario}  (dim={dim}D)")
    print(f"{'='*64}")

    for n in resolutions:
        total = n ** dim
        inputs = make_inputs(n, nsteps, scenario, dim)

        if dry_run:
            print(f"\n--- n={n} inputs ---")
            print(inputs)
            continue

        print(f"  n={n:4d}  ({total:>10,} cells)  ...", end="", flush=True)
        wall, rc, stdout, stderr = run_solver(exe, inputs, timeout)

        msg = ""
        if rc != 0:
            # Try to extract the first error line from stdout/stderr
            for line in (stdout + stderr).splitlines():
                if any(kw in line.lower() for kw in ("error", "abort", "failed")):
                    msg = line.strip()[:100]
                    break
            if not msg:
                msg = (stderr or stdout).strip()[:100]

        r = RunResult(
            scenario=scenario,
            n_cells=n,
            total_cells=total,
            wall_time_s=wall,
            nsteps=nsteps,
            returncode=rc,
            stdout=stdout,
            stderr=stderr,
            error_msg=msg,
        )
        results.append(r)

        status = "✓" if r.ok else f"✗  rc={rc}"
        print(f"  {wall:6.2f} s  {r.steps_per_second:5.1f} steps/s"
              f"  {r.cells_per_step_per_s:>12,.0f} cells/step/s  {status}")
        if not r.ok:
            print(f"       {msg}")

    return results


# ---------------------------------------------------------------------------
# Analysis and validation
# ---------------------------------------------------------------------------

def check_monotone(results: list[RunResult], tol: float = SCALING_TOLERANCE) -> list[str]:
    """
    Check that wall-time is monotonically non-decreasing with grid size.

    A small tolerance *tol* (fraction of the larger time) is allowed to
    forgive timing jitter when the two runs are very close.
    Returns a list of error strings (empty = pass).
    """
    errors = []
    ok_results = [r for r in results if r.ok]
    for i in range(1, len(ok_results)):
        prev, curr = ok_results[i - 1], ok_results[i]
        # Allow curr.wall_time_s < prev.wall_time_s only if difference < tol×max
        if curr.wall_time_s < prev.wall_time_s:
            delta = prev.wall_time_s - curr.wall_time_s
            ref   = max(prev.wall_time_s, curr.wall_time_s)
            if ref > 0 and delta / ref > tol:
                errors.append(
                    f"  Non-monotone: n={prev.n_cells} → n={curr.n_cells}  "
                    f"({prev.wall_time_s:.2f} s → {curr.wall_time_s:.2f} s)"
                )
    return errors


def estimate_scaling_exponent(results: list[RunResult]) -> Optional[float]:
    """
    Least-squares estimate of scaling exponent α in T ∝ N^α from all
    successful runs.  Returns None if fewer than 2 successful results.
    """
    import math
    ok = [r for r in results if r.ok and r.n_cells > 0 and r.wall_time_s > 0]
    if len(ok) < 2:
        return None
    log_n = [math.log(r.n_cells)       for r in ok]
    log_t = [math.log(r.wall_time_s)   for r in ok]
    n_pts  = len(ok)
    mean_n = sum(log_n) / n_pts
    mean_t = sum(log_t) / n_pts
    num    = sum((x - mean_n) * (y - mean_t) for x, y in zip(log_n, log_t))
    den    = sum((x - mean_n) ** 2           for x in log_n)
    return num / den if den > 0 else None


# ---------------------------------------------------------------------------
# CSV writer
# ---------------------------------------------------------------------------

def write_csv(results: list[RunResult], out_path: str) -> None:
    fieldnames = [
        "scenario", "n_cells", "total_cells",
        "wall_time_s", "nsteps", "steps_per_second",
        "cells_per_step_per_s", "returncode",
    ]
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in results:
            w.writerow({
                "scenario":            r.scenario,
                "n_cells":             r.n_cells,
                "total_cells":         r.total_cells,
                "wall_time_s":         f"{r.wall_time_s:.4f}",
                "nsteps":              r.nsteps,
                "steps_per_second":    f"{r.steps_per_second:.3f}",
                "cells_per_step_per_s":f"{r.cells_per_step_per_s:.1f}",
                "returncode":          r.returncode,
            })
    print(f"\nTiming results written to: {out_path}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--exe", required=True, metavar="PATH",
                   help="Path to the levelset executable")
    p.add_argument("--dim", type=int, choices=[2, 3], default=3, metavar="{2,3}",
                   help="Spatial dimension of the build (default: 3)")
    p.add_argument("--nsteps", type=int, default=DEFAULT_NSTEPS, metavar="N",
                   help=f"Number of solver steps per run (default: {DEFAULT_NSTEPS})")
    p.add_argument("--resolutions", nargs="+", type=int, metavar="N",
                   help="Custom list of cell counts per side (overrides defaults)")
    p.add_argument("--timeout", type=float, default=DEFAULT_TIMEOUT, metavar="T",
                   help=f"Per-run timeout in seconds (default: {DEFAULT_TIMEOUT})")
    p.add_argument("--output", default="timing_results.csv", metavar="FILE",
                   help="Output CSV file (default: timing_results.csv)")
    p.add_argument("--dry-run", action="store_true",
                   help="Print generated inputs and exit without running solver")
    p.add_argument("--skip-scaling", action="store_true",
                   help="Skip the scaling-law check (still checks for run failures)")
    p.add_argument("--scenarios", nargs="+",
                   choices=["levelset", "farsite", "balbi", "cruz_crown", 
                            "cheney_gould", "fbp_o1a", "lautenberger"],
                   default=["levelset", "farsite"],
                   help="Scenarios to benchmark (default: levelset and farsite)")
    return p


def main() -> None:
    args = _build_parser().parse_args()

    exe = Path(args.exe)
    if not args.dry_run and not exe.exists():
        sys.exit(f"ERROR: solver executable not found: {exe}")

    resolutions = (args.resolutions or
                   (DEFAULT_RESOLUTIONS_2D if args.dim == 2
                    else DEFAULT_RESOLUTIONS_3D))

    print(f"Timing benchmark: dim={args.dim}D  nsteps={args.nsteps}")
    print(f"Resolutions: {resolutions}")
    print(f"Scenarios:   {args.scenarios}")
    print(f"Executable:  {exe}")

    all_results: list[RunResult] = []
    all_errors: list[str] = []

    for scenario in args.scenarios:
        results = run_scenario(
            exe=str(exe),
            scenario=scenario,
            resolutions=resolutions,
            nsteps=args.nsteps,
            dim=args.dim,
            timeout=args.timeout,
            dry_run=args.dry_run,
        )
        all_results.extend(results)

        if args.dry_run:
            continue

        # --- Check for solver failures ---
        for r in results:
            if not r.ok:
                all_errors.append(
                    f"Solver failed: scenario={scenario} n={r.n_cells} "
                    f"rc={r.returncode}  {r.error_msg}"
                )

        # --- Scaling analysis ---
        ok = [r for r in results if r.ok]
        if len(ok) >= 2:
            alpha = estimate_scaling_exponent(ok)
            if alpha is not None:
                dim = args.dim
                expected_lo = float(dim) * 0.8
                expected_hi = float(dim) * 1.6
                in_range = expected_lo <= alpha <= expected_hi
                flag = "✓" if in_range else "⚠ (outside expected range)"
                print(f"\n  Estimated scaling exponent α = {alpha:.2f}  "
                      f"(expected [{expected_lo:.1f}, {expected_hi:.1f}])  {flag}")
                if not in_range and not args.skip_scaling:
                    all_errors.append(
                        f"Scaling exponent α={alpha:.2f} is outside the expected "
                        f"range [{expected_lo:.1f}, {expected_hi:.1f}] for "
                        f"scenario={scenario} dim={dim}D"
                    )

        # --- Monotonicity check ---
        if not args.skip_scaling:
            mono_errors = check_monotone(ok)
            if mono_errors:
                print("\n  WARNING: non-monotone wall-time detected:")
                for e in mono_errors:
                    print(f"   {e}")
                all_errors.extend(mono_errors)

    if args.dry_run:
        return

    # --- Write CSV ---
    if all_results:
        write_csv(all_results, args.output)

    # --- Summary table ---
    print(f"\n{'='*64}")
    print("Summary")
    print(f"{'='*64}")
    print(f"{'Scenario':<12} {'N':>6} {'Cells':>10} "
          f"{'Time(s)':>9} {'Steps/s':>9} {'Cells/step/s':>14}")
    print("-" * 64)
    for r in all_results:
        status = "✓" if r.ok else "✗"
        print(f"{r.scenario:<12} {r.n_cells:>6} {r.total_cells:>10,} "
              f"{r.wall_time_s:>9.2f} {r.steps_per_second:>9.1f} "
              f"{r.cells_per_step_per_s:>14,.0f}  {status}")

    # --- Pass / fail ---
    if all_errors:
        print(f"\nFAILED ({len(all_errors)} error(s)):")
        for e in all_errors:
            print(f"  ERROR: {e}")
        sys.exit(1)

    print("\nAll timing benchmark checks PASSED.")


if __name__ == "__main__":
    main()
