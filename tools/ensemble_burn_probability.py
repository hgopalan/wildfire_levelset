#!/usr/bin/env python3
"""
ensemble_burn_probability.py – FSPro-style ensemble burn probability driver
for wildfire_levelset.

Runs the solver N times with perturbed wind speed, wind direction, and dead
fuel moisture, then accumulates per-cell burn counts from the final-step
``phi_negative_NNNN.dat`` output files and writes a probability map.

Workflow
--------
1. Read a template ``inputs.i`` file.
2. Generate N parameter samples (Latin hypercube by default, or random).
3. For each sample:
   a. Create a scratch run directory.
   b. Write a modified ``inputs.i`` with the perturbed parameters.
   c. Run the solver executable.
   d. Read the final ``phi_negative_NNNN.dat`` (all phi < 0 cell centres).
4. Aggregate per-cell burn counts on a regular grid.
5. Write ``burn_probability.csv`` (X, Y, P_burn columns) and optional GeoJSON.

Perturbation model
------------------
Three independent parameters are perturbed:
  • wind speed factor  (multiplicative, lognormal)
  • wind direction     (additive offset, normal)
  • 1-hr dead fuel moisture offset (additive, normal)

The distributions are parameterised by their mean (= 1.0 for speed factor,
0.0 for direction and moisture offset) and a standard deviation.  The samples
are generated using a Latin hypercube or, if ``scipy`` is unavailable, simple
uniform random sampling.

Usage examples
--------------
  # Basic: 50 runs, ±20% wind speed, ±15° direction, ±2% moisture
  python3 tools/ensemble_burn_probability.py \\
      --exe ./wildfire_levelset \\
      --inputs inputs.i \\
      --n-runs 50 \\
      --wind-speed-sigma 0.20 \\
      --wind-dir-sigma 15.0 \\
      --moisture-sigma 0.02

  # More runs, georeferenced probability map
  python3 tools/ensemble_burn_probability.py \\
      --exe ./wildfire_levelset \\
      --inputs inputs.i \\
      --n-runs 200 \\
      --wind-speed-sigma 0.30 \\
      --wind-dir-sigma 20.0 \\
      --moisture-sigma 0.03 \\
      --resolution 30 \\
      --out burn_probability.csv \\
      --geojson burn_probability.geojson

  # Deterministic seed for reproducibility
  python3 tools/ensemble_burn_probability.py \\
      --exe ./wildfire_levelset \\
      --inputs inputs.i \\
      --n-runs 50 \\
      --seed 42

Options
-------
  --exe FILE           Solver executable path (default: ./wildfire_levelset)
  --inputs FILE        Template inputs.i file (default: inputs.i)
  --n-runs N           Number of ensemble members (default: 50)
  --wind-speed-sigma S Std dev of multiplicative wind-speed factor (default: 0.20)
  --wind-dir-sigma D   Std dev of additive wind-direction offset [deg] (default: 15.0)
  --moisture-sigma M   Std dev of additive M_d1 offset [fraction] (default: 0.02)
  --resolution R       Grid spacing for probability accumulation [m] (default: auto)
  --out FILE           Output burn probability CSV (default: burn_probability.csv)
  --geojson FILE       Optional GeoJSON probability raster output
  --work-dir DIR       Base working directory for run scratch dirs (default: /tmp/ensemble)
  --seed N             Random seed (default: 0 = system clock)
  --sampler TYPE       Sampling method: "lhs" or "random" (default: lhs)
  --jobs J             Number of parallel solver runs (default: 1 = sequential)
  --keep-runs          Keep individual run directories after aggregation

References
----------
  Finney, M.A. et al. (2011). A method for ensemble wildland fire simulation.
    Environmental Modelling & Software, 26(10), 1352-1359.
  Iman, R.L. & Conover, W.J. (1980). Small sample sensitivity analysis
    techniques for computer models. Communications in Statistics, A9(17).
"""

from __future__ import annotations

import argparse
import glob
import json
import math
import os
import random
import re
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Optional, Tuple


# ---------------------------------------------------------------------------
# Latin hypercube sampling (pure Python, no scipy required)
# ---------------------------------------------------------------------------

def lhs_uniform(n_samples: int, n_dims: int, rng: random.Random) -> List[List[float]]:
    """Return n_samples × n_dims matrix of Latin-hypercube uniform [0,1] samples."""
    result = [[0.0] * n_dims for _ in range(n_samples)]
    for d in range(n_dims):
        order = list(range(n_samples))
        rng.shuffle(order)
        for i, pos in enumerate(order):
            result[pos][d] = (i + rng.random()) / n_samples
    return result


_MIN_Q_THRESHOLD = 1e-300  # minimum argument for log to avoid -inf
_MIN_MOISTURE    = 0.01    # minimum dead fuel moisture fraction (physical lower bound)
_MAX_MOISTURE    = 0.60    # maximum dead fuel moisture fraction (physical upper bound)


def norm_ppf(p: float) -> float:
    """Rational approximation of the normal quantile (Beasley-Springer-Moro)."""
    # Coefficients from Hart (1968)
    a = [2.515517, 0.802853, 0.010328]
    b = [1.432788, 0.189269, 0.001308]
    if p < 0.5:
        sign, q = -1, p
    else:
        sign, q = 1, 1.0 - p
    t = math.sqrt(-2.0 * math.log(max(q, _MIN_Q_THRESHOLD)))
    num = a[0] + a[1]*t + a[2]*t*t
    den = 1.0 + b[0]*t + b[1]*t*t + b[2]*t*t*t
    return sign * (t - num/den)


def lhs_normal(n_samples: int, n_dims: int, rng: random.Random) -> List[List[float]]:
    """Return n_samples × n_dims matrix of LHS normal(0,1) samples."""
    u = lhs_uniform(n_samples, n_dims, rng)
    return [[norm_ppf(max(1e-6, min(1 - 1e-6, u[i][d]))) for d in range(n_dims)]
            for i in range(n_samples)]


# ---------------------------------------------------------------------------
# Parameter generation
# ---------------------------------------------------------------------------

def generate_samples(
    n: int,
    speed_sigma: float,
    dir_sigma: float,
    moisture_sigma: float,
    sampler: str,
    rng: random.Random,
) -> List[Dict[str, float]]:
    """Generate n parameter samples as a list of dicts."""
    if sampler == "lhs":
        raw = lhs_normal(n, 3, rng)
    else:
        raw = [[rng.gauss(0, 1) for _ in range(3)] for _ in range(n)]

    samples = []
    for row in raw:
        # wind_speed_factor: lognormal – exp(sigma * z) so mean ≈ 1.0
        speed_factor = math.exp(speed_sigma * row[0])
        dir_offset   = dir_sigma * row[1]          # [deg]
        moist_offset = moisture_sigma * row[2]      # [fraction]
        samples.append({
            "speed_factor": speed_factor,
            "dir_offset": dir_offset,
            "moist_offset": moist_offset,
        })
    return samples


# ---------------------------------------------------------------------------
# inputs.i manipulation
# ---------------------------------------------------------------------------

_FLOAT_RE = re.compile(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?")


def _set_param(lines: List[str], key: str, value: str) -> List[str]:
    """Replace or append a key = value line in an inputs file."""
    new_line = f"{key} = {value}\n"
    prefix   = key + " "
    for i, line in enumerate(lines):
        stripped = line.lstrip()
        if stripped.startswith(prefix) or stripped.startswith(key + "="):
            lines[i] = new_line
            return lines
    lines.append(new_line)
    return lines


def _get_float(lines: List[str], key: str) -> Optional[float]:
    """Return float value of key from lines, or None."""
    prefix = key + " "
    for line in lines:
        stripped = line.lstrip()
        if stripped.startswith(prefix) or stripped.startswith(key + "="):
            m = _FLOAT_RE.search(stripped.split("=", 1)[-1])
            if m:
                return float(m.group())
    return None


def make_perturbed_inputs(
    template_lines: List[str],
    sample: Dict[str, float],
) -> List[str]:
    """Return a modified copy of template_lines with perturbed parameters."""
    lines = list(template_lines)

    # Read current wind
    ux = _get_float(lines, "u_x") or 0.0
    uy = _get_float(lines, "u_y") or 0.0
    speed = math.sqrt(ux**2 + uy**2)
    angle = math.atan2(uy, ux)

    # Apply speed factor and direction offset
    new_speed = speed * sample["speed_factor"]
    new_angle = angle + math.radians(sample["dir_offset"])
    new_ux = new_speed * math.cos(new_angle)
    new_uy = new_speed * math.sin(new_angle)

    lines = _set_param(lines, "u_x", f"{new_ux:.6f}")
    lines = _set_param(lines, "u_y", f"{new_uy:.6f}")

    # Apply moisture offset (clamp to [0.01, 0.60])
    m_d1 = (_get_float(lines, "rothermel.M_d1") or 0.08) + sample["moist_offset"]
    m_d1 = max(_MIN_MOISTURE, min(_MAX_MOISTURE, m_d1))
    lines = _set_param(lines, "rothermel.M_d1", f"{m_d1:.4f}")

    return lines


# ---------------------------------------------------------------------------
# Read phi_negative_NNNN.dat output
# ---------------------------------------------------------------------------

def read_phi_negative(run_dir: str) -> List[Tuple[float, float]]:
    """Return list of (x, y) burned cell centres from the last phi_negative file."""
    pattern = os.path.join(run_dir, "phi_negative_*.dat")
    files   = sorted(glob.glob(pattern))
    if not files:
        return []
    points = []
    with open(files[-1]) as f:
        for line in f:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            parts = stripped.split()
            if len(parts) >= 2:
                try:
                    points.append((float(parts[0]), float(parts[1])))
                except ValueError:
                    pass
    return points


# ---------------------------------------------------------------------------
# Grid accumulation
# ---------------------------------------------------------------------------

def accumulate_on_grid(
    all_points: List[List[Tuple[float, float]]],
    resolution: float,
) -> Tuple[Dict[Tuple[int, int], int], float, float, float]:
    """
    Snap all burn points to a regular grid of size ``resolution`` and
    count how many ensemble members burned each cell.

    Returns (counts_dict, x_origin, y_origin, resolution).
    counts_dict: (ix, iy) -> count
    """
    xs = [p[0] for run in all_points for p in run]
    ys = [p[1] for run in all_points for p in run]
    if not xs:
        return {}, 0.0, 0.0, resolution

    x_min, y_min = min(xs), min(ys)
    counts: Dict[Tuple[int, int], int] = {}

    for run_pts in all_points:
        visited: set = set()
        for x, y in run_pts:
            ix = int((x - x_min) / resolution)
            iy = int((y - y_min) / resolution)
            key = (ix, iy)
            if key not in visited:
                visited.add(key)
                counts[key] = counts.get(key, 0) + 1

    return counts, x_min, y_min, resolution


# ---------------------------------------------------------------------------
# Write outputs
# ---------------------------------------------------------------------------

def write_burn_probability_csv(
    out_path: str,
    counts: Dict[Tuple[int, int], int],
    n_runs: int,
    x_origin: float,
    y_origin: float,
    resolution: float,
) -> None:
    """Write burn probability CSV: X, Y, P_burn."""
    with open(out_path, "w") as f:
        f.write("# Ensemble burn probability map\n")
        f.write(f"# n_runs = {n_runs}  resolution = {resolution} m\n")
        f.write("X,Y,P_burn\n")
        for (ix, iy), cnt in sorted(counts.items()):
            x = x_origin + (ix + 0.5) * resolution
            y = y_origin + (iy + 0.5) * resolution
            p = cnt / n_runs
            f.write(f"{x:.2f},{y:.2f},{p:.4f}\n")
    print(f"Wrote burn probability CSV: {out_path}  ({len(counts)} cells)")


def write_burn_probability_geojson(
    out_path: str,
    counts: Dict[Tuple[int, int], int],
    n_runs: int,
    x_origin: float,
    y_origin: float,
    resolution: float,
) -> None:
    """Write burn probability as GeoJSON FeatureCollection of square polygons."""
    features = []
    for (ix, iy), cnt in sorted(counts.items()):
        x0 = x_origin + ix * resolution
        y0 = y_origin + iy * resolution
        x1, y1 = x0 + resolution, y0 + resolution
        p = cnt / n_runs
        features.append({
            "type": "Feature",
            "geometry": {
                "type": "Polygon",
                "coordinates": [[[x0, y0], [x1, y0], [x1, y1], [x0, y1], [x0, y0]]],
            },
            "properties": {"P_burn": round(p, 4), "count": cnt},
        })
    gj = {"type": "FeatureCollection", "features": features}
    with open(out_path, "w") as f:
        json.dump(gj, f)
    print(f"Wrote burn probability GeoJSON: {out_path}  ({len(features)} cells)")


# ---------------------------------------------------------------------------
# Single run execution
# ---------------------------------------------------------------------------

def run_member(
    run_id: int,
    exe: str,
    template_lines: List[str],
    sample: Dict[str, float],
    work_dir: str,
    keep: bool,
) -> Optional[List[Tuple[float, float]]]:
    """Execute one ensemble member. Returns list of burned (x,y) or None on failure."""
    run_dir = os.path.join(work_dir, f"run_{run_id:04d}")
    os.makedirs(run_dir, exist_ok=True)

    perturbed = make_perturbed_inputs(template_lines, sample)
    inputs_path = os.path.join(run_dir, "inputs.i")
    with open(inputs_path, "w") as f:
        f.writelines(perturbed)

    log_path = os.path.join(run_dir, "run.log")
    try:
        with open(log_path, "w") as log:
            result = subprocess.run(
                [exe, inputs_path],
                cwd=run_dir,
                stdout=log,
                stderr=subprocess.STDOUT,
                timeout=3600,
            )
        if result.returncode != 0:
            print(f"  WARNING: run {run_id} exited with code {result.returncode}", flush=True)
        pts = read_phi_negative(run_dir)
        print(f"  Run {run_id:4d}: speed×{sample['speed_factor']:.2f}  "
              f"dir+{sample['dir_offset']:.1f}°  M_d1+{sample['moist_offset']:.3f}  "
              f"→ {len(pts)} burned cells", flush=True)
        return pts
    except subprocess.TimeoutExpired:
        print(f"  WARNING: run {run_id} timed out", flush=True)
        return None
    except Exception as exc:
        print(f"  WARNING: run {run_id} failed: {exc}", flush=True)
        return None
    finally:
        if not keep:
            # Remove large AMReX plotfiles but keep the phi_negative .dat files
            for subdir in glob.glob(os.path.join(run_dir, "plt*")):
                shutil.rmtree(subdir, ignore_errors=True)
            for chkdir in glob.glob(os.path.join(run_dir, "chk*")):
                shutil.rmtree(chkdir, ignore_errors=True)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv=None):
    parser = argparse.ArgumentParser(
        description="FSPro-style ensemble burn probability driver for wildfire_levelset.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--exe",              default="./wildfire_levelset",
                        help="Solver executable (default: ./wildfire_levelset)")
    parser.add_argument("--inputs",           default="inputs.i",
                        help="Template inputs.i file (default: inputs.i)")
    parser.add_argument("--n-runs",           type=int, default=50,
                        help="Number of ensemble members (default: 50)")
    parser.add_argument("--wind-speed-sigma", type=float, default=0.20,
                        help="Std dev of lognormal wind-speed multiplier (default: 0.20)")
    parser.add_argument("--wind-dir-sigma",   type=float, default=15.0,
                        help="Std dev of wind-direction offset [deg] (default: 15.0)")
    parser.add_argument("--moisture-sigma",   type=float, default=0.02,
                        help="Std dev of M_d1 offset [fraction] (default: 0.02)")
    parser.add_argument("--resolution",       type=float, default=None,
                        help="Grid resolution for probability map [m] (default: auto)")
    parser.add_argument("--out",              default="burn_probability.csv",
                        help="Output CSV path (default: burn_probability.csv)")
    parser.add_argument("--geojson",          default=None,
                        help="Optional output GeoJSON path")
    parser.add_argument("--work-dir",         default=None,
                        help="Scratch directory for run dirs (default: /tmp/ensemble_XXXXX)")
    parser.add_argument("--seed",             type=int, default=0,
                        help="Random seed (0 = system clock) (default: 0)")
    parser.add_argument("--sampler",          default="lhs",
                        choices=["lhs", "random"],
                        help="Sampling method (default: lhs)")
    parser.add_argument("--jobs",             type=int, default=1,
                        help="Number of parallel runs (default: 1)")
    parser.add_argument("--keep-runs",        action="store_true",
                        help="Keep individual run directories after aggregation")

    args = parser.parse_args(argv)

    # Validate
    if not os.path.isfile(args.inputs):
        print(f"ERROR: template inputs file not found: {args.inputs}", file=sys.stderr)
        sys.exit(1)
    if not os.path.isfile(args.exe):
        print(f"ERROR: solver executable not found: {args.exe}", file=sys.stderr)
        sys.exit(1)

    with open(args.inputs) as f:
        template_lines = f.readlines()

    rng = random.Random(args.seed if args.seed != 0 else None)

    print(f"Generating {args.n_runs} ensemble samples ({args.sampler}) ...")
    samples = generate_samples(
        args.n_runs,
        args.wind_speed_sigma,
        args.wind_dir_sigma,
        args.moisture_sigma,
        args.sampler,
        rng,
    )

    work_dir = args.work_dir or tempfile.mkdtemp(prefix="ensemble_")
    os.makedirs(work_dir, exist_ok=True)
    print(f"Working directory: {work_dir}")
    print(f"Running {args.n_runs} members (jobs={args.jobs}) ...")

    all_points: List[List[Tuple[float, float]]] = []
    failed = 0

    if args.jobs > 1:
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            futures = {
                pool.submit(run_member, i, args.exe, template_lines,
                            samples[i], work_dir, args.keep_runs): i
                for i in range(args.n_runs)
            }
            for fut in as_completed(futures):
                pts = fut.result()
                if pts is not None:
                    all_points.append(pts)
                else:
                    failed += 1
    else:
        for i, sample in enumerate(samples):
            pts = run_member(i, args.exe, template_lines, sample,
                             work_dir, args.keep_runs)
            if pts is not None:
                all_points.append(pts)
            else:
                failed += 1

    n_success = len(all_points)
    print(f"\nCompleted {n_success}/{args.n_runs} runs ({failed} failed).")

    if n_success == 0:
        print("ERROR: all runs failed; cannot compute probability map.", file=sys.stderr)
        sys.exit(1)

    # Auto-detect resolution
    all_flat = [p for run in all_points for p in run]
    if args.resolution:
        resolution = args.resolution
    elif len(all_flat) > 1:
        # Estimate typical cell spacing from the bounding box and point count
        xs = [p[0] for p in all_flat]
        ys = [p[1] for p in all_flat]
        area = (max(xs) - min(xs)) * (max(ys) - min(ys))
        resolution = math.sqrt(area / max(len(all_flat), 1)) if area > 0 else 30.0
        resolution = max(1.0, round(resolution, 0))
        print(f"Auto-detected grid resolution: {resolution:.0f} m")
    else:
        resolution = 30.0
        print(f"Using default grid resolution: {resolution:.0f} m")

    counts, x0, y0, res = accumulate_on_grid(all_points, resolution)

    write_burn_probability_csv(args.out, counts, n_success, x0, y0, res)

    if args.geojson:
        write_burn_probability_geojson(args.geojson, counts, n_success, x0, y0, res)

    if not args.keep_runs and args.work_dir is None:
        # Clean up the auto-created temp directory
        shutil.rmtree(work_dir, ignore_errors=True)
        print(f"Cleaned up working directory: {work_dir}")


if __name__ == "__main__":
    main()
