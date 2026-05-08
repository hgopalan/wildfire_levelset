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
   c. Run the solver executable (serial or MPI).
   d. Read the final ``phi_negative_NNNN.dat`` (all phi < 0 cell centres).
4. Aggregate per-cell burn counts on a regular grid.
5. Write ``burn_probability.csv`` (X, Y, P_burn columns) and optional GeoJSON.

Additional outputs (Cap 7 – Exceedance probability curves)
-----------------------------------------------------------
``--area-exceedance``
    After all ensemble runs, computes and writes a fire-area exceedance curve
    (CCDF): the fraction of runs in which the total burned area exceeds each
    area value.  The curve is saved as ``area_exceedance.csv`` (columns:
    ``area_ha``, ``P_exceed``) and, optionally, as a PNG plot with
    ``--area-exceedance-plot``.

``--fl-thresholds M1 M2 …``
    Conditional flame-length exceedance probability maps: P(FL > M) at each
    grid cell, saved as CSV files.  Requires ``plot_int > 0`` in ``inputs.i``
    so that plotfiles are written.

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

MPI launch
----------
Pass ``--mpi-ranks N`` (N ≥ 1) to launch each solver member via MPI::

    mpirun -n N ./wildfire_levelset inputs.i amrex.use_gpu_aware_mpi=0

The MPI launcher command can be changed with ``--mpirun`` (default: ``mpirun``).
When ``--mpi-ranks`` is 0 (the default) the solver is run serially without
any MPI wrapper.

Usage examples
--------------
  # Basic serial: 50 runs, ±20% wind speed, ±15° direction, ±2% moisture
  python3 tools/ensemble_burn_probability.py \\
      --exe ./wildfire_levelset \\
      --inputs inputs.i \\
      --n-runs 50 \\
      --wind-speed-sigma 0.20 \\
      --wind-dir-sigma 15.0 \\
      --moisture-sigma 0.02

  # MPI: each ensemble member uses 4 MPI ranks
  python3 tools/ensemble_burn_probability.py \\
      --exe ./wildfire_levelset \\
      --inputs inputs.i \\
      --n-runs 50 \\
      --mpi-ranks 4

  # MPI with mpiexec launcher and 8 ranks per member, 2 members in parallel
  python3 tools/ensemble_burn_probability.py \\
      --exe ./wildfire_levelset \\
      --inputs inputs.i \\
      --n-runs 100 \\
      --mpirun mpiexec \\
      --mpi-ranks 8 \\
      --jobs 2

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

  # Area exceedance CCDF (Cap 7)
  python3 tools/ensemble_burn_probability.py \\
      --exe ./wildfire_levelset \\
      --inputs inputs.i \\
      --n-runs 50 \\
      --area-exceedance \\
      --area-exceedance-plot

  # Deterministic seed for reproducibility
  python3 tools/ensemble_burn_probability.py \\
      --exe ./wildfire_levelset \\
      --inputs inputs.i \\
      --n-runs 50 \\
      --seed 42

Options
-------
  --exe FILE              Solver executable path (default: ./wildfire_levelset)
  --inputs FILE           Template inputs.i file (default: inputs.i)
  --n-runs N              Number of ensemble members (default: 50)
  --mpi-ranks N           MPI ranks per ensemble member; 0 = serial (default: 0)
  --mpirun CMD            MPI launcher command (default: mpirun)
  --mpi-args ARGS         Extra space-separated arguments inserted after the
                          launcher but before "-n N" (default: "")
  --wind-speed-sigma S    Std dev of multiplicative wind-speed factor (default: 0.20)
  --wind-dir-sigma D      Std dev of additive wind-direction offset [deg] (default: 15.0)
  --moisture-sigma M      Std dev of additive M_d1 offset [fraction] (default: 0.02)
  --resolution R          Grid spacing for probability accumulation [m] (default: auto)
  --out FILE              Output burn probability CSV (default: burn_probability.csv)
  --geojson FILE          Optional GeoJSON probability raster output
  --work-dir DIR          Base working directory for run scratch dirs (default: /tmp/ensemble)
  --seed N                Random seed (default: 0 = system clock)
  --sampler TYPE          Sampling method: "lhs" or "random" (default: lhs)
  --jobs J                Number of parallel solver runs (default: 1 = sequential)
  --keep-runs             Keep individual run directories after aggregation
  --area-exceedance       Write area exceedance CCDF (Cap 7)
  --area-exceedance-out F  Path for area exceedance CSV (default: area_exceedance.csv)
  --area-exceedance-plot  Save area exceedance CCDF plot PNG
  --area-exceedance-plot-out F  Path for exceedance PNG (default: area_exceedance.png)

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
from typing import Dict, List, Optional, Set, Tuple


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
# Flame-length reader (reads from the final AMReX plotfile in a run dir)
# ---------------------------------------------------------------------------

def _read_plotfile_field(run_dir: str, field: str) -> Dict[Tuple[float, float], float]:
    """Read a named scalar field from the last plt#### directory in *run_dir*.

    Returns a dict mapping (x, y) → value for all cells in the plotfile.
    Returns an empty dict if no plotfile is found or the field is absent.
    """
    plt_dirs = sorted(glob.glob(os.path.join(run_dir, "plt[0-9][0-9][0-9][0-9]")))
    if not plt_dirs:
        return {}
    plt_path = plt_dirs[-1]

    # ---- Parse Header ----
    header_path = os.path.join(plt_path, "Header")
    if not os.path.isfile(header_path):
        return {}
    with open(header_path) as fh:
        lines = [l.rstrip("\n") for l in fh]

    try:
        idx = 1
        ncomp = int(lines[idx]); idx += 1
        varnames = []
        for _ in range(ncomp):
            varnames.append(lines[idx].strip()); idx += 1
        _spacedim = int(lines[idx]); idx += 1
        _time = float(lines[idx]); idx += 1
        _finest = int(lines[idx]); idx += 1
        problo = list(map(float, lines[idx].split())); idx += 1
        probhi = list(map(float, lines[idx].split())); idx += 1
        idx += 1  # ref ratios
        idx += 1  # box count
        box_line = lines[idx]; idx += 1
        nums = [int(t) for t in box_line.replace("(","").replace(")","").replace(","," ").split()
                if t.lstrip("-").isdigit()]
        if len(nums) < 4:
            return {}
        nx = nums[2] - nums[0] + 1
        ny = nums[3] - nums[1] + 1
    except (IndexError, ValueError):
        return {}

    if field not in varnames:
        return {}

    cell_h = os.path.join(plt_path, "Level_0", "Cell_H")
    if not os.path.isfile(cell_h):
        return {}

    # ---- Parse Cell_H to find FAB files ----
    with open(cell_h) as fh:
        content = fh.read()

    patches = []
    for line in content.splitlines():
        line = line.strip()
        if not line or line.startswith("FAB") or (line and line[0].isdigit()):
            continue
        parts = line.split()
        if len(parts) >= 2 and "Cell" in parts[0]:
            try:
                patches.append((parts[0], int(parts[1])))
            except (ValueError, IndexError):
                pass

    import struct as _struct
    import numpy as _np

    ic = varnames.index(field)
    result: Dict[Tuple[float, float], float] = {}
    dx = (probhi[0] - problo[0]) / nx
    dy = (probhi[1] - problo[1]) / ny

    fab_dir = os.path.join(plt_path, "Level_0")
    for file_rel, offset in patches:
        fab_path = os.path.join(fab_dir, os.path.basename(file_rel))
        if not os.path.isfile(fab_path):
            fab_path = os.path.join(plt_path, file_rel)
        if not os.path.isfile(fab_path):
            continue
        try:
            with open(fab_path, "rb") as fb:
                fb.seek(offset)
                hdr = b""
                while True:
                    byte = fb.read(1)
                    if not byte:
                        break
                    hdr += byte
                    if byte == b"\n" and b")" in hdr:
                        break
                hdr_str = hdr.decode("ascii", errors="replace")
                hnums = [int(t) for t in hdr_str.replace("(","").replace(")","").replace(","," ").split()
                         if t.lstrip("-").isdigit()]
                if len(hnums) < 6:
                    continue
                fab_ncomp = hnums[0]
                ixlo, iylo, ixhi, iyhi = hnums[2], hnums[3], hnums[4], hnums[5]
                fab_nx = ixhi - ixlo + 1
                fab_ny = iyhi - iylo + 1
                n_vals = fab_ncomp * fab_nx * fab_ny
                raw = fb.read(n_vals * 8)
                if len(raw) < n_vals * 8:
                    fb.seek(offset)
                    while True:
                        b2 = fb.read(1)
                        if not b2:
                            break
                        hdr2 = b""
                        hdr2 += b2
                        if b2 == b"\n" and b")" in hdr2:
                            break
                    raw = fb.read(n_vals * 4)
                    dtype = _np.float32
                else:
                    dtype = _np.float64
                arr = _np.frombuffer(raw, dtype=dtype)
                if arr.size != fab_ncomp * fab_nx * fab_ny:
                    continue
                arr = arr.reshape((fab_ncomp, fab_nx, fab_ny), order="F").transpose(0, 2, 1)
                # Extract component ic for all cells in this FAB
                for j in range(fab_ny):
                    for i in range(fab_nx):
                        x = problo[0] + (ixlo + i + 0.5) * dx
                        y = problo[1] + (iylo + j + 0.5) * dy
                        result[(round(x, 2), round(y, 2))] = float(arr[ic, j, i])
        except Exception:
            continue

    return result


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


def accumulate_flame_length(
    all_fl_data: List[Optional[Dict[Tuple[float, float], float]]],
    all_points:  List[List[Tuple[float, float]]],
    thresholds:  List[float],
    resolution:  float,
) -> Tuple[Dict[float, Dict[Tuple[int, int], int]], float, float]:
    """Accumulate per-cell counts where flame length exceeds each threshold.

    Parameters
    ----------
    all_fl_data : list of per-run flame_length dicts {(x,y)->fl_m} or None.
    all_points  : list of per-run burned (x,y) point lists.
    thresholds  : flame length thresholds [m].
    resolution  : grid cell size [m].

    Returns
    -------
    exceedance_counts : {threshold: {(ix,iy): count}}
    x_min, y_min : grid origin
    """
    xs = [p[0] for run in all_points for p in run]
    ys = [p[1] for run in all_points for p in run]
    if not xs:
        return {t: {} for t in thresholds}, 0.0, 0.0

    x_min, y_min = min(xs), min(ys)
    exceedance: Dict[float, Dict[Tuple[int, int], int]] = {t: {} for t in thresholds}

    for run_pts, fl_dict in zip(all_points, all_fl_data):
        if fl_dict is None:
            continue
        visited: set = set()
        for x, y in run_pts:
            ix = int((x - x_min) / resolution)
            iy = int((y - y_min) / resolution)
            key = (ix, iy)
            if key in visited:
                continue
            visited.add(key)
            # Look up flame length in the plotfile dict (nearest key)
            fl_val = fl_dict.get((round(x, 2), round(y, 2)), 0.0)
            for t in thresholds:
                if fl_val > t:
                    exceedance[t][key] = exceedance[t].get(key, 0) + 1

    return exceedance, x_min, y_min


def write_flame_length_exceedance_csv(
    out_path: str,
    exceedance: Dict[Tuple[int, int], int],
    n_runs: int,
    threshold: float,
    x_origin: float,
    y_origin: float,
    resolution: float,
) -> None:
    """Write P(FL > threshold) CSV: X, Y, P_fl_exceed."""
    with open(out_path, "w") as f:
        f.write(f"# Conditional flame length exceedance  FL > {threshold:.2f} m\n")
        f.write(f"# n_runs = {n_runs}  resolution = {resolution} m\n")
        f.write("X,Y,P_fl_exceed\n")
        for (ix, iy), cnt in sorted(exceedance.items()):
            x = x_origin + (ix + 0.5) * resolution
            y = y_origin + (iy + 0.5) * resolution
            p = cnt / n_runs
            f.write(f"{x:.2f},{y:.2f},{p:.4f}\n")
    print(f"  Wrote P(FL>{threshold:.1f}m) CSV: {out_path}  ({len(exceedance)} cells)")


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
# Cap 7: Fire-area exceedance CCDF
# ---------------------------------------------------------------------------

def compute_run_areas(
    all_points: List[List[Tuple[float, float]]],
    resolution: float,
) -> List[float]:
    """Compute total burned area [ha] for each ensemble run.

    Each burned cell is counted as ``resolution² m²`` and converted to
    hectares.  A cell is counted once per run even if multiple burned points
    map to the same grid cell.

    Parameters
    ----------
    all_points : list[list[tuple[float, float]]]
        Per-run lists of burned (x, y) point coordinates.
    resolution : float
        Grid cell size [m] used to snap burned points.

    Returns
    -------
    list[float]
        Total burned area [ha] for each successful ensemble run,
        in the same order as ``all_points``.
    """
    cell_area_ha = resolution * resolution / 1.0e4
    areas = []
    for run_pts in all_points:
        cells: Set[Tuple[int, int]] = set()
        for x, y in run_pts:
            ix = int(x / resolution)
            iy = int(y / resolution)
            cells.add((ix, iy))
        areas.append(len(cells) * cell_area_ha)
    return areas


def compute_area_exceedance(areas: List[float], n_points: int = 200) -> List[Tuple[float, float]]:
    """Compute the complementary CDF (exceedance curve) of fire area.

    Returns a list of ``(area_ha, P_exceed)`` pairs where
    ``P_exceed = fraction of runs with burned area ≥ area_ha``.
    The curve is sampled at *n_points* linearly spaced area values from
    0 to the maximum observed area.

    Parameters
    ----------
    areas : list[float]
        Burned area [ha] for each ensemble run.
    n_points : int
        Number of points on the exceedance curve (default: 200).

    Returns
    -------
    list[tuple[float, float]]
        List of ``(area_ha, P_exceed)`` pairs, sorted by ascending area.
    """
    if not areas:
        return []
    n = len(areas)
    max_area = max(areas)
    # Sample the CCDF at n_points linearly spaced threshold values
    thresholds = [max_area * i / max(n_points - 1, 1) for i in range(n_points)]
    curve = []
    for thr in thresholds:
        p_exceed = sum(1 for a in areas if a >= thr) / n
        curve.append((round(thr, 4), round(p_exceed, 6)))
    return curve


def write_area_exceedance_csv(
    curve: List[Tuple[float, float]],
    out_path: str,
    n_runs: int,
) -> None:
    """Write the area exceedance curve to a CSV file.

    The output CSV has two columns: ``area_ha`` and ``P_exceed``.

    Parameters
    ----------
    curve : list[tuple[float, float]]
        Output of :func:`compute_area_exceedance`.
    out_path : str
        Destination file path.
    n_runs : int
        Number of ensemble runs (written as a header comment).
    """
    with open(out_path, "w") as f:
        f.write(f"# Fire area exceedance curve  n_runs={n_runs}\n")
        f.write("area_ha,P_exceed\n")
        for area, p in curve:
            f.write(f"{area:.4f},{p:.6f}\n")
    print(f"Wrote area exceedance CSV: {out_path}  ({len(curve)} points)")


def plot_area_exceedance(
    curve: List[Tuple[float, float]],
    out_path: str,
    n_runs: int,
) -> None:
    """Save the area exceedance curve as a PNG plot.

    Requires ``matplotlib``.

    Parameters
    ----------
    curve : list[tuple[float, float]]
        Output of :func:`compute_area_exceedance`.
    out_path : str
        Destination PNG file path.
    n_runs : int
        Number of ensemble runs (used in plot title).
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("WARNING: matplotlib not available – skipping exceedance plot.",
              file=sys.stderr)
        return

    areas = [a for a, _ in curve]
    probs = [p for _, p in curve]

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(areas, probs, linewidth=2, color="firebrick")
    ax.fill_between(areas, probs, alpha=0.15, color="firebrick")
    ax.set_xlabel("Burned Area [ha]")
    ax.set_ylabel("P(burned area ≥ A)")
    ax.set_title(f"Fire Area Exceedance Curve  (n_runs={n_runs})")
    ax.set_xlim(left=0)
    ax.set_ylim(0, 1)
    ax.grid(True, linestyle="--", alpha=0.5)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"Saved area exceedance plot → {out_path}")


# ---------------------------------------------------------------------------
# Single run execution
# ---------------------------------------------------------------------------

def build_launch_cmd(
    exe: str,
    inputs_path: str,
    mpi_ranks: int,
    mpirun: str,
    mpi_extra_args: List[str],
) -> List[str]:
    """Build the full command list for launching the solver.

    Serial (mpi_ranks == 0):
        [exe, inputs_path]

    MPI (mpi_ranks >= 1):
        [mpirun, *mpi_extra_args, "-n", str(mpi_ranks), exe, inputs_path]
    """
    if mpi_ranks <= 0:
        return [exe, inputs_path]
    return [mpirun, *mpi_extra_args, "-n", str(mpi_ranks), exe, inputs_path]


def run_member(
    run_id: int,
    exe: str,
    template_lines: List[str],
    sample: Dict[str, float],
    work_dir: str,
    keep: bool,
    mpi_ranks: int = 0,
    mpirun: str = "mpirun",
    mpi_extra_args: Optional[List[str]] = None,
    csv_files: Optional[List[str]] = None,
    read_flame_length: bool = False,
) -> Optional[Tuple[List[Tuple[float, float]], Optional[Dict[Tuple[float, float], float]]]]:
    """Execute one ensemble member.

    Returns (burned_points, flame_length_dict) or None on failure.

    Parameters
    ----------
    run_id : int
        Ensemble member index.
    exe : str
        Path to the solver executable.
    template_lines : list[str]
        Lines of the template inputs.i file.
    sample : dict
        Perturbed parameter values for this member.
    work_dir : str
        Base directory for per-run scratch directories.
    keep : bool
        If False, plotfile and checkpoint sub-directories are removed after
        the run to save disk space (phi_negative .dat files are kept).
    mpi_ranks : int, optional
        Number of MPI ranks per member.  0 (default) → serial.
    mpirun : str, optional
        MPI launcher executable name or path (default: ``"mpirun"``).
    mpi_extra_args : list[str] or None, optional
        Extra arguments inserted between the launcher and ``-n N``.
    csv_files : list[str] or None, optional
        Pre-collected list of CSV file paths to copy into the run directory.
    read_flame_length : bool, optional
        If True, read the ``flame_length`` field from the final plotfile and
        return it as a dict {(x, y): fl_m}.  Requires plot_int > 0 in
        the template inputs file (default: False).
    """
    if mpi_extra_args is None:
        mpi_extra_args = []
    if csv_files is None:
        csv_files = []

    run_dir = os.path.join(work_dir, f"run_{run_id:04d}")
    os.makedirs(run_dir, exist_ok=True)

    # Copy CSV files from the original working directory into the run directory
    # so that the solver can find them when it runs with cwd=run_dir.
    for csv_file in csv_files:
        shutil.copy2(csv_file, run_dir)

    perturbed = make_perturbed_inputs(template_lines, sample)
    inputs_path = os.path.join(run_dir, "inputs.i")
    with open(inputs_path, "w") as f:
        f.writelines(perturbed)

    cmd = build_launch_cmd(exe, inputs_path, mpi_ranks, mpirun, mpi_extra_args)

    log_path = os.path.join(run_dir, "run.log")
    try:
        with open(log_path, "w") as log:
            result = subprocess.run(
                cmd,
                cwd=run_dir,
                stdout=log,
                stderr=subprocess.STDOUT,
                timeout=3600,
            )
        if result.returncode != 0:
            print(f"  WARNING: run {run_id} exited with code {result.returncode}", flush=True)
        pts = read_phi_negative(run_dir)
        fl_dict: Optional[Dict[Tuple[float, float], float]] = None
        if read_flame_length:
            fl_dict = _read_plotfile_field(run_dir, "flame_length")
        launch_info = (f"mpi×{mpi_ranks}" if mpi_ranks > 0 else "serial")
        print(f"  Run {run_id:4d} [{launch_info}]: speed×{sample['speed_factor']:.2f}  "
              f"dir+{sample['dir_offset']:.1f}°  M_d1+{sample['moist_offset']:.3f}  "
              f"→ {len(pts)} burned cells"
              + (f"  FL cells: {len(fl_dict)}" if fl_dict is not None else ""),
              flush=True)
        return pts, fl_dict
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
    # ---- MPI options ----
    parser.add_argument("--mpi-ranks",        type=int, default=0,
                        help="MPI ranks per member; 0 = serial (default: 0)")
    parser.add_argument("--mpirun",           default="mpirun",
                        help="MPI launcher command (default: mpirun)")
    parser.add_argument("--mpi-args",         default="",
                        help="Extra space-separated args inserted after the launcher "
                             "and before '-n N', e.g. '--bind-to core' (default: '')")
    # ---- perturbation ----
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
    # ---- flame length exceedance ----
    parser.add_argument("--fl-thresholds",    nargs="+", type=float, default=None,
                        metavar="M",
                        help="Compute P(FL > M) maps for each flame-length threshold [m]. "
                             "Requires plot_int > 0 in inputs.i so plotfiles are written. "
                             "Example: --fl-thresholds 0.5 1.0 2.0 4.0")
    parser.add_argument("--fl-out-prefix",   default="fl_exceedance",
                        help="Output file prefix for P(FL>M) CSVs "
                             "(default: fl_exceedance; produces fl_exceedance_0.5m.csv etc.)")
    # ---- Cap 7: area exceedance CCDF ----
    parser.add_argument("--area-exceedance",  action="store_true",
                        help="Compute and write a fire-area exceedance curve "
                             "P(burned_area >= A) as a CCDF across ensemble runs.")
    parser.add_argument("--area-exceedance-out", default="area_exceedance.csv",
                        help="Output CSV path for area exceedance curve "
                             "(default: area_exceedance.csv)")
    parser.add_argument("--area-exceedance-plot", action="store_true",
                        help="Save area exceedance curve as a PNG plot "
                             "(requires matplotlib)")
    parser.add_argument("--area-exceedance-plot-out", default="area_exceedance.png",
                        help="Output PNG path for area exceedance plot "
                             "(default: area_exceedance.png)")

    args = parser.parse_args(argv)

    fl_thresholds: List[float] = sorted(args.fl_thresholds) if args.fl_thresholds else []
    do_fl = len(fl_thresholds) > 0

    # Validate
    if not os.path.isfile(args.inputs):
        print(f"ERROR: template inputs file not found: {args.inputs}", file=sys.stderr)
        sys.exit(1)
    if not os.path.isfile(args.exe):
        print(f"ERROR: solver executable not found: {args.exe}", file=sys.stderr)
        sys.exit(1)

    # Validate and resolve MPI launcher when MPI is requested
    mpi_ranks = args.mpi_ranks
    mpirun_cmd = args.mpirun
    mpi_extra: List[str] = args.mpi_args.split() if args.mpi_args.strip() else []

    if mpi_ranks > 0:
        resolved = shutil.which(mpirun_cmd)
        if resolved is None:
            print(f"ERROR: MPI launcher '{mpirun_cmd}' not found on PATH. "
                  "Install an MPI implementation or use --mpirun to specify the path.",
                  file=sys.stderr)
            sys.exit(1)
        mpirun_cmd = resolved
        print(f"MPI launcher: {mpirun_cmd}  ranks per member: {mpi_ranks}")
        if mpi_extra:
            print(f"MPI extra args: {mpi_extra}")
    else:
        print("Launch mode: serial (use --mpi-ranks N for MPI)")

    if do_fl:
        print(f"Flame length exceedance thresholds [m]: {fl_thresholds}")
    if args.area_exceedance:
        print(f"Area exceedance CCDF: enabled → {args.area_exceedance_out}")

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

    # Collect CSV files from the current directory once; they will be copied
    # into each per-member run directory before the solver is launched.
    csv_files: List[str] = glob.glob(os.path.join(os.getcwd(), "*.csv"))

    all_points: List[List[Tuple[float, float]]] = []
    all_fl_data: List[Optional[Dict[Tuple[float, float], float]]] = []
    failed = 0

    if args.jobs > 1:
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            futures = {
                pool.submit(run_member, i, args.exe, template_lines,
                            samples[i], work_dir, args.keep_runs,
                            mpi_ranks, mpirun_cmd, mpi_extra, csv_files, do_fl): i
                for i in range(args.n_runs)
            }
            for fut in as_completed(futures):
                result = fut.result()
                if result is not None:
                    pts, fl_dict = result
                    all_points.append(pts)
                    all_fl_data.append(fl_dict)
                else:
                    failed += 1
    else:
        for i, sample in enumerate(samples):
            result = run_member(i, args.exe, template_lines, sample,
                                work_dir, args.keep_runs,
                                mpi_ranks, mpirun_cmd, mpi_extra, csv_files, do_fl)
            if result is not None:
                pts, fl_dict = result
                all_points.append(pts)
                all_fl_data.append(fl_dict)
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

    # ---- Flame length exceedance maps ----
    if do_fl and any(fl is not None for fl in all_fl_data):
        print(f"\nComputing flame length exceedance maps for {len(fl_thresholds)} threshold(s) ...")
        exceedance, xe0, ye0 = accumulate_flame_length(
            all_fl_data, all_points, fl_thresholds, resolution
        )
        for t in fl_thresholds:
            exc_counts = exceedance[t]
            if not exc_counts:
                print(f"  No cells exceed FL > {t:.1f} m")
                continue
            suffix = f"{t:.1f}m".replace(".", "p")
            out_fl = f"{args.fl_out_prefix}_{suffix}.csv"
            write_flame_length_exceedance_csv(
                out_fl, exc_counts, n_success, t, xe0, ye0, resolution
            )
    elif do_fl:
        print("WARNING: --fl-thresholds requested but no flame_length data was read "
              "(ensure plot_int > 0 in inputs.i).", file=sys.stderr)

    # ---- Cap 7: Area exceedance CCDF ----
    if args.area_exceedance:
        print("\nComputing fire-area exceedance curve ...")
        run_areas = compute_run_areas(all_points, resolution)
        print(f"  Area range: {min(run_areas):.2f} – {max(run_areas):.2f} ha"
              f"  (mean={sum(run_areas)/len(run_areas):.2f} ha)")
        curve = compute_area_exceedance(run_areas)
        write_area_exceedance_csv(curve, args.area_exceedance_out, n_success)
        if args.area_exceedance_plot:
            plot_area_exceedance(curve, args.area_exceedance_plot_out, n_success)

    if not args.keep_runs and args.work_dir is None:
        # Clean up the auto-created temp directory
        shutil.rmtree(work_dir, ignore_errors=True)
        print(f"Cleaned up working directory: {work_dir}")


if __name__ == "__main__":
    main()
