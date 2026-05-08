#!/usr/bin/env python3
"""
topographic_horizon_analysis.py
================================
Standalone Python utility for analyzing and visualizing the FARSITE
topographic horizon scan implemented in ``src/solar_radiation.H``.

This script serves both as **executable documentation** of the horizon-scan
algorithm and as a practical analysis tool.  It can:

1. **Compute horizon angles** from an XYZ terrain CSV file (same format as
   ``rothermel.terrain_file``) or a synthetic canyon profile.
2. **Diagnose canyon shading** — compare the solar elevation with the per-cell
   topographic horizon to determine which cells are in ridge shadow.
3. **Plot results** — optional heatmaps of horizon angles and shade fraction
   (requires *matplotlib*; gracefully skipped when unavailable).

Algorithm (matches ``compute_topographic_horizon_angles`` in C++)
-----------------------------------------------------------------
For each cell (i, j) and each of the 8 FARSITE compass directions
(N, NE, E, SE, S, SW, W, NW at 45° intervals, clockwise from North):

    θ_horizon[d](i, j) = max_{s=1..S}  atan2( z(i+s·di, j+s·dj) − z(i, j),
                                               s · step_dist[d] )

where step_dist is the physical cell spacing in that direction (√(dx²+dy²) for
diagonals).  The scan stops at the domain boundary or at *max_scan_dist_m*.

For shading at solar azimuth φ_s:
    d_lo  = floor(φ_s_deg / 45)  mod 8
    d_hi  = (d_lo + 1) mod 8
    frac  = (φ_s_deg − d_lo · 45) / 45
    θ_hz  = θ_horizon[d_lo] · (1 − frac) + θ_horizon[d_hi] · frac
    shaded = (solar_elevation < θ_hz)

Usage
-----
::

    # Analyze a real terrain file
    python topographic_horizon_analysis.py terrain.csv \\
        --nx 64 --ny 64 \\
        --solar-elevation 10.0 --solar-azimuth 90.0 \\
        --max-scan-dist 800

    # Generate and analyze a synthetic canyon (no terrain file needed)
    python topographic_horizon_analysis.py --synthetic \\
        --solar-elevation 10.0 --solar-azimuth 90.0

    # Print only the fraction of shaded cells (useful in scripts)
    python topographic_horizon_analysis.py --synthetic --quiet

References
----------
- Finney, M.A. (2004). FARSITE: Fire Area Simulator. USDA FS RMRS-RP-4 (Rev.).
- Andrews, P.L. (2018). Rothermel surface fire spread model. USDA FS RMRS-GTR-371.
"""

from __future__ import annotations

import argparse
import math
import sys
from typing import List, Tuple

# ---------------------------------------------------------------------------
# Type aliases
# ---------------------------------------------------------------------------
Grid2D = List[List[float]]   # [row][col] = elevation [m]


# ---------------------------------------------------------------------------
# Constants: 8 FARSITE compass directions (CW from N)
# ---------------------------------------------------------------------------
_DIR_LABELS = ["N", "NE", "E", "SE", "S", "SW", "W", "NW"]
_DIR_DI     = [ 0,   1,   1,   1,   0,  -1,  -1,  -1]   # column step (+E)
_DIR_DJ     = [ 1,   1,   0,  -1,  -1,  -1,   0,   1]   # row step    (+N)
_DIR_AZ_DEG = [  0,  45,  90, 135, 180, 225, 270, 315]


# ---------------------------------------------------------------------------
# Terrain helpers
# ---------------------------------------------------------------------------

def _read_terrain_csv(path: str) -> Tuple[List[float], List[float], List[float]]:
    """Read an XYZ terrain CSV (comment lines start with '#')."""
    xs, ys, zs = [], [], []
    with open(path) as fh:
        for line in fh:
            line = line.split("#", 1)[0].strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) >= 3:
                xs.append(float(parts[0]))
                ys.append(float(parts[1]))
                zs.append(float(parts[2]))
    return xs, ys, zs


def _xyz_to_grid(xs: List[float], ys: List[float], zs: List[float],
                 nx: int, ny: int) -> Tuple[Grid2D, float, float, float, float]:
    """Bin scattered XYZ points onto a regular nx×ny grid by nearest-point lookup.

    Returns (grid, x_lo, x_hi, y_lo, y_hi) where grid[j][i] is elevation at
    column i, row j.
    """
    x_lo, x_hi = min(xs), max(xs)
    y_lo, y_hi = min(ys), max(ys)
    dx = (x_hi - x_lo) / max(nx - 1, 1)
    dy = (y_hi - y_lo) / max(ny - 1, 1)

    grid: Grid2D = [[0.0] * nx for _ in range(ny)]
    counts       = [[0]   * nx for _ in range(ny)]

    for (x, y, z) in zip(xs, ys, zs):
        col = min(int(round((x - x_lo) / dx)), nx - 1)
        row = min(int(round((y - y_lo) / dy)), ny - 1)
        grid[row][col] += z
        counts[row][col] += 1

    for j in range(ny):
        for i in range(nx):
            if counts[j][i] > 0:
                grid[j][i] /= counts[j][i]

    return grid, x_lo, x_hi, y_lo, y_hi


def _make_synthetic_canyon(nx: int = 32, ny: int = 32,
                            x_lo: float = 0.0, x_hi: float = 1000.0,
                            y_lo: float = 0.0, y_hi: float = 1000.0,
                            ridge_elev: float = 200.0,
                            floor_elev: float = 50.0,
                            ns_gradient: float = 20.0) -> Grid2D:
    """Return a synthetic ridge-canyon-ridge elevation grid.

    The canyon runs N-S; ridges flank the E and W walls.
    Nominal E/W horizon angle from floor ≈ atan2(ridge-floor, domain_width/2).
    """
    cx    = 0.5 * (x_lo + x_hi)
    half_w = 0.5 * (x_hi - x_lo)
    sigma  = 0.35
    grid: Grid2D = []
    for j in range(ny):
        row = []
        norm_y = j / max(ny - 1, 1)
        for i in range(nx):
            x = x_lo + i * (x_hi - x_lo) / max(nx - 1, 1)
            norm_x = abs(x - cx) / half_w
            ridge_frac = math.exp(-((1.0 - norm_x) ** 2) / (2.0 * sigma ** 2))
            z = floor_elev + (ridge_elev - floor_elev) * ridge_frac
            z += ns_gradient * norm_y
            row.append(z)
        grid.append(row)
    return grid


# ---------------------------------------------------------------------------
# Horizon-scan core (mirrors C++ compute_topographic_horizon_angles)
# ---------------------------------------------------------------------------

def compute_horizon_angles(
        grid: Grid2D,
        dx: float,
        dy: float,
        max_scan_dist_m: float = 0.0,
) -> List[List[List[float]]]:
    """Return per-cell 8-direction horizon angles [rad].

    Returns horizon[d][j][i] for direction d, row j, column i.
    Angles in radians; initial floor is −π/2.
    """
    ny = len(grid)
    nx = len(grid[0]) if ny > 0 else 0
    diag = math.sqrt(dx * dx + dy * dy)
    step_dist = [dy, diag, dx, diag, dy, diag, dx, diag]

    max_steps: int
    if max_scan_dist_m > 0.0:
        max_steps = int(math.ceil(max_scan_dist_m / min(dx, dy)))
    else:
        max_steps = max(nx, ny)

    horizon = [[[- math.pi / 2.0] * nx for _ in range(ny)] for _ in range(8)]

    for j in range(ny):
        for i in range(nx):
            z_src = grid[j][i]
            for d in range(8):
                di = _DIR_DI[d]
                dj = _DIR_DJ[d]
                sd = step_dist[d]
                max_angle = -math.pi / 2.0
                for s in range(1, max_steps + 1):
                    it = i + s * di
                    jt = j + s * dj
                    if it < 0 or it >= nx or jt < 0 or jt >= ny:
                        break
                    z_tgt  = grid[jt][it]
                    angle  = math.atan2(z_tgt - z_src, s * sd)
                    if angle > max_angle:
                        max_angle = angle
                horizon[d][j][i] = max_angle

    return horizon


def interpolated_horizon(
        horizon: List[List[List[float]]],
        solar_azimuth_deg: float,
        j: int, i: int,
) -> float:
    """Return the interpolated horizon angle [rad] for the given solar azimuth."""
    az = solar_azimuth_deg % 360.0
    dir_f = az / 45.0
    d_lo  = int(dir_f) % 8
    d_hi  = (d_lo + 1) % 8
    frac  = dir_f - math.floor(dir_f)
    return horizon[d_lo][j][i] * (1.0 - frac) + horizon[d_hi][j][i] * frac


def compute_shade_fraction(
        horizon: List[List[List[float]]],
        solar_elevation_deg: float,
        solar_azimuth_deg: float,
) -> List[List[float]]:
    """Return per-cell shade fraction (0 = insolated, 1 = shaded).

    A cell is shaded when the solar elevation is below the interpolated
    topographic horizon angle in the solar-azimuth direction.
    """
    ny = len(horizon[0])
    nx = len(horizon[0][0]) if ny > 0 else 0
    elev_rad = math.radians(solar_elevation_deg)

    if elev_rad <= 0.0:
        return [[1.0] * nx for _ in range(ny)]   # below horizon globally

    shade: List[List[float]] = []
    for j in range(ny):
        row = []
        for i in range(nx):
            hz = interpolated_horizon(horizon, solar_azimuth_deg, j, i)
            row.append(1.0 if elev_rad < hz else 0.0)
        shade.append(row)
    return shade


# ---------------------------------------------------------------------------
# Statistics helpers
# ---------------------------------------------------------------------------

def _grid_stats(grid_2d: List[List[float]]) -> dict:
    flat = [v for row in grid_2d for v in row]
    return {
        "min": min(flat),
        "max": max(flat),
        "mean": sum(flat) / len(flat),
        "count": len(flat),
    }


# ---------------------------------------------------------------------------
# Optional matplotlib visualization
# ---------------------------------------------------------------------------

def _try_plot(grid: Grid2D, horizon: List[List[List[float]]],
              shade: List[List[float]],
              solar_elevation_deg: float, solar_azimuth_deg: float,
              output_prefix: str = "horizon_analysis") -> None:
    """Produce heatmap figures if matplotlib is available."""
    try:
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        print("INFO: matplotlib not installed – skipping plots.")
        return

    ny = len(grid)
    nx = len(grid[0])

    elev_arr  = np.array(grid)
    shade_arr = np.array(shade)

    # Horizon angle for the solar-azimuth direction
    hz_arr = np.array([
        [interpolated_horizon(horizon, solar_azimuth_deg, j, i)
         for i in range(nx)]
        for j in range(ny)
    ])

    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    fig.suptitle(
        f"Topographic Horizon Analysis\n"
        f"Solar elevation {solar_elevation_deg:.1f}°, "
        f"azimuth {solar_azimuth_deg:.1f}°",
        fontsize=11,
    )

    im0 = axes[0].imshow(elev_arr, origin="lower", cmap="terrain")
    axes[0].set_title("Terrain elevation [m]")
    fig.colorbar(im0, ax=axes[0])

    im1 = axes[1].imshow(np.degrees(hz_arr), origin="lower",
                          cmap="RdYlGn_r", vmin=-10, vmax=40)
    axes[1].set_title(f"Horizon angle [°] (az={solar_azimuth_deg:.0f}°)")
    fig.colorbar(im1, ax=axes[1])

    im2 = axes[2].imshow(shade_arr, origin="lower", cmap="gray_r",
                          vmin=0, vmax=1)
    axes[2].set_title("Shade fraction (1 = ridge-shaded)")
    fig.colorbar(im2, ax=axes[2])

    out_path = f"{output_prefix}.png"
    fig.tight_layout()
    fig.savefig(out_path, dpi=120)
    plt.close(fig)
    print(f"Saved plot → '{out_path}'")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("terrain_csv", nargs="?",
                   help="XYZ terrain CSV file (X Y Z, one point per line). "
                        "Required unless --synthetic is given.")
    p.add_argument("--synthetic", action="store_true",
                   help="Generate and use a synthetic ridge-canyon-ridge terrain.")
    p.add_argument("--nx", type=int, default=32,
                   help="Grid columns when gridding a scattered CSV (default: 32).")
    p.add_argument("--ny", type=int, default=32,
                   help="Grid rows when gridding a scattered CSV (default: 32).")
    p.add_argument("--solar-elevation", type=float, default=10.0, metavar="DEG",
                   help="Solar elevation angle [°] (default: 10.0).")
    p.add_argument("--solar-azimuth", type=float, default=90.0, metavar="DEG",
                   help="Solar azimuth [°, FARSITE CW from N] (default: 90.0 = E).")
    p.add_argument("--max-scan-dist", type=float, default=0.0, metavar="M",
                   help="Maximum ray-march distance [m]; 0 = full domain (default: 0).")
    p.add_argument("--plot", action="store_true",
                   help="Save heatmap PNGs (requires matplotlib).")
    p.add_argument("--plot-prefix", default="horizon_analysis", metavar="PREFIX",
                   help="Output filename prefix for plots (default: horizon_analysis).")
    p.add_argument("--quiet", action="store_true",
                   help="Suppress verbose output; print only the shaded-cell fraction.")
    return p


def main(argv=None):
    args = _build_parser().parse_args(argv)

    # ------------------------------------------------------------------
    # Build or load terrain grid
    # ------------------------------------------------------------------
    if args.synthetic:
        if not args.quiet:
            print("Using synthetic ridge-canyon-ridge terrain …")
        grid = _make_synthetic_canyon(nx=args.nx, ny=args.ny)
        x_lo, x_hi = 0.0, 1000.0
        y_lo, y_hi = 0.0, 1000.0
    elif args.terrain_csv:
        if not args.quiet:
            print(f"Reading terrain from '{args.terrain_csv}' …")
        xs, ys, zs = _read_terrain_csv(args.terrain_csv)
        if not args.quiet:
            print(f"  Loaded {len(xs)} terrain points.")
        grid, x_lo, x_hi, y_lo, y_hi = _xyz_to_grid(xs, ys, zs,
                                                      args.nx, args.ny)
    else:
        print("ERROR: provide a terrain CSV file or use --synthetic.",
              file=sys.stderr)
        sys.exit(1)

    ny = len(grid)
    nx = len(grid[0]) if ny > 0 else 0
    dx = (x_hi - x_lo) / max(nx - 1, 1)
    dy = (y_hi - y_lo) / max(ny - 1, 1)

    elev_stats = _grid_stats(grid)
    if not args.quiet:
        print(f"Grid: {nx}×{ny}  dx={dx:.1f} m  dy={dy:.1f} m")
        print(f"Elevation: min={elev_stats['min']:.1f} m  "
              f"max={elev_stats['max']:.1f} m  "
              f"mean={elev_stats['mean']:.1f} m")

    # ------------------------------------------------------------------
    # Compute horizon angles
    # ------------------------------------------------------------------
    if not args.quiet:
        print(f"\nComputing 8-direction horizon angles "
              f"(max_scan_dist={args.max_scan_dist:.0f} m) …")
    horizon = compute_horizon_angles(grid, dx, dy, args.max_scan_dist)

    if not args.quiet:
        for d, label in enumerate(_DIR_LABELS):
            stats = _grid_stats(horizon[d])
            print(f"  {label:2s} (az={_DIR_AZ_DEG[d]:3d}°): "
                  f"min={math.degrees(stats['min']):+6.1f}°  "
                  f"max={math.degrees(stats['max']):+6.1f}°  "
                  f"mean={math.degrees(stats['mean']):+5.1f}°")

    # ------------------------------------------------------------------
    # Shade fraction at specified solar position
    # ------------------------------------------------------------------
    shade = compute_shade_fraction(
        horizon, args.solar_elevation, args.solar_azimuth
    )
    shade_stats = _grid_stats(shade)
    pct_shaded  = 100.0 * shade_stats["mean"]

    if args.quiet:
        print(f"{pct_shaded:.1f}")
    else:
        print(f"\nShade fraction summary "
              f"(solar el={args.solar_elevation:.1f}°, "
              f"az={args.solar_azimuth:.1f}°):")
        print(f"  Shaded cells: {pct_shaded:.1f}%  "
              f"(1 = ridge-shaded, 0 = insolated)")

        # Horizon angle at domain centre
        jc, ic = ny // 2, nx // 2
        hz_centre = interpolated_horizon(
            horizon, args.solar_azimuth, jc, ic)
        print(f"  Horizon angle at domain centre: "
              f"{math.degrees(hz_centre):.1f}°")
        if args.solar_elevation < math.degrees(hz_centre):
            print(f"  → Centre cell IS in ridge shadow "
                  f"(sun {args.solar_elevation:.1f}° < hz {math.degrees(hz_centre):.1f}°)")
        else:
            print(f"  → Centre cell is NOT in ridge shadow "
                  f"(sun {args.solar_elevation:.1f}° ≥ hz {math.degrees(hz_centre):.1f}°)")

    # ------------------------------------------------------------------
    # Optional plot
    # ------------------------------------------------------------------
    if args.plot:
        _try_plot(grid, horizon, shade,
                  args.solar_elevation, args.solar_azimuth,
                  args.plot_prefix)


if __name__ == "__main__":
    main()
