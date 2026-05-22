#!/usr/bin/env python3
"""
minimum_travel_path.py – Extract minimum travel time paths from wildfire simulations.

Reads the ``arrival_time`` field from AMReX plotfiles and computes the minimum
travel time (MTT) path from the ignition point(s) to a specified destination.
The path follows the steepest descent gradient of the arrival_time field.

This tool is useful for:
  - Understanding fire spread paths
  - Identifying critical control points
  - Analyzing firebreak placement effectiveness
  - Visualizing dominant fire progression routes

Requirements
------------
  pip install numpy matplotlib

Usage
-----
  # Extract path to a destination point (x, y) in meters
  python3 tools/minimum_travel_path.py plt0100 --dest 5000 3000 --output path.csv

  # Extract path with visualization
  python3 tools/minimum_travel_path.py plt0100 --dest 5000 3000 \\
      --output path.csv --plot path.png

  # Multiple paths to different destinations
  python3 tools/minimum_travel_path.py plt0100 --dest 5000 3000 6000 3500 \\
      --output paths.csv --plot paths.png

Outputs
-------
  CSV file with columns: path_id, x, y, arrival_time_s, arrival_time_min
  Optional PNG plot showing the path overlaid on arrival time contours

References
----------
  Finney, M.A. (2002). Fire growth using minimum travel time methods.
  Canadian Journal of Forest Research, 32(8), 1420–1424.
"""

from __future__ import annotations

import argparse
import csv
import os
import struct
import sys
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle
    _HAS_MPL = True
except ImportError:
    _HAS_MPL = False


# ---------------------------------------------------------------------------
# AMReX plotfile reader (minimal inline version)
# ---------------------------------------------------------------------------

def read_plotfile_header(plotfile_dir: Path):
    """Parse AMReX Header to get grid info and variable names."""
    header_path = plotfile_dir / "Header"
    if not header_path.exists():
        raise FileNotFoundError(f"No Header in {plotfile_dir}")
    
    with open(header_path) as fh:
        lines = [l.strip() for l in fh]
    
    # Line 1: version
    # Line 2: number of variables
    n_vars = int(lines[1])
    # Lines 3..3+n_vars-1: variable names
    varnames = [lines[2 + i] for i in range(n_vars)]
    
    # Find coordinate lines
    dim_line_idx = 2 + n_vars
    ndim = int(lines[dim_line_idx])
    
    # Problo / probhi
    problo_line = dim_line_idx + 2
    probhi_line = problo_line + 1
    problo = [float(x) for x in lines[problo_line].split()]
    probhi = [float(x) for x in lines[probhi_line].split()]
    
    # Grid dimensions (find "((0,0) (nx-1, ny-1) (0))" pattern)
    # Look for the box definition
    for i, line in enumerate(lines):
        if line.startswith("(("):
            parts = line.replace("(", "").replace(")", "").split()
            if len(parts) >= 2 * ndim:
                hi_idx = [int(parts[ndim + j]) for j in range(ndim)]
                nx, ny = hi_idx[0] + 1, hi_idx[1] + 1
                break
    else:
        raise ValueError("Could not parse grid dimensions from Header")
    
    return varnames, problo, probhi, nx, ny


def read_arrival_time_field(plotfile_dir: Path):
    """Read arrival_time field from plotfile. Returns (data, problo, probhi, nx, ny)."""
    varnames, problo, probhi, nx, ny = read_plotfile_header(plotfile_dir)
    
    if "arrival_time" not in varnames:
        raise ValueError(f"arrival_time field not found in {plotfile_dir}")
    
    var_idx = varnames.index("arrival_time")
    
    # Read binary data from Level_0/Cell_D_<fab_number>
    level_dir = plotfile_dir / "Level_0"
    if not level_dir.exists():
        raise FileNotFoundError(f"Level_0 directory not found in {plotfile_dir}")
    
    # Assuming single FAB for simplicity (works for single-box grids)
    fab_files = sorted(level_dir.glob("Cell_D_*"))
    if not fab_files:
        raise FileNotFoundError(f"No Cell_D_* files in {level_dir}")
    
    fab_path = fab_files[0]
    
    # Read FAB data (IEEE double precision, Fortran column-major order)
    with open(fab_path, "rb") as fh:
        # FABs start with a header; skip to data
        # Simple approach: read all and reshape
        # Each FAB stores all variables for all cells
        n_cells = nx * ny
        n_vars = len(varnames)
        
        # FAB format: Fortran column-major, all variables interleaved per cell
        # Skip the FAB header (varies; we'll read raw data)
        fh.seek(0, 2)  # Go to end
        file_size = fh.tell()
        fh.seek(0)
        
        # Read all doubles
        n_doubles = file_size // 8
        data_all = np.fromfile(fh, dtype=np.float64, count=n_doubles)
    
    # Reshape: AMReX stores as (var, k, j, i) in Fortran order
    # For 2D: (var, j, i) with i fastest
    # Reshape to (n_vars, ny, nx)
    try:
        data_all = data_all.reshape((n_vars, ny, nx), order='F')
    except ValueError:
        # Try without n_vars separation (single-var FAB)
        data_all = data_all.reshape((ny, nx), order='F')[np.newaxis, ...]
    
    arrival = data_all[var_idx, :, :]
    
    # Replace sentinel values (e.g., -1, 1e30) with NaN
    arrival = np.where((arrival < 0) | (arrival > 1e20), np.nan, arrival)
    
    return arrival, problo, probhi, nx, ny


# ---------------------------------------------------------------------------
# MTT path extraction
# ---------------------------------------------------------------------------

def extract_mtt_path(
    arrival: np.ndarray,
    problo: List[float],
    probhi: List[float],
    nx: int,
    ny: int,
    dest_x: float,
    dest_y: float,
    max_steps: int = 10000
) -> List[Tuple[float, float, float]]:
    """
    Extract MTT path from destination point back to ignition.
    
    Returns list of (x, y, arrival_time) tuples along the path.
    Path is traced backward from destination following steepest descent.
    """
    dx = (probhi[0] - problo[0]) / nx
    dy = (probhi[1] - problo[1]) / ny
    
    # Convert destination to grid indices
    i = int((dest_x - problo[0]) / dx)
    j = int((dest_y - problo[1]) / dy)
    
    # Clamp to grid bounds
    i = max(0, min(nx - 1, i))
    j = max(0, min(ny - 1, j))
    
    path = []
    visited = set()
    
    for step in range(max_steps):
        if (i, j) in visited:
            break  # Prevent infinite loops
        visited.add((i, j))
        
        # Current arrival time
        t = arrival[j, i]
        if np.isnan(t):
            break  # Reached unburned area
        
        # Physical coordinates
        x = problo[0] + (i + 0.5) * dx
        y = problo[1] + (j + 0.5) * dy
        path.append((x, y, float(t)))
        
        # Check if we've reached the ignition (arrival_time ≈ 0)
        if t < 1.0:  # Within 1 second of ignition
            break
        
        # Find steepest descent neighbor (minimum arrival time)
        min_t = t
        next_i, next_j = i, j
        
        for di in [-1, 0, 1]:
            for dj in [-1, 0, 1]:
                if di == 0 and dj == 0:
                    continue
                ni, nj = i + di, j + dj
                if 0 <= ni < nx and 0 <= nj < ny:
                    nt = arrival[nj, ni]
                    if not np.isnan(nt) and nt < min_t:
                        min_t = nt
                        next_i, next_j = ni, nj
        
        if next_i == i and next_j == j:
            # No downhill neighbor found
            break
        
        i, j = next_i, next_j
    
    return path


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def write_paths_csv(
    paths: List[List[Tuple[float, float, float]]],
    output_path: str
) -> None:
    """Write MTT paths to CSV file."""
    with open(output_path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["path_id", "x", "y", "arrival_time_s", "arrival_time_min"])
        
        for path_id, path in enumerate(paths):
            for x, y, t in path:
                writer.writerow([path_id, x, y, t, t / 60.0])
    
    print(f"Wrote {sum(len(p) for p in paths)} path points to {output_path}")


def plot_paths(
    arrival: np.ndarray,
    problo: List[float],
    probhi: List[float],
    nx: int,
    ny: int,
    paths: List[List[Tuple[float, float, float]]],
    destinations: List[Tuple[float, float]],
    output_path: str
) -> None:
    """Create a plot of MTT paths overlaid on arrival time contours."""
    if not _HAS_MPL:
        print("WARNING: matplotlib not available – skipping plot.", file=sys.stderr)
        return
    
    # Create coordinate arrays
    x = np.linspace(problo[0], probhi[0], nx + 1)
    y = np.linspace(problo[1], probhi[1], ny + 1)
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Plot arrival time contours
    valid_mask = ~np.isnan(arrival)
    if np.any(valid_mask):
        levels = np.linspace(np.nanmin(arrival), np.nanmax(arrival), 20)
        contour = ax.contourf(
            x[:-1] + (x[1] - x[0]) / 2,
            y[:-1] + (y[1] - y[0]) / 2,
            arrival / 60.0,  # Convert to minutes
            levels=levels / 60.0,
            cmap='YlOrRd',
            alpha=0.7
        )
        cbar = fig.colorbar(contour, ax=ax, label='Arrival time [min]')
    
    # Plot paths
    colors = plt.cm.tab10(np.linspace(0, 1, len(paths)))
    for path_id, (path, color) in enumerate(zip(paths, colors)):
        if not path:
            continue
        xs, ys, ts = zip(*path)
        ax.plot(xs, ys, 'o-', color=color, linewidth=2, markersize=3,
                label=f'Path {path_id}', zorder=10)
        
        # Mark start (destination) and end (ignition)
        ax.plot(xs[0], ys[0], 'o', color=color, markersize=10,
                markeredgecolor='white', markeredgewidth=2, zorder=11)
        ax.plot(xs[-1], ys[-1], '*', color='yellow', markersize=15,
                markeredgecolor='black', markeredgewidth=1.5, zorder=11)
    
    # Mark destination points
    for dest_x, dest_y in destinations:
        ax.add_patch(Circle((dest_x, dest_y), radius=(probhi[0]-problo[0])/100,
                           color='blue', alpha=0.5, zorder=9))
    
    ax.set_xlabel('X [m]')
    ax.set_ylabel('Y [m]')
    ax.set_title('Minimum Travel Time Paths')
    ax.legend(loc='best')
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)
    print(f"Saved plot → {output_path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Extract minimum travel time paths from wildfire simulation.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "plotfile",
        help="AMReX plotfile directory (e.g., plt0100) containing arrival_time field.",
    )
    parser.add_argument(
        "--dest",
        nargs="+",
        type=float,
        required=True,
        metavar="X Y",
        help="Destination point(s) in physical coordinates [m]. "
             "Provide pairs: --dest x1 y1 x2 y2 ...",
    )
    parser.add_argument(
        "--output",
        "-o",
        default="mtt_path.csv",
        help="Output CSV file path (default: mtt_path.csv).",
    )
    parser.add_argument(
        "--plot",
        default=None,
        metavar="FILE",
        help="Save a plot of the path(s) to FILE (e.g., path.png).",
    )
    parser.add_argument(
        "--max-steps",
        type=int,
        default=10000,
        help="Maximum number of steps when tracing path (default: 10000).",
    )
    
    args = parser.parse_args(argv)
    
    # Parse destination points
    if len(args.dest) % 2 != 0:
        print("ERROR: --dest requires pairs of coordinates (x y x y ...)", file=sys.stderr)
        sys.exit(1)
    
    destinations = [(args.dest[i], args.dest[i+1]) for i in range(0, len(args.dest), 2)]
    
    # Read plotfile
    plotfile_dir = Path(args.plotfile)
    if not plotfile_dir.exists():
        print(f"ERROR: plotfile not found: {plotfile_dir}", file=sys.stderr)
        sys.exit(1)
    
    print(f"Reading arrival_time from {plotfile_dir} ...")
    try:
        arrival, problo, probhi, nx, ny = read_arrival_time_field(plotfile_dir)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
    
    print(f"Grid: {nx} × {ny}")
    print(f"Domain: [{problo[0]:.1f}, {probhi[0]:.1f}] × [{problo[1]:.1f}, {probhi[1]:.1f}] m")
    
    # Extract paths
    paths = []
    for dest_idx, (dest_x, dest_y) in enumerate(destinations):
        print(f"\nExtracting path {dest_idx} to destination ({dest_x:.1f}, {dest_y:.1f}) ...")
        path = extract_mtt_path(arrival, problo, probhi, nx, ny, dest_x, dest_y, args.max_steps)
        paths.append(path)
        
        if path:
            t_start = path[0][2]
            t_end = path[-1][2]
            print(f"  Path length: {len(path)} points")
            print(f"  Arrival time at destination: {t_start/60:.2f} min")
            print(f"  Arrival time at ignition: {t_end/60:.2f} min")
            print(f"  Travel time: {(t_start - t_end)/60:.2f} min")
        else:
            print(f"  WARNING: No path found (destination may be in unburned area)")
    
    # Write CSV
    write_paths_csv(paths, args.output)
    
    # Optional plot
    if args.plot:
        plot_paths(arrival, problo, probhi, nx, ny, paths, destinations, args.plot)


if __name__ == "__main__":
    main()
