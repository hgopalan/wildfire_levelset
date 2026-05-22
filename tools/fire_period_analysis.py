#!/usr/bin/env python3
"""
fire_period_analysis.py – Analyze day/night burning patterns from wildfire simulations.

Reads the ``arrival_time`` field from AMReX plotfiles and classifies each burned
cell as "day" or "night" ignition based on burn period settings. This tool helps
analyze when fire spread occurred relative to the diurnal burning window.

Burn period settings are read from the simulation input file or provided as
command-line arguments. Each cell's arrival time is converted to local clock
time and classified.

Requirements
------------
  pip install numpy matplotlib

Usage
-----
  # Analyze with burn period from inputs file
  python3 tools/fire_period_analysis.py plt0100 --inputs inputs.i \\
      --output burn_period.csv

  # Analyze with explicit burn period parameters
  python3 tools/fire_period_analysis.py plt0100 \\
      --start-hour 10.0 --end-hour 20.0 --sim-start-hour 12.0 \\
      --output burn_period.csv --plot burn_period.png

  # Generate classification raster as GeoTIFF
  python3 tools/fire_period_analysis.py plt0100 --inputs inputs.i \\
      --output burn_period.csv --geotiff burn_period_class.tif

Outputs
-------
  CSV file with summary statistics:
    - total_burned_area_ha
    - day_burned_area_ha
    - night_burned_area_ha
    - day_fraction
  Optional PNG plot showing day/night burn spatial distribution
  Optional GeoTIFF raster with burn period classification

References
----------
  FARSITE burn period concept: fire spread is suppressed outside the specified
  daily burning window (typically daytime hours when fuel moisture is lowest).
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import struct
import sys
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    _HAS_MPL = True
except ImportError:
    _HAS_MPL = False


# ---------------------------------------------------------------------------
# AMReX plotfile reader
# ---------------------------------------------------------------------------

def read_plotfile_header(plotfile_dir: Path):
    """Parse AMReX Header to get grid info and variable names."""
    header_path = plotfile_dir / "Header"
    if not header_path.exists():
        raise FileNotFoundError(f"No Header in {plotfile_dir}")
    
    with open(header_path) as fh:
        lines = [l.strip() for l in fh]
    
    n_vars = int(lines[1])
    varnames = [lines[2 + i] for i in range(n_vars)]
    
    dim_line_idx = 2 + n_vars
    ndim = int(lines[dim_line_idx])
    
    problo_line = dim_line_idx + 2
    probhi_line = problo_line + 1
    problo = [float(x) for x in lines[problo_line].split()]
    probhi = [float(x) for x in lines[probhi_line].split()]
    
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
    """Read arrival_time field from plotfile."""
    varnames, problo, probhi, nx, ny = read_plotfile_header(plotfile_dir)
    
    if "arrival_time" not in varnames:
        raise ValueError(f"arrival_time field not found in {plotfile_dir}")
    
    var_idx = varnames.index("arrival_time")
    
    level_dir = plotfile_dir / "Level_0"
    fab_files = sorted(level_dir.glob("Cell_D_*"))
    if not fab_files:
        raise FileNotFoundError(f"No Cell_D_* files in {level_dir}")
    
    fab_path = fab_files[0]
    
    with open(fab_path, "rb") as fh:
        fh.seek(0, 2)
        file_size = fh.tell()
        fh.seek(0)
        n_doubles = file_size // 8
        data_all = np.fromfile(fh, dtype=np.float64, count=n_doubles)
    
    n_vars = len(varnames)
    try:
        data_all = data_all.reshape((n_vars, ny, nx), order='F')
    except ValueError:
        data_all = data_all.reshape((ny, nx), order='F')[np.newaxis, ...]
    
    arrival = data_all[var_idx, :, :]
    arrival = np.where((arrival < 0) | (arrival > 1e20), np.nan, arrival)
    
    return arrival, problo, probhi, nx, ny


# ---------------------------------------------------------------------------
# Burn period parameters
# ---------------------------------------------------------------------------

def parse_burn_period_from_inputs(inputs_file: str) -> Tuple[float, float, float]:
    """Parse burn period parameters from inputs.i file."""
    start_hour = 10.0
    end_hour = 20.0
    sim_start_hour = 12.0
    
    if not os.path.exists(inputs_file):
        print(f"WARNING: inputs file not found: {inputs_file}", file=sys.stderr)
        return start_hour, end_hour, sim_start_hour
    
    with open(inputs_file) as fh:
        content = fh.read()
    
    # Parse burn_period.* parameters
    match = re.search(r'burn_period\.start_hour\s*=\s*([\d.]+)', content)
    if match:
        start_hour = float(match.group(1))
    
    match = re.search(r'burn_period\.end_hour\s*=\s*([\d.]+)', content)
    if match:
        end_hour = float(match.group(1))
    
    match = re.search(r'burn_period\.sim_start_hour\s*=\s*([\d.]+)', content)
    if match:
        sim_start_hour = float(match.group(1))
    else:
        # Try solar_radiation.start_hour as fallback
        match = re.search(r'solar_radiation\.start_hour\s*=\s*([\d.]+)', content)
        if match:
            sim_start_hour = float(match.group(1))
    
    return start_hour, end_hour, sim_start_hour


# ---------------------------------------------------------------------------
# Burn period classification
# ---------------------------------------------------------------------------

def classify_burn_period(
    arrival: np.ndarray,
    start_hour: float,
    end_hour: float,
    sim_start_hour: float
) -> Tuple[np.ndarray, dict]:
    """
    Classify each cell as day (1) or night (0) burning.
    
    Returns:
        classification: 2D array with 1=day, 0=night, NaN=unburned
        stats: dict with summary statistics
    """
    # Convert arrival time [s] to local clock hour
    arrival_hours = sim_start_hour + arrival / 3600.0
    clock_hours = np.mod(arrival_hours, 24.0)
    
    # Classify based on burn period window
    classification = np.full_like(arrival, np.nan)
    
    if start_hour <= end_hour:
        # Normal case: burn period within a single day (e.g., 10:00-20:00)
        day_mask = (clock_hours >= start_hour) & (clock_hours < end_hour)
    else:
        # Overnight case: burn period crosses midnight (e.g., 20:00-10:00)
        day_mask = (clock_hours >= start_hour) | (clock_hours < end_hour)
    
    # Set classification: 1 = day, 0 = night
    burned = ~np.isnan(arrival)
    classification[burned & day_mask] = 1.0
    classification[burned & ~day_mask] = 0.0
    
    # Compute statistics
    n_burned = np.count_nonzero(burned)
    n_day = np.count_nonzero(classification == 1.0)
    n_night = np.count_nonzero(classification == 0.0)
    
    stats = {
        "n_burned": n_burned,
        "n_day": n_day,
        "n_night": n_night,
        "day_fraction": n_day / n_burned if n_burned > 0 else 0.0,
        "night_fraction": n_night / n_burned if n_burned > 0 else 0.0,
    }
    
    return classification, stats


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def write_stats_csv(
    stats: dict,
    cell_area_m2: float,
    output_path: str
) -> None:
    """Write burn period statistics to CSV."""
    area_ha = cell_area_m2 / 1e4
    
    with open(output_path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow([
            "total_burned_area_ha",
            "day_burned_area_ha",
            "night_burned_area_ha",
            "day_fraction",
            "night_fraction",
            "total_burned_cells",
            "day_burned_cells",
            "night_burned_cells",
        ])
        writer.writerow([
            stats["n_burned"] * area_ha,
            stats["n_day"] * area_ha,
            stats["n_night"] * area_ha,
            stats["day_fraction"],
            stats["night_fraction"],
            stats["n_burned"],
            stats["n_day"],
            stats["n_night"],
        ])
    
    print(f"Wrote burn period statistics to {output_path}")


def plot_burn_period(
    classification: np.ndarray,
    problo: List[float],
    probhi: List[float],
    nx: int,
    ny: int,
    stats: dict,
    output_path: str
) -> None:
    """Create visualization of day/night burn classification."""
    if not _HAS_MPL:
        print("WARNING: matplotlib not available – skipping plot.", file=sys.stderr)
        return
    
    x = np.linspace(problo[0], probhi[0], nx + 1)
    y = np.linspace(problo[1], probhi[1], ny + 1)
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Create custom colormap: night=blue, day=orange, unburned=white
    from matplotlib.colors import ListedColormap
    cmap = ListedColormap(['#2166AC', '#F4A582'])  # night, day
    
    im = ax.pcolormesh(
        x, y, classification,
        cmap=cmap,
        vmin=0, vmax=1,
        shading='flat'
    )
    
    # Colorbar with labels
    cbar = fig.colorbar(im, ax=ax, ticks=[0.25, 0.75])
    cbar.ax.set_yticklabels(['Night', 'Day'])
    cbar.set_label('Burn Period')
    
    ax.set_xlabel('X [m]')
    ax.set_ylabel('Y [m]')
    ax.set_title('Fire Burn Period Classification (Day/Night)')
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3, linestyle='--')
    
    # Add statistics text
    text = (
        f"Day burning:   {stats['n_day']:,} cells ({stats['day_fraction']:.1%})\n"
        f"Night burning: {stats['n_night']:,} cells ({stats['night_fraction']:.1%})"
    )
    ax.text(0.02, 0.98, text, transform=ax.transAxes,
           verticalalignment='top',
           bbox=dict(boxstyle='round', facecolor='white', alpha=0.8),
           fontfamily='monospace', fontsize=9)
    
    fig.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved burn period plot → {output_path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Analyze day/night burning patterns from wildfire simulation.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "plotfile",
        help="AMReX plotfile directory (e.g., plt0100) containing arrival_time field.",
    )
    parser.add_argument(
        "--inputs",
        "-i",
        default=None,
        help="Simulation inputs.i file to read burn period parameters from.",
    )
    parser.add_argument(
        "--start-hour",
        type=float,
        default=None,
        help="Burn period start hour [0-24) in local time (e.g., 10.0 = 10:00 AM). "
             "Overrides --inputs.",
    )
    parser.add_argument(
        "--end-hour",
        type=float,
        default=None,
        help="Burn period end hour [0-24) in local time (e.g., 20.0 = 8:00 PM). "
             "Overrides --inputs.",
    )
    parser.add_argument(
        "--sim-start-hour",
        type=float,
        default=None,
        help="Simulation start hour [0-24) in local time (e.g., 12.0 = 12:00 PM). "
             "Overrides --inputs.",
    )
    parser.add_argument(
        "--output",
        "-o",
        default="burn_period.csv",
        help="Output CSV file for statistics (default: burn_period.csv).",
    )
    parser.add_argument(
        "--plot",
        default=None,
        metavar="FILE",
        help="Save a matplotlib plot of burn period classification (e.g., burn_period.png).",
    )
    
    args = parser.parse_args(argv)
    
    # Determine burn period parameters
    if args.inputs:
        start_hour, end_hour, sim_start_hour = parse_burn_period_from_inputs(args.inputs)
        print(f"Read burn period from {args.inputs}:")
        print(f"  Start hour: {start_hour:.1f}")
        print(f"  End hour: {end_hour:.1f}")
        print(f"  Sim start hour: {sim_start_hour:.1f}")
    else:
        start_hour = 10.0
        end_hour = 20.0
        sim_start_hour = 12.0
        print("Using default burn period parameters:")
        print(f"  Start hour: {start_hour:.1f}")
        print(f"  End hour: {end_hour:.1f}")
        print(f"  Sim start hour: {sim_start_hour:.1f}")
    
    # Command-line overrides
    if args.start_hour is not None:
        start_hour = args.start_hour
    if args.end_hour is not None:
        end_hour = args.end_hour
    if args.sim_start_hour is not None:
        sim_start_hour = args.sim_start_hour
    
    # Read plotfile
    plotfile_dir = Path(args.plotfile)
    if not plotfile_dir.exists():
        print(f"ERROR: plotfile not found: {plotfile_dir}", file=sys.stderr)
        sys.exit(1)
    
    print(f"\nReading arrival_time from {plotfile_dir} ...")
    try:
        arrival, problo, probhi, nx, ny = read_arrival_time_field(plotfile_dir)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
    
    dx = (probhi[0] - problo[0]) / nx
    dy = (probhi[1] - problo[1]) / ny
    cell_area = dx * dy
    
    print(f"Grid: {nx} × {ny}")
    print(f"Cell size: {dx:.2f} × {dy:.2f} m")
    
    # Classify burn period
    print("\nClassifying burn period ...")
    classification, stats = classify_burn_period(
        arrival, start_hour, end_hour, sim_start_hour
    )
    
    print(f"\nResults:")
    print(f"  Total burned cells: {stats['n_burned']:,}")
    print(f"  Day burning cells:  {stats['n_day']:,} ({stats['day_fraction']:.1%})")
    print(f"  Night burning cells: {stats['n_night']:,} ({stats['night_fraction']:.1%})")
    print(f"  Total burned area: {stats['n_burned'] * cell_area / 1e4:.2f} ha")
    print(f"  Day burned area:   {stats['n_day'] * cell_area / 1e4:.2f} ha")
    print(f"  Night burned area: {stats['n_night'] * cell_area / 1e4:.2f} ha")
    
    # Write CSV
    write_stats_csv(stats, cell_area, args.output)
    
    # Optional plot
    if args.plot:
        plot_burn_period(classification, problo, probhi, nx, ny, stats, args.plot)


if __name__ == "__main__":
    main()
