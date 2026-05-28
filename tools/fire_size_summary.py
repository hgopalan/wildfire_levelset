#!/usr/bin/env python3
"""
fire_size_summary.py – Tabulate and plot fire size statistics over time.

Reads a ``fire_stats.csv`` file written by the wildfire_levelset solver
(enabled with ``fire_stats_file = fire_stats.csv`` in inputs.i) and produces:

  1. An ASCII summary table printed to stdout.
  2. Optional PNG/PDF plots of burned area, perimeter length, and cumulative
     emissions vs. simulation time (requires matplotlib).

The CSV has columns::

  step, time_s, burned_area_ha, perimeter_km,
  active_front_cells, total_co2_kg_m2, total_co_kg_m2, total_pm25_kg_m2

Usage
-----
  # Print table only
  python3 tools/fire_size_summary.py fire_stats.csv

  # Print table and save PNG plots
  python3 tools/fire_size_summary.py fire_stats.csv --plot --outdir fire_plots

  # Save to a different plot format
  python3 tools/fire_size_summary.py fire_stats.csv --plot --fmt pdf \\
      --outdir fire_plots

  # Output a reformatted CSV with time in minutes
  python3 tools/fire_size_summary.py fire_stats.csv --csv summary_min.csv

References
----------
  FARSITE (Finney 1998) produces a similar text summary of burned area and
  perimeter at each time-step output.
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from pathlib import Path
from typing import Dict, List, Optional

try:
    import numpy as np
    _HAS_NUMPY = True
except ImportError:
    _HAS_NUMPY = False

# ---------------------------------------------------------------------------
# CSV reader
# ---------------------------------------------------------------------------

def read_fire_stats(csv_path: str) -> List[Dict[str, float]]:
    """Read fire_stats.csv and return a list of row dicts."""
    rows = []
    with open(csv_path, newline="") as fh:
        # Skip comment lines before handing to DictReader
        lines = [l for l in fh if not l.strip().startswith("#")]
    reader = csv.DictReader(lines)
    for row in reader:
        parsed = {}
        for k, v in row.items():
            if k is None:
                continue
            try:
                parsed[k.strip()] = float(v.strip()) if isinstance(v, str) else float(v)
            except (ValueError, AttributeError):
                parsed[k.strip()] = v
        rows.append(parsed)
    return rows


# ---------------------------------------------------------------------------
# ASCII table
# ---------------------------------------------------------------------------

def print_summary_table(rows: List[Dict[str, float]]) -> None:
    """Print a formatted ASCII summary table to stdout."""
    if not rows:
        print("(no data)")
        return

    cols = [
        ("step",             "Step",           "6.0f"),
        ("time_s",           "Time [s]",       "9.1f"),
        ("time_min",         "Time [min]",     "10.2f"),
        ("burned_area_ha",   "Area [ha]",      "10.2f"),
        ("perimeter_km",     "Perim [km]",     "11.3f"),
        ("active_front_cells","Active cells",  "12.0f"),
        ("growth_rate_ha_min","Growth [ha/m]", "12.3f"),
        ("total_co2_kg_m2",  "CO2 [kg/m2]",   "12.4e"),
        ("total_pm25_kg_m2", "PM2.5[kg/m2]",  "12.4e"),
    ]

    # Build header
    header = "  ".join(f"{label:>{int(fmt.split('.')[0].lstrip('0') or 6)}s}"
                       for _, label, fmt in cols)
    sep = "  ".join("-" * max(len(label), int(fmt.split('.')[0].lstrip('0') or 6))
                    for _, label, fmt in cols)

    print(header)
    print(sep)

    for row in rows:
        parts = []
        for key, _label, fmt in cols:
            if key == "time_min":
                val = row.get("time_s", 0.0) / 60.0
            else:
                val = row.get(key, 0.0)
            try:
                parts.append(f"{val:{fmt}}")
            except (ValueError, TypeError):
                parts.append(f"{str(val):>12s}")
        print("  ".join(parts))

    # Final summary line
    if rows:
        last = rows[-1]
        t_h = last.get("time_s", 0.0) / 3600.0
        area = last.get("burned_area_ha", 0.0)
        perim = last.get("perimeter_km", 0.0)
        major = last.get("major_axis_m", 0.0)
        minor = last.get("minor_axis_m", 0.0)
        print()
        print(f"Simulation duration : {t_h:.2f} h")
        print(f"Final burned area   : {area:.2f} ha  ({area/100.0:.4f} km²)")
        print(f"Final perimeter     : {perim:.3f} km")
        if major > 0 and minor > 0:
            print(f"Fire ellipse axes   : major={major:.1f} m, minor={minor:.1f} m")
            print(f"Ellipse eccentricity: {np.sqrt(1 - (minor/major)**2):.3f}" if _HAS_NUMPY else "")
        print(f"Total time steps    : {len(rows)}")
        
        # Print percentile statistics if numpy is available
        if _HAS_NUMPY:
            print()
            print_percentile_stats(rows)


# ---------------------------------------------------------------------------
# Percentile Statistics
# ---------------------------------------------------------------------------

def print_percentile_stats(rows: List[Dict[str, float]]) -> None:
    """Print percentile statistics for fire growth metrics."""
    if not _HAS_NUMPY or len(rows) < 2:
        return
    
    print("=" * 78)
    print("PERCENTILE STATISTICS (Fire Growth Distribution)")
    print("=" * 78)
    
    # Extract key metrics
    areas = np.array([r.get("burned_area_ha", 0.0) for r in rows])
    perims = np.array([r.get("perimeter_km", 0.0) for r in rows])
    active = np.array([r.get("active_front_cells", 0.0) for r in rows])
    head_ros = np.array([r.get("head_ros_ms", 0.0) for r in rows])
    growth_rate = np.array([r.get("growth_rate_ha_min", 0.0) for r in rows])
    major_axis = np.array([r.get("major_axis_m", 0.0) for r in rows])
    minor_axis = np.array([r.get("minor_axis_m", 0.0) for r in rows])
    
    metrics = [
        ("Burned Area [ha]", areas),
        ("Perimeter [km]", perims),
        ("Active Front Cells", active),
        ("Head ROS [m/s]", head_ros),
        ("Growth Rate [ha/min]", growth_rate),
        ("Major Axis [m]", major_axis),
        ("Minor Axis [m]", minor_axis),
    ]
    
    print(f"{'Metric':<30} {'10th %':>12} {'50th % (Median)':>18} {'90th %':>12}")
    print("-" * 80)
    
    for label, data in metrics:
        # Skip if all data is zero
        if np.all(data == 0.0):
            continue
            
        # For growth rate, handle negative values (shrinking phase)
        if label == "Growth Rate [ha/min]" and np.any(data != 0):
            # Show all percentiles including negative ones
            p10 = np.percentile(data, 10)
            p50 = np.percentile(data, 50)
            p90 = np.percentile(data, 90)
            print(f"{label:<30} {p10:>12.4f} {p50:>18.4f} {p90:>12.4f}")
        elif np.any(data > 0):
            p10 = np.percentile(data, 10)
            p50 = np.percentile(data, 50)
            p90 = np.percentile(data, 90)
            print(f"{label:<30} {p10:>12.4f} {p50:>18.4f} {p90:>12.4f}")
    
    print()


# ---------------------------------------------------------------------------
# Optional matplotlib plots
# ---------------------------------------------------------------------------

def make_plots(rows: List[Dict[str, float]], outdir: Path, fmt: str) -> None:
    """Save matplotlib plots of fire size vs. time."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("WARNING: matplotlib not available – skipping plots.", file=sys.stderr)
        return

    outdir.mkdir(parents=True, exist_ok=True)

    times_min = [r.get("time_s", 0.0) / 60.0 for r in rows]

    # ---- Plot 1: Burned area and perimeter ----
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6), sharex=True)

    areas = [r.get("burned_area_ha", 0.0) for r in rows]
    perims = [r.get("perimeter_km", 0.0) for r in rows]

    ax1.plot(times_min, areas, "r-o", ms=3, label="Burned area")
    ax1.set_ylabel("Burned area [ha]")
    ax1.legend(loc="upper left")
    ax1.grid(True, alpha=0.3)

    ax2.plot(times_min, perims, "b-s", ms=3, label="Perimeter")
    ax2.set_ylabel("Perimeter [km]")
    ax2.set_xlabel("Simulation time [min]")
    ax2.legend(loc="upper left")
    ax2.grid(True, alpha=0.3)

    fig.suptitle("Fire size evolution")
    fig.tight_layout()
    out1 = outdir / f"fire_area_perimeter.{fmt}"
    fig.savefig(out1, dpi=150)
    plt.close(fig)
    print(f"Saved → {out1}")

    # ---- Plot 2: Emissions ----
    co2   = [r.get("total_co2_kg_m2", 0.0) for r in rows]
    co    = [r.get("total_co_kg_m2", 0.0)  for r in rows]
    pm25  = [r.get("total_pm25_kg_m2", 0.0) for r in rows]

    if any(v > 0 for v in co2):
        fig2, ax = plt.subplots(figsize=(8, 4))
        ax.plot(times_min, co2,  "r-",  label="CO₂ [kg/m²]")
        ax.plot(times_min, co,   "b--", label="CO [kg/m²]")
        ax.plot(times_min, pm25, "g:",  label="PM₂.₅ [kg/m²]")
        ax.set_xlabel("Simulation time [min]")
        ax.set_ylabel("Mean domain emissions [kg/m²]")
        ax.set_title("Fire emissions over time")
        ax.legend()
        ax.grid(True, alpha=0.3)
        fig2.tight_layout()
        out2 = outdir / f"fire_emissions.{fmt}"
        fig2.savefig(out2, dpi=150)
        plt.close(fig2)
        print(f"Saved → {out2}")

    # ---- Plot 3: Active front cells ----
    active = [r.get("active_front_cells", 0.0) for r in rows]
    if any(v > 0 for v in active):
        fig3, ax = plt.subplots(figsize=(8, 3))
        ax.fill_between(times_min, active, color="orange", alpha=0.5, label="Active front cells")
        ax.plot(times_min, active, "k-", lw=0.8)
        ax.set_xlabel("Simulation time [min]")
        ax.set_ylabel("Cell count")
        ax.set_title("Active fire front cells")
        ax.legend()
        ax.grid(True, alpha=0.3)
        fig3.tight_layout()
        out3 = outdir / f"fire_front_cells.{fmt}"
        fig3.savefig(out3, dpi=150)
        plt.close(fig3)
        print(f"Saved → {out3}")
    
    # ---- Plot 1b: Fire ellipse axes (if available) ----
    major_axes = [r.get("major_axis_m", 0.0) for r in rows]
    minor_axes = [r.get("minor_axis_m", 0.0) for r in rows]
    
    if any(v > 0 for v in major_axes):
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.plot(times_min, major_axes, "r-o", ms=3, label="Major axis")
        ax.plot(times_min, minor_axes, "b-s", ms=3, label="Minor axis")
        ax.set_xlabel("Simulation time [min]")
        ax.set_ylabel("Axis length [m]")
        ax.set_title("Fire ellipse axes evolution")
        ax.legend()
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        out1b = outdir / f"fire_ellipse_axes.{fmt}"
        fig.savefig(out1b, dpi=150)
        plt.close(fig)
        print(f"Saved → {out1b}")
    
    # ---- Plot 4: Fire growth rate ----
    growth_rate = [r.get("growth_rate_ha_min", 0.0) for r in rows]
    if any(v != 0 for v in growth_rate):
        fig4, ax = plt.subplots(figsize=(8, 4))
        ax.plot(times_min, growth_rate, "g-^", ms=3, label="Growth rate")
        ax.axhline(y=0, color='k', linestyle='--', lw=0.5, alpha=0.5)
        ax.set_xlabel("Simulation time [min]")
        ax.set_ylabel("Growth rate [ha/min]")
        ax.set_title("Fire growth rate over time")
        ax.legend()
        ax.grid(True, alpha=0.3)
        fig4.tight_layout()
        out4 = outdir / f"fire_growth_rate.{fmt}"
        fig4.savefig(out4, dpi=150)
        plt.close(fig4)
        print(f"Saved → {out4}")


# ---------------------------------------------------------------------------
# Optional CSV export
# ---------------------------------------------------------------------------

def export_csv(rows: List[Dict[str, float]], out_path: str) -> None:
    """Write enriched CSV with time in minutes added."""
    fieldnames = list(rows[0].keys()) + ["time_min"] if rows else []
    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            row_out = dict(row)
            row_out["time_min"] = round(row.get("time_s", 0.0) / 60.0, 4)
            writer.writerow(row_out)
    print(f"Wrote enriched CSV → {out_path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Tabulate and plot wildfire_levelset fire size statistics.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "csv",
        help="fire_stats.csv file written by the solver.",
    )
    parser.add_argument(
        "--plot",
        action="store_true",
        help="Save matplotlib PNG/PDF plots of fire size and emissions.",
    )
    parser.add_argument(
        "--outdir",
        default="fire_plots",
        metavar="DIR",
        help="Output directory for plots (default: fire_plots).",
    )
    parser.add_argument(
        "--fmt",
        default="png",
        choices=["png", "pdf", "svg"],
        help="Plot file format (default: png).",
    )
    parser.add_argument(
        "--csv",
        dest="out_csv",
        default=None,
        metavar="FILE",
        help="Save an enriched CSV with time_min column (optional).",
    )

    args = parser.parse_args(argv)

    if not os.path.isfile(args.csv):
        print(f"ERROR: file not found: {args.csv}", file=sys.stderr)
        sys.exit(1)

    rows = read_fire_stats(args.csv)
    if not rows:
        print("No data rows found in CSV.", file=sys.stderr)
        sys.exit(1)

    print_summary_table(rows)

    if args.plot:
        make_plots(rows, Path(args.outdir), args.fmt)

    if args.out_csv:
        export_csv(rows, args.out_csv)


if __name__ == "__main__":
    main()
