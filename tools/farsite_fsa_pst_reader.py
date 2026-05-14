#!/usr/bin/env python3
"""
farsite_fsa_pst_reader.py – Read and inspect FARSITE Fire Spread Atlas (.fsa)
and Post-processing Statistics (.pst) files written by wildfire_levelset.

The solver writes these files when the following inputs are set::

    fsa_file = fire.fsa    # cumulative perimeter archive
    pst_file = fire.pst    # per-step statistics with max intensity/flame-length

FSA file format
---------------
Each snapshot block looks like::

    SNAPSHOT <step> <time_s> <n_points>
    <x1> <y1>
    <x2> <y2>
    ...
    END_SNAPSHOT

PST file format
---------------
CSV with header::

    step,time_h,burned_ha,perimeter_km,active_cells,
    max_fireline_intensity_kW_m,max_flame_length_m,max_spotting_dist_m

This script can:

  1. **Inspect** an FSA file: list all snapshots with step, time, point count.
  2. **Extract** a specific snapshot by step number (print X Y pairs).
  3. **Plot** fire perimeter polygons from the FSA archive (requires matplotlib).
  4. **Summarise** a PST file as an ASCII table.
  5. **Plot** PST time-series charts of intensity and flame length.

Usage examples
--------------
List all FSA snapshots::

    python3 tools/farsite_fsa_pst_reader.py --fsa fire.fsa

Extract snapshot at step 30::

    python3 tools/farsite_fsa_pst_reader.py --fsa fire.fsa --extract-step 30

Plot all perimeter snapshots::

    python3 tools/farsite_fsa_pst_reader.py --fsa fire.fsa --plot-perimeters

Print PST statistics summary::

    python3 tools/farsite_fsa_pst_reader.py --pst fire.pst

Plot PST time-series::

    python3 tools/farsite_fsa_pst_reader.py --pst fire.pst --plot-pst

References
----------
Finney, M.A. (1998). FARSITE: Fire Area Simulator – Model Development and
  Evaluation. USDA Forest Service Research Paper RMRS-RP-4.
"""

from __future__ import annotations

import argparse
import csv
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional, Tuple


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class FsaSnapshot:
    """One perimeter snapshot from a .fsa archive."""
    step: int
    time_s: float
    points: List[Tuple[float, float]] = field(default_factory=list)


@dataclass
class PstRow:
    """One row from a .pst statistics file."""
    step: int
    time_h: float
    burned_ha: float
    perimeter_km: float
    active_cells: int
    max_intensity_kW_m: float
    max_flame_length_m: float
    max_spotting_dist_m: float


# ---------------------------------------------------------------------------
# FSA reader
# ---------------------------------------------------------------------------

def read_fsa(path: str) -> List[FsaSnapshot]:
    """Parse a .fsa Fire Spread Atlas file; returns list of FsaSnapshot."""
    snapshots: List[FsaSnapshot] = []
    p = Path(path)
    if not p.exists():
        print(f"ERROR: file not found: {path}", file=sys.stderr)
        return snapshots

    with open(p) as f:
        current: Optional[FsaSnapshot] = None
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if line.startswith("SNAPSHOT"):
                parts = line.split()
                step = int(parts[1])
                time_s = float(parts[2])
                current = FsaSnapshot(step=step, time_s=time_s)
                snapshots.append(current)
            elif line == "END_SNAPSHOT":
                current = None
            elif current is not None:
                x, y = map(float, line.split())
                current.points.append((x, y))
    return snapshots


# ---------------------------------------------------------------------------
# PST reader
# ---------------------------------------------------------------------------

def read_pst(path: str) -> List[PstRow]:
    """Parse a .pst post-processing statistics file; returns list of PstRow."""
    rows: List[PstRow] = []
    p = Path(path)
    if not p.exists():
        print(f"ERROR: file not found: {path}", file=sys.stderr)
        return rows

    with open(p, newline="") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            # Skip header row (contains non-numeric first token)
            try:
                int(line.split(",")[0])
            except ValueError:
                continue
            parts = [s.strip() for s in line.split(",")]
            if len(parts) < 8:
                continue
            try:
                rows.append(PstRow(
                    step=int(parts[0]),
                    time_h=float(parts[1]),
                    burned_ha=float(parts[2]),
                    perimeter_km=float(parts[3]),
                    active_cells=int(parts[4]),
                    max_intensity_kW_m=float(parts[5]),
                    max_flame_length_m=float(parts[6]),
                    max_spotting_dist_m=float(parts[7]),
                ))
            except (ValueError, IndexError):
                continue
    return rows


# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------

def print_fsa_summary(snapshots: List[FsaSnapshot]) -> None:
    """Print an ASCII table of all FSA snapshots."""
    print(f"{'Step':>6}  {'Time (h)':>10}  {'N points':>10}")
    print("-" * 32)
    for s in snapshots:
        print(f"{s.step:>6}  {s.time_s / 3600.0:>10.3f}  {len(s.points):>10}")
    print(f"\n{len(snapshots)} snapshot(s) total.")


def print_pst_summary(rows: List[PstRow]) -> None:
    """Print an ASCII table of PST statistics."""
    hdr = (f"{'Step':>6}  {'Time(h)':>8}  {'Burned(ha)':>11}  {'Perim(km)':>10}"
           f"  {'MaxFI(kW/m)':>12}  {'MaxFL(m)':>9}  {'MaxSpot(m)':>11}")
    print(hdr)
    print("-" * len(hdr))
    for r in rows:
        print(f"{r.step:>6}  {r.time_h:>8.3f}  {r.burned_ha:>11.2f}  "
              f"{r.perimeter_km:>10.3f}  {r.max_intensity_kW_m:>12.1f}  "
              f"{r.max_flame_length_m:>9.2f}  {r.max_spotting_dist_m:>11.1f}")
    if rows:
        print(f"\n{len(rows)} row(s) | final burned area: {rows[-1].burned_ha:.2f} ha")


def extract_snapshot(snapshots: List[FsaSnapshot], step: int) -> None:
    """Print X Y pairs for the snapshot matching `step`."""
    matches = [s for s in snapshots if s.step == step]
    if not matches:
        print(f"No snapshot found for step {step}.", file=sys.stderr)
        sys.exit(1)
    snap = matches[-1]  # last match if duplicate
    print(f"# Snapshot step={snap.step}  time={snap.time_s:.1f} s")
    print("X,Y")
    for x, y in snap.points:
        print(f"{x:.2f},{y:.2f}")


def plot_perimeters(snapshots: List[FsaSnapshot], outfile: Optional[str] = None) -> None:
    """Plot all perimeter polygons colour-coded by time."""
    try:
        import matplotlib.pyplot as plt
        import matplotlib.cm as cm
        import numpy as np
    except ImportError:
        print("matplotlib is required for --plot-perimeters", file=sys.stderr)
        sys.exit(1)

    fig, ax = plt.subplots(figsize=(8, 8))
    times = [s.time_s / 3600.0 for s in snapshots]
    t_min, t_max = min(times, default=0), max(times, default=1)
    if t_max == t_min:
        t_max = t_min + 1

    cmap = cm.YlOrRd
    for snap in snapshots:
        if not snap.points:
            continue
        xs = [p[0] for p in snap.points] + [snap.points[0][0]]
        ys = [p[1] for p in snap.points] + [snap.points[0][1]]
        t_norm = (snap.time_s / 3600.0 - t_min) / (t_max - t_min)
        ax.plot(xs, ys, color=cmap(t_norm), linewidth=1.0, alpha=0.8)

    sm = plt.cm.ScalarMappable(cmap=cmap,
                                norm=plt.Normalize(vmin=t_min, vmax=t_max))
    sm.set_array([])
    fig.colorbar(sm, ax=ax, label="Time (h)")
    ax.set_xlabel("X [m]")
    ax.set_ylabel("Y [m]")
    ax.set_title(f"Fire Spread Atlas – {len(snapshots)} perimeter snapshots")
    ax.set_aspect("equal")
    plt.tight_layout()
    if outfile:
        fig.savefig(outfile, dpi=150)
        print(f"Saved perimeter plot to {outfile}")
    else:
        plt.show()


def plot_pst(rows: List[PstRow], outfile: Optional[str] = None) -> None:
    """Plot PST time-series for intensity and flame length."""
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib is required for --plot-pst", file=sys.stderr)
        sys.exit(1)

    if not rows:
        print("No PST data to plot.", file=sys.stderr)
        return

    t = [r.time_h for r in rows]
    fi = [r.max_intensity_kW_m for r in rows]
    fl = [r.max_flame_length_m for r in rows]
    ha = [r.burned_ha for r in rows]

    fig, axes = plt.subplots(3, 1, figsize=(8, 9), sharex=True)
    axes[0].plot(t, ha,  color="#e63946", linewidth=1.5)
    axes[0].set_ylabel("Burned area (ha)")
    axes[0].grid(True, alpha=0.4)

    axes[1].plot(t, fi, color="#f4a261", linewidth=1.5)
    axes[1].set_ylabel("Max fireline intensity (kW/m)")
    axes[1].grid(True, alpha=0.4)

    axes[2].plot(t, fl, color="#2a9d8f", linewidth=1.5)
    axes[2].set_ylabel("Max flame length (m)")
    axes[2].set_xlabel("Time (h)")
    axes[2].grid(True, alpha=0.4)

    fig.suptitle("FARSITE Post-processing Statistics")
    plt.tight_layout()
    if outfile:
        fig.savefig(outfile, dpi=150)
        print(f"Saved PST plot to {outfile}")
    else:
        plt.show()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Read and inspect wildfire_levelset .fsa and .pst files.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__)
    p.add_argument("--fsa",  metavar="FILE", help=".fsa Fire Spread Atlas file")
    p.add_argument("--pst",  metavar="FILE", help=".pst statistics file")
    p.add_argument("--extract-step", metavar="N", type=int,
                   help="Extract FSA perimeter snapshot for step N (CSV to stdout)")
    p.add_argument("--plot-perimeters", action="store_true",
                   help="Plot all FSA perimeter polygons (requires matplotlib)")
    p.add_argument("--plot-pst", action="store_true",
                   help="Plot PST time-series charts (requires matplotlib)")
    p.add_argument("--output", "-o", metavar="FILE",
                   help="Output file for plots (PNG/PDF); default: show interactively")
    return p


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    if not args.fsa and not args.pst:
        parser.print_help()
        sys.exit(0)

    if args.fsa:
        snapshots = read_fsa(args.fsa)
        if args.extract_step is not None:
            extract_snapshot(snapshots, args.extract_step)
        elif args.plot_perimeters:
            plot_perimeters(snapshots, outfile=args.output)
        else:
            print_fsa_summary(snapshots)

    if args.pst:
        rows = read_pst(args.pst)
        if args.plot_pst:
            plot_pst(rows, outfile=args.output)
        else:
            print_pst_summary(rows)


if __name__ == "__main__":
    main()
