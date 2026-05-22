#!/usr/bin/env python3
"""
isochrone_extractor.py – Extract fire arrival-time isochrones from AMReX plotfiles.

Reads the ``arrival_time`` field (simulation time [s] at which each cell first
ignited) from one or more wildfire_levelset AMReX 2-D plotfiles and writes
time-of-arrival isochrone polygons as GeoJSON.

Each isochrone is the boundary of the region that had ignited by time T.  By
default, isochrones are extracted at equally-spaced intervals; custom times can
also be specified.

Requirements
------------
  pip install numpy matplotlib

Usage
-----
  # Isochrones at 600-s intervals from the final plotfile
  python3 tools/isochrone_extractor.py plt0100 --interval 600 --outdir iso_out

  # Custom isochrone times [s]
  python3 tools/isochrone_extractor.py plt0100 --times 300 600 900 1200 \\
      --outdir iso_out

  # Batch: process all plt#### directories, 10-minute isochrones
  python3 tools/isochrone_extractor.py --all --interval 600 --outdir iso_out

  # Georeference with a UTM origin (easting, northing)
  python3 tools/isochrone_extractor.py plt0100 --interval 600 \\
      --utm-origin 450000 4200000 --outdir iso_out

Outputs
-------
  <outdir>/<plotfile>_isochrones.geojson
    GeoJSON FeatureCollection of Polygon features, one per isochrone time.
    Feature properties include:
      time_s      – isochrone time [s]
      time_min    – isochrone time [min]
      label       – human-readable label (e.g. "t=10 min")
      level_index – index in the sorted list of isochrone times

References
----------
  FARSITE (Finney 1998) time-of-arrival isochrones use the same concept:
  each isochrone marks the fire perimeter at a specific elapsed time.
"""

from __future__ import annotations

import argparse
import json
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
    _HAS_MPL = True
except ImportError:
    _HAS_MPL = False


# ---------------------------------------------------------------------------
# Reuse the AMReX plotfile parser from plotfile_to_geotiff.py if available;
# otherwise define a minimal inline version.
# ---------------------------------------------------------------------------

_SCRIPT_DIR = Path(__file__).parent


def _try_import_geotiff_parser():
    """Try to import helpers from plotfile_to_geotiff in the same directory."""
    try:
        sys.path.insert(0, str(_SCRIPT_DIR))
        from plotfile_to_geotiff import _parse_header, _read_fab_data
        return _parse_header, _read_fab_data
    except ImportError:
        return None, None


_parse_header_ext, _read_fab_data_ext = _try_import_geotiff_parser()


# ---------------------------------------------------------------------------
# Minimal inline plotfile reader (fallback if geotiff module unavailable)
# ---------------------------------------------------------------------------

def _parse_header_inline(plotfile_dir: Path):
    """Minimal Header parser.  Returns (varnames, problo, probhi, nx, ny)."""
    header_path = plotfile_dir / "Header"
    if not header_path.exists():
        raise FileNotFoundError(f"No Header in {plotfile_dir}")
    with open(header_path) as fh:
        lines = [l.rstrip("\n") for l in fh]
    idx = 1
    ncomp = int(lines[idx]); idx += 1
    varnames = []
    for _ in range(ncomp):
        varnames.append(lines[idx].strip()); idx += 1
    # spacedim
    spacedim = int(lines[idx]); idx += 1
    if spacedim != 2:
        raise ValueError(f"Only 2-D plotfiles supported (spacedim={spacedim})")
    _time = float(lines[idx]); idx += 1
    finest_level = int(lines[idx]); idx += 1
    problo = list(map(float, lines[idx].split())); idx += 1
    probhi = list(map(float, lines[idx].split())); idx += 1
    idx += 1  # ref ratios
    # domain box for Level_0
    idx += 1  # box-array count
    box_line = lines[idx]; idx += 1
    nums = [int(t) for t in box_line.replace("(","").replace(")","").replace(","," ").split()
            if t.lstrip("-").isdigit()]
    if len(nums) >= 4:
        nx = nums[2] - nums[0] + 1
        ny = nums[3] - nums[1] + 1
    else:
        nx, ny = 64, 64
    return varnames, problo, probhi, nx, ny


def _read_fab_data_inline(cell_h_path: Path, ncomp: int, nx: int, ny: int) -> np.ndarray:
    """Minimal FAB reader.  Returns array (ncomp, ny, nx)."""
    with open(cell_h_path) as fh:
        content = fh.read()

    patches = []
    for line in content.splitlines():
        line = line.strip()
        if not line or line.startswith("FAB") or line[0].isdigit():
            continue
        parts = line.split()
        if len(parts) >= 2 and "Cell" in parts[0]:
            try:
                patches.append((parts[0], int(parts[1])))
            except (ValueError, IndexError):
                pass

    data = np.full((ncomp, ny, nx), np.nan, dtype=np.float64)
    fab_dir = cell_h_path.parent

    for file_rel, offset in patches:
        fab_path = fab_dir / Path(file_rel).name
        if not fab_path.exists():
            fab_path = fab_dir.parent / file_rel
        if not fab_path.exists():
            continue
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
            nums = [int(t) for t in hdr_str.replace("(","").replace(")","").replace(","," ").split()
                    if t.lstrip("-").isdigit()]
            if len(nums) < 6:
                continue
            fab_ncomp = nums[0]
            ixlo, iylo, ixhi, iyhi = nums[2], nums[3], nums[4], nums[5]
            fab_nx = ixhi - ixlo + 1
            fab_ny = iyhi - iylo + 1
            prec = np.float64
            n_vals = fab_ncomp * fab_nx * fab_ny
            raw = fb.read(n_vals * 8)
            if len(raw) < n_vals * 8:
                prec = np.float32
                fb.seek(offset)
                hdr2 = b""
                while True:
                    b2 = fb.read(1)
                    if not b2:
                        break
                    hdr2 += b2
                    if b2 == b"\n" and b")" in hdr2:
                        break
                raw = fb.read(n_vals * 4)
            arr = np.frombuffer(raw, dtype=prec)
            if arr.size == fab_ncomp * fab_nx * fab_ny:
                arr = arr.reshape((fab_ncomp, fab_nx, fab_ny), order="F").transpose(0, 2, 1)
                for ic in range(min(fab_ncomp, ncomp)):
                    data[ic, iylo:iyhi+1, ixlo:ixhi+1] = arr[ic]
    return data


# ---------------------------------------------------------------------------
# Plotfile reader (using imported helpers when available)
# ---------------------------------------------------------------------------

def read_arrival_time(plotfile_dir: Path) -> Tuple[np.ndarray, list, list, int, int]:
    """Read the arrival_time field from an AMReX plotfile.

    Returns
    -------
    arrival : ndarray, shape (ny, nx)  – arrival time [s]; NaN for unburned.
    problo  : [x0, y0]
    probhi  : [x1, y1]
    nx, ny  : domain dimensions
    """
    if _parse_header_ext is not None:
        varnames, problo, probhi, nx, ny, finest_level, ref_ratios, level_domains = \
            _parse_header_ext(plotfile_dir)
    else:
        varnames, problo, probhi, nx, ny = _parse_header_inline(plotfile_dir)

    if "arrival_time" not in varnames:
        raise ValueError(
            f"'arrival_time' not found in {plotfile_dir}.  "
            "Run the solver with plot_int > 0 so the full plotfile is written."
        )

    cell_h = plotfile_dir / "Level_0" / "Cell_H"
    if not cell_h.exists():
        raise FileNotFoundError(f"Level_0/Cell_H not found in {plotfile_dir}")

    ncomp = len(varnames)
    if _read_fab_data_ext is not None:
        data = _read_fab_data_ext(cell_h, ncomp, nx, ny)
    else:
        data = _read_fab_data_inline(cell_h, ncomp, nx, ny)

    ic = varnames.index("arrival_time")
    arrival = data[ic].copy()

    # Sentinel -1.0 (unburned cells) → NaN
    arrival[arrival < 0.0] = np.nan

    return arrival, problo, probhi, nx, ny


# ---------------------------------------------------------------------------
# Isochrone extraction
# ---------------------------------------------------------------------------

def extract_isochrones(
    arrival: np.ndarray,
    problo: list,
    probhi: list,
    nx: int,
    ny: int,
    iso_times: List[float],
    utm_origin: Optional[Tuple[float, float]] = None,
) -> dict:
    """Extract isochrone polygons as a GeoJSON FeatureCollection.

    For each time T in *iso_times*, the isochrone is the boundary of
    ``arrival_time <= T`` — i.e., the fire perimeter at elapsed time T.
    """
    if not _HAS_MPL:
        raise ImportError(
            "matplotlib is required for isochrone extraction.  "
            "Install it with:  pip install matplotlib"
        )

    dx = (probhi[0] - problo[0]) / nx
    dy = (probhi[1] - problo[1]) / ny

    xs = np.array([problo[0] + (i + 0.5) * dx for i in range(nx)])
    ys = np.array([problo[1] + (j + 0.5) * dy for j in range(ny)])

    easting_offset  = utm_origin[0] if utm_origin else 0.0
    northing_offset = utm_origin[1] if utm_origin else 0.0

    features = []
    for level_idx, t in enumerate(sorted(iso_times)):
        # Create a binary mask: 1 where ignited by time t, 0 elsewhere
        # arrival NaN = unburned → treat as "not yet ignited"
        mask = np.where(np.isnan(arrival), 0.0, np.where(arrival <= t, 1.0, 0.0))

        if mask.max() == 0.0:
            continue  # fire has not reached this isochrone time

        fig, ax = plt.subplots()
        cs = ax.contour(xs, ys, mask, levels=[0.5])
        plt.close(fig)

        for path in cs.get_paths():
            verts = path.vertices
            if len(verts) < 3:
                continue
            coords = [
                [v[0] + easting_offset, v[1] + northing_offset]
                for v in verts
            ]
            # Close the ring
            if coords[0] != coords[-1]:
                coords.append(coords[0])
            features.append({
                "type": "Feature",
                "geometry": {
                    "type": "Polygon",
                    "coordinates": [coords],
                },
                "properties": {
                    "time_s":      t,
                    "time_min":    round(t / 60.0, 2),
                    "label":       f"t={t/60:.1f} min",
                    "level_index": level_idx,
                },
            })

    return {"type": "FeatureCollection", "features": features}


# ---------------------------------------------------------------------------
# Isochrone time generation
# ---------------------------------------------------------------------------

def make_iso_times(
    arrival: np.ndarray,
    interval: Optional[float],
    custom_times: Optional[List[float]],
    n_iso: int,
) -> List[float]:
    """Return a sorted list of isochrone times [s]."""
    if custom_times:
        return sorted(custom_times)

    valid = arrival[~np.isnan(arrival)]
    if valid.size == 0:
        return []

    t_min = float(valid.min())
    t_max = float(valid.max())

    if interval is not None:
        # Build times from the first multiple of interval ≥ t_min up to t_max
        start = max(interval, interval * int(t_min / interval + 1))
        times = []
        t = start
        while t <= t_max + 1e-6:
            times.append(t)
            t += interval
        return times

    # Default: n_iso equally-spaced times
    return list(np.linspace(t_min + (t_max - t_min) / (n_iso + 1),
                            t_max,
                            n_iso))


# ---------------------------------------------------------------------------
# Visualization
# ---------------------------------------------------------------------------

def plot_isochrones(
    arrival: np.ndarray,
    problo: List[float],
    probhi: List[float],
    nx: int,
    ny: int,
    geojson: dict,
    output_path: str
) -> None:
    """Create a visualization of isochrones with time labels."""
    if not _HAS_MPL:
        return
    
    # Create coordinate arrays for cell centers
    x = np.linspace(problo[0], probhi[0], nx + 1)
    y = np.linspace(problo[1], probhi[1], ny + 1)
    
    fig, ax = plt.subplots(figsize=(12, 10))
    
    # Plot arrival time as filled contours
    valid_mask = ~np.isnan(arrival)
    if np.any(valid_mask):
        levels = np.linspace(np.nanmin(arrival), np.nanmax(arrival), 30)
        contourf = ax.contourf(
            x[:-1] + (x[1] - x[0]) / 2,
            y[:-1] + (y[1] - y[0]) / 2,
            arrival / 60.0,  # Convert to minutes
            levels=levels / 60.0,
            cmap='YlOrRd',
            alpha=0.8
        )
        cbar = fig.colorbar(contourf, ax=ax, label='Arrival time [min]')
    
    # Plot isochrone contours with labels
    colors = plt.cm.viridis(np.linspace(0, 1, len(geojson["features"])))
    
    for idx, feature in enumerate(geojson["features"]):
        props = feature["properties"]
        coords_list = feature["geometry"]["coordinates"]
        
        # Draw each polygon ring
        for ring in coords_list:
            if len(ring) < 3:
                continue
            xs, ys = zip(*ring)
            ax.plot(xs, ys, '-', color=colors[idx], linewidth=2, 
                   label=props["label"])
            
            # Add time label at the first point
            ax.text(xs[0], ys[0], f" {props['label']}", 
                   fontsize=9, color=colors[idx],
                   bbox=dict(boxstyle='round,pad=0.3', facecolor='white', 
                            alpha=0.7, edgecolor=colors[idx]))
    
    ax.set_xlabel('X [m]')
    ax.set_ylabel('Y [m]')
    ax.set_title('Fire Arrival Time Isochrones')
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3, linestyle='--')
    
    # Add legend if not too many isochrones
    if len(geojson["features"]) <= 15:
        ax.legend(loc='best', fontsize=8)
    
    fig.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved isochrone plot → {output_path}")


# ---------------------------------------------------------------------------
# Per-plotfile processing
# ---------------------------------------------------------------------------

def process_plotfile(
    plotfile_dir: Path,
    outdir: Path,
    interval: Optional[float],
    custom_times: Optional[List[float]],
    n_iso: int,
    utm_origin: Optional[Tuple[float, float]],
    plot_output: Optional[str] = None,
):
    print(f"\nProcessing {plotfile_dir} …")
    try:
        arrival, problo, probhi, nx, ny = read_arrival_time(plotfile_dir)
    except (FileNotFoundError, ValueError) as exc:
        print(f"  ERROR: {exc}")
        return

    valid = arrival[~np.isnan(arrival)]
    if valid.size == 0:
        print("  No burned cells (arrival_time all NaN) – skipping.")
        return

    t_max_min = float(valid.max()) / 60.0
    print(f"  Grid {nx}×{ny}  burned cells: {valid.size}  "
          f"max arrival time: {t_max_min:.1f} min")

    iso_times = make_iso_times(arrival, interval, custom_times, n_iso)
    if not iso_times:
        print("  WARNING: no isochrone times in range – skipping.")
        return
    print(f"  Extracting {len(iso_times)} isochrones at: "
          + ", ".join(f"{t/60:.1f} min" for t in iso_times))

    try:
        geojson = extract_isochrones(
            arrival, problo, probhi, nx, ny, iso_times, utm_origin=utm_origin
        )
    except ImportError as exc:
        print(f"  ERROR: {exc}")
        return

    n_feat = len(geojson["features"])
    if n_feat == 0:
        print("  WARNING: no isochrone contours found (fire may not span any interval).")
        return

    outdir.mkdir(parents=True, exist_ok=True)
    stem = plotfile_dir.name
    out_path = outdir / f"{stem}_isochrones.geojson"
    with open(out_path, "w") as fh:
        json.dump(geojson, fh, indent=2)
    print(f"  Wrote {n_feat} isochrone polygon(s) → {out_path}")
    
    # Optional plot
    if plot_output and _HAS_MPL:
        plot_isochrones(arrival, problo, probhi, nx, ny, geojson, plot_output)



# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Extract fire arrival-time isochrones from AMReX plotfiles.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "plotfile",
        nargs="?",
        help="Plotfile directory (e.g. plt0100).  Omit with --all.",
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="Process all plt#### directories in the current directory.",
    )
    parser.add_argument(
        "--interval",
        type=float,
        metavar="SECONDS",
        help="Isochrone interval [s].  Isochrones are written at multiples of "
             "this interval from the first ignition to the final arrival time.",
    )
    parser.add_argument(
        "--times",
        nargs="+",
        type=float,
        metavar="T",
        help="Explicit isochrone times [s] (overrides --interval).",
    )
    parser.add_argument(
        "--n-iso",
        type=int,
        default=10,
        metavar="N",
        help="Number of equally-spaced isochrones when neither --interval nor "
             "--times is given (default: 10).",
    )
    parser.add_argument(
        "--outdir",
        default="iso_out",
        metavar="DIR",
        help="Output directory for GeoJSON files (default: iso_out).",
    )
    parser.add_argument(
        "--utm-origin",
        nargs=2,
        type=float,
        metavar=("EASTING", "NORTHING"),
        help="UTM origin of the simulation domain [m].  Added to simulation "
             "coordinates to produce absolute UTM coordinates.",
    )
    parser.add_argument(
        "--plot",
        default=None,
        metavar="FILE",
        help="Save a matplotlib plot of isochrones with labels (e.g., isochrones.png).",
    )

    args = parser.parse_args(argv)
    outdir = Path(args.outdir)
    utm_origin = tuple(args.utm_origin) if args.utm_origin else None

    if args.all:
        dirs = sorted(Path(".").glob("plt[0-9][0-9][0-9][0-9]"))
        if not dirs:
            sys.exit("No plt#### directories found in the current directory.")
        for d in dirs:
            # For --all mode, create per-plotfile plot names
            plot_name = str(outdir / f"{d.name}_isochrones.png") if args.plot else None
            process_plotfile(d, outdir, args.interval, args.times, args.n_iso, 
                           utm_origin, plot_name)
    elif args.plotfile:
        process_plotfile(
            Path(args.plotfile), outdir,
            args.interval, args.times, args.n_iso, utm_origin, args.plot,
        )
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
