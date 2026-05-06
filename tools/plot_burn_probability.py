#!/usr/bin/env python3
"""
plot_burn_probability.py – Visualise ensemble burn probability maps.

Reads one or more ``burn_probability.csv`` files written by
``ensemble_burn_probability.py`` and produces publication-quality raster
plots using matplotlib.

The CSV format is::

  # Ensemble burn probability map
  # n_runs = 50  resolution = 30 m
  X,Y,P_burn

Each row represents a grid cell centre (X, Y) in simulation physical
coordinates [m] (UTM easting/northing when the domain is georeferenced).
``P_burn`` is the fraction of ensemble members [0, 1] in which that cell
burned.

Usage
-----
  # Basic heatmap, auto-colourbar, output PNG
  python3 tools/plot_burn_probability.py burn_probability.csv

  # Custom colour map and thresholds
  python3 tools/plot_burn_probability.py burn_probability.csv \\
      --cmap YlOrRd --vmin 0.1 --vmax 1.0 --out bp_map.png

  # Compare two runs side by side
  python3 tools/plot_burn_probability.py run_A/burn_probability.csv \\
      run_B/burn_probability.csv \\
      --labels "Run A (dry)" "Run B (wet)" --out comparison.png

  # Overlay probability contour lines at 10 %, 25 %, 50 %, 75 %
  python3 tools/plot_burn_probability.py burn_probability.csv \\
      --contours 0.10 0.25 0.50 0.75

  # Export as GeoTIFF (requires rasterio)
  python3 tools/plot_burn_probability.py burn_probability.csv \\
      --geotiff burn_prob.tif --epsg 32613

  # Difference map between two runs
  python3 tools/plot_burn_probability.py \\
      --diff run_A/burn_probability.csv run_B/burn_probability.csv \\
      --out diff_map.png

Requirements
------------
  pip install matplotlib numpy
  pip install rasterio   # optional, for --geotiff export

References
----------
  Finney, M.A. et al. (2011). A method for ensemble wildland fire simulation.
    Environmental Modelling & Software, 26(10), 1352-1359.
  FSPro (Fire Spread Probability) produces analogous burn probability maps.
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple


# ---------------------------------------------------------------------------
# CSV reader
# ---------------------------------------------------------------------------

def read_burn_probability_csv(
    csv_path: str,
) -> Tuple[Dict[Tuple[float, float], float], float, int]:
    """Read burn_probability.csv.

    Returns
    -------
    points : dict  (x, y) → P_burn
    resolution : float  grid cell size [m] (parsed from comment, else 0)
    n_runs : int  number of ensemble members (parsed from comment, else 0)
    """
    points: Dict[Tuple[float, float], float] = {}
    resolution = 0.0
    n_runs = 0

    with open(csv_path, newline="") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                # Try to parse metadata from comment lines
                # e.g. "# n_runs = 50  resolution = 30 m"
                if "n_runs" in line:
                    import re
                    m = re.search(r"n_runs\s*=\s*(\d+)", line)
                    if m:
                        n_runs = int(m.group(1))
                    m = re.search(r"resolution\s*=\s*([\d.]+)", line)
                    if m:
                        resolution = float(m.group(1))
                continue
            if line.lower().startswith("x") or line.lower().startswith('"x"'):
                continue  # header row
            parts = line.split(",")
            if len(parts) < 3:
                continue
            try:
                x, y, p = float(parts[0]), float(parts[1]), float(parts[2])
                points[(x, y)] = p
            except ValueError:
                continue

    return points, resolution, n_runs


# ---------------------------------------------------------------------------
# Grid builder
# ---------------------------------------------------------------------------

def build_grid(
    points: Dict[Tuple[float, float], float],
    resolution: float,
) -> Tuple["np.ndarray", "np.ndarray", "np.ndarray", float]:
    """Convert scattered (x, y, p) points to a regular 2-D numpy grid.

    Returns (X_1d, Y_1d, grid_2d, resolution).
    grid_2d has shape (len(Y_1d), len(X_1d)).
    Cells not present in *points* are set to NaN.
    """
    import numpy as np

    if not points:
        raise ValueError("No data points in CSV.")

    xs = sorted({xy[0] for xy in points})
    ys = sorted({xy[1] for xy in points})

    # Auto-detect resolution if not provided
    if resolution <= 0.0:
        if len(xs) > 1:
            diffs = [xs[i+1] - xs[i] for i in range(len(xs)-1) if xs[i+1] > xs[i]]
            resolution = min(diffs) if diffs else 1.0
        else:
            resolution = 1.0

    x_arr = np.array(xs)
    y_arr = np.array(ys)
    grid = np.full((len(y_arr), len(x_arr)), np.nan)

    # Build index maps
    x_idx = {x: i for i, x in enumerate(xs)}
    y_idx = {y: j for j, y in enumerate(ys)}

    for (x, y), p in points.items():
        xi = x_idx.get(x)
        yi = y_idx.get(y)
        if xi is not None and yi is not None:
            grid[yi, xi] = p

    return x_arr, y_arr, grid, resolution


# ---------------------------------------------------------------------------
# Single-file plot
# ---------------------------------------------------------------------------

def plot_single(
    ax,
    x_arr, y_arr, grid,
    title: str = "",
    cmap: str = "YlOrRd",
    vmin: float = 0.0,
    vmax: float = 1.0,
    contours: Optional[List[float]] = None,
    show_colorbar: bool = True,
    resolution: float = 0.0,
    n_runs: int = 0,
):
    """Plot one burn probability grid on *ax*."""
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.colors import BoundaryNorm, LinearSegmentedColormap

    # Mask NaN
    masked = np.ma.masked_invalid(grid)

    # Use pcolormesh for correct cell placement
    dx = (x_arr[1] - x_arr[0]) / 2 if len(x_arr) > 1 else resolution / 2
    dy = (y_arr[1] - y_arr[0]) / 2 if len(y_arr) > 1 else resolution / 2
    # Cell edges
    xe = np.concatenate([[x_arr[0] - dx], x_arr[:-1] + np.diff(x_arr)/2,
                          [x_arr[-1] + dx]])
    ye = np.concatenate([[y_arr[0] - dy], y_arr[:-1] + np.diff(y_arr)/2,
                          [y_arr[-1] + dy]])

    im = ax.pcolormesh(xe, ye, masked, cmap=cmap, vmin=vmin, vmax=vmax,
                       shading="flat", rasterized=True)

    if contours:
        X, Y = np.meshgrid(x_arr, y_arr)
        cs = ax.contour(X, Y, masked.filled(np.nan), levels=sorted(contours),
                        colors="k", linewidths=0.8, linestyles="--")
        ax.clabel(cs, fmt=lambda v: f"{v*100:.0f}%", fontsize=7, inline=True)

    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("X [m]")
    ax.set_ylabel("Y [m]")

    subtitle_parts = []
    if n_runs > 0:
        subtitle_parts.append(f"n_runs = {n_runs}")
    if resolution > 0:
        subtitle_parts.append(f"res = {resolution:.0f} m")
    subtitle = "  |  ".join(subtitle_parts)

    ax.set_title(f"{title}\n{subtitle}" if subtitle else title, fontsize=10)

    if show_colorbar:
        cbar = ax.get_figure().colorbar(im, ax=ax, fraction=0.04, pad=0.02)
        cbar.set_label("Burn probability [-]")

    return im


# ---------------------------------------------------------------------------
# Difference map
# ---------------------------------------------------------------------------

def plot_diff(
    ax,
    x_arr_a, y_arr_a, grid_a,
    x_arr_b, y_arr_b, grid_b,
    label_a: str = "A",
    label_b: str = "B",
):
    """Plot P_burn(A) – P_burn(B) as a diverging map."""
    import numpy as np

    # Align on common grid
    xs = sorted(set(x_arr_a.tolist()) | set(x_arr_b.tolist()))
    ys = sorted(set(y_arr_a.tolist()) | set(y_arr_b.tolist()))
    x_union = np.array(xs)
    y_union = np.array(ys)

    def _expand(x_src, y_src, g_src):
        out = np.full((len(y_union), len(x_union)), np.nan)
        xi_map = {x: i for i, x in enumerate(xs)}
        yi_map = {y: j for j, y in enumerate(ys)}
        for r, y in enumerate(y_src):
            for c, x in enumerate(x_src):
                if yi_map.get(y) is not None and xi_map.get(x) is not None:
                    out[yi_map[y], xi_map[x]] = g_src[r, c]
        return out

    ga = _expand(x_arr_a, y_arr_a, grid_a)
    gb = _expand(x_arr_b, y_arr_b, grid_b)
    diff = ga - gb

    masked = np.ma.masked_invalid(diff)
    vmax_abs = max(abs(float(np.nanmin(diff))), abs(float(np.nanmax(diff))), 0.01)

    dx = (x_union[1] - x_union[0]) / 2 if len(x_union) > 1 else 1
    dy = (y_union[1] - y_union[0]) / 2 if len(y_union) > 1 else 1
    xe = np.concatenate([[x_union[0]-dx], x_union[:-1]+np.diff(x_union)/2, [x_union[-1]+dx]])
    ye = np.concatenate([[y_union[0]-dy], y_union[:-1]+np.diff(y_union)/2, [y_union[-1]+dy]])

    im = ax.pcolormesh(xe, ye, masked, cmap="RdBu_r",
                       vmin=-vmax_abs, vmax=vmax_abs,
                       shading="flat", rasterized=True)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("X [m]")
    ax.set_ylabel("Y [m]")
    ax.set_title(f"ΔP = {label_a} − {label_b}", fontsize=10)
    cbar = ax.get_figure().colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    cbar.set_label("ΔP burn [-]")
    return im


# ---------------------------------------------------------------------------
# GeoTIFF export
# ---------------------------------------------------------------------------

def export_geotiff(
    x_arr, y_arr, grid,
    out_path: str,
    resolution: float,
    epsg: Optional[int] = None,
):
    """Write the burn probability grid as a GeoTIFF (requires rasterio)."""
    try:
        import rasterio
        from rasterio.transform import from_origin
        from rasterio.crs import CRS
    except ImportError:
        print("ERROR: rasterio not installed.  Run:  pip install rasterio",
              file=sys.stderr)
        return

    import numpy as np

    dx = (x_arr[1] - x_arr[0]) if len(x_arr) > 1 else resolution
    dy = (y_arr[1] - y_arr[0]) if len(y_arr) > 1 else resolution

    west  = float(x_arr[0]) - dx / 2
    north = float(y_arr[-1]) + dy / 2

    transform = from_origin(west, north, dx, dy)

    # flip so row 0 = northernmost
    arr = np.flipud(grid).astype(np.float32)

    profile = {
        "driver":    "GTiff",
        "dtype":     "float32",
        "width":     len(x_arr),
        "height":    len(y_arr),
        "count":     1,
        "transform": transform,
        "compress":  "deflate",
        "nodata":    float("nan"),
    }
    if epsg is not None:
        profile["crs"] = CRS.from_epsg(epsg)

    with rasterio.open(out_path, "w", **profile) as dst:
        dst.write(arr, 1)
        dst.update_tags(variable="burn_probability")
    print(f"Wrote GeoTIFF → {out_path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Plot ensemble burn probability maps from burn_probability.csv.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "csvfiles",
        nargs="*",
        metavar="CSV",
        help="One or more burn_probability.csv files to plot.",
    )
    parser.add_argument(
        "--diff",
        nargs=2,
        metavar=("CSV_A", "CSV_B"),
        help="Produce a difference map: P_burn(A) – P_burn(B).",
    )
    parser.add_argument(
        "--labels",
        nargs="+",
        metavar="LABEL",
        help="Panel titles matching the order of CSV inputs (default: file names).",
    )
    parser.add_argument(
        "--cmap",
        default="YlOrRd",
        help="Matplotlib colormap for the probability raster (default: YlOrRd).",
    )
    parser.add_argument(
        "--vmin",
        type=float,
        default=0.0,
        help="Colourbar lower bound (default: 0.0).",
    )
    parser.add_argument(
        "--vmax",
        type=float,
        default=1.0,
        help="Colourbar upper bound (default: 1.0).",
    )
    parser.add_argument(
        "--contours",
        nargs="+",
        type=float,
        metavar="P",
        help="Overlay probability contour lines at these values, e.g. 0.1 0.5 0.9.",
    )
    parser.add_argument(
        "--out",
        default="burn_probability_map.png",
        help="Output figure path (default: burn_probability_map.png).  "
             "Format is inferred from the extension: .png, .pdf, .svg.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=150,
        help="Figure DPI (default: 150).",
    )
    parser.add_argument(
        "--geotiff",
        metavar="FILE",
        help="Also export the first CSV as a GeoTIFF (requires rasterio).",
    )
    parser.add_argument(
        "--epsg",
        type=int,
        metavar="CODE",
        help="EPSG code for GeoTIFF CRS (e.g. 32613 for UTM zone 13N).",
    )

    args = parser.parse_args(argv)

    # Validate inputs
    if not args.csvfiles and not args.diff:
        parser.print_help()
        sys.exit(1)

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        sys.exit("ERROR: matplotlib and numpy are required.  "
                 "Run:  pip install matplotlib numpy")

    # ---- Difference map mode ----
    if args.diff:
        csv_a, csv_b = args.diff
        for f in [csv_a, csv_b]:
            if not os.path.isfile(f):
                sys.exit(f"ERROR: file not found: {f}")

        pts_a, res_a, n_a = read_burn_probability_csv(csv_a)
        pts_b, res_b, n_b = read_burn_probability_csv(csv_b)
        xa, ya, ga, resa = build_grid(pts_a, res_a)
        xb, yb, gb, resb = build_grid(pts_b, res_b)

        labels = args.labels or [Path(csv_a).stem, Path(csv_b).stem]
        la = labels[0] if len(labels) > 0 else Path(csv_a).stem
        lb = labels[1] if len(labels) > 1 else Path(csv_b).stem

        fig, axes = plt.subplots(1, 3, figsize=(16, 5))
        plot_single(axes[0], xa, ya, ga, title=la, cmap=args.cmap,
                    vmin=args.vmin, vmax=args.vmax, contours=args.contours,
                    n_runs=n_a, resolution=resa)
        plot_single(axes[1], xb, yb, gb, title=lb, cmap=args.cmap,
                    vmin=args.vmin, vmax=args.vmax, contours=args.contours,
                    n_runs=n_b, resolution=resb)
        plot_diff(axes[2], xa, ya, ga, xb, yb, gb, la, lb)

        fig.suptitle("Ensemble burn probability comparison", fontsize=12)
        fig.tight_layout()
        fig.savefig(args.out, dpi=args.dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved → {args.out}")

        if args.geotiff:
            export_geotiff(xa, ya, ga, args.geotiff, resa, args.epsg)
        return

    # ---- Single or multi-panel mode ----
    csv_files = args.csvfiles
    for f in csv_files:
        if not os.path.isfile(f):
            sys.exit(f"ERROR: file not found: {f}")

    n_panels = len(csv_files)
    labels = args.labels or [Path(f).stem for f in csv_files]
    # Pad labels if fewer provided than files
    while len(labels) < n_panels:
        labels.append(Path(csv_files[len(labels)]).stem)

    fig_w = max(6, 6 * n_panels)
    fig, axes = plt.subplots(1, n_panels, figsize=(fig_w, 5),
                             squeeze=False)
    axes = axes[0]  # flatten to 1-D list

    first_pts = first_res = first_xa = first_ya = first_ga = None

    for ax, csv_path, label in zip(axes, csv_files, labels):
        pts, res, n_runs = read_burn_probability_csv(csv_path)
        if not pts:
            ax.text(0.5, 0.5, f"No data in\n{csv_path}",
                    ha="center", va="center", transform=ax.transAxes)
            continue

        xa, ya, ga, resa = build_grid(pts, res)
        if first_xa is None:
            first_xa, first_ya, first_ga, first_res = xa, ya, ga, resa

        plot_single(ax, xa, ya, ga, title=label, cmap=args.cmap,
                    vmin=args.vmin, vmax=args.vmax, contours=args.contours,
                    show_colorbar=(n_panels == 1 or ax is axes[-1]),
                    n_runs=n_runs, resolution=resa)

        # Print summary stats
        valid = ga[~np.isnan(ga)]
        if valid.size > 0:
            high_burn = (valid >= 0.5).sum()
            print(f"{label}:  {valid.size} cells  "
                  f"mean P={float(valid.mean()):.3f}  "
                  f"max P={float(valid.max()):.3f}  "
                  f"cells ≥50%: {high_burn} ({100*high_burn/valid.size:.1f}%)")

    title = "Ensemble burn probability"
    if n_panels == 1:
        title = f"{title} — {labels[0]}"
    fig.suptitle(title, fontsize=12)
    fig.tight_layout()
    fig.savefig(args.out, dpi=args.dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved → {args.out}")

    if args.geotiff and first_xa is not None:
        export_geotiff(first_xa, first_ya, first_ga, args.geotiff,
                       first_res, args.epsg)


if __name__ == "__main__":
    main()
