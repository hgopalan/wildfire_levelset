#!/usr/bin/env python3
"""ember_cascade_analysis.py — Post-processing tool for the flux-based ember cascade model.

Reads the ``ember_cascade_flux`` and ``ember_cascade_ignition`` fields from
AMReX plotfiles produced by the ``ember_cascade.enable = 1`` model and
generates summary statistics, time-series tables, and optional visualisations.

Physical background
-------------------
The ember cascade model computes a continuous landing-flux density
[embers/m²/s] at each grid cell by summing Gaussian kernels from all
active fire-front source cells:

.. math::

   F_{\\text{land}}(x, y) = \\sum_{s} F_{\\text{src},s}
       \\cdot \\frac{A_{\\text{cell}}}{2\\pi\\sigma_s^2}
       \\cdot \\exp\\!\\left(-\\frac{\\Delta x_s^2 + \\Delta y_s^2}{2\\sigma_s^2}\\right)

where :math:`\\Delta x_s = x - (x_s + u_s t_f)`, :math:`t_f = H_z/v_t`, and
:math:`\\sigma_s = \\sigma_{\\text{base}} + k_\\sigma H_z`.

This tool extracts those fields for analysis.

Usage
-----
Analyse all plotfiles in a directory::

    python3 tools/ember_cascade_analysis.py --plt-dir ./plt_dir

Analyse specific plotfile steps and write a CSV summary::

    python3 tools/ember_cascade_analysis.py --plt-dir ./plt_dir \\
        --steps 0010 0020 0050 --csv ember_summary.csv

Visualise the landing flux map at step 50::

    python3 tools/ember_cascade_analysis.py --plt-dir ./plt_dir \\
        --plot 0050 --output flux_map.png

Dependencies
------------
``numpy`` (required), ``matplotlib`` (optional, for ``--plot``)

The tool uses only AMReX plotfile Header parsing and raw binary FAB reading —
it does NOT require the ``amrex`` Python bindings or any C extension.

Examples
--------
Typical workflow::

    # 1. Run simulation
    ./build/levelset regtest/spotting/ember_cascade_flux/inputs.i

    # 2. Summarise ember cascade activity across all steps
    python3 tools/ember_cascade_analysis.py --plt-dir . --csv ember_cascade.csv

    # 3. Plot the landing flux map at the final step
    python3 tools/ember_cascade_analysis.py --plt-dir . --plot 0100 \\
        --output final_flux.png
"""

import argparse
import os
import re
import struct
import sys
from pathlib import Path


# ---------------------------------------------------------------------------
# AMReX plotfile reader (Header + VisMF FAB, no external dependencies)
# ---------------------------------------------------------------------------

def _parse_plt_header(plt_dir):
    """Parse an AMReX plotfile Header file.

    Returns
    -------
    dict with keys:
        var_names : list[str]
        spacedim  : int
        time      : float
        prob_lo   : list[float]
        prob_hi   : list[float]
        nx, ny    : int  (2-D cell counts)
        level_dir : str  (path to the Level_0 data)
    """
    header_path = os.path.join(plt_dir, "Header")
    if not os.path.isfile(header_path):
        raise FileNotFoundError(f"No Header file in {plt_dir}")

    with open(header_path) as f:
        lines = [l.rstrip("\n") for l in f.readlines()]

    idx = 0
    _version = lines[idx]; idx += 1                  # HyperCLaw-V1.1
    nvar = int(lines[idx]); idx += 1
    var_names = [lines[idx + i] for i in range(nvar)]; idx += nvar
    spacedim = int(lines[idx]); idx += 1
    time = float(lines[idx]); idx += 1
    _finest = int(lines[idx]); idx += 1              # finest_level
    prob_lo = list(map(float, lines[idx].split())); idx += 1
    prob_hi = list(map(float, lines[idx].split())); idx += 1

    # Parse bbox: "((ix0,iy0,...) (ix1,iy1,...) (0,0,...))"
    bbox_line = lines[idx]; idx += 1
    coords = re.findall(r'\(([^)]+)\)', bbox_line)
    lo_idx = list(map(int, coords[0].split(',')))
    hi_idx = list(map(int, coords[1].split(',')))
    nx = hi_idx[0] - lo_idx[0] + 1
    ny = hi_idx[1] - lo_idx[1] + 1

    _step = lines[idx]; idx += 1
    _cell_sizes = lines[idx]; idx += 1
    _coord = lines[idx]; idx += 1
    _reserved = lines[idx]; idx += 1
    level_data = lines[idx].strip(); idx += 1        # "Level_0/Cell"
    level_dir = os.path.join(plt_dir, os.path.dirname(level_data))

    return {
        "var_names": var_names,
        "spacedim": spacedim,
        "time": time,
        "prob_lo": prob_lo,
        "prob_hi": prob_hi,
        "nx": nx,
        "ny": ny,
        "level_dir": level_dir,
    }


def _read_fab_field(level_dir, var_names, field_name, nx, ny):
    """Read a single field from an AMReX Level_0 FAB data file.

    Parameters
    ----------
    level_dir : str
    var_names : list[str]
    field_name : str   name of the variable to extract
    nx, ny    : int    2-D grid dimensions

    Returns
    -------
    numpy.ndarray of shape (ny, nx)  or None if field not found
    """
    import numpy as np

    if field_name not in var_names:
        return None

    comp = var_names.index(field_name)
    nvar = len(var_names)

    # Read Cell_H to locate the data file and byte offset
    cell_h = os.path.join(level_dir, "Cell_H")
    with open(cell_h) as f:
        lines = [l.rstrip() for l in f.readlines()]

    _vmf_version = lines[0]   # VisMF-V1
    _nfab = int(lines[1])
    _ncomp = int(lines[2])
    _nghost = int(lines[3])
    data_entry = lines[4]     # "((lo) (hi) (0)) nvar  filename:offset"

    parts = data_entry.split()
    file_ref = parts[-1]
    data_filename, offset_str = file_ref.rsplit(":", 1)
    offset = int(offset_str)
    data_path = os.path.join(level_dir, data_filename)

    # Read the ASCII FAB header at `offset` to find the binary data start
    with open(data_path, "rb") as f:
        f.seek(offset)
        fab_header_line = f.readline()
        binary_start = f.tell()
        n_cells = nx * ny                          # 2-D (nz == 1 in 2-D builds)
        # Each component is n_cells doubles
        f.seek(binary_start + comp * n_cells * 8)
        raw = f.read(n_cells * 8)

    values = np.array(struct.unpack(f"{n_cells}d", raw)).reshape((ny, nx))
    return values


# ---------------------------------------------------------------------------
# Plotfile discovery
# ---------------------------------------------------------------------------

def discover_plotfiles(plt_dir):
    """Return sorted list of (step, path) for all plt* directories."""
    entries = []
    for name in sorted(os.listdir(plt_dir)):
        m = re.fullmatch(r"plt(\d+)", name)
        if m:
            full = os.path.join(plt_dir, name)
            if os.path.isdir(full) and os.path.isfile(os.path.join(full, "Header")):
                entries.append((int(m.group(1)), full))
    return entries


# ---------------------------------------------------------------------------
# Analysis helpers
# ---------------------------------------------------------------------------

def analyse_plotfile(plt_path):
    """Extract ember cascade summary statistics from a single plotfile.

    Returns
    -------
    dict or None if neither ``ember_cascade_flux`` nor
    ``ember_cascade_ignition`` is present.
    """
    try:
        import numpy as np
    except ImportError:
        sys.exit("numpy is required: pip install numpy")

    hdr = _parse_plt_header(plt_path)
    var_names = hdr["var_names"]
    nx, ny = hdr["nx"], hdr["ny"]
    level_dir = hdr["level_dir"]

    if ("ember_cascade_flux" not in var_names and
            "ember_cascade_ignition" not in var_names):
        return None

    result = {"step": int(os.path.basename(plt_path).lstrip("plt")),
              "time_s": hdr["time"]}

    flux = _read_fab_field(level_dir, var_names, "ember_cascade_flux", nx, ny)
    if flux is not None:
        result["max_flux"] = float(np.max(flux))
        result["mean_flux"] = float(np.mean(flux))
        result["n_active_cells"] = int(np.count_nonzero(flux > 0))
        result["total_flux"] = float(np.sum(flux))

    ign = _read_fab_field(level_dir, var_names, "ember_cascade_ignition", nx, ny)
    if ign is not None:
        result["n_ignitions"] = int(np.sum(ign > 0.5))

    return result


# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------

def write_csv(rows, csv_path):
    """Write analysis rows to a CSV file."""
    if not rows:
        print("No ember cascade data found.")
        return
    fieldnames = list(rows[0].keys())
    with open(csv_path, "w") as f:
        f.write(",".join(fieldnames) + "\n")
        for row in rows:
            f.write(",".join(str(row.get(k, "")) for k in fieldnames) + "\n")
    print(f"Wrote {len(rows)} rows to {csv_path}")


def plot_flux_map(plt_path, output_path=None):
    """Visualise the ember landing-flux density map from a single plotfile."""
    try:
        import numpy as np
        import matplotlib.pyplot as plt
        import matplotlib.colors as mcolors
    except ImportError:
        sys.exit("matplotlib and numpy are required for --plot: pip install matplotlib numpy")

    hdr = _parse_plt_header(plt_path)
    var_names = hdr["var_names"]
    nx, ny = hdr["nx"], hdr["ny"]
    level_dir = hdr["level_dir"]
    prob_lo = hdr["prob_lo"]
    prob_hi = hdr["prob_hi"]

    flux = _read_fab_field(level_dir, var_names, "ember_cascade_flux", nx, ny)
    ign = _read_fab_field(level_dir, var_names, "ember_cascade_ignition", nx, ny)
    phi = _read_fab_field(level_dir, var_names, "phi", nx, ny)

    if flux is None:
        sys.exit(f"ember_cascade_flux not found in {plt_path}")

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Extent for imshow
    extent = [prob_lo[0], prob_hi[0], prob_lo[1], prob_hi[1]]

    # Left panel: landing flux density
    ax = axes[0]
    flux_plot = np.where(flux > 0, flux, np.nan)
    vmax = np.nanmax(flux_plot) if np.any(np.isfinite(flux_plot)) else 1.0
    im = ax.imshow(flux_plot, origin="lower", extent=extent,
                   norm=mcolors.LogNorm(vmin=max(vmax * 1e-4, 1e-10), vmax=vmax),
                   cmap="hot_r", interpolation="nearest")
    if phi is not None:
        ax.contour(np.linspace(prob_lo[0], prob_hi[0], nx),
                   np.linspace(prob_lo[1], prob_hi[1], ny),
                   phi, levels=[0.0], colors="red", linewidths=1.0)
    plt.colorbar(im, ax=ax, label="Landing flux [embers/m²/s]")
    ax.set_xlabel("Easting [m]")
    ax.set_ylabel("Northing [m]")
    step_str = os.path.basename(plt_path)
    ax.set_title(f"Ember cascade flux density — {step_str}")

    # Right panel: ignition flag
    ax = axes[1]
    if ign is not None:
        ign_plot = np.where(ign > 0.5, 1.0, np.nan)
        ax.imshow(ign_plot, origin="lower", extent=extent,
                  cmap="Reds", vmin=0, vmax=1, interpolation="nearest")
    if phi is not None:
        ax.contour(np.linspace(prob_lo[0], prob_hi[0], nx),
                   np.linspace(prob_lo[1], prob_hi[1], ny),
                   phi, levels=[0.0], colors="red", linewidths=1.0)
    ax.set_xlabel("Easting [m]")
    ax.set_ylabel("Northing [m]")
    n_ign = int(np.sum(ign > 0.5)) if ign is not None else 0
    ax.set_title(f"Ember ignitions ({n_ign} cells) — {step_str}")

    plt.tight_layout()
    if output_path:
        plt.savefig(output_path, dpi=150)
        print(f"Saved plot to {output_path}")
    else:
        plt.show()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--plt-dir", default=".",
                        help="Directory containing plt* plotfile folders (default: current dir)")
    parser.add_argument("--steps", nargs="*",
                        help="Zero-padded step numbers to analyse (default: all)")
    parser.add_argument("--csv", metavar="FILE",
                        help="Write summary statistics to a CSV file")
    parser.add_argument("--plot", metavar="STEP",
                        help="Plot flux map for the given zero-padded step (e.g. 0050)")
    parser.add_argument("--output", metavar="PNG",
                        help="Save plot to this file instead of showing interactively")
    args = parser.parse_args()

    plt_dir = args.plt_dir
    plotfiles = discover_plotfiles(plt_dir)

    if not plotfiles:
        print(f"No plt* plotfiles found in {plt_dir}")
        sys.exit(0)

    # Filter to requested steps
    if args.steps:
        requested = {int(s) for s in args.steps}
        plotfiles = [(s, p) for s, p in plotfiles if s in requested]

    if args.plot:
        target_step = int(args.plot)
        matches = [p for s, p in plotfiles if s == target_step]
        if not matches:
            sys.exit(f"Step {target_step} not found in {plt_dir}")
        plot_flux_map(matches[0], args.output)
        return

    # Analyse all (or selected) plotfiles
    rows = []
    for step, plt_path in plotfiles:
        result = analyse_plotfile(plt_path)
        if result is not None:
            rows.append(result)

    if not rows:
        print("No plotfiles with ember_cascade fields found.")
        sys.exit(0)

    # Print table
    header = (f"{'Step':>6}  {'Time(s)':>10}  {'MaxFlux':>12}  "
              f"{'MeanFlux':>12}  {'ActiveCells':>12}  {'Ignitions':>10}")
    print(header)
    print("-" * len(header))
    for row in rows:
        print(f"{row.get('step', ''):>6}  "
              f"{row.get('time_s', 0):>10.1f}  "
              f"{row.get('max_flux', 0):>12.4e}  "
              f"{row.get('mean_flux', 0):>12.4e}  "
              f"{row.get('n_active_cells', 0):>12d}  "
              f"{row.get('n_ignitions', 0):>10d}")

    if args.csv:
        write_csv(rows, args.csv)


if __name__ == "__main__":
    main()
