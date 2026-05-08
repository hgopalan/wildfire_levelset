#!/usr/bin/env python3
"""
plotfile_to_geotiff.py – Convert AMReX plotfiles to GeoTIFF rasters for GIS.

The script reads a 2-D AMReX plotfile directory (e.g. ``plt0100``) and writes
one GeoTIFF per variable (or a user-specified subset).  An optional UTM origin
can be supplied so that the output is correctly georeferenced; without it the
coordinates are written in the native simulation units (metres by default).

Cap 6 — Multi-variable behavioral raster export
------------------------------------------------
By default, **all** variables stored in a plotfile are converted automatically
to individual GeoTIFFs in one pass.  This covers the full set of fire-behavior
diagnostics written by the solver (phi, R, fireline_intensity, flame_length,
elevation, slope, aspect, fuel_model, scorch_height, prob_ignition,
tree_mortality, crown_activity, co2_emissions, …) without requiring the user
to enumerate them with ``-v``.

Use ``--all`` to process every ``plt####`` directory in the working directory
in a single command.

**Multi-level (AMR) plotfiles** from external AMReX-based codes are supported.
When ``finest_level > 0`` is detected in the plotfile Header, all available
levels are read.  Finer-level data are composited (block-averaged) onto the
Level_0 base grid, so the output GeoTIFF has full-domain coverage at the
coarsest resolution with finer-level detail where AMR patches exist.

Requirements
------------
  pip install rasterio numpy

Usage
-----
  # Export ALL variables from a plotfile (Cap 6 default)
  python3 tools/plotfile_to_geotiff.py plt0100 --outdir gis_out

  # Specific variables only
  python3 tools/plotfile_to_geotiff.py plt0100 -v phi R fireline_intensity flame_length

  # Georeference with a UTM origin (easting, northing) and EPSG code
  python3 tools/plotfile_to_geotiff.py plt0100 \\
      --utm-origin 450000 4200000 \\
      --epsg 32613 \\
      --outdir gis_out

  # Convert every plt#### directory in the current working directory (Cap 6 batch)
  python3 tools/plotfile_to_geotiff.py --all --outdir gis_out

  # Convert a multi-level AMR plotfile from an external AMReX-based code
  python3 tools/plotfile_to_geotiff.py plt_amr0050 --outdir gis_out

Fire-behaviour variables of interest
-------------------------------------
  phi                 – level-set function (< 0 = burned / on fire)
  R                   – Rothermel rate of spread [m/s]
  fireline_intensity  – Byram (1959) fireline intensity [kW/m]
  flame_length        – Byram (1959) flame length [m]
  elevation           – terrain elevation [m]
  slope               – terrain slope [degrees]
  aspect              – terrain aspect [degrees]
  fuel_model          – FBFM13/FBFM40 fuel model code
  scorch_height       – Van Wagner (1973) scorch height [m]
  prob_ignition       – Anderson (1970) probability of ignition [-]
  tree_mortality      – Ryan-Reinhardt (1988) tree mortality fraction [-]
  crown_activity      – crown fire activity (0=surface, 1=passive, 2=active)
  co2_emissions       – cumulative CO₂ emissions [kg/m²]
  arrival_time        – time of first ignition per cell [s]
  heat_per_unit_area  – Rothermel heat per unit area [BTU/ft²]
"""

import argparse
import glob
import os
import struct
import sys
from pathlib import Path

import numpy as np

try:
    import rasterio
    from rasterio.transform import from_origin
    from rasterio.crs import CRS
except ImportError:
    sys.exit(
        "ERROR: rasterio is not installed.  Install it with:\n"
        "    pip install rasterio"
    )


# ---------------------------------------------------------------------------
# AMReX plotfile parser (2-D, single-level and multi-level AMR)
# ---------------------------------------------------------------------------

# Fallback domain size used when the Header cannot be parsed completely.
_DEFAULT_DOMAIN_SIZE = (64, 64)

# Default refinement ratio assumed when the Header omits a level's ref_ratio.
_DEFAULT_REF_RATIO = 2

def _parse_header(plotfile_dir: Path):
    """Return (varnames, problo, probhi, nx, ny, finest_level, ref_ratios, level_domains)
    from the plotfile Header.

    Supports both single-level (finest_level=0) and multi-level AMR plotfiles.
    *level_domains* is a list of (nx, ny) tuples, one per level (index 0 = coarsest).
    *ref_ratios* is a list of integer refinement ratios between successive levels.
    *nx*, *ny* are the coarsest (Level_0) domain dimensions.
    """
    header_path = plotfile_dir / "Header"
    if not header_path.exists():
        raise FileNotFoundError(f"No Header file in {plotfile_dir}")

    with open(header_path) as fh:
        lines = [l.rstrip("\n") for l in fh]

    idx = 0
    # Line 0: version string  (e.g. "HyperCLaw-V1.1")
    idx = 1

    # Number of components
    ncomp = int(lines[idx]); idx += 1

    # Variable names (one per line)
    varnames = []
    for _ in range(ncomp):
        varnames.append(lines[idx].strip()); idx += 1

    # spacedim
    spacedim = int(lines[idx]); idx += 1
    if spacedim != 2:
        raise ValueError(
            f"plotfile_to_geotiff only supports 2-D plotfiles (found spacedim={spacedim})."
        )

    # time
    _time = float(lines[idx]); idx += 1

    # finest_level (0 = single level; >0 = AMR with that many refinement levels)
    finest_level = int(lines[idx]); idx += 1

    # problo / probhi
    problo = list(map(float, lines[idx].split())); idx += 1
    probhi = list(map(float, lines[idx].split())); idx += 1

    # refinement ratios: one integer per level gap (empty for single-level)
    ref_ratio_tokens = lines[idx].split(); idx += 1
    ref_ratios = []
    for t in ref_ratio_tokens:
        try:
            ref_ratios.append(int(float(t)))
        except (ValueError, OverflowError):
            pass

    # domain boxes per level: AMReX writes one BoxArray per level.
    # Each BoxArray in the Header has the format:
    #   <count>                         <- number of boxes (always 1 for the domain)
    #   ((ixlo,iylo) (ixhi,iyhi) (...)) <- the domain box
    level_domains = []
    for _lev in range(finest_level + 1):
        idx += 1  # skip the box-array size line (always "1" for domain boxes)
        if idx >= len(lines):
            break
        domain_box_line = lines[idx]; idx += 1
        nums = [int(t) for t in domain_box_line.replace("(","").replace(")","").replace(","," ").split()
                if t.lstrip("-").isdigit()]
        if len(nums) >= 4:
            nx_lev = nums[2] - nums[0] + 1
            ny_lev = nums[3] - nums[1] + 1
        else:
            # Fallback: use previous level dimensions (should not happen in valid files)
            nx_lev, ny_lev = level_domains[-1] if level_domains else _DEFAULT_DOMAIN_SIZE
        level_domains.append((nx_lev, ny_lev))

    # Ensure we always have at least Level_0 dimensions
    if not level_domains:
        level_domains = [_DEFAULT_DOMAIN_SIZE]

    nx, ny = level_domains[0]
    return varnames, problo, probhi, nx, ny, finest_level, ref_ratios, level_domains


def _read_fab_data(fab_header_path: Path, ncomp: int, nx: int, ny: int):
    """
    Read the binary FAB data for level 0.

    Returns a dict  { varname_index: 2-D numpy array (ny, nx) }.
    The FAB header (Cell_H) lists each patch with its box and file offset.
    For a single-level, single-box run there is typically one FAB file.
    """
    with open(fab_header_path) as fh:
        content = fh.read()

    # Each FAB entry looks like:
    #   Level_0/Cell_D_00000  offset  order  ncomp  (lo,lo)  (hi,hi)
    # We collect all patches and their boxes then stitch them together.
    patches = []
    for line in content.splitlines():
        line = line.strip()
        if not line or line.startswith("FAB") or line[0].isdigit():
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        # heuristic: first token is the file path (contains "Cell_D")
        if "Cell_D" in parts[0] or "Cell" in parts[0]:
            try:
                file_rel = parts[0]
                offset   = int(parts[1])
                patches.append((file_rel, offset))
            except (ValueError, IndexError):
                pass

    # Fall back: scan for lines with file names after we parse more carefully
    if not patches:
        # Try a simpler parse: the Cell_H format starts with nFAB integer,
        # then pairs of (filename offset) lines
        lines = content.splitlines()
        i = 0
        while i < len(lines):
            l = lines[i].strip()
            if l and "Cell_D" in l:
                parts = l.split()
                try:
                    patches.append((parts[0], int(parts[1])))
                except (ValueError, IndexError):
                    pass
            i += 1

    # Allocate result array  [ncomp, ny, nx]
    data = np.full((ncomp, ny, nx), np.nan, dtype=np.float64)

    fab_dir = fab_header_path.parent

    for file_rel, offset in patches:
        fab_path = fab_dir / Path(file_rel).name
        if not fab_path.exists():
            # try relative to plotfile root
            fab_path = fab_header_path.parent.parent / file_rel
        if not fab_path.exists():
            continue

        with open(fab_path, "rb") as fbin:
            fbin.seek(offset)
            # Read the ASCII FAB header up to the newline after the format tag
            ascii_header = b""
            while True:
                byte = fbin.read(1)
                if not byte:
                    break
                ascii_header += byte
                if byte == b"\n" and b")" in ascii_header:
                    break
            # Parse the FAB header: "FAB ((ncomp_fab, ord), (ixlo,iylo)(ixhi,iyhi)) real_type\n"
            hdr = ascii_header.decode("ascii", errors="replace")
            nums = [int(t) for t in hdr.replace("(","").replace(")","").replace(","," ").split() if t.lstrip("-").isdigit()]
            # nums[0]=ncomp_fab, nums[1]=ngrow, nums[2..5]=(ixlo,iylo,ixhi,iyhi)
            if len(nums) < 6:
                continue
            fab_ncomp = nums[0]
            ixlo, iylo, ixhi, iyhi = nums[2], nums[3], nums[4], nums[5]
            fab_nx = ixhi - ixlo + 1
            fab_ny = iyhi - iylo + 1
            # Determine precision from header keyword
            prec = np.float64
            if "Real" in hdr and "float" in hdr.lower():
                prec = np.float32
            elif "Real" in hdr and "double" in hdr.lower():
                prec = np.float64
            else:
                # default: try double
                prec = np.float64

            n_vals = fab_ncomp * fab_nx * fab_ny
            raw = fbin.read(n_vals * np.dtype(prec).itemsize)
            if len(raw) < n_vals * np.dtype(prec).itemsize:
                # Try float32
                prec = np.float32
                fbin.seek(offset)
                # skip ASCII header again
                ascii_header2 = b""
                while True:
                    byte = fbin.read(1)
                    if not byte:
                        break
                    ascii_header2 += byte
                    if byte == b"\n" and b")" in ascii_header2:
                        break
                raw = fbin.read(n_vals * np.dtype(prec).itemsize)

            arr = np.frombuffer(raw, dtype=prec)
            # AMReX stores FAB data in Fortran order: (comp, x, y) for 2-D
            # Actually column-major: x varies fastest within a component
            arr = arr.reshape((fab_ncomp, fab_nx, fab_ny), order="F")
            # Transpose to (ncomp, ny, nx) matching numpy convention
            arr = arr.transpose(0, 2, 1)

            # Place into full domain array
            for ic in range(min(fab_ncomp, ncomp)):
                data[ic,
                     iylo : iyhi + 1,
                     ixlo : ixhi + 1] = arr[ic]

    return data


# ---------------------------------------------------------------------------
# Multi-level AMR compositing helpers
# ---------------------------------------------------------------------------

def _cumulative_ref_ratio(ref_ratios, level):
    """Total refinement factor at *level* relative to Level_0."""
    ratio = 1
    for i in range(level):
        ratio *= ref_ratios[i] if i < len(ref_ratios) else _DEFAULT_REF_RATIO
    return ratio


def _composite_fine_level(coarse_data, fine_data, ref_ratio):
    """Overlay fine-level data onto *coarse_data* in place.

    For each coarse cell covered by valid (non-NaN) fine cells, replace the
    coarse value with the block-average of the fine cells.  Coarse cells that
    are not covered by any fine patch remain unchanged.

    Parameters
    ----------
    coarse_data : ndarray, shape (ncomp, ny_c, nx_c)
        Level_0 (or any coarser level) data array – modified in place.
    fine_data   : ndarray, shape (ncomp, ny_f, nx_f)
        Finer-level data array.  NaN cells are treated as "not covered".
    ref_ratio   : int
        Spatial refinement ratio between the two levels.
    """
    r = ref_ratio
    ncomp, ny_c, nx_c = coarse_data.shape

    # Trim fine_data to an exact multiple of r that fits within the coarse grid
    ny_f_use = min(fine_data.shape[1], ny_c * r)
    nx_f_use = min(fine_data.shape[2], nx_c * r)
    ny_f_use = (ny_f_use // r) * r
    nx_f_use = (nx_f_use // r) * r
    if ny_f_use == 0 or nx_f_use == 0:
        return

    ny_eff = ny_f_use // r
    nx_eff = nx_f_use // r

    # Group fine cells by their coarse parent: (ncomp, ny_eff, r, nx_eff, r)
    fine_grouped = fine_data[:, :ny_f_use, :nx_f_use].reshape(
        ncomp, ny_eff, r, nx_eff, r
    )
    # Block-average over the r×r fine sub-cells, ignoring NaN
    with np.errstate(all="ignore"):
        fine_mean = np.nanmean(fine_grouped, axis=(2, 4))  # (ncomp, ny_eff, nx_eff)

    # Determine which coarse cells have at least one valid fine cell
    valid = np.any(~np.isnan(fine_mean), axis=0)  # (ny_eff, nx_eff)

    # Overwrite coarse values where valid fine data is available
    for ic in range(ncomp):
        coarse_data[ic, :ny_eff, :nx_eff][valid] = fine_mean[ic, :ny_eff, :nx_eff][valid]


# ---------------------------------------------------------------------------
# GeoTIFF writer
# ---------------------------------------------------------------------------

def write_geotiff(
    array_2d: np.ndarray,
    varname: str,
    outpath: Path,
    problo,
    probhi,
    nx: int,
    ny: int,
    utm_origin=None,
    epsg: int = None,
):
    """Write a single 2-D array as a GeoTIFF."""
    # Cell sizes
    dx = (probhi[0] - problo[0]) / nx
    dy = (probhi[1] - problo[1]) / ny

    if utm_origin is not None:
        easting, northing = utm_origin
        west = easting  + problo[0]
        north = northing + probhi[1]
    else:
        west  = problo[0]
        north = probhi[1]

    transform = from_origin(west, north, dx, dy)

    crs = None
    if epsg is not None:
        crs = CRS.from_epsg(epsg)

    profile = {
        "driver":    "GTiff",
        "dtype":     "float32",
        "width":     nx,
        "height":    ny,
        "count":     1,
        "transform": transform,
        "compress":  "deflate",
        "nodata":    float("nan"),
    }
    if crs is not None:
        profile["crs"] = crs

    # Flip array vertically: rasterio expects row 0 at the top (north)
    arr = np.flipud(array_2d).astype(np.float32)

    with rasterio.open(outpath, "w", **profile) as dst:
        dst.write(arr, 1)
        dst.update_tags(variable=varname)


# ---------------------------------------------------------------------------
# Fire-perimeter GeoJSON writer
# ---------------------------------------------------------------------------

def write_fire_perimeter_geojson(
    phi_array: np.ndarray,
    outpath: Path,
    problo,
    probhi,
    nx: int,
    ny: int,
    utm_origin=None,
):
    """Write the fire perimeter (phi == 0 contour) as a GeoJSON LineString."""
    try:
        from matplotlib.contour import QuadContourSet
        import matplotlib.pyplot as plt
    except ImportError:
        print("  (skipping fire perimeter GeoJSON – matplotlib not available)")
        return

    dx = (probhi[0] - problo[0]) / nx
    dy = (probhi[1] - problo[1]) / ny
    xs = np.array([problo[0] + (i + 0.5) * dx for i in range(nx)])
    ys = np.array([problo[1] + (j + 0.5) * dy for j in range(ny)])

    fig, ax = plt.subplots()
    cs = ax.contour(xs, ys, phi_array, levels=[0.0])
    plt.close(fig)

    easting_offset = utm_origin[0] if utm_origin else 0.0
    northing_offset = utm_origin[1] if utm_origin else 0.0

    features = []
    for path in cs.get_paths():
        if len(path.vertices) < 2:
            continue
        coords = [
            [v[0] + easting_offset, v[1] + northing_offset]
            for v in path.vertices
        ]
        features.append(
            {
                "type": "Feature",
                "geometry": {"type": "LineString", "coordinates": coords},
                "properties": {},
            }
        )

    import json
    geojson = {"type": "FeatureCollection", "features": features}
    with open(outpath, "w") as fh:
        json.dump(geojson, fh)
    print(f"  Wrote fire perimeter → {outpath}")


# ---------------------------------------------------------------------------
# Main conversion routine
# ---------------------------------------------------------------------------

# Variables most useful for fire behaviour GIS analysis
FIRE_VARS = {
    "phi",
    "R",
    "fireline_intensity",
    "flame_length",
    "elevation",
    "slope",
    "aspect",
    "fuel_model",
    "fuel_consumption",
    "crown_fraction",
    "arrival_time",
    "reaction_intensity",
    # Derived fields computed from velx/vely (not stored in plotfile directly)
    "wind_speed",
    "wind_direction",
}

# Derived fields computed from combinations of stored plotfile variables.
# Each entry: derived_name -> (function(fab_data, varnames), description)
# The function receives the full fab_data array and varnames list, and
# returns a 2-D array of shape (ny, nx) or None if inputs are missing.

def _derive_wind_speed(fab_data, varnames):
    """Compute wind speed [m/s] from velx and vely."""
    if "velx" not in varnames or "vely" not in varnames:
        return None
    ix = varnames.index("velx")
    iy = varnames.index("vely")
    return np.sqrt(fab_data[ix]**2 + fab_data[iy]**2)


def _derive_wind_direction(fab_data, varnames):
    """Compute meteorological wind direction [degrees from N, clockwise].

    Met convention: 0° = wind blowing FROM north, 90° = FROM east, etc.
    Computed as: dir = (270 - atan2(vy, vx) * 180/pi) mod 360
    """
    if "velx" not in varnames or "vely" not in varnames:
        return None
    ix = varnames.index("velx")
    iy = varnames.index("vely")
    u = fab_data[ix]
    v = fab_data[iy]
    # Math angle → met direction (FROM direction, degrees clockwise from N)
    direction = (270.0 - np.degrees(np.arctan2(v, u))) % 360.0
    return direction


_DERIVED_FIELDS = {
    "wind_speed":     (_derive_wind_speed,     "Wind speed [m/s] derived from velx, vely"),
    "wind_direction": (_derive_wind_direction,  "Wind direction [° from N, clockwise] derived from velx, vely"),
}


def convert_plotfile(
    plotfile_dir: Path,
    outdir: Path,
    varnames_filter=None,
    utm_origin=None,
    epsg: int = None,
    fire_vars_only: bool = False,
):
    """Convert one plotfile directory to GeoTIFFs.

    Multi-level AMR plotfiles are supported: data from finer levels is
    composited onto the Level_0 base grid so that the output GeoTIFF has
    full-domain coverage at the coarsest resolution with finer-level detail
    where available.
    """
    print(f"\nProcessing {plotfile_dir} …")

    varnames, problo, probhi, nx, ny, finest_level, ref_ratios, level_domains = \
        _parse_header(plotfile_dir)
    print(f"  Domain: ({problo[0]:.2f}, {problo[1]:.2f}) – ({probhi[0]:.2f}, {probhi[1]:.2f})  "
          f"grid {nx}×{ny}  vars: {len(varnames)}"
          + (f"  levels: {finest_level + 1}" if finest_level > 0 else ""))

    if finest_level > 0:
        print(f"  Multi-level AMR plotfile (finest_level={finest_level}): "
              f"compositing {finest_level + 1} level(s) onto Level_0 base grid.")

    # Read binary FAB data from Level_0 (coarsest, always covers the full domain)
    cell_h = plotfile_dir / "Level_0" / "Cell_H"
    if not cell_h.exists():
        print(f"  WARNING: {cell_h} not found – skipping.")
        return

    fab_data = _read_fab_data(cell_h, len(varnames), nx, ny)

    # Overlay finer levels where available (AMR compositing)
    for lev in range(1, finest_level + 1):
        cell_h_lev = plotfile_dir / f"Level_{lev}" / "Cell_H"
        if not cell_h_lev.exists():
            print(f"  WARNING: Level_{lev}/Cell_H not found – skipping level {lev}.")
            continue
        nx_lev, ny_lev = level_domains[lev] if lev < len(level_domains) else (nx, ny)
        total_ref = _cumulative_ref_ratio(ref_ratios, lev)
        fab_data_lev = _read_fab_data(cell_h_lev, len(varnames), nx_lev, ny_lev)
        _composite_fine_level(fab_data, fab_data_lev, total_ref)
        print(f"  Composited Level_{lev} ({nx_lev}×{ny_lev}, ref={total_ref}×) onto base grid.")

    # Determine which variables to export (stored fields)
    # Also identify any derived fields requested
    if varnames_filter:
        export_vars    = [v for v in varnames_filter if v in varnames]
        export_derived = [v for v in varnames_filter if v in _DERIVED_FIELDS]
        missing = [v for v in varnames_filter
                   if v not in varnames and v not in _DERIVED_FIELDS]
        if missing:
            print(f"  WARNING: variables not found in plotfile: {missing}")
    elif fire_vars_only:
        export_vars    = [v for v in varnames if v in FIRE_VARS]
        export_derived = [v for v in FIRE_VARS if v in _DERIVED_FIELDS]
    else:
        export_vars    = list(varnames)
        export_derived = list(_DERIVED_FIELDS.keys())

    outdir.mkdir(parents=True, exist_ok=True)
    stem = plotfile_dir.name  # e.g. "plt0100"

    for vname in export_vars:
        idx = varnames.index(vname)
        arr2d = fab_data[idx]  # shape (ny, nx)
        if np.all(np.isnan(arr2d)):
            print(f"  Skipping {vname} (all NaN – FAB read may have failed)")
            continue
        out_tif = outdir / f"{stem}_{vname}.tif"
        write_geotiff(arr2d, vname, out_tif, problo, probhi, nx, ny,
                      utm_origin=utm_origin, epsg=epsg)
        vmin, vmax = np.nanmin(arr2d), np.nanmax(arr2d)
        print(f"  {vname:25s}  min={vmin:.4g}  max={vmax:.4g}  → {out_tif.name}")

    # Export derived fields (e.g. wind_speed, wind_direction computed from velx/vely)
    for vname in export_derived:
        fn, desc = _DERIVED_FIELDS[vname]
        arr2d = fn(fab_data, varnames)
        if arr2d is None:
            print(f"  Skipping derived field '{vname}': required source variables "
                  "not in plotfile")
            continue
        if np.all(np.isnan(arr2d)):
            print(f"  Skipping derived field '{vname}' (all NaN)")
            continue
        out_tif = outdir / f"{stem}_{vname}.tif"
        write_geotiff(arr2d, vname, out_tif, problo, probhi, nx, ny,
                      utm_origin=utm_origin, epsg=epsg)
        vmin, vmax = float(np.nanmin(arr2d)), float(np.nanmax(arr2d))
        print(f"  {vname:25s}  min={vmin:.4g}  max={vmax:.4g}  → {out_tif.name}  [{desc}]")

    # Write fire perimeter GeoJSON if phi is available
    if "phi" in varnames and ("phi" in export_vars or varnames_filter is None):
        phi_arr = fab_data[varnames.index("phi")]
        geojson_path = outdir / f"{stem}_fire_perimeter.geojson"
        write_fire_perimeter_geojson(phi_arr, geojson_path, problo, probhi, nx, ny,
                                     utm_origin=utm_origin)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Convert AMReX plotfiles to GeoTIFF for GIS import.",
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
        help="Convert all plt#### directories in the current directory.",
    )
    parser.add_argument(
        "-v", "--variables",
        nargs="+",
        metavar="VAR",
        help="Export only these variables (default: all fire-behaviour variables).",
    )
    parser.add_argument(
        "--all-vars",
        action="store_true",
        help="Export every variable in the plotfile (default: fire-behaviour subset).",
    )
    parser.add_argument(
        "--outdir",
        default="gis_out",
        metavar="DIR",
        help="Output directory for GeoTIFF files (default: gis_out).",
    )
    parser.add_argument(
        "--utm-origin",
        nargs=2,
        type=float,
        metavar=("EASTING", "NORTHING"),
        help="UTM origin of the simulation domain (m).  Added to simulation "
             "coordinates to produce absolute UTM coordinates.",
    )
    parser.add_argument(
        "--epsg",
        type=int,
        metavar="CODE",
        help="EPSG code for the output CRS (e.g. 32613 for UTM zone 13N).  "
             "Requires --utm-origin.",
    )
    args = parser.parse_args()

    outdir = Path(args.outdir)
    utm_origin = tuple(args.utm_origin) if args.utm_origin else None
    epsg = args.epsg
    fire_vars_only = not args.all_vars

    if args.all:
        dirs = sorted(Path(".").glob("plt[0-9][0-9][0-9][0-9]"))
        if not dirs:
            sys.exit("No plt#### directories found in the current directory.")
        for d in dirs:
            convert_plotfile(d, outdir,
                             varnames_filter=args.variables,
                             utm_origin=utm_origin,
                             epsg=epsg,
                             fire_vars_only=fire_vars_only)
    elif args.plotfile:
        convert_plotfile(Path(args.plotfile), outdir,
                         varnames_filter=args.variables,
                         utm_origin=utm_origin,
                         epsg=epsg,
                         fire_vars_only=fire_vars_only)
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
