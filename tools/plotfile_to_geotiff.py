#!/usr/bin/env python3
"""
plotfile_to_geotiff.py – Convert AMReX plotfiles to GeoTIFF rasters for GIS.

The script reads a 2-D AMReX plotfile directory (e.g. ``plt0100``) and writes
one GeoTIFF per variable (or a user-specified subset).  An optional UTM origin
can be supplied so that the output is correctly georeferenced; without it the
coordinates are written in the native simulation units (metres by default).

Requirements
------------
  pip install rasterio numpy

Usage
-----
  # Convert all variables in a plotfile, simulation coordinates only
  python3 tools/plotfile_to_geotiff.py plt0100 --outdir gis_out

  # Specific variables only
  python3 tools/plotfile_to_geotiff.py plt0100 -v phi R fireline_intensity flame_length

  # Georeference with a UTM origin (easting, northing) and EPSG code
  python3 tools/plotfile_to_geotiff.py plt0100 \\
      --utm-origin 450000 4200000 \\
      --epsg 32613 \\
      --outdir gis_out

  # Convert every plt#### directory in the current working directory
  python3 tools/plotfile_to_geotiff.py --all --outdir gis_out

Fire-behaviour variables of interest
-------------------------------------
  phi                 – level-set function (< 0 = burned / on fire)
  R                   – Rothermel rate of spread [m/s]
  fireline_intensity  – Byram (1959) fireline intensity [kW/m]
  flame_length        – Byram (1959) flame length [m]
  elevation           – terrain elevation [m]
  slope               – terrain slope [degrees]
  fuel_model          – FBFM13/FBFM40 fuel model code
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
# AMReX plotfile parser (2-D, single level, no AMR)
# ---------------------------------------------------------------------------

def _parse_header(plotfile_dir: Path):
    """Return (varnames, problo, probhi, nx, ny) from the plotfile Header."""
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

    # finest_level
    finest_level = int(lines[idx]); idx += 1
    if finest_level != 0:
        raise ValueError(
            "plotfile_to_geotiff does not support AMR plotfiles (finest_level > 0)."
        )

    # problo / probhi
    problo = list(map(float, lines[idx].split())); idx += 1
    probhi = list(map(float, lines[idx].split())); idx += 1

    # refinement ratios (one per level gap; skip for single level)
    idx += 1  # skip refinement ratio line

    # domain boxes per level: format is ((lo),(hi))
    # for level 0 only
    idx += 1  # skip the enclosing box-array size line
    domain_box_line = lines[idx]; idx += 1
    # parse "(( x0,y0) (x1,y1) (...))" – extract integers
    nums = [int(t) for t in domain_box_line.replace("(","").replace(")","").replace(","," ").split() if t.lstrip("-").isdigit()]
    nx = nums[2] - nums[0] + 1
    ny = nums[3] - nums[1] + 1

    return varnames, problo, probhi, nx, ny


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
}


def convert_plotfile(
    plotfile_dir: Path,
    outdir: Path,
    varnames_filter=None,
    utm_origin=None,
    epsg: int = None,
    fire_vars_only: bool = False,
):
    """Convert one plotfile directory to GeoTIFFs."""
    print(f"\nProcessing {plotfile_dir} …")

    varnames, problo, probhi, nx, ny = _parse_header(plotfile_dir)
    print(f"  Domain: ({problo[0]:.2f}, {problo[1]:.2f}) – ({probhi[0]:.2f}, {probhi[1]:.2f})  "
          f"grid {nx}×{ny}  vars: {len(varnames)}")

    # Read binary FAB data
    cell_h = plotfile_dir / "Level_0" / "Cell_H"
    if not cell_h.exists():
        print(f"  WARNING: {cell_h} not found – skipping.")
        return

    fab_data = _read_fab_data(cell_h, len(varnames), nx, ny)

    # Determine which variables to export
    if varnames_filter:
        export_vars = [v for v in varnames_filter if v in varnames]
        missing = [v for v in varnames_filter if v not in varnames]
        if missing:
            print(f"  WARNING: variables not found in plotfile: {missing}")
    elif fire_vars_only:
        export_vars = [v for v in varnames if v in FIRE_VARS]
    else:
        export_vars = list(varnames)

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
