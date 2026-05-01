#!/usr/bin/env python3
"""
dem_to_xyz.py - Convert Digital Elevation Map (DEM) files to X Y Z format.

Supported input formats:
  - Arc/Info ASCII Grid (.asc)         - GDAL/rasterio readable rasters (GeoTIFF, .tif/.tiff)
  - SRTM HGT binary files (.hgt)

Output format (matches wildfire_levelset terrain file convention):
  # X Y Z (meters)
  x1  y1  z1
  x2  y2  z2
  ...

Coordinates are output in the native projection of the input file (metres
for projected rasters such as UTM, or decimal degrees for geographic rasters
such as WGS-84 SRTM files).  For geographic rasters an approximate
metre-scaled coordinate system can be requested with --project-utm.

Usage:
  python3 dem_to_xyz.py input.asc output.csv
  python3 dem_to_xyz.py input.tif  output.csv --nodata -9999
  python3 dem_to_xyz.py input.hgt  output.csv --project-utm
  python3 dem_to_xyz.py input.asc  output.csv --subsample 4

Options:
  --nodata VALUE      Override the no-data sentinel value in the file.
  --project-utm       Reproject geographic (lon/lat) coordinates to UTM.
  --subsample N       Keep every N-th point in each dimension (default: 1).
  --help              Show this message and exit.
"""

import sys
import os
import struct
import argparse
import math


# ---------------------------------------------------------------------------
# Arc/Info ASCII Grid reader
# ---------------------------------------------------------------------------

def _read_asc(path):
    """Return (x_flat, y_flat, z_flat) numpy arrays from an Arc/Info ASCII Grid."""
    import numpy as np

    header = {}
    data_lines = []
    header_keys = {"ncols", "nrows", "xllcorner", "xllcenter",
                   "yllcorner", "yllcenter", "cellsize", "nodata_value"}

    with open(path, "r") as fh:
        for line in fh:
            stripped = line.strip()
            if not stripped:
                continue
            parts = stripped.split()
            if parts[0].lower() in header_keys:
                header[parts[0].lower()] = parts[1]
            else:
                # First non-header line — start collecting grid rows
                data_lines.append(stripped)
                break
        for line in fh:
            stripped = line.strip()
            if stripped:
                data_lines.append(stripped)

    ncols = int(header["ncols"])
    nrows = int(header["nrows"])
    cellsize = float(header["cellsize"])
    nodata = float(header.get("nodata_value", -9999))

    # X origin: xllcorner is the left edge; xllcenter is the centre of the SW cell
    if "xllcenter" in header:
        x0 = float(header["xllcenter"])
    else:
        x0 = float(header["xllcorner"]) + 0.5 * cellsize

    if "yllcenter" in header:
        y0 = float(header["yllcenter"])
    else:
        y0 = float(header["yllcorner"]) + 0.5 * cellsize

    grid = np.array([list(map(float, row.split())) for row in data_lines],
                    dtype=np.float64)
    # Arc/Info grids are stored top-to-bottom (first row = northernmost)
    grid = grid[::-1, :]  # flip to bottom-up so row 0 = southernmost

    cols = np.arange(ncols, dtype=np.float64)
    rows = np.arange(nrows, dtype=np.float64)
    xx, yy = np.meshgrid(x0 + cols * cellsize, y0 + rows * cellsize)

    mask = grid != nodata
    return xx[mask], yy[mask], grid[mask], nodata


# ---------------------------------------------------------------------------
# SRTM HGT reader
# ---------------------------------------------------------------------------

def _read_hgt(path):
    """Return (lon_flat, lat_flat, z_flat) from an SRTM1 or SRTM3 HGT file.

    Filename must encode the SW corner, e.g. N37W120.hgt.
    """
    import numpy as np

    basename = os.path.splitext(os.path.basename(path))[0].upper()
    lat_sign = 1 if basename[0] == "N" else -1
    lon_sign = 1 if "E" in basename else -1
    lat_str = basename[1:3]
    lon_str = basename[4:7] if "E" in basename else basename[4:7]
    lat_sw = lat_sign * int(lat_str)
    lon_sw = lon_sign * int(lon_str)

    file_size = os.path.getsize(path)
    if file_size == 2 * 3601 * 3601:
        samples = 3601  # SRTM1 (1 arc-second)
    elif file_size == 2 * 1201 * 1201:
        samples = 1201  # SRTM3 (3 arc-second)
    else:
        raise ValueError(
            f"Unexpected HGT file size {file_size}. "
            "Expected SRTM1 (3601×3601) or SRTM3 (1201×1201)."
        )

    data = np.frombuffer(open(path, "rb").read(), dtype=">i2").astype(np.float64)
    data = data.reshape((samples, samples))
    nodata = -32768.0
    data[data == nodata] = np.nan

    # HGT files are stored top-to-bottom (first row = northernmost)
    data = data[::-1, :]  # flip so row 0 = southernmost

    arc_sec = 1.0 / (samples - 1)  # 1/3600 for SRTM1, 3/3600 for SRTM3
    lons = lon_sw + np.arange(samples) * arc_sec
    lats = lat_sw + np.arange(samples) * arc_sec
    lon2d, lat2d = np.meshgrid(lons, lats)

    mask = ~np.isnan(data)
    return lon2d[mask], lat2d[mask], data[mask], nodata


# ---------------------------------------------------------------------------
# GeoTIFF / generic raster reader (rasterio)
# ---------------------------------------------------------------------------

def _read_rasterio(path, nodata_override=None):
    """Return (x_flat, y_flat, z_flat) from any rasterio-readable raster."""
    import numpy as np
    try:
        import rasterio
        from rasterio.transform import xy as rasterio_xy
    except ImportError:
        raise ImportError(
            "rasterio is required to read GeoTIFF files. "
            "Install with: pip install rasterio"
        )

    with rasterio.open(path) as ds:
        band = ds.read(1).astype(np.float64)
        nodata = nodata_override if nodata_override is not None else ds.nodata
        if nodata is None:
            nodata = -9999.0

        # Build row/col index arrays
        rows_idx, cols_idx = np.where(band != nodata)
        xs, ys = rasterio_xy(ds.transform, rows_idx, cols_idx)
        xs = np.array(xs, dtype=np.float64)
        ys = np.array(ys, dtype=np.float64)
        zs = band[rows_idx, cols_idx]

    return xs, ys, zs, nodata


# ---------------------------------------------------------------------------
# UTM projection helper
# ---------------------------------------------------------------------------

def _project_to_utm(lons, lats):
    """Project geographic lon/lat (WGS-84) arrays to UTM (metres).

    The UTM zone is chosen from the centre of the data extent.
    Returns (x_utm, y_utm) in metres.
    """
    try:
        from pyproj import Transformer
    except ImportError:
        raise ImportError(
            "pyproj is required for UTM projection. "
            "Install with: pip install pyproj"
        )

    import numpy as np
    center_lon = float(np.nanmean(lons))
    center_lat = float(np.nanmean(lats))

    zone = int((center_lon + 180) / 6) + 1
    hemisphere = "north" if center_lat >= 0 else "south"
    epsg = 32600 + zone if hemisphere == "north" else 32700 + zone

    transformer = Transformer.from_crs("EPSG:4326", f"EPSG:{epsg}", always_xy=True)
    x_utm, y_utm = transformer.transform(lons, lats)
    return np.array(x_utm), np.array(y_utm)


# ---------------------------------------------------------------------------
# Subsampling helper
# ---------------------------------------------------------------------------

def _subsample(xs, ys, zs, factor, ncols=None):
    """Keep every *factor*-th point.  If ncols is known the 2-D structure is
    preserved; otherwise a simple stride is applied to the flattened arrays."""
    if factor <= 1:
        return xs, ys, zs
    import numpy as np
    if ncols is not None:
        nrows = len(xs) // ncols
        xs2 = xs.reshape(nrows, ncols)[::factor, ::factor].ravel()
        ys2 = ys.reshape(nrows, ncols)[::factor, ::factor].ravel()
        zs2 = zs.reshape(nrows, ncols)[::factor, ::factor].ravel()
    else:
        xs2 = xs[::factor]
        ys2 = ys[::factor]
        zs2 = zs[::factor]
    return xs2, ys2, zs2


# ---------------------------------------------------------------------------
# Writer
# ---------------------------------------------------------------------------

def _write_xyz(path, xs, ys, zs):
    """Write X Y Z space-separated file with comment header."""
    with open(path, "w") as fh:
        fh.write("# X Y Z (meters)\n")
        for x, y, z in zip(xs, ys, zs):
            fh.write(f"{x:.2f} {y:.2f} {z:.6f}\n")


# ---------------------------------------------------------------------------
# Main conversion routine
# ---------------------------------------------------------------------------

def convert_dem(input_path, output_path,
                nodata_override=None, project_utm=False, subsample=1):
    """Convert a DEM file to the X Y Z format used by wildfire_levelset.

    Parameters
    ----------
    input_path : str
        Path to the input DEM file (.asc, .hgt, .tif, .tiff).
    output_path : str
        Path for the output CSV/text file.
    nodata_override : float or None
        Override the no-data value detected from the file.
    project_utm : bool
        If True, reproject geographic lon/lat to UTM metres before writing.
    subsample : int
        Keep every *subsample*-th point in each dimension (≥1).
    """
    ext = os.path.splitext(input_path)[1].lower()

    if ext == ".asc":
        xs, ys, zs, _ = _read_asc(input_path)
        if nodata_override is not None:
            import numpy as np
            mask = zs != nodata_override
            xs, ys, zs = xs[mask], ys[mask], zs[mask]
    elif ext == ".hgt":
        xs, ys, zs, _ = _read_hgt(input_path)
        import numpy as np
        xs, ys, zs = xs[~np.isnan(zs)], ys[~np.isnan(zs)], zs[~np.isnan(zs)]
    elif ext in (".tif", ".tiff", ".img", ".dem"):
        xs, ys, zs, _ = _read_rasterio(input_path, nodata_override)
    else:
        # Try rasterio as a fallback for unknown extensions
        print(f"Unknown extension '{ext}'; attempting to read with rasterio.")
        xs, ys, zs, _ = _read_rasterio(input_path, nodata_override)

    if project_utm:
        print("Projecting lon/lat coordinates to UTM …")
        xs, ys = _project_to_utm(xs, ys)

    if subsample > 1:
        xs, ys, zs = _subsample(xs, ys, zs, subsample)
        print(f"Subsampled to every {subsample}-th point: {len(xs)} points retained.")

    _write_xyz(output_path, xs, ys, zs)
    print(f"Wrote {len(xs)} terrain points to '{output_path}'.")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _build_parser():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("input",  help="Input DEM file (.asc, .hgt, .tif, .tiff)")
    parser.add_argument("output", help="Output X Y Z text file")
    parser.add_argument(
        "--nodata", type=float, default=None,
        metavar="VALUE",
        help="Override the no-data sentinel value in the input file",
    )
    parser.add_argument(
        "--project-utm", action="store_true",
        help="Reproject geographic lon/lat to UTM metres (requires pyproj)",
    )
    parser.add_argument(
        "--subsample", type=int, default=1, metavar="N",
        help="Keep every N-th point in each dimension (default: 1 = all points)",
    )
    return parser


def main(argv=None):
    parser = _build_parser()
    args = parser.parse_args(argv)

    if not os.path.isfile(args.input):
        print(f"ERROR: input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)
    if args.subsample < 1:
        print("ERROR: --subsample must be >= 1", file=sys.stderr)
        sys.exit(1)

    convert_dem(
        args.input,
        args.output,
        nodata_override=args.nodata,
        project_utm=args.project_utm,
        subsample=args.subsample,
    )


if __name__ == "__main__":
    main()
