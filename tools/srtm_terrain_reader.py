#!/usr/bin/env python3
"""
srtm_terrain_reader.py - Download SRTM elevation data and write a UTM terrain XYZ file.

Downloads SRTM1 elevation data for a lat/lon bounding box using the ``elevation``
Python package and writes a UTM terrain XYZ file suitable for use as
``rothermel.terrain_file`` in the wildfire_levelset solver.

Usage examples
--------------
  # Download SRTM and write terrain.xyz
  python3 srtm_terrain_reader.py \\
      --lat-min 40 --lat-max 40.5 \\
      --lon-min -105 --lon-max -104.5

  # Save to a specific file with subsampling
  python3 srtm_terrain_reader.py \\
      --lat-min 40 --lat-max 40.5 \\
      --lon-min -105 --lon-max -104.5 \\
      --terrain terrain_out.xyz \\
      --subsample 2

  # Keep intermediate SRTM GeoTIFF
  python3 srtm_terrain_reader.py \\
      --lat-min 40 --lat-max 40.5 \\
      --lon-min -105 --lon-max -104.5 \\
      --tif srtm_raw.tif

Options
-------
  --lat-min / --lat-max   Latitude  bounds (WGS-84 degrees; required).
  --lon-min / --lon-max   Longitude bounds (WGS-84 degrees; required).
  --terrain FILE          Output terrain XYZ file (default: terrain.xyz).
  --tif FILE              Intermediate SRTM GeoTIFF path (temp file if not given).
  --subsample N           Keep every N-th point (default: 1).
  --help                  Show this message and exit.

Requires: pip install elevation rasterio numpy pyproj scipy
"""

import argparse
import math
import os
import sys
import tempfile


# ===========================================================================
# Shared projection helpers
# ===========================================================================

def _latlon_to_utm(lats, lons):
    """Project (lat, lon) arrays (WGS-84) to UTM metres.

    The UTM zone is determined from the centre of the domain.
    Returns (x_utm, y_utm) both in metres.
    """
    try:
        from pyproj import Transformer
    except ImportError:
        raise ImportError(
            "pyproj is required for UTM projection. "
            "Install with: pip install pyproj"
        )
    import numpy as np

    lats = np.asarray(lats, dtype=np.float64)
    lons = np.asarray(lons, dtype=np.float64)

    center_lat = float(np.mean(lats))
    center_lon = float(np.mean(lons))
    zone = int((center_lon + 180.0) / 6.0) + 1
    hemisphere = "north" if center_lat >= 0.0 else "south"
    epsg = 32600 + zone if hemisphere == "north" else 32700 + zone

    transformer = Transformer.from_crs("EPSG:4326", f"EPSG:{epsg}", always_xy=True)
    x_utm, y_utm = transformer.transform(lons.ravel(), lats.ravel())
    return (
        np.array(x_utm, dtype=np.float64).reshape(lats.shape),
        np.array(y_utm, dtype=np.float64).reshape(lats.shape),
    )


# ===========================================================================
# SRTM terrain helpers
# ===========================================================================

def download_srtm(lat_min, lat_max, lon_min, lon_max, out_tif):
    """Download SRTM1 elevation data clipped to the given bounding box."""
    try:
        import elevation
    except ImportError:
        raise ImportError(
            "elevation is required to download SRTM data. "
            "Install with: pip install elevation"
        )
    elevation.clean()
    elevation.clip(
        bounds=(lon_min, lat_min, lon_max, lat_max),
        output=out_tif,
        product="SRTM1",
    )


def _fix_srtm_zeros(z):
    """Replace spurious zero-elevation pixels with interpolated or mean values.

    SRTM tiles sometimes contain random zero values at tile borders that are
    not genuine sea-level points.

    * Isolated spurious zeros (at least one immediate non-zero neighbour) are
      filled by linear interpolation from the nearest valid neighbours using
      ``scipy.interpolate.griddata``.
    * Clusters of nearby zeros where the surrounding window at a larger scale
      contains non-zero values (i.e. the zeros are abrupt relative to the rest
      of the grid) are replaced with the mean elevation of the entire grid.
      This handles the case where multiple adjacent zero pixels prevent the
      local-neighbour check from identifying them as spurious.

    Parameters
    ----------
    z : numpy.ndarray (2-D, float)
        Elevation grid, potentially containing spurious zeros.

    Returns
    -------
    numpy.ndarray
        Copy of *z* with spurious zeros replaced.
    """
    import numpy as np

    zero_mask = z == 0.0
    if not np.any(zero_mask):
        return z

    from scipy.ndimage import generic_filter

    spurious = generic_filter(
        z,
        lambda p: float((p[len(p) // 2] == 0.0) and np.any(p > 0.0)),
        size=3,
        mode="nearest",
    ).astype(bool)

    z_out = z.copy()
    h, w = z.shape
    rows, cols = np.mgrid[0:h, 0:w]

    if np.any(spurious):
        print(f"  Fixing {int(np.sum(spurious))} spurious zero elevation pixels …")
        valid = ~spurious
        try:
            from scipy.interpolate import griddata

            pts_valid = np.column_stack([rows[valid].ravel(), cols[valid].ravel()])
            pts_spurious = np.column_stack(
                [rows[spurious].ravel(), cols[spurious].ravel()]
            )
            z_out[spurious] = griddata(
                pts_valid, z[valid].ravel(), pts_spurious, method="linear",
                fill_value=0.0,
            )
        except ImportError:
            for r, c in zip(rows[spurious].ravel(), cols[spurious].ravel()):
                neighbours = []
                for dr in (-1, 0, 1):
                    for dc in (-1, 0, 1):
                        if dr == 0 and dc == 0:
                            continue
                        rr, cc = r + dr, c + dc
                        if 0 <= rr < h and 0 <= cc < w and z[rr, cc] > 0.0:
                            neighbours.append(z[rr, cc])
                if neighbours:
                    z_out[r, c] = float(np.mean(neighbours))

    remaining_zeros = z_out == 0.0
    if np.any(remaining_zeros) and np.any(z_out > 0.0):
        grid_mean = float(np.mean(z_out[z_out > 0.0]))
        abrupt = generic_filter(
            z_out,
            lambda p: float((p[len(p) // 2] == 0.0) and np.any(p > 0.0)),
            size=11,
            mode="nearest",
        ).astype(bool)
        if np.any(abrupt):
            print(
                f"  Fixing {int(np.sum(abrupt))} abrupt zero-cluster pixels "
                f"with grid mean ({grid_mean:.1f} m) …"
            )
            z_out[abrupt] = grid_mean

    return z_out


def _smooth_terrain_border(z, fraction=0.2):
    """Smooth the outer *fraction* of the terrain grid to reduce tile-seam artefacts.

    A Gaussian-smoothed version of the full grid is blended into the border
    region.  The blend weight transitions linearly from 1 (fully smoothed) at
    the grid edges to 0 (original values) at the interior boundary of the
    border band, so there is no hard discontinuity.

    Parameters
    ----------
    z : numpy.ndarray (2-D, float)
        Elevation grid.
    fraction : float
        Width of the border band expressed as a fraction of the grid dimension
        (default 0.2 → outer 20 % on each side).

    Returns
    -------
    numpy.ndarray
        Smoothed copy of *z*.
    """
    import numpy as np

    try:
        from scipy.ndimage import gaussian_filter
    except ImportError:
        print(
            "WARNING: scipy not available; skipping terrain border smoothing.",
            file=sys.stderr,
        )
        return z

    h, w = z.shape
    border_h = max(1, int(h * fraction))
    border_w = max(1, int(w * fraction))

    sigma = max(1.0, min(border_h, border_w) / 4.0)
    z_smooth = gaussian_filter(z, sigma=sigma)

    weight = np.zeros((h, w), dtype=np.float64)

    ramp_v = np.linspace(1.0, 0.0, border_h, endpoint=False)
    weight[:border_h, :] = np.maximum(weight[:border_h, :],
                                       ramp_v[:, np.newaxis])
    weight[-border_h:, :] = np.maximum(weight[-border_h:, :],
                                        ramp_v[::-1, np.newaxis])

    ramp_h = np.linspace(1.0, 0.0, border_w, endpoint=False)
    weight[:, :border_w] = np.maximum(weight[:, :border_w],
                                       ramp_h[np.newaxis, :])
    weight[:, -border_w:] = np.maximum(weight[:, -border_w:],
                                        ramp_h[np.newaxis, ::-1])

    return z * (1.0 - weight) + z_smooth * weight


def tiff_to_xyz_utm(tif_path, subsample=1):
    """Read a GeoTIFF and return (utm_x, utm_y, z) 2-D arrays in UTM metres."""
    import numpy as np
    try:
        import rasterio
        from rasterio.transform import xy as rasterio_xy
    except ImportError:
        raise ImportError(
            "rasterio is required to read GeoTIFF files. "
            "Install with: pip install rasterio"
        )
    try:
        from pyproj import Transformer
    except ImportError:
        raise ImportError(
            "pyproj is required for UTM projection. "
            "Install with: pip install pyproj"
        )

    with rasterio.open(tif_path) as src:
        z = src.read(1).astype(float)
        z = np.where(z < 0, 0.0, z)

        h, w = z.shape
        rows_idx, cols_idx = np.meshgrid(np.arange(h), np.arange(w), indexing="ij")
        lon, lat = rasterio_xy(src.transform, rows_idx.ravel(), cols_idx.ravel())
        lon = np.array(lon, dtype=np.float64).reshape(h, w)
        lat = np.array(lat, dtype=np.float64).reshape(h, w)

    z = _fix_srtm_zeros(z)
    z = _smooth_terrain_border(z, fraction=0.2)

    center_lon = float(np.nanmean(lon))
    center_lat = float(np.nanmean(lat))
    zone = int((center_lon + 180.0) / 6.0) + 1
    epsg = (32600 + zone) if center_lat >= 0.0 else (32700 + zone)

    transformer = Transformer.from_crs("EPSG:4326", f"EPSG:{epsg}", always_xy=True)
    utm_x_flat, utm_y_flat = transformer.transform(lon.ravel(), lat.ravel())
    utm_x = np.array(utm_x_flat, dtype=np.float64).reshape(h, w)
    utm_y = np.array(utm_y_flat, dtype=np.float64).reshape(h, w)

    if subsample > 1:
        utm_x = utm_x[::subsample, ::subsample]
        utm_y = utm_y[::subsample, ::subsample]
        z     = z    [::subsample, ::subsample]

    return utm_x, utm_y, z


def write_terrain_xyz(utm_x, utm_y, z, output_path):
    """Write a terrain XYZ file compatible with ``rothermel.terrain_file``."""
    import numpy as np

    data = np.column_stack([
        np.asarray(utm_x).ravel(),
        np.asarray(utm_y).ravel(),
        np.asarray(z).ravel(),
    ])
    with open(output_path, "w") as fh:
        fh.write("# X Y Z (meters)\n")
        for row in data:
            fh.write(f"{row[0]:.3f} {row[1]:.3f} {row[2]:.3f}\n")
    print(f"Wrote {len(data)} terrain points to '{output_path}'.")


def create_terrain_xyz(output_path, lat_min, lat_max, lon_min, lon_max,
                       tif_path=None, subsample=1):
    """Download SRTM data and write a UTM terrain XYZ file."""
    _tmp_tif = tif_path is None
    if _tmp_tif:
        tmp = tempfile.NamedTemporaryFile(suffix=".tif", delete=False)
        tif_path = tmp.name
        tmp.close()

    try:
        print(f"Downloading SRTM data for bbox "
              f"({lat_min},{lon_min}) – ({lat_max},{lon_max}) …")
        download_srtm(lat_min, lat_max, lon_min, lon_max, tif_path)

        print("Converting SRTM GeoTIFF → UTM XYZ …")
        utm_x, utm_y, z = tiff_to_xyz_utm(tif_path, subsample=subsample)

        write_terrain_xyz(utm_x, utm_y, z, output_path)
    finally:
        if _tmp_tif and os.path.isfile(tif_path):
            os.remove(tif_path)


def create_terrain_xyz_return_grid(lat_min, lat_max, lon_min, lon_max,
                                   output_path, tif_path=None, subsample=1):
    """Download SRTM, write terrain XYZ, and return (utm_x, utm_y, z) arrays.

    Like ``create_terrain_xyz`` but also returns the 2-D grid arrays so they
    can be used for wind interpolation.
    """
    _tmp_tif = tif_path is None
    if _tmp_tif:
        tmp = tempfile.NamedTemporaryFile(suffix=".tif", delete=False)
        tif_path = tmp.name
        tmp.close()

    try:
        print(f"Downloading SRTM data for bbox "
              f"({lat_min},{lon_min}) – ({lat_max},{lon_max}) …")
        download_srtm(lat_min, lat_max, lon_min, lon_max, tif_path)

        print("Converting SRTM GeoTIFF → UTM XYZ …")
        utm_x, utm_y, z = tiff_to_xyz_utm(tif_path, subsample=subsample)

        write_terrain_xyz(utm_x, utm_y, z, output_path)
    finally:
        if _tmp_tif and os.path.isfile(tif_path):
            os.remove(tif_path)

    return utm_x, utm_y, z


def _compute_slope_aspect_from_srtm_tif(tif_path):
    """Compute slope and aspect from a downloaded SRTM GeoTIFF.

    Derives slope (degrees from horizontal) and aspect (degrees clockwise
    from North) using finite differences on the elevation array.  Cell sizes
    are converted from geographic degrees to approximate metres at the
    raster's centre latitude so that the gradient has consistent units.

    Spurious zero-elevation pixels are corrected with :func:`_fix_srtm_zeros`
    before computing the gradient (matching the treatment in
    :func:`tiff_to_xyz_utm`).

    Parameters
    ----------
    tif_path : str
        Path to an SRTM GeoTIFF (expected to be in EPSG:4326).

    Returns
    -------
    tuple
        ``(elev_data, slope_data, aspect_data, transform, crs, nodata)``
        where all three data arrays share the same *transform* and *crs* as
        the input raster so they can be passed directly to
        :func:`landscape_writer.assemble_landscape`.
    """
    import numpy as np
    try:
        import rasterio
    except ImportError:
        raise ImportError(
            "rasterio is required to read SRTM GeoTIFF files. "
            "Install with: pip install rasterio"
        )

    with rasterio.open(tif_path) as src:
        elev_data = src.read(1).astype(np.float64)
        transform = src.transform
        crs = src.crs
        nodata = src.nodata

    elev_data = np.where(elev_data < 0, 0.0, elev_data)
    elev_data = _fix_srtm_zeros(elev_data)

    h, w = elev_data.shape

    src_epsg = None
    try:
        src_epsg = crs.to_epsg()
    except (AttributeError, Exception):
        pass

    if src_epsg == 4326 or (src_epsg is None and abs(transform.a) < 1.0):
        center_col = w / 2.0
        center_row = h / 2.0
        lon_c, lat_c = transform * (center_col, center_row)
        lat_rad = math.radians(lat_c)
        m_per_deg_lat = 111320.0
        m_per_deg_lon = 111320.0 * math.cos(lat_rad)
        cell_size_x = abs(transform.a) * m_per_deg_lon
        cell_size_y = abs(transform.e) * m_per_deg_lat
    else:
        cell_size_x = abs(transform.a)
        cell_size_y = abs(transform.e)

    dz_row, dz_col = np.gradient(elev_data, cell_size_y, cell_size_x)

    dz_east = dz_col
    dz_north = -dz_row

    slope_rad = np.arctan(np.sqrt(dz_east ** 2 + dz_north ** 2))
    slope_data = np.degrees(slope_rad)

    aspect_data = (
        np.degrees(np.arctan2(dz_east, dz_north)) + 180.0
    ) % 360.0

    return elev_data, slope_data, aspect_data, transform, crs, nodata


# ===========================================================================
# CLI
# ===========================================================================

def _build_parser():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument("--lat-min", type=float, required=True,
                        help="Minimum latitude (WGS-84 decimal degrees)")
    parser.add_argument("--lat-max", type=float, required=True,
                        help="Maximum latitude (WGS-84 decimal degrees)")
    parser.add_argument("--lon-min", type=float, required=True,
                        help="Minimum longitude (WGS-84 decimal degrees)")
    parser.add_argument("--lon-max", type=float, required=True,
                        help="Maximum longitude (WGS-84 decimal degrees)")

    parser.add_argument("--terrain", default="terrain.xyz",
                        help="Output terrain XYZ file (default: terrain.xyz)")
    parser.add_argument("--tif", default=None, metavar="FILE",
                        help="Path for the intermediate SRTM GeoTIFF "
                             "(temporary file used and removed if not specified)")
    parser.add_argument("--subsample", type=int, default=1, metavar="N",
                        help="Keep every N-th point in each dimension (default: 1)")

    return parser


def main(argv=None):
    parser = _build_parser()
    args = parser.parse_args(argv)

    if args.lat_min >= args.lat_max or args.lon_min >= args.lon_max:
        print(
            "ERROR: bounding box must satisfy "
            "lon-min < lon-max and lat-min < lat-max",
            file=sys.stderr,
        )
        sys.exit(1)

    if args.subsample < 1:
        print("ERROR: --subsample must be >= 1", file=sys.stderr)
        sys.exit(1)

    tif_path = os.path.abspath(args.tif) if args.tif is not None else None
    out_terrain = os.path.abspath(args.terrain)

    create_terrain_xyz(
        out_terrain,
        lat_min=args.lat_min, lat_max=args.lat_max,
        lon_min=args.lon_min, lon_max=args.lon_max,
        tif_path=tif_path,
        subsample=args.subsample,
    )

    # Report UTM extents for use with prob_lo / prob_hi
    try:
        import numpy as np
        data = []
        with open(out_terrain) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.split()
                if len(parts) >= 2:
                    data.append((float(parts[0]), float(parts[1])))
        if data:
            xs = [d[0] for d in data]
            ys = [d[1] for d in data]
            print(
                f"\nUTM extents of terrain data (use for prob_lo / prob_hi):\n"
                f"  prob_lo_x = {min(xs):.0f}\n"
                f"  prob_lo_y = {min(ys):.0f}\n"
                f"  prob_hi_x = {max(xs):.0f}\n"
                f"  prob_hi_y = {max(ys):.0f}\n"
            )
    except Exception:
        pass


if __name__ == "__main__":
    main()
