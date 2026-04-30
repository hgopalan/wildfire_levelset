#!/usr/bin/env python3
"""
wrf_to_terrain_wind.py - Extract terrain and wind data from a WRF-style netCDF file.

WRF uses a staggered Arakawa C-grid:
  - Mass-point variables (XLAT, XLONG, HGT_M) are defined at the centre of
    each cell: shape (Time, south_north, west_east).
  - U wind component is staggered in the west-east direction:
    shape (Time, bottom_top, south_north, west_east_stag), i.e. nx+1 columns.
  - V wind component is staggered in the south-north direction:
    shape (Time, bottom_top, south_north_stag, west_east), i.e. ny+1 rows.

This script:
  1. Reads the first time step (index 0) unless --time-index is supplied.
  2. Reads the lowest model level (index 0 for bottom_top) for U and V.
  3. Destaggered U and V to mass-point locations (cell centres):
       u_mass[j, i] = 0.5 * (U[j, i] + U[j, i+1])
       v_mass[j, i] = 0.5 * (V[j, i] + V[j+1, i])
  4. Converts XLAT/XLONG to UTM metres using pyproj.
  5. Writes two output files:
       terrain file : utm_x  utm_y  z          (one row per mass point)
       wind file    : utm_x  utm_y  u  v        (one row per mass point)

Usage:
  python3 wrf_to_terrain_wind.py wrfout_d01 terrain.csv wind.csv
  python3 wrf_to_terrain_wind.py wrfout_d01 terrain.csv wind.csv --time-index 2
  python3 wrf_to_terrain_wind.py wrfout_d01 terrain.csv wind.csv --level 1
  python3 wrf_to_terrain_wind.py wrfout_d01 terrain.csv wind.csv --subsample 2
  python3 wrf_to_terrain_wind.py wrfout_d01 terrain.csv wind.csv --help

Options:
  --time-index N   Time index to extract (default: 0).
  --level N        Vertical level index for U and V (default: 0 = lowest).
  --subsample N    Keep every N-th point in each dimension (default: 1).
  --help           Show this message and exit.
"""

import argparse
import os
import sys


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _open_nc(path):
    try:
        import netCDF4 as nc
    except ImportError:
        raise ImportError(
            "netCDF4 is required to read WRF output files. "
            "Install with: pip install netCDF4"
        )
    return nc.Dataset(path, "r")


def _get_var(ds, *candidates):
    """Return the first variable found in *candidates*; raise if none found."""
    for name in candidates:
        if name in ds.variables:
            return ds.variables[name]
    raise KeyError(
        f"None of the expected variables {candidates} found in file. "
        f"Available variables: {list(ds.variables.keys())}"
    )


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


def _destagger_u(U_stag):
    """Average staggered U (shape ny, nx+1) to mass points (shape ny, nx)."""
    return 0.5 * (U_stag[:, :-1] + U_stag[:, 1:])


def _destagger_v(V_stag):
    """Average staggered V (shape ny+1, nx) to mass points (shape ny, nx)."""
    return 0.5 * (V_stag[:-1, :] + V_stag[1:, :])


def _write_terrain(path, x_utm, y_utm, z):
    """Write terrain file: utm_x  utm_y  z (one row per mass point)."""
    import numpy as np
    with open(path, "w") as fh:
        fh.write("# utm_x utm_y z (meters)\n")
        for x, y, zv in zip(x_utm.ravel(), y_utm.ravel(), z.ravel()):
            fh.write(f"{x:.2f} {y:.2f} {float(zv):.6f}\n")


def _write_wind(path, x_utm, y_utm, u, v):
    """Write wind file: utm_x  utm_y  u  v (m/s, one row per mass point)."""
    with open(path, "w") as fh:
        fh.write("# utm_x utm_y u v (m/s)\n")
        for x, y, uv, vv in zip(
            x_utm.ravel(), y_utm.ravel(), u.ravel(), v.ravel()
        ):
            fh.write(f"{x:.2f} {y:.2f} {float(uv):.6f} {float(vv):.6f}\n")


# ---------------------------------------------------------------------------
# Main conversion routine
# ---------------------------------------------------------------------------

def convert_wrf(wrf_path, terrain_out, wind_out,
                time_index=0, level=0, subsample=1):
    """Extract terrain and wind data from a WRF-style netCDF file.

    Parameters
    ----------
    wrf_path : str
        Path to the WRF output netCDF file.
    terrain_out : str
        Output path for the terrain file (utm_x, utm_y, z).
    wind_out : str
        Output path for the wind file (utm_x, utm_y, u, v).
    time_index : int
        Time snapshot index to extract (default 0).
    level : int
        Vertical model level for wind (default 0 = lowest level).
    subsample : int
        Keep every *subsample*-th point in each dimension (default 1).
    """
    import numpy as np

    ds = _open_nc(wrf_path)

    # --- lat/lon at mass points ---
    lat_var = _get_var(ds, "XLAT", "lat", "latitude", "LAT")
    lon_var = _get_var(ds, "XLONG", "lon", "longitude", "LON")

    lat_arr = np.array(lat_var[:], dtype=np.float64)
    lon_arr = np.array(lon_var[:], dtype=np.float64)

    # Remove time dimension if present (WRF stores XLAT as (Time, ny, nx))
    if lat_arr.ndim == 3:
        lat_2d = lat_arr[time_index]
        lon_2d = lon_arr[time_index]
    else:
        lat_2d = lat_arr
        lon_2d = lon_arr

    # --- terrain height ---
    hgt_var = _get_var(ds, "HGT_M", "HGT", "hgt", "TERRAIN", "terrain", "elevation")
    hgt_arr = np.array(hgt_var[:], dtype=np.float64)
    if hgt_arr.ndim == 3:
        hgt_2d = hgt_arr[time_index]
    else:
        hgt_2d = hgt_arr

    # --- U wind (staggered west-east) ---
    u_var = _get_var(ds, "U", "u", "U10", "u10")
    u_arr = np.array(u_var[:], dtype=np.float64)

    # Handle different dimensionalities
    if u_arr.ndim == 4:
        # (Time, bottom_top, south_north, west_east_stag)
        U_2d_stag = u_arr[time_index, level]
    elif u_arr.ndim == 3:
        # (Time, south_north, west_east_stag) — e.g. U10 after surface extraction
        U_2d_stag = u_arr[time_index]
    else:
        U_2d_stag = u_arr

    # --- V wind (staggered south-north) ---
    v_var = _get_var(ds, "V", "v", "V10", "v10")
    v_arr = np.array(v_var[:], dtype=np.float64)

    if v_arr.ndim == 4:
        V_2d_stag = v_arr[time_index, level]
    elif v_arr.ndim == 3:
        V_2d_stag = v_arr[time_index]
    else:
        V_2d_stag = v_arr

    ds.close()

    ny, nx = lat_2d.shape

    # --- Destagger U and V to mass-point locations ---
    # U: shape (ny, nx+1) → (ny, nx)
    if U_2d_stag.shape == (ny, nx + 1):
        u_mass = _destagger_u(U_2d_stag)
    elif U_2d_stag.shape == (ny, nx):
        # Already at mass points (e.g. U10 from some WRF configurations)
        u_mass = U_2d_stag
    else:
        raise ValueError(
            f"Unexpected U shape {U_2d_stag.shape} for domain ({ny}, {nx})"
        )

    # V: shape (ny+1, nx) → (ny, nx)
    if V_2d_stag.shape == (ny + 1, nx):
        v_mass = _destagger_v(V_2d_stag)
    elif V_2d_stag.shape == (ny, nx):
        v_mass = V_2d_stag
    else:
        raise ValueError(
            f"Unexpected V shape {V_2d_stag.shape} for domain ({ny}, {nx})"
        )

    # --- Project to UTM ---
    print(f"Projecting {ny}×{nx} grid to UTM …")
    x_utm, y_utm = _latlon_to_utm(lat_2d, lon_2d)

    # --- Subsample ---
    if subsample > 1:
        sl = slice(None, None, subsample)
        lat_2d  = lat_2d[sl, sl]
        lon_2d  = lon_2d[sl, sl]
        hgt_2d  = hgt_2d[sl, sl]
        u_mass  = u_mass[sl, sl]
        v_mass  = v_mass[sl, sl]
        x_utm   = x_utm[sl, sl]
        y_utm   = y_utm[sl, sl]
        ny_out, nx_out = lat_2d.shape
        print(f"Subsampled to every {subsample}-th point: {ny_out}×{nx_out} = "
              f"{ny_out * nx_out} points retained.")

    # --- Write output files ---
    _write_terrain(terrain_out, x_utm, y_utm, hgt_2d)
    print(f"Wrote terrain file: '{terrain_out}' "
          f"({x_utm.size} points, utm_x utm_y z)")

    _write_wind(wind_out, x_utm, y_utm, u_mass, v_mass)
    print(f"Wrote wind file:    '{wind_out}' "
          f"({x_utm.size} points, utm_x utm_y u v)")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _build_parser():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("wrf_file",    help="Input WRF netCDF file")
    parser.add_argument("terrain_out", help="Output terrain file (utm_x utm_y z)")
    parser.add_argument("wind_out",    help="Output wind file (utm_x utm_y u v)")
    parser.add_argument(
        "--time-index", type=int, default=0, metavar="N",
        help="Time snapshot index to extract (default: 0)",
    )
    parser.add_argument(
        "--level", type=int, default=0, metavar="N",
        help="Vertical level index for U/V extraction (default: 0 = lowest)",
    )
    parser.add_argument(
        "--subsample", type=int, default=1, metavar="N",
        help="Keep every N-th point in each dimension (default: 1 = all points)",
    )
    return parser


def main(argv=None):
    parser = _build_parser()
    args = parser.parse_args(argv)

    if not os.path.isfile(args.wrf_file):
        print(f"ERROR: WRF file not found: {args.wrf_file}", file=sys.stderr)
        sys.exit(1)
    if args.subsample < 1:
        print("ERROR: --subsample must be >= 1", file=sys.stderr)
        sys.exit(1)

    convert_wrf(
        args.wrf_file,
        args.terrain_out,
        args.wind_out,
        time_index=args.time_index,
        level=args.level,
        subsample=args.subsample,
    )


if __name__ == "__main__":
    main()
