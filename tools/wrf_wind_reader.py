#!/usr/bin/env python3
"""
wrf_wind_reader.py - Extract wind data from a WRF netCDF file and write CSV wind files.

Reads U and V wind components from a WRF-style netCDF file and writes
one or more CSV wind files (one per requested time step) suitable for use
as ``velocity_file`` in the wildfire_levelset solver.

Key behaviours
--------------
* The bounding box (lat/lon min/max) is derived automatically from the WRF
  netCDF file: the domain centre is computed and the bounding box is clipped
  to **±0.45 degrees** in each direction.  Explicit ``--lat-min/max/lon-min/max``
  values can optionally be passed to override the auto-derived bounds.
* ``--time-range T1:TN`` extracts WRF time steps T1 through TN **inclusive**
  (0-based indices).  When multiple time steps are requested, output filenames
  are derived from ``--wind`` by inserting ``_N`` before the extension, e.g.
  ``wind.csv``, ``wind_1.csv``, ``wind_2.csv``, …
* ``--time-index N`` (single time step, default 0) is the legacy flag and
  takes precedence when ``--time-range`` is not supplied.
* ``--terrain-file PATH`` or ``--grid-resolution M`` can be specified to
  interpolate the extracted wind onto a specific spatial grid.
* Optionally generates an ``inputs.i`` stub for a FARSITE run.

Usage examples
--------------
  # Extract wind at t=0 (bbox from WRF)
  python3 wrf_wind_reader.py --wrf-file wrfout_d01 --wind wind.csv

  # Extract a range of time steps
  python3 wrf_wind_reader.py --wrf-file wrfout_d01 --wind wind.csv --time-range 0:5

  # Extract at specific lat/lon bounds
  python3 wrf_wind_reader.py \\
      --wrf-file wrfout_d01 --wind wind.csv \\
      --lat-min 40 --lat-max 40.5 \\
      --lon-min -105 --lon-max -104.5

  # Interpolate wind to a terrain XYZ grid
  python3 wrf_wind_reader.py \\
      --wrf-file wrfout_d01 --wind wind.csv \\
      --terrain-file terrain.xyz

  # Interpolate wind to a regular grid at 30 m resolution
  python3 wrf_wind_reader.py \\
      --wrf-file wrfout_d01 --wind wind.csv \\
      --grid-resolution 30

Options
-------
  --wrf-file FILE       WRF netCDF file (required).
  --wind FILE           Output wind CSV file (default: wind.csv).
  --lat-min / --lat-max Latitude  bounds (WGS-84 degrees); optional override.
  --lon-min / --lon-max Longitude bounds (WGS-84 degrees); optional override.
  --time-range T1:TN    Inclusive range of WRF time indices (e.g. '0:4').
  --time-index N        Single WRF time index (default: 0).
  --level N             WRF vertical level for wind (default: 0 = lowest).
  --subsample N         Keep every N-th grid point (default: 1).
  --terrain-file PATH   Existing terrain XYZ file to interpolate wind onto.
  --grid-resolution M   Target grid resolution in metres for regular interpolation.
  --no-inputs           Skip automatic generation of inputs.i.
  --inputs FILE         Output inputs.i filename (default: inputs.i).
  --help                Show this message and exit.

Requires: pip install netCDF4 numpy pyproj
Optional: pip install scipy (for --terrain-file or --grid-resolution interpolation)
"""

import argparse
import os
import sys

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
# WRF helpers
# ===========================================================================

def _open_nc(path):
    """Open a netCDF4 file and return the dataset."""
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


def _destagger_u(U_stag):
    """Average staggered U (shape ny, nx+1) to mass points (shape ny, nx)."""
    return 0.5 * (U_stag[:, :-1] + U_stag[:, 1:])


def _destagger_v(V_stag):
    """Average staggered V (shape ny+1, nx) to mass points (shape ny, nx)."""
    return 0.5 * (V_stag[:-1, :] + V_stag[1:, :])


def _write_wind(path, x_utm, y_utm, u, v):
    """Write wind file: utm_x  utm_y  u  v (m/s, one row per point)."""
    with open(path, "w") as fh:
        fh.write("# utm_x utm_y u v (m/s)\n")
        for x, y, uv, vv in zip(
            x_utm.ravel(), y_utm.ravel(), u.ravel(), v.ravel()
        ):
            fh.write(f"{x:.2f} {y:.2f} {float(uv):.6f} {float(vv):.6f}\n")


_WRF_BBOX_HALF_SPAN = 0.45  # degrees; SRTM download limit is 1°×1°


def read_wrf_bbox(wrf_path):
    """Read a clipped lat/lon bounding box from a WRF netCDF file.

    The bounding box is centred on the domain centre and extends
    ``_WRF_BBOX_HALF_SPAN`` degrees (0.45°) in each direction, keeping the
    total span below the SRTM 1-degree-tile download limit.

    Parameters
    ----------
    wrf_path : str
        Path to the WRF output netCDF file.

    Returns
    -------
    tuple of float
        (lat_min, lat_max, lon_min, lon_max) in WGS-84 decimal degrees,
        each being ``center ± 0.45``.
    """
    import numpy as np

    ds = _open_nc(wrf_path)
    try:
        lat_var = _get_var(ds, "XLAT", "lat", "latitude", "LAT")
        lon_var = _get_var(ds, "XLONG", "XLON", "lon", "longitude", "LON")

        lat_arr = np.array(lat_var[:], dtype=np.float64)
        lon_arr = np.array(lon_var[:], dtype=np.float64)
    finally:
        ds.close()

    if lat_arr.ndim == 3:
        lat_arr = lat_arr[0]
        lon_arr = lon_arr[0]

    center_lat = 0.5 * (float(np.min(lat_arr)) + float(np.max(lat_arr)))
    center_lon = 0.5 * (float(np.min(lon_arr)) + float(np.max(lon_arr)))

    return (
        center_lat - _WRF_BBOX_HALF_SPAN,
        center_lat + _WRF_BBOX_HALF_SPAN,
        center_lon - _WRF_BBOX_HALF_SPAN,
        center_lon + _WRF_BBOX_HALF_SPAN,
    )


def _parse_time_range(time_range_str):
    """Parse a 'T1:TN' string and return a list of integer indices.

    The range is **inclusive** on both ends: '0:3' → [0, 1, 2, 3].
    """
    parts = time_range_str.split(":")
    if len(parts) != 2:
        raise ValueError(
            f"--time-range must be in the form 'T1:TN', got: {time_range_str!r}"
        )
    try:
        t1, tn = int(parts[0]), int(parts[1])
    except ValueError:
        raise ValueError(
            f"--time-range indices must be integers, got: {time_range_str!r}"
        )
    if t1 > tn:
        raise ValueError(
            f"--time-range T1 must be <= TN, got {t1}:{tn}"
        )
    return list(range(t1, tn + 1))


def _wind_output_path(base_path, position):
    """Return the wind output path for *position* (0-based) within a time range.

    position=0 returns *base_path* unchanged (e.g. ``wind.csv``).
    position>0 inserts ``_<position>`` before the file extension
    (e.g. position=1 → ``wind_1.csv``, position=2 → ``wind_2.csv``).
    """
    if position == 0:
        return base_path
    root, ext = os.path.splitext(base_path)
    return f"{root}_{position}{ext}"


def read_wrf_time_spacing(wrf_path):
    """Return the time spacing (seconds) between consecutive WRF snapshots.

    Tries ``XTIME`` (float, minutes from start) first, then ``Times``
    (character array in ``YYYY-MM-DD_HH:MM:SS`` format).  Falls back to
    3600 s (1 hour) if neither variable is present or only one snapshot
    exists.
    """
    import numpy as np

    ds = _open_nc(wrf_path)
    try:
        if "XTIME" in ds.variables:
            xtime = np.array(ds.variables["XTIME"][:], dtype=np.float64)
            if xtime.size > 1:
                return float((xtime[1] - xtime[0]) * 60.0)

        if "Times" in ds.variables:
            times_var = ds.variables["Times"]
            if times_var.shape[0] > 1:
                try:
                    from datetime import datetime

                    def _parse_wrf_time(chars):
                        s = "".join(c.decode("utf-8") if isinstance(c, bytes)
                                    else c for c in chars).strip().replace("_", "T")
                        for fmt in ("%Y-%m-%dT%H:%M:%S", "%Y-%m-%d %H:%M:%S"):
                            try:
                                return datetime.strptime(s, fmt)
                            except ValueError:
                                continue
                        return None

                    t0 = _parse_wrf_time(times_var[0])
                    t1 = _parse_wrf_time(times_var[1])
                    if t0 is not None and t1 is not None:
                        return float((t1 - t0).total_seconds())
                except Exception:
                    pass
    finally:
        ds.close()

    print("WARNING: Could not determine WRF time spacing; defaulting to 3600 s",
          file=sys.stderr)
    return 3600.0


def _wrf_bbox_to_utm_domain_bounds(lat_min, lat_max, lon_min, lon_max):
    """Convert a WRF lat/lon bounding box to UTM domain bounds."""
    import numpy as np

    corner_lats = np.array([lat_min, lat_min, lat_max, lat_max])
    corner_lons = np.array([lon_min, lon_max, lon_min, lon_max])
    x_utm, y_utm = _latlon_to_utm(corner_lats, corner_lons)
    return (float(np.min(x_utm)), float(np.min(y_utm)),
            float(np.max(x_utm)), float(np.max(y_utm)))


def _read_domain_bounds(file_path):
    """Read x/y min/max from a terrain XYZ or landscape LCP file.

    Returns
    -------
    tuple of float or None
        (x_min, y_min, x_max, y_max), or ``None`` if the file cannot be read.
    """
    xs, ys = [], []
    try:
        with open(file_path) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        xs.append(float(parts[0]))
                        ys.append(float(parts[1]))
                    except ValueError:
                        pass
    except OSError:
        return None
    if not xs:
        return None
    return min(xs), min(ys), max(xs), max(ys)


def interpolate_wind_to_grid(wrf_x, wrf_y, u, v, target_x, target_y):
    """Interpolate WRF wind components (u, v) onto a target (x, y) grid.

    Uses scipy's linear ``griddata`` interpolation.  Points outside the WRF
    domain are filled with 0.0.
    """
    import numpy as np
    try:
        from scipy.interpolate import griddata
    except ImportError:
        raise ImportError(
            "scipy is required for wind interpolation. "
            "Install with: pip install scipy"
        )

    wrf_pts = np.column_stack([
        np.asarray(wrf_x).ravel(),
        np.asarray(wrf_y).ravel(),
    ])
    target_pts = np.column_stack([
        np.asarray(target_x).ravel(),
        np.asarray(target_y).ravel(),
    ])
    target_shape = np.asarray(target_x).shape

    u_interp = griddata(wrf_pts, np.asarray(u).ravel(),
                        target_pts, method="linear", fill_value=0.0)
    v_interp = griddata(wrf_pts, np.asarray(v).ravel(),
                        target_pts, method="linear", fill_value=0.0)

    return u_interp.reshape(target_shape), v_interp.reshape(target_shape)


def extract_wrf_wind(wrf_path, time_index=0, level=0, subsample=1,
                     lat_min=None, lat_max=None,
                     lon_min=None, lon_max=None):
    """Extract wind and coordinate data from a WRF netCDF file.

    Parameters
    ----------
    wrf_path : str
        Path to the WRF output netCDF file.
    time_index : int
        Time snapshot index to extract (default 0).
    level : int
        Vertical model level for wind (default 0 = lowest level).
    subsample : int
        Keep every *subsample*-th point in each dimension (default 1).
    lat_min, lat_max, lon_min, lon_max : float or None
        Optional lat/lon bounding box (WGS-84 degrees).  When all four are
        provided the WRF grid is clipped to the tightest rectangular sub-grid
        whose grid points all fall within the specified bounds before UTM
        projection and subsampling.

    Returns
    -------
    x_utm, y_utm, u_mass, v_mass : numpy.ndarray (2-D)
        UTM coordinates and destaggered wind components.
    """
    import numpy as np

    ds = _open_nc(wrf_path)

    lat_var = _get_var(ds, "XLAT", "lat", "latitude", "LAT")
    lon_var = _get_var(ds, "XLONG", "XLON", "lon", "longitude", "LON")

    lat_arr = np.array(lat_var[:], dtype=np.float64)
    lon_arr = np.array(lon_var[:], dtype=np.float64)

    if lat_arr.ndim == 3:
        lat_2d = lat_arr[0]
        lon_2d = lon_arr[0]
    else:
        lat_2d = lat_arr
        lon_2d = lon_arr

    u_var = _get_var(ds, "U", "u", "U10", "u10")
    u_arr = np.array(u_var[:], dtype=np.float64)
    if u_arr.ndim == 4:
        U_2d_stag = u_arr[time_index, level]
    elif u_arr.ndim == 3:
        U_2d_stag = u_arr[time_index]
    else:
        U_2d_stag = u_arr

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

    if U_2d_stag.shape == (ny, nx + 1):
        u_mass = _destagger_u(U_2d_stag)
    elif U_2d_stag.shape == (ny, nx):
        u_mass = U_2d_stag
    else:
        raise ValueError(
            f"Unexpected U shape {U_2d_stag.shape} for domain ({ny}, {nx})"
        )

    if V_2d_stag.shape == (ny + 1, nx):
        v_mass = _destagger_v(V_2d_stag)
    elif V_2d_stag.shape == (ny, nx):
        v_mass = V_2d_stag
    else:
        raise ValueError(
            f"Unexpected V shape {V_2d_stag.shape} for domain ({ny}, {nx})"
        )

    clip_bounds = (lat_min is not None and lat_max is not None
                   and lon_min is not None and lon_max is not None)
    if clip_bounds:
        mask = ((lat_2d >= lat_min) & (lat_2d <= lat_max) &
                (lon_2d >= lon_min) & (lon_2d <= lon_max))
        rows, cols = np.where(mask)
        if rows.size == 0:
            raise ValueError(
                f"No WRF grid points fall within the lat/lon bounds "
                f"lat=[{lat_min}, {lat_max}] lon=[{lon_min}, {lon_max}]."
            )
        r0, r1 = int(rows.min()), int(rows.max()) + 1
        c0, c1 = int(cols.min()), int(cols.max()) + 1
        lat_2d = lat_2d[r0:r1, c0:c1]
        lon_2d = lon_2d[r0:r1, c0:c1]
        u_mass = u_mass[r0:r1, c0:c1]
        v_mass = v_mass[r0:r1, c0:c1]
        ny, nx = lat_2d.shape
        print(f"Clipped WRF wind grid to lat=[{lat_min:.4f}, {lat_max:.4f}] "
              f"lon=[{lon_min:.4f}, {lon_max:.4f}]: {ny}×{nx} points.")

    print(f"Projecting {ny}×{nx} WRF grid to UTM …")
    x_utm, y_utm = _latlon_to_utm(lat_2d, lon_2d)

    if subsample > 1:
        sl = slice(None, None, subsample)
        x_utm  = x_utm [sl, sl]
        y_utm  = y_utm [sl, sl]
        u_mass = u_mass[sl, sl]
        v_mass = v_mass[sl, sl]
        ny_out, nx_out = x_utm.shape
        print(f"Subsampled to {ny_out}×{nx_out} = {ny_out * nx_out} WRF points.")

    return x_utm, y_utm, u_mass, v_mass


def extract_wrf_terrain(wrf_path, output_path, subsample=1,
                        lat_min=None, lat_max=None,
                        lon_min=None, lon_max=None):
    """Extract terrain height from a WRF netCDF file and write an XYZ file.

    Returns
    -------
    x_utm, y_utm : numpy.ndarray (2-D)
        UTM coordinates of the terrain grid points.
    """
    import numpy as np

    ds = _open_nc(wrf_path)
    try:
        lat_var = _get_var(ds, "XLAT", "lat", "latitude", "LAT")
        lon_var = _get_var(ds, "XLONG", "XLON", "lon", "longitude", "LON")
        hgt_var = _get_var(ds, "HGT_M", "HGT", "hgt", "TERRAIN",
                           "terrain", "elevation")
        lat_arr = np.array(lat_var[:], dtype=np.float64)
        lon_arr = np.array(lon_var[:], dtype=np.float64)
        hgt_arr = np.array(hgt_var[:], dtype=np.float64)
    finally:
        ds.close()

    lat_2d = lat_arr[0] if lat_arr.ndim == 3 else lat_arr
    lon_2d = lon_arr[0] if lon_arr.ndim == 3 else lon_arr
    hgt_2d = hgt_arr[0] if hgt_arr.ndim == 3 else hgt_arr

    clip_bounds = (lat_min is not None and lat_max is not None
                   and lon_min is not None and lon_max is not None)
    if clip_bounds:
        mask = ((lat_2d >= lat_min) & (lat_2d <= lat_max) &
                (lon_2d >= lon_min) & (lon_2d <= lon_max))
        rows, cols = np.where(mask)
        if rows.size == 0:
            raise ValueError(
                f"No WRF grid points fall within the lat/lon bounds "
                f"lat=[{lat_min}, {lat_max}] lon=[{lon_min}, {lon_max}]."
            )
        r0, r1 = int(rows.min()), int(rows.max()) + 1
        c0, c1 = int(cols.min()), int(cols.max()) + 1
        lat_2d = lat_2d[r0:r1, c0:c1]
        lon_2d = lon_2d[r0:r1, c0:c1]
        hgt_2d = hgt_2d[r0:r1, c0:c1]
        ny, nx = lat_2d.shape
        print(f"Clipped WRF terrain grid to lat=[{lat_min:.4f}, {lat_max:.4f}] "
              f"lon=[{lon_min:.4f}, {lon_max:.4f}]: {ny}×{nx} points.")

    ny, nx = lat_2d.shape
    print(f"Projecting {ny}×{nx} WRF terrain grid to UTM …")
    x_utm, y_utm = _latlon_to_utm(lat_2d, lon_2d)

    if subsample > 1:
        sl = slice(None, None, subsample)
        x_utm = x_utm[sl, sl]
        y_utm = y_utm[sl, sl]
        hgt_2d = hgt_2d[sl, sl]
        ny_out, nx_out = x_utm.shape
        print(f"Subsampled to {ny_out}×{nx_out} = {ny_out * nx_out} points.")

    with open(output_path, "w") as fh:
        fh.write("# utm_x utm_y z (meters)\n")
        for xv, yv, zv in zip(x_utm.ravel(), y_utm.ravel(), hgt_2d.ravel()):
            fh.write(f"{xv:.2f} {yv:.2f} {float(zv):.6f}\n")
    print(f"Wrote WRF terrain file: '{output_path}' ({x_utm.size} points)")

    return x_utm, y_utm


def write_inputs_file(output_path,
                      terrain_file=None,
                      landscape_file=None,
                      wind_base_file=None,
                      multi_time=False,
                      wind_time_spacing=None,
                      final_time=None,
                      domain_bounds=None):
    """Write an ``inputs.i`` file configured for a pure FARSITE run.

    The generated file enables FARSITE elliptical fire spread and disables
    firebrand spotting, crown-fire initiation, and level-set advection.
    """
    bounds = domain_bounds

    with open(output_path, "w") as fh:
        fh.write(
            "# inputs.i – FARSITE run (no spotting, no crown fire)\n"
            "# Auto-generated by wrf_wind_reader.py\n"
            "#\n"
            "# Edit the ignition box extents and fuel model before running.\n\n"
        )

        fh.write("# Grid & domain\n")
        if bounds is not None:
            x_lo, y_lo, x_hi, y_hi = bounds
            n_cell_x = max(1, round((x_hi - x_lo) / 30.0))
            n_cell_y = max(1, round((y_hi - y_lo) / 30.0))
            fh.write(f"n_cell_x = {n_cell_x}\n")
            fh.write(f"n_cell_y = {n_cell_y}\n")
            fh.write(f"max_grid_size = {max(1, max(n_cell_x, n_cell_y) // 2)}\n")
            fh.write(f"prob_lo_x = {x_lo:.2f}\n")
            fh.write(f"prob_lo_y = {y_lo:.2f}\n")
            fh.write(f"prob_hi_x = {x_hi:.2f}\n")
            fh.write(f"prob_hi_y = {y_hi:.2f}\n")
        else:
            fh.write("# n_cell_x = <round((x_max - x_min) / 30)>  # 30 m resolution\n")
            fh.write("# n_cell_y = <round((y_max - y_min) / 30)>  # 30 m resolution\n")
            fh.write("# max_grid_size = <n_cell_x or n_cell_y / 2>\n")
            fh.write("# prob_lo_x = <x_min from terrain/landscape>\n")
            fh.write("# prob_lo_y = <y_min from terrain/landscape>\n")
            fh.write("# prob_hi_x = <x_max from terrain/landscape>\n")
            fh.write("# prob_hi_y = <y_max from terrain/landscape>\n")
        fh.write("\n")

        fh.write("# Time & output\n")
        if final_time is not None and final_time > 0.0:
            fh.write(f"final_time = {final_time:.1f}\n")
        else:
            fh.write("# final_time = <simulation_duration_seconds>\n")
        fh.write("cfl = 0.5\n")
        fh.write("plot_int = 10\n")
        fh.write("reinit_int = -1\n\n")

        fh.write("# Wind\n")
        if wind_base_file is not None:
            fh.write(f"velocity_file = {os.path.basename(wind_base_file)}\n")
        else:
            fh.write("# velocity_file = wind.csv\n")
        if multi_time:
            fh.write("use_time_dependent_wind = 1\n")
            if wind_time_spacing is not None:
                fh.write(f"wind_time_spacing = {wind_time_spacing:.1f}\n")
        else:
            fh.write("use_time_dependent_wind = 0\n")
        fh.write("\n")

        fh.write("# Initial ignition source (adjust extents to your fire location)\n")
        fh.write("source_type = box\n")
        if bounds is not None:
            x_lo, y_lo, x_hi, y_hi = bounds
            w = x_hi - x_lo
            h = y_hi - y_lo
            fh.write(f"box_xmin = {x_lo + 0.45 * w:.2f}\n")
            fh.write(f"box_xmax = {x_lo + 0.55 * w:.2f}\n")
            fh.write(f"box_ymin = {y_lo + 0.45 * h:.2f}\n")
            fh.write(f"box_ymax = {y_lo + 0.55 * h:.2f}\n")
        else:
            fh.write("# box_xmin / box_xmax / box_ymin / box_ymax = <UTM metres>\n")
        fh.write("\n")

        fh.write("# Rothermel fuel model (edit as needed)\n")
        if landscape_file is None:
            fh.write(
                "# No landscape file: defaulting to Southern California chaparral "
                "(NFFL FM4)\n"
            )
        fh.write("rothermel.fuel_model = FM4\n")
        fh.write("rothermel.M_f = 0.08\n")
        if terrain_file is not None:
            fh.write(f"rothermel.terrain_file = {os.path.basename(terrain_file)}\n")
        if landscape_file is not None:
            fh.write(
                f"rothermel.landscape_file = {os.path.basename(landscape_file)}\n"
            )
        fh.write("\n")

        fh.write("# FARSITE ellipse model (Richards 1990)\n")
        fh.write("skip_levelset = 1\n")
        fh.write("farsite.enable = 1\n")
        fh.write("farsite.use_anderson_LW = 1\n")
        fh.write("farsite.phi_threshold = 0.1\n\n")

        fh.write("# Disabled sub-models\n")
        fh.write("spotting.enable = 0\n")
        fh.write("crown.enable = 0\n")
        fh.write("albini_spotting.enable = 0\n")

    print(f"Wrote inputs file: '{output_path}'.")


# ===========================================================================
# Target grid construction for interpolation
# ===========================================================================

def _read_terrain_xyz_grid(terrain_path):
    """Read UTM x/y coordinates from a terrain XYZ file as 1-D arrays."""
    xs, ys = [], []
    with open(terrain_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    xs.append(float(parts[0]))
                    ys.append(float(parts[1]))
                except ValueError:
                    pass
    import numpy as np
    return np.array(xs, dtype=np.float64), np.array(ys, dtype=np.float64)


def _build_regular_grid(x_lo, y_lo, x_hi, y_hi, resolution):
    """Build a regular UTM grid with the given cell size in metres."""
    import numpy as np
    xs = np.arange(x_lo, x_hi + resolution, resolution)
    ys = np.arange(y_lo, y_hi + resolution, resolution)
    xg, yg = np.meshgrid(xs, ys)
    return xg, yg


# ===========================================================================
# CLI
# ===========================================================================

def _build_parser():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "--wrf-file", required=True, metavar="FILE",
        help="WRF netCDF file (required).",
    )
    parser.add_argument("--wind", default="wind.csv",
                        help="Output wind CSV file (default: wind.csv)")

    parser.add_argument("--lat-min", type=float, default=None,
                        help="Minimum latitude  (WGS-84 decimal degrees); "
                             "optional override of WRF-derived bbox")
    parser.add_argument("--lat-max", type=float, default=None,
                        help="Maximum latitude  (WGS-84 decimal degrees)")
    parser.add_argument("--lon-min", type=float, default=None,
                        help="Minimum longitude (WGS-84 decimal degrees)")
    parser.add_argument("--lon-max", type=float, default=None,
                        help="Maximum longitude (WGS-84 decimal degrees)")

    parser.add_argument(
        "--time-range", default=None, metavar="T1:TN",
        help=(
            "Inclusive range of WRF time indices to extract, e.g. '0:4' "
            "extracts indices 0,1,2,3,4.  Overrides --time-index when given."
        ),
    )
    parser.add_argument(
        "--time-index", type=int, default=0, metavar="N",
        help="Single WRF time index to extract (default: 0)",
    )
    parser.add_argument(
        "--level", type=int, default=0, metavar="N",
        help="WRF vertical level for U/V (default: 0 = lowest)",
    )
    parser.add_argument("--subsample", type=int, default=1, metavar="N",
                        help="Keep every N-th grid point (default: 1)")

    # Interpolation target
    interp = parser.add_argument_group(
        "interpolation target (optional; requires scipy)"
    )
    interp.add_argument(
        "--terrain-file", default=None, metavar="PATH",
        help="Existing terrain XYZ file; wind is interpolated onto its (x,y) grid.",
    )
    interp.add_argument(
        "--grid-resolution", type=float, default=None, metavar="M",
        help="Target grid resolution in metres for regular-grid interpolation.",
    )

    # inputs.i generation
    parser.add_argument(
        "--inputs", default="inputs.i", metavar="FILE",
        help="Output inputs.i file for a FARSITE run (default: inputs.i)",
    )
    parser.add_argument(
        "--no-inputs", action="store_true",
        help="Skip automatic generation of the inputs.i file",
    )

    return parser


def main(argv=None):
    parser = _build_parser()
    args = parser.parse_args(argv)

    if not os.path.isfile(args.wrf_file):
        print(f"ERROR: WRF file not found: {args.wrf_file}", file=sys.stderr)
        sys.exit(1)

    # Resolve bounding box
    wrf_lat_min, wrf_lat_max, wrf_lon_min, wrf_lon_max = read_wrf_bbox(args.wrf_file)

    # Allow user to override with explicit lat/lon bounds
    have_explicit_bbox = all(v is not None for v in [
        args.lat_min, args.lat_max, args.lon_min, args.lon_max
    ])
    if have_explicit_bbox:
        lat_min = args.lat_min
        lat_max = args.lat_max
        lon_min = args.lon_min
        lon_max = args.lon_max
        print(f"Using user-specified bounding box: "
              f"lat=[{lat_min:.4f}, {lat_max:.4f}] "
              f"lon=[{lon_min:.4f}, {lon_max:.4f}]")
    else:
        lat_min, lat_max, lon_min, lon_max = (
            wrf_lat_min, wrf_lat_max, wrf_lon_min, wrf_lon_max
        )
        print(f"Bounding box from WRF file (centre ±0.45°): "
              f"lat=[{lat_min:.4f}, {lat_max:.4f}] "
              f"lon=[{lon_min:.4f}, {lon_max:.4f}]")

    if args.subsample < 1:
        print("ERROR: --subsample must be >= 1", file=sys.stderr)
        sys.exit(1)

    if args.time_range is not None:
        try:
            time_indices = _parse_time_range(args.time_range)
        except ValueError as exc:
            print(f"ERROR: {exc}", file=sys.stderr)
            sys.exit(1)
    else:
        time_indices = [args.time_index]

    multi_time = len(time_indices) > 1
    wind_base_out = os.path.abspath(args.wind)

    # Determine interpolation target grid (if any)
    target_x = target_y = None

    if args.terrain_file is not None:
        print(f"Loading terrain grid from '{args.terrain_file}' for interpolation …")
        tx, ty = _read_terrain_xyz_grid(args.terrain_file)
        target_x = tx.reshape(-1, 1)
        target_y = ty.reshape(-1, 1)
    elif args.grid_resolution is not None:
        # Build a regular grid spanning the bbox at the requested resolution
        x_lo, y_lo, x_hi, y_hi = _wrf_bbox_to_utm_domain_bounds(
            lat_min, lat_max, lon_min, lon_max
        )
        print(
            f"Building regular {args.grid_resolution:.0f} m target grid "
            f"over UTM extents x=[{x_lo:.0f}, {x_hi:.0f}] "
            f"y=[{y_lo:.0f}, {y_hi:.0f}] …"
        )
        target_x, target_y = _build_regular_grid(
            x_lo, y_lo, x_hi, y_hi, args.grid_resolution
        )

    # Extract wind for each requested time step
    for pos, t_idx in enumerate(time_indices):
        print(f"\nExtracting WRF wind for time index {t_idx} …")
        wrf_x, wrf_y, u_mass, v_mass = extract_wrf_wind(
            args.wrf_file,
            time_index=t_idx,
            level=args.level,
            subsample=args.subsample,
            lat_min=lat_min,
            lat_max=lat_max,
            lon_min=lon_min,
            lon_max=lon_max,
        )

        wind_out = (_wind_output_path(wind_base_out, pos)
                    if multi_time else wind_base_out)

        if target_x is not None:
            print(
                f"Interpolating wind from WRF grid "
                f"({wrf_x.size} pts) to target grid "
                f"({target_x.size} pts) …"
            )
            u_interp, v_interp = interpolate_wind_to_grid(
                wrf_x, wrf_y, u_mass, v_mass, target_x, target_y
            )
            _write_wind(wind_out, target_x, target_y, u_interp, v_interp)
            print(f"Wrote wind file: '{wind_out}' "
                  f"({target_x.size} points on target grid, "
                  f"utm_x utm_y u v, time index {t_idx})")
        else:
            _write_wind(wind_out, wrf_x, wrf_y, u_mass, v_mass)
            print(f"Wrote wind file: '{wind_out}' ({wrf_x.size} points, "
                  f"utm_x utm_y u v, time index {t_idx})")

    # Auto-generate inputs.i
    if not args.no_inputs:
        domain_bounds = None
        if args.terrain_file is not None and os.path.isfile(args.terrain_file):
            domain_bounds = _read_domain_bounds(args.terrain_file)

        if domain_bounds is None:
            domain_bounds = _wrf_bbox_to_utm_domain_bounds(
                lat_min, lat_max, lon_min, lon_max
            )

        wind_time_spacing = read_wrf_time_spacing(args.wrf_file)

        final_time = None
        if args.time_range is not None and wind_time_spacing is not None:
            t_indices = _parse_time_range(args.time_range)
            if len(t_indices) > 1:
                t1 = t_indices[0]
                tn = t_indices[-1]
                final_time = (tn - t1) * wind_time_spacing

        write_inputs_file(
            output_path=os.path.abspath(args.inputs),
            wind_base_file=wind_base_out,
            multi_time=multi_time,
            wind_time_spacing=wind_time_spacing,
            final_time=final_time,
            domain_bounds=domain_bounds,
        )


if __name__ == "__main__":
    main()
