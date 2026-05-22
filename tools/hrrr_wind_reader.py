#!/usr/bin/env python3
"""
hrrr_wind_reader.py - Extract wind data from HRRR data and write CSV wind files.

Reads U and V wind components from High Resolution Rapid Refresh (HRRR) data
and writes a CSV wind file (one per requested time step) suitable for use
as ``velocity_file`` in the wildfire_levelset solver.

Key behaviours
--------------
* Fetches 10-meter U and V wind components from NOAA's HRRR archive.
* A regular grid is created within the specified lat/lon bounding box.
* Wind values are bilinearly interpolated onto the regular grid.
* Coordinates are projected to UTM for output.
* Multiple time indices can be requested; output filenames are adjusted accordingly.
* Optionally generates an ``inputs.i`` stub for a FARSITE run.

Usage examples
--------------
  # Extract wind at a specific date/time with lat/lon bounds
  python3 hrrr_wind_reader.py \\
      --lat-min 40 --lat-max 40.5 \\
      --lon-min -105 --lon-max -104.5 \\
      --date-time "2023-08-15 14:00" \\
      --wind wind.csv

  # Use grid resolution (default ~100m)
  python3 hrrr_wind_reader.py \\
      --lat-min 40 --lat-max 40.5 \\
      --lon-min -105 --lon-max -104.5 \\
      --date-time "2023-08-15 14:00" \\
      --wind wind.csv \\
      --grid-resolution 30

  # Multiple analysis times (every 6 hours)
  python3 hrrr_wind_reader.py \\
      --lat-min 40 --lat-max 40.5 \\
      --lon-min -105 --lon-max -104.5 \\
      --date-time "2023-08-15 00:00" \\
      --wind wind.csv \\
      --time-indices 0 6 12 18

Options
-------
  --lat-min / --lat-max Latitude  bounds (WGS-84 degrees; required).
  --lon-min / --lon-max Longitude bounds (WGS-84 degrees; required).
  --date-time DT        ISO-format datetime (YYYY-MM-DD HH:MM, required).
                        Uses HRRR analysis data (f00 forecast).
  --wind FILE           Output wind CSV file (default: wind.csv).
  --grid-resolution M   Target grid resolution in metres (default: 100).
  --time-indices T1 T2  Forecast lead times in hours (default: 0).
                        E.g. --time-indices 0 1 2 3 extracts f00, f01, f02, f03.
  --no-inputs           Skip automatic generation of inputs.i.
  --inputs FILE         Output inputs.i filename (default: inputs.i).
  --help                Show this message and exit.

Requires: pip install herbie-data numpy pyproj
Optional: pip install scipy (for smoother interpolation)
"""

import argparse
import os
import sys
from datetime import datetime, timedelta

import numpy as np


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
# HRRR helpers
# ===========================================================================

def _fetch_hrrr_data(date_time_str, time_index, lat_min, lat_max, lon_min, lon_max):
    """Fetch HRRR wind data for a specific time.

    Parameters
    ----------
    date_time_str : str
        ISO-format datetime string (YYYY-MM-DD HH:MM).
    time_index : int
        Forecast lead time in hours (0 for analysis, 1 for f01, etc.).
    lat_min, lat_max, lon_min, lon_max : float
        Bounding box in WGS-84 degrees.

    Returns
    -------
    u : numpy.ndarray
        U wind component (m/s) on 2D grid.
    v : numpy.ndarray
        V wind component (m/s) on 2D grid.
    lat : numpy.ndarray
        Latitude coordinates of grid points.
    lon : numpy.ndarray
        Longitude coordinates of grid points.
    """
    try:
        from herbie import Herbie
    except ImportError:
        raise ImportError(
            "herbie-data is required to fetch HRRR data. "
            "Install with: pip install herbie-data"
        )

    # Parse the date-time string
    try:
        dt = datetime.strptime(date_time_str, "%Y-%m-%d %H:%M")
    except ValueError:
        raise ValueError(
            f"Invalid date-time format: {date_time_str!r}. "
            f"Expected format: YYYY-MM-DD HH:MM"
        )

    try:
        # Fetch HRRR data using Herbie
        hrrr = Herbie(
            dt,
            model="hrrr",
            product="sfc",  # Surface level data (includes 10m wind)
            fxx=time_index,  # Forecast hour
        )

        # Download and parse wind components
        # The variable search strings find the 10m U and V components.
        # Use cfgrib engine for reading GRIB2 format HRRR data files.
        backend_kwargs = {'engine': 'cfgrib'}
        u_data = hrrr.xarray(":10 m U wind component:", backend_kwargs=backend_kwargs)
        v_data = hrrr.xarray(":10 m V wind component:", backend_kwargs=backend_kwargs)

        # Extract the data variables (first available variable in each)
        u_var_name = list(u_data.data_vars)[0] if u_data.data_vars else None
        v_var_name = list(v_data.data_vars)[0] if v_data.data_vars else None

        if u_var_name is None or v_var_name is None:
            raise ValueError(
                f"Could not find U or V wind variables in HRRR data. "
                f"Available U vars: {list(u_data.data_vars)}, "
                f"Available V vars: {list(v_data.data_vars)}"
            )

        # Extract 2D fields
        u_field = u_data[u_var_name].values
        v_field = v_data[v_var_name].values

        # Handle dimensions - remove leading time/level dimensions if present
        if u_field.ndim > 2:
            # Typically (time, lat, lon) -> take first time
            u_field = u_field[0, :, :] if u_field.ndim == 3 else u_field[0, 0, :, :]
            v_field = v_field[0, :, :] if v_field.ndim == 3 else v_field[0, 0, :, :]

        # Get coordinate arrays
        lat_coord_name = 'latitude' if 'latitude' in u_data.coords else 'lat'
        lon_coord_name = 'longitude' if 'longitude' in u_data.coords else 'lon'

        lat_coords = u_data.coords[lat_coord_name].values
        lon_coords = u_data.coords[lon_coord_name].values

        # Create 2D coordinate grids if needed
        if lat_coords.ndim == 1 and lon_coords.ndim == 1:
            lons_2d, lats_2d = np.meshgrid(lon_coords, lat_coords)
        else:
            lats_2d, lons_2d = lat_coords, lon_coords

        # Clip to bounding box
        mask = ((lons_2d >= lon_min) & (lons_2d <= lon_max) &
                (lats_2d >= lat_min) & (lats_2d <= lat_max))

        u_clipped = u_field[mask]
        v_clipped = v_field[mask]
        lat_clipped = lats_2d[mask]
        lon_clipped = lons_2d[mask]

        # Ensure we have data
        if len(u_clipped) == 0:
            raise ValueError(
                f"No HRRR data found within bounds: "
                f"lat [{lat_min}, {lat_max}], lon [{lon_min}, {lon_max}]"
            )

        return u_clipped, v_clipped, lat_clipped, lon_clipped

    except Exception as e:
        raise RuntimeError(
            f"Failed to fetch HRRR data for {date_time_str} (f{time_index:02d}): {e}"
        )


def _create_regular_grid(lat_min, lat_max, lon_min, lon_max, resolution_m):
    """Create a regular grid in the specified lat/lon bounds.

    Parameters
    ----------
    lat_min, lat_max, lon_min, lon_max : float
        Bounding box in WGS-84 degrees.
    resolution_m : float
        Desired grid resolution in metres.

    Returns
    -------
    lat : numpy.ndarray
        1D array of latitude coordinates.
    lon : numpy.ndarray
        1D array of longitude coordinates.
    """
    # Approximate metres per degree at mid-latitude
    center_lat = 0.5 * (lat_min + lat_max)
    m_per_deg_lat = 111320.0  # metres per degree latitude
    m_per_deg_lon = 111320.0 * np.cos(np.radians(center_lat))  # adjusted for latitude

    # Compute grid spacing in degrees
    dlat = resolution_m / m_per_deg_lat
    dlon = resolution_m / m_per_deg_lon

    # Create 1D coordinate arrays
    lat_arr = np.arange(lat_min, lat_max + dlat, dlat)
    lon_arr = np.arange(lon_min, lon_max + dlon, dlon)

    return lat_arr, lon_arr


def _interpolate_to_grid(u_hrrr, v_hrrr, lat_hrrr, lon_hrrr,
                         lat_grid, lon_grid):
    """Interpolate HRRR wind data onto a regular grid.

    Uses bilinear interpolation to project HRRR data onto the target grid.
    Falls back to nearest-neighbor if scipy is not available.

    Parameters
    ----------
    u_hrrr, v_hrrr : numpy.ndarray
        Wind components from HRRR data.
    lat_hrrr, lon_hrrr : numpy.ndarray
        Latitude/longitude of HRRR data points.
    lat_grid, lon_grid : numpy.ndarray
        1D arrays defining the target grid (will be converted to 2D).

    Returns
    -------
    u_interp, v_interp : numpy.ndarray
        Interpolated wind components on the regular grid (2D).
    lats_2d, lons_2d : numpy.ndarray
        2D coordinate arrays of the regular grid.
    """
    # Create 2D grid from 1D coordinates
    lons_2d, lats_2d = np.meshgrid(lon_grid, lat_grid)
    points_grid = np.column_stack([lons_2d.ravel(), lats_2d.ravel()])

    # Stack HRRR coordinates
    points_hrrr = np.column_stack([lon_hrrr.ravel(), lat_hrrr.ravel()])

    try:
        from scipy.interpolate import griddata
        # Use scipy's griddata for smooth interpolation
        u_interp = griddata(
            points_hrrr, u_hrrr.ravel(), points_grid, method='linear'
        )
        v_interp = griddata(
            points_hrrr, v_hrrr.ravel(), points_grid, method='linear'
        )
    except ImportError:
        # Fallback to nearest neighbor using pure numpy
        print(
            "Warning: scipy not installed; using nearest-neighbor interpolation. "
            "Install scipy for smoother results: pip install scipy"
        )
        # Compute distances between all grid points and HRRR points
        # For each grid point, find the nearest HRRR point
        u_interp = np.zeros(len(points_grid))
        v_interp = np.zeros(len(points_grid))
        for i, grid_pt in enumerate(points_grid):
            dists = np.sqrt(
                (points_hrrr[:, 0] - grid_pt[0])**2 +
                (points_hrrr[:, 1] - grid_pt[1])**2
            )
            nearest_idx = np.argmin(dists)
            u_interp[i] = u_hrrr.ravel()[nearest_idx]
            v_interp[i] = v_hrrr.ravel()[nearest_idx]

    u_interp = u_interp.reshape(lats_2d.shape)
    v_interp = v_interp.reshape(lats_2d.shape)

    return u_interp, v_interp, lats_2d, lons_2d


def _write_wind_csv(path, x_utm, y_utm, u, v):
    """Write wind CSV file in solver format.

    Format: # utm_x utm_y u v (m/s)
    One row per grid point (flattened from 2D grid).
    """
    with open(path, "w") as fh:
        fh.write("# utm_x utm_y u v (m/s)\n")
        for x, y, uv, vv in zip(
            x_utm.ravel(), y_utm.ravel(), u.ravel(), v.ravel()
        ):
            fh.write(f"{x:.2f} {y:.2f} {float(uv):.6f} {float(vv):.6f}\n")


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


def _generate_inputs_stub(output_path, wind_path, time_s=3600):
    """Generate a minimal inputs.i stub file for the solver.

    Parameters
    ----------
    output_path : str
        Path to write the inputs.i file.
    wind_path : str
        Path to the wind CSV file (relative or absolute).
    time_s : float
        Simulation time in seconds (default: 3600 s = 1 hour).
    """
    stub = f"""# Minimal inputs.i generated by hrrr_wind_reader.py

# Domain and time
domain.mesh_type    = "xyz"
domain.end_time     = {time_s}

# Wind input
rothermel.velocity_file = {os.path.basename(wind_path)}

# Other parameters (stub - fill in as needed)
rothermel.terrain_file  = terrain.xyz
rothermel.moisture_file = moisture.csv
"""
    with open(output_path, "w") as fh:
        fh.write(stub)


# ===========================================================================
# Main entry point
# ===========================================================================

def main():
    """Parse arguments and extract HRRR wind data."""
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "--lat-min", type=float, required=True,
        help="Minimum latitude (WGS-84 degrees)"
    )
    parser.add_argument(
        "--lat-max", type=float, required=True,
        help="Maximum latitude (WGS-84 degrees)"
    )
    parser.add_argument(
        "--lon-min", type=float, required=True,
        help="Minimum longitude (WGS-84 degrees)"
    )
    parser.add_argument(
        "--lon-max", type=float, required=True,
        help="Maximum longitude (WGS-84 degrees)"
    )
    parser.add_argument(
        "--date-time", type=str, required=True,
        help="ISO-format datetime (YYYY-MM-DD HH:MM)"
    )
    parser.add_argument(
        "--wind", type=str, default="wind.csv",
        help="Output wind CSV file (default: wind.csv)"
    )
    parser.add_argument(
        "--grid-resolution", type=float, default=100.0,
        help="Target grid resolution in metres (default: 100)"
    )
    parser.add_argument(
        "--time-indices", type=int, nargs="+", default=[0],
        help="Forecast lead times in hours (default: 0)"
    )
    parser.add_argument(
        "--no-inputs", action="store_true",
        help="Skip generation of inputs.i stub"
    )
    parser.add_argument(
        "--inputs", type=str, default="inputs.i",
        help="Output inputs.i filename (default: inputs.i)"
    )

    args = parser.parse_args()

    # Validate bounds
    if args.lat_min >= args.lat_max:
        print("Error: --lat-min must be < --lat-max", file=sys.stderr)
        sys.exit(1)
    if args.lon_min >= args.lon_max:
        print("Error: --lon-min must be < --lon-max", file=sys.stderr)
        sys.exit(1)

    print(
        f"Extracting HRRR wind data for bounds:\n"
        f"  lat: [{args.lat_min}, {args.lat_max}]°\n"
        f"  lon: [{args.lon_min}, {args.lon_max}]°\n"
        f"Grid resolution: {args.grid_resolution} m"
    )

    # Create regular grid
    lat_grid, lon_grid = _create_regular_grid(
        args.lat_min, args.lat_max, args.lon_min, args.lon_max,
        args.grid_resolution
    )
    print(f"Created {len(lon_grid)} × {len(lat_grid)} grid points")

    # Process each time index
    for pos, time_idx in enumerate(args.time_indices):
        try:
            print(f"\nProcessing f{time_idx:02d} (position {pos})...")

            # Fetch HRRR data
            u_hrrr, v_hrrr, lat_hrrr, lon_hrrr = _fetch_hrrr_data(
                args.date_time, time_idx,
                args.lat_min, args.lat_max, args.lon_min, args.lon_max
            )
            print(f"  Fetched {len(u_hrrr)} HRRR points")

            # Interpolate to regular grid
            u_grid, v_grid, lats_2d, lons_2d = _interpolate_to_grid(
                u_hrrr, v_hrrr, lat_hrrr, lon_hrrr,
                lat_grid, lon_grid
            )

            # Project to UTM
            x_utm, y_utm = _latlon_to_utm(lats_2d, lons_2d)

            # Write wind CSV
            wind_output = _wind_output_path(args.wind, pos)
            _write_wind_csv(wind_output, x_utm, y_utm, u_grid, v_grid)
            print(f"  Wrote wind data to {wind_output}")

        except Exception as e:
            print(f"Error processing f{time_idx:02d}: {e}", file=sys.stderr)
            sys.exit(1)

    # Generate inputs.i stub if requested
    if not args.no_inputs:
        _generate_inputs_stub(args.inputs, args.wind)
        print(f"\nGenerated inputs.i stub: {args.inputs}")

    print("\nDone!")


if __name__ == "__main__":
    main()
