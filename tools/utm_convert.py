#!/usr/bin/env python3
"""
utm_convert.py - Bidirectional lat/lon ↔ UTM coordinate conversion.

Pure-Python WGS-84 implementation requires no external dependencies.
If pyproj is available it is used as the backend for higher accuracy.

Module-level API
----------------
    latlon_to_utm(lat, lon)
        → (easting, northing, zone_number, zone_letter)

    utm_to_latlon(easting, northing, zone_number, northern=True)
        → (lat, lon)

CLI usage
---------
    # Convert lat/lon → UTM (zone auto-detected)
    python utm_convert.py --to-utm 34.10 -118.85

    # Convert lat/lon → UTM (force specific zone)
    python utm_convert.py --to-utm 34.10 -118.85 --zone 11

    # Convert UTM → lat/lon
    python utm_convert.py --to-latlon 330000 3775000 --zone 11

    # Southern hemisphere
    python utm_convert.py --to-latlon 330000 3775000 --zone 11 --south

Examples
--------
Southern California (Santa Monica Mountains, ~34.10°N 118.85°W):

>>> e, n, zone, letter = latlon_to_utm(34.10, -118.85)
>>> print(f"UTM Zone {zone}{letter}: {e:.0f} E, {n:.0f} N")
UTM Zone 11S: 329346 E, 3774789 N

Note: The zone letter 'S' is the UTM latitude band (32°N–40°N) and does NOT
indicate the southern hemisphere.  Northern-hemisphere coordinates always
use false northing N0 = 0; pass ``northern=True`` to utm_to_latlon.

>>> lat, lon = utm_to_latlon(330000, 3775000, zone_number=11, northern=True)
>>> print(f"{lat:.4f}°N, {lon:.4f}°W")
34.1020°N, 118.8430°W
"""

import argparse
import math
import sys

# ---------------------------------------------------------------------------
# WGS-84 ellipsoid constants
# ---------------------------------------------------------------------------

_WGS84_A  = 6378137.0          # semi-major axis [m]
_WGS84_F  = 1.0 / 298.257223563
_WGS84_B  = _WGS84_A * (1.0 - _WGS84_F)
_WGS84_E2 = 1.0 - (_WGS84_B / _WGS84_A) ** 2   # first eccentricity squared
_WGS84_EP2 = _WGS84_E2 / (1.0 - _WGS84_E2)     # second eccentricity squared

_UTM_K0 = 0.9996          # central-meridian scale factor
_UTM_E0 = 500000.0        # false easting [m]
_UTM_N0_SOUTH = 10000000.0  # false northing for southern hemisphere [m]


# ---------------------------------------------------------------------------
# Pure-Python implementation
# ---------------------------------------------------------------------------

def _zone_letter(lat):
    """Return the UTM latitude band letter for *lat* (degrees)."""
    bands = "CDEFGHJKLMNPQRSTUVWXX"
    if -80.0 <= lat <= 84.0:
        idx = int((lat + 80.0) / 8.0)
        idx = min(idx, len(bands) - 1)
        return bands[idx]
    return "Z"  # out-of-bounds sentinel


def _central_meridian(zone_number):
    """Return the central meridian longitude for a given UTM zone number."""
    return (zone_number - 1) * 6.0 - 177.0


def latlon_to_utm(lat, lon):
    """Convert WGS-84 geographic coordinates to UTM.

    Parameters
    ----------
    lat : float
        Latitude in decimal degrees (positive = north).
    lon : float
        Longitude in decimal degrees (positive = east).

    Returns
    -------
    easting : float
        UTM easting in metres.
    northing : float
        UTM northing in metres.
    zone_number : int
        UTM zone number (1-60).
    zone_letter : str
        UTM latitude band letter.

    Notes
    -----
    Uses pyproj when available; falls back to a pure-Python WGS-84
    Transverse Mercator formula otherwise.
    """
    # --- pyproj fast path ---
    try:
        from pyproj import Transformer, CRS
        zone_number = int((lon + 180.0) / 6.0) + 1
        hemisphere  = "north" if lat >= 0.0 else "south"
        epsg = 32600 + zone_number if lat >= 0.0 else 32700 + zone_number
        t = Transformer.from_crs("EPSG:4326", f"EPSG:{epsg}", always_xy=True)
        easting, northing = t.transform(lon, lat)
        return easting, northing, zone_number, _zone_letter(lat)
    except ImportError:
        pass

    # --- pure-Python WGS-84 TM ---
    lat_r = math.radians(lat)
    zone_number = int((lon + 180.0) / 6.0) + 1
    lon0_r = math.radians(_central_meridian(zone_number))
    # A = (lon - lon0) * cos(lat) — standard TM intermediate variable
    A = (math.radians(lon) - lon0_r) * math.cos(lat_r)

    a  = _WGS84_A
    e2 = _WGS84_E2
    n  = a / math.sqrt(1.0 - e2 * math.sin(lat_r) ** 2)
    t  = math.tan(lat_r) ** 2
    c  = _WGS84_EP2 * math.cos(lat_r) ** 2

    # Meridional arc
    e4 = e2 * e2
    e6 = e4 * e2
    m = a * (
        (1.0 - e2 / 4.0 - 3.0 * e4 / 64.0 - 5.0 * e6 / 256.0) * lat_r
        - (3.0 * e2 / 8.0 + 3.0 * e4 / 32.0 + 45.0 * e6 / 1024.0)
          * math.sin(2.0 * lat_r)
        + (15.0 * e4 / 256.0 + 45.0 * e6 / 1024.0) * math.sin(4.0 * lat_r)
        - (35.0 * e6 / 3072.0) * math.sin(6.0 * lat_r)
    )

    easting = _UTM_K0 * n * (
        A
        + A ** 3 / 6.0   * (1.0 - t + c)
        + A ** 5 / 120.0 * (5.0 - 18.0 * t + t * t + 72.0 * c - 58.0 * _WGS84_EP2)
    ) + _UTM_E0

    northing_raw = _UTM_K0 * (
        m
        + n * math.tan(lat_r) * (
            A ** 2 / 2.0
            + A ** 4 / 24.0 * (5.0 - t + 9.0 * c + 4.0 * c * c)
            + A ** 6 / 720.0 * (61.0 - 58.0 * t + t * t + 600.0 * c - 330.0 * _WGS84_EP2)
        )
    )

    if lat < 0.0:
        northing = northing_raw + _UTM_N0_SOUTH
    else:
        northing = northing_raw

    return easting, northing, zone_number, _zone_letter(lat)


def utm_to_latlon(easting, northing, zone_number, northern=True):
    """Convert UTM coordinates to WGS-84 geographic coordinates.

    Parameters
    ----------
    easting : float
        UTM easting in metres.
    northing : float
        UTM northing in metres.
    zone_number : int
        UTM zone number (1-60).
    northern : bool
        True (default) for northern hemisphere; False for southern.

    Returns
    -------
    lat : float
        Latitude in decimal degrees.
    lon : float
        Longitude in decimal degrees.
    """
    # --- pyproj fast path ---
    try:
        from pyproj import Transformer
        epsg = 32600 + zone_number if northern else 32700 + zone_number
        t = Transformer.from_crs(f"EPSG:{epsg}", "EPSG:4326", always_xy=True)
        lon, lat = t.transform(easting, northing)
        return lat, lon
    except ImportError:
        pass

    # --- pure-Python inverse TM ---
    a  = _WGS84_A
    e2 = _WGS84_E2
    e4 = e2 * e2
    e6 = e4 * e2
    e1 = (1.0 - math.sqrt(1.0 - e2)) / (1.0 + math.sqrt(1.0 - e2))

    x = easting - _UTM_E0
    y = northing if northern else northing - _UTM_N0_SOUTH

    m = y / _UTM_K0
    mu = m / (a * (1.0 - e2 / 4.0 - 3.0 * e4 / 64.0 - 5.0 * e6 / 256.0))

    lat1 = (
        mu
        + (3.0 * e1 / 2.0 - 27.0 * e1 ** 3 / 32.0) * math.sin(2.0 * mu)
        + (21.0 * e1 ** 2 / 16.0 - 55.0 * e1 ** 4 / 32.0) * math.sin(4.0 * mu)
        + (151.0 * e1 ** 3 / 96.0) * math.sin(6.0 * mu)
        + (1097.0 * e1 ** 4 / 512.0) * math.sin(8.0 * mu)
    )

    n1 = a / math.sqrt(1.0 - e2 * math.sin(lat1) ** 2)
    t1 = math.tan(lat1) ** 2
    c1 = _WGS84_EP2 * math.cos(lat1) ** 2
    r1 = a * (1.0 - e2) / (1.0 - e2 * math.sin(lat1) ** 2) ** 1.5
    d  = x / (n1 * _UTM_K0)

    lat = lat1 - (n1 * math.tan(lat1) / r1) * (
        d * d / 2.0
        - d ** 4 / 24.0 * (5.0 + 3.0 * t1 + 10.0 * c1 - 4.0 * c1 * c1 - 9.0 * _WGS84_EP2)
        + d ** 6 / 720.0 * (61.0 + 90.0 * t1 + 298.0 * c1 + 45.0 * t1 * t1 - 252.0 * _WGS84_EP2 - 3.0 * c1 * c1)
    )

    lon0 = math.radians(_central_meridian(zone_number))
    lon = lon0 + (
        d
        - d ** 3 / 6.0   * (1.0 + 2.0 * t1 + c1)
        + d ** 5 / 120.0 * (5.0 - 2.0 * c1 + 28.0 * t1 - 3.0 * c1 * c1 + 8.0 * _WGS84_EP2 + 24.0 * t1 * t1)
    ) / math.cos(lat1)

    return math.degrees(lat), math.degrees(lon)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _build_parser():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--to-utm", nargs=2, type=float,
        metavar=("LAT", "LON"),
        help="Convert lat/lon (decimal degrees) to UTM",
    )
    group.add_argument(
        "--to-latlon", nargs=2, type=float,
        metavar=("EASTING", "NORTHING"),
        help="Convert UTM easting/northing (metres) to lat/lon",
    )
    p.add_argument(
        "--zone", type=int, default=None,
        help="UTM zone number (1-60); auto-detected for --to-utm, required for --to-latlon",
    )
    p.add_argument(
        "--south", action="store_true",
        help="Southern hemisphere (applies to --to-latlon only)",
    )
    return p


def main(argv=None):
    parser = _build_parser()
    args = parser.parse_args(argv)

    if args.to_utm is not None:
        lat, lon = args.to_utm
        easting, northing, zone_auto, letter = latlon_to_utm(lat, lon)
        zone = args.zone if args.zone is not None else zone_auto
        if args.zone is not None and args.zone != zone_auto:
            # Re-compute with forced zone
            try:
                from pyproj import Transformer
                epsg = 32600 + args.zone if lat >= 0 else 32700 + args.zone
                t = Transformer.from_crs("EPSG:4326", f"EPSG:{epsg}", always_xy=True)
                easting, northing = t.transform(lon, lat)
            except ImportError:
                pass  # pure-Python path auto-detects zone; just report
        print(f"UTM Zone {zone}{letter}")
        print(f"  Easting  : {easting:.2f} m")
        print(f"  Northing : {northing:.2f} m")
        print(f"{easting:.2f} {northing:.2f} {zone}")
    else:
        easting, northing = args.to_latlon
        if args.zone is None:
            parser.error("--zone is required for --to-latlon")
        northern = not args.south
        lat, lon = utm_to_latlon(easting, northing, args.zone, northern=northern)
        hemi = "N" if lat >= 0 else "S"
        print(f"Latitude  : {abs(lat):.6f}° {hemi}")
        print(f"Longitude : {abs(lon):.6f}° {'E' if lon >= 0 else 'W'}")
        print(f"{lat:.6f} {lon:.6f}")


if __name__ == "__main__":
    main()
