#!/usr/bin/env python3
"""
perimeter_to_shapefile.py – Convert wildfire_levelset fire perimeter files
to Esri Shapefile (.shp) format.

Reads one or more fire perimeter files in any of the supported formats
produced by wildfire_levelset:

  • ``perimeter_NNNN.geojson``  – GeoJSON FeatureCollection with a Polygon
  • ``perimeter_NNNN.csv``      – two-column CSV with X,Y coordinates
  • ``phi_negative_NNNN.dat``   – whitespace-delimited X Y point cloud (phi < 0)
  • ``phi_envelope_NNNN.dat``   – whitespace-delimited X Y convex hull vertices

and writes each to a corresponding ``.shp`` file (with the standard sidecar
files ``.dbf``, ``.shx``, ``.prj``).

When a UTM EPSG code is provided the ``.prj`` file is written with the correct
coordinate reference system so GIS applications can georeference the output.

Requirements
------------
  pip install shapely pyshp

Optional (for richer CRS support):
  pip install pyproj

Usage examples
--------------
  # Convert a single GeoJSON perimeter to shapefile
  python3 tools/perimeter_to_shapefile.py perimeter_0100.geojson

  # Convert all GeoJSON perimeters, assign UTM zone 13N (EPSG:32613)
  python3 tools/perimeter_to_shapefile.py perimeter_*.geojson --epsg 32613

  # Convert CSV perimeters, write to a different directory
  python3 tools/perimeter_to_shapefile.py perimeter_*.csv --outdir shapefiles/

  # Convert all phi_negative point clouds to Multipoint shapefiles
  python3 tools/perimeter_to_shapefile.py phi_negative_*.dat --point-cloud

  # Batch convert everything in the current directory
  python3 tools/perimeter_to_shapefile.py --all

Options
-------
  INPUT [INPUT ...]    Input perimeter file(s). Glob patterns are supported on
                       systems where the shell does not expand them.
  --outdir DIR         Output directory (default: same as input file)
  --epsg N             EPSG code for the CRS written to the .prj file
                       (e.g. 32613 for UTM zone 13N); omit for simulation coords
  --point-cloud        Write phi_negative .dat files as Multipoint rather than
                       attempting to form a polygon ring
  --all                Convert all supported perimeter files in the current directory
  --overwrite          Overwrite existing .shp files (default: skip)

References
----------
  ESRI Shapefile Technical Description (1998).
  https://support.esri.com/en/white-paper/279
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import sys
from pathlib import Path
from typing import List, Optional, Tuple


# ---------------------------------------------------------------------------
# Optional imports
# ---------------------------------------------------------------------------

try:
    import shapefile  # pyshp
    HAS_PYSHP = True
except ImportError:
    HAS_PYSHP = False

try:
    from shapely.geometry import Polygon, MultiPoint, mapping
    from shapely.ops import unary_union
    HAS_SHAPELY = True
except ImportError:
    HAS_SHAPELY = False

try:
    import pyproj
    HAS_PYPROJ = True
except ImportError:
    HAS_PYPROJ = False


# ---------------------------------------------------------------------------
# WKT CRS definitions for common UTM zones (subset; EPSG:326xx north)
# ---------------------------------------------------------------------------

def _utm_wkt(epsg: int) -> str:
    """Return a minimal WKT string for a UTM zone given its EPSG code."""
    if HAS_PYPROJ:
        try:
            crs = pyproj.CRS.from_epsg(epsg)
            return crs.to_wkt()
        except Exception:
            pass
    # Fallback: minimal WKT for EPSG:326xx (UTM north) and 327xx (UTM south)
    if 32601 <= epsg <= 32660:
        zone = epsg - 32600
        return (
            f'PROJCS["WGS 84 / UTM zone {zone}N",'
            'GEOGCS["WGS 84",DATUM["WGS_1984",'
            'SPHEROID["WGS 84",6378137,298.257223563]],'
            'PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]],'
            f'PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],'
            f'PARAMETER["central_meridian",{-183 + zone*6}],'
            'PARAMETER["scale_factor",0.9996],'
            'PARAMETER["false_easting",500000],PARAMETER["false_northing",0],'
            'UNIT["metre",1]]'
        )
    if 32701 <= epsg <= 32760:
        zone = epsg - 32700
        return (
            f'PROJCS["WGS 84 / UTM zone {zone}S",'
            'GEOGCS["WGS 84",DATUM["WGS_1984",'
            'SPHEROID["WGS 84",6378137,298.257223563]],'
            'PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]],'
            f'PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],'
            f'PARAMETER["central_meridian",{-183 + zone*6}],'
            'PARAMETER["scale_factor",0.9996],'
            'PARAMETER["false_easting",500000],'
            'PARAMETER["false_northing",10000000],'
            'UNIT["metre",1]]'
        )
    return ""


# ---------------------------------------------------------------------------
# Readers
# ---------------------------------------------------------------------------

def _read_geojson(path: str) -> List[List[Tuple[float, float]]]:
    """Read GeoJSON and return list of rings [(x,y), ...]."""
    with open(path) as f:
        gj = json.load(f)
    rings = []
    features = gj.get("features", []) if gj.get("type") == "FeatureCollection" else [gj]
    for feat in features:
        geom = feat.get("geometry", feat)
        gtype = geom.get("type", "")
        coords = geom.get("coordinates", [])
        if gtype == "Polygon":
            rings.extend([list(map(tuple, r)) for r in coords])
        elif gtype == "MultiPolygon":
            for poly in coords:
                rings.extend([list(map(tuple, r)) for r in poly])
    return rings


def _read_csv(path: str) -> List[Tuple[float, float]]:
    """Read X,Y CSV and return list of points."""
    pts = []
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#") or s.lower().startswith("x"):
                continue
            parts = s.replace(",", " ").split()
            if len(parts) >= 2:
                try:
                    pts.append((float(parts[0]), float(parts[1])))
                except ValueError:
                    pass
    return pts


def _read_dat(path: str) -> List[Tuple[float, float]]:
    """Read whitespace-delimited X Y .dat file and return list of points."""
    pts = []
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.split()
            if len(parts) >= 2:
                try:
                    pts.append((float(parts[0]), float(parts[1])))
                except ValueError:
                    pass
    return pts


# ---------------------------------------------------------------------------
# Point cloud → polygon (convex hull fallback)
# ---------------------------------------------------------------------------

def _points_to_ring(pts: List[Tuple[float, float]]) -> List[Tuple[float, float]]:
    """
    Convert a set of cell-centre points to a polygon ring.
    Uses shapely convex_hull when available; falls back to a simple
    bounding-box outline otherwise.
    """
    if not pts:
        return []
    if HAS_SHAPELY and len(pts) >= 3:
        try:
            hull = MultiPoint(pts).convex_hull
            if hull.geom_type == "Polygon":
                return list(hull.exterior.coords)
            if hull.geom_type == "LineString":
                return list(hull.coords)
        except Exception:
            pass
    # Bounding box fallback
    xs = [p[0] for p in pts]
    ys = [p[1] for p in pts]
    x0, x1 = min(xs), max(xs)
    y0, y1 = min(ys), max(ys)
    return [(x0, y0), (x1, y0), (x1, y1), (x0, y1), (x0, y0)]


# ---------------------------------------------------------------------------
# Shapefile writing (pyshp)
# ---------------------------------------------------------------------------

def _write_shapefile(
    out_stem: str,
    rings: List[List[Tuple[float, float]]],
    step: Optional[int],
    time_s: Optional[float],
    epsg: Optional[int],
    as_multipoint: bool = False,
    pts: Optional[List[Tuple[float, float]]] = None,
) -> None:
    """Write a Polygon or Multipoint shapefile from rings/points."""
    if not HAS_PYSHP:
        print("  ERROR: pyshp not installed; run:  pip install pyshp", file=sys.stderr)
        return

    if as_multipoint and pts:
        shp = shapefile.Writer(out_stem, shapefile.MULTIPOINT)
        shp.field("step",   "N")
        shp.field("time_s", "F", decimal=2)
        shp.multipoint(pts)
        shp.record(step or 0, time_s or 0.0)
    else:
        if not rings:
            print(f"  WARNING: no rings found; skipping {out_stem}.shp")
            return
        shp = shapefile.Writer(out_stem, shapefile.POLYGON)
        shp.field("step",   "N")
        shp.field("time_s", "F", decimal=2)
        shp.poly([list(ring) for ring in rings])
        shp.record(step or 0, time_s or 0.0)

    shp.close()

    # Write .prj file
    if epsg:
        wkt = _utm_wkt(epsg)
        if wkt:
            with open(out_stem + ".prj", "w") as f:
                f.write(wkt)

    print(f"  Wrote {out_stem}.shp")


# ---------------------------------------------------------------------------
# Detect step and time from filename
# ---------------------------------------------------------------------------

def _parse_step(name: str) -> Optional[int]:
    """Extract step number from filename like perimeter_0100.csv."""
    import re
    m = re.search(r"_(\d+)\.", name)
    return int(m.group(1)) if m else None


def _read_and_validate(path: str, suffix: str) -> Optional[object]:
    """Read a perimeter file and return the parsed data, or None on empty/error."""
    if suffix == ".geojson":
        rings = _read_geojson(path)
        if not rings:
            print(f"  WARNING: no rings in {path}")
            return None
        return rings
    elif suffix == ".csv":
        pts = _read_csv(path)
        if not pts:
            print(f"  WARNING: no points in {path}")
            return None
        return pts
    elif suffix == ".dat":
        pts = _read_dat(path)
        if not pts:
            print(f"  WARNING: no points in {path}")
            return None
        return pts
    return None


# ---------------------------------------------------------------------------
# Per-file conversion
# ---------------------------------------------------------------------------

def convert_file(
    path: str,
    outdir: Optional[str],
    epsg: Optional[int],
    overwrite: bool,
    point_cloud: bool,
) -> bool:
    """Convert one input file to shapefile. Returns True on success."""
    p = Path(path)
    suffix = p.suffix.lower()
    stem   = p.stem
    out_parent = Path(outdir) if outdir else p.parent
    out_parent.mkdir(parents=True, exist_ok=True)
    out_stem = str(out_parent / stem)

    if not overwrite and os.path.isfile(out_stem + ".shp"):
        print(f"  Skipping (exists): {out_stem}.shp  (use --overwrite to replace)")
        return True

    step = _parse_step(p.name)
    data = _read_and_validate(path, suffix)
    if data is None:
        return False

    if suffix == ".geojson":
        _write_shapefile(out_stem, data, step, None, epsg)
    elif suffix == ".csv":
        pts = data
        ring = pts if (pts[0] == pts[-1]) else pts + [pts[0]]
        _write_shapefile(out_stem, [ring], step, None, epsg)
    elif suffix == ".dat":
        pts = data
        if point_cloud:
            _write_shapefile(out_stem, [], step, None, epsg,
                             as_multipoint=True, pts=pts)
        else:
            ring = _points_to_ring(pts)
            _write_shapefile(out_stem, [ring], step, None, epsg)
    else:
        print(f"  WARNING: unsupported file type {suffix} for {path}")
        return False

    return True


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Convert wildfire_levelset perimeter files to Esri Shapefile.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("inputs", nargs="*", metavar="INPUT",
                        help="Input perimeter file(s) (.geojson, .csv, .dat)")
    parser.add_argument("--outdir",      default=None,
                        help="Output directory (default: same as input file)")
    parser.add_argument("--epsg",        type=int, default=None,
                        help="EPSG code for .prj CRS (e.g. 32613 for UTM zone 13N)")
    parser.add_argument("--point-cloud", action="store_true",
                        help="Write phi_negative .dat as Multipoint (default: polygon)")
    parser.add_argument("--all",         action="store_true",
                        help="Convert all supported perimeter files in current directory")
    parser.add_argument("--overwrite",   action="store_true",
                        help="Overwrite existing .shp files")

    args = parser.parse_args(argv)

    if not HAS_PYSHP:
        print("ERROR: pyshp is required.  Install it with:  pip install pyshp",
              file=sys.stderr)
        sys.exit(1)

    files = list(args.inputs)
    if args.all:
        for pattern in [
            "perimeter_*.geojson",
            "perimeter_*.csv",
            "phi_envelope_*.dat",
            "phi_negative_*.dat",
        ]:
            files.extend(sorted(glob.glob(pattern)))

    # Expand any remaining glob patterns (for shells that don't auto-expand)
    expanded = []
    for f in files:
        matches = sorted(glob.glob(f))
        expanded.extend(matches if matches else [f])
    files = expanded

    if not files:
        parser.print_help()
        print("\nERROR: no input files specified.", file=sys.stderr)
        sys.exit(1)

    ok = 0
    for path in files:
        if not os.path.isfile(path):
            print(f"  WARNING: file not found: {path}", file=sys.stderr)
            continue
        print(f"Converting {path} ...")
        success = convert_file(path, args.outdir, args.epsg,
                               args.overwrite, args.point_cloud)
        if success:
            ok += 1

    print(f"\nConverted {ok}/{len(files)} file(s).")


if __name__ == "__main__":
    main()
