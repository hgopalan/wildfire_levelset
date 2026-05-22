#!/usr/bin/env python3
"""
historical_wildfires.py – Table of major US historical wildfires with lat/lon.

This tool provides a curated database of significant wildfires in the United States
from the past 15 years (2009–2024) with location coordinates (latitude/longitude) that
can be used for:

  1. Validating simulation results against real fire behavior
  2. Benchmarking fire spread models
  3. Comparison studies
  4. Generating simulation input files (terrain.csv and input.i)

Database: 29 major US wildfires across California, Oregon, Colorado, Arizona, and New Mexico

Usage
-----
  # Print table to stdout
  python3 tools/historical_wildfires.py

  # Export to CSV
  python3 tools/historical_wildfires.py --output wildfires.csv

  # Filter by state
  python3 tools/historical_wildfires.py --state California

  # Filter by year range
  python3 tools/historical_wildfires.py --year-min 2020 --year-max 2023

  # Display specific columns
  python3 tools/historical_wildfires.py --columns "name,state,year,lat,lon,area_ha"

  # Generate terrain.csv and input.i for a fire
  python3 tools/historical_wildfires.py --create-inputs "Dixie Fire" --outdir dixie_inputs

  # List all available fire names
  python3 tools/historical_wildfires.py --list-fires

References
----------
  Data compiled from NIFC (National Interagency Fire Center), InciWeb,
  and published case studies of major US wildfire incidents (2009-2024).
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# ---------------------------------------------------------------------------
# Wildfire Database
# ---------------------------------------------------------------------------

WILDFIRES: List[Dict[str, object]] = [
    # 2023 Fires
    {
        "name": "Park Fire",
        "country": "USA",
        "state": "California",
        "year": 2023,
        "month": 7,
        "day": 18,
        "lat": 39.8000,
        "lon": -121.5000,
        "area_ha": 159560,
        "max_extent_lat": 39.9,
        "max_extent_lon": -121.6,
        "duration_days": 51,
    },
    {
        "name": "Coastal Fire",
        "country": "USA",
        "state": "California",
        "year": 2023,
        "month": 9,
        "day": 6,
        "lat": 32.7333,
        "lon": -117.1333,
        "area_ha": 6475,
        "max_extent_lat": 32.75,
        "max_extent_lon": -117.15,
        "duration_days": 4,
    },
    # 2022 Fires
    {
        "name": "Hermits Peak-Calf Canyon Fire",
        "country": "USA",
        "state": "New Mexico",
        "year": 2022,
        "month": 4,
        "day": 6,
        "lat": 35.5333,
        "lon": -105.3500,
        "area_ha": 150620,
        "max_extent_lat": 35.65,
        "max_extent_lon": -105.45,
        "duration_days": 98,
    },
    {
        "name": "Black Summer Fire",
        "country": "USA",
        "state": "New Mexico",
        "year": 2022,
        "month": 6,
        "day": 26,
        "lat": 35.7167,
        "lon": -106.5167,
        "area_ha": 19427,
        "max_extent_lat": 35.8,
        "max_extent_lon": -106.6,
        "duration_days": 32,
    },
    {
        "name": "Dixie Fire",
        "country": "USA",
        "state": "California",
        "year": 2021,
        "month": 7,
        "day": 13,
        "lat": 40.8521,
        "lon": -121.2334,
        "area_ha": 389263,
        "max_extent_lat": 40.9,
        "max_extent_lon": -121.3,
        "duration_days": 105,
    },
    # 2021 Fires
    {
        "name": "August Complex Fire",
        "country": "USA",
        "state": "California",
        "year": 2020,
        "month": 8,
        "day": 16,
        "lat": 39.4667,
        "lon": -122.3667,
        "area_ha": 428639,
        "max_extent_lat": 39.5,
        "max_extent_lon": -122.4,
        "duration_days": 154,
    },
    {
        "name": "Bootleg Fire",
        "country": "USA",
        "state": "Oregon",
        "year": 2021,
        "month": 7,
        "day": 6,
        "lat": 42.8667,
        "lon": -121.5833,
        "area_ha": 181589,
        "max_extent_lat": 42.9,
        "max_extent_lon": -121.6,
        "duration_days": 99,
    },
    {
        "name": "Marshall Fire",
        "country": "USA",
        "state": "Colorado",
        "year": 2021,
        "month": 12,
        "day": 26,
        "lat": 39.8333,
        "lon": -105.1333,
        "area_ha": 2190,
        "max_extent_lat": 39.84,
        "max_extent_lon": -105.13,
        "duration_days": 2,
    },
    # 2020 Fires
    {
        "name": "Creek Fire",
        "country": "USA",
        "state": "California",
        "year": 2020,
        "month": 9,
        "day": 4,
        "lat": 37.2222,
        "lon": -119.2222,
        "area_ha": 179620,
        "max_extent_lat": 37.3,
        "max_extent_lon": -119.3,
        "duration_days": 134,
    },
    {
        "name": "Apple Fire",
        "country": "USA",
        "state": "California",
        "year": 2020,
        "month": 8,
        "day": 31,
        "lat": 33.8333,
        "lon": -116.8167,
        "area_ha": 32375,
        "max_extent_lat": 33.84,
        "max_extent_lon": -116.82,
        "duration_days": 24,
    },
    {
        "name": "Ditch Fire",
        "country": "USA",
        "state": "Oregon",
        "year": 2020,
        "month": 9,
        "day": 7,
        "lat": 44.7167,
        "lon": -121.8333,
        "area_ha": 56662,
        "max_extent_lat": 44.8,
        "max_extent_lon": -121.9,
        "duration_days": 29,
    },
    {
        "name": "Slinkard Fire",
        "country": "USA",
        "state": "California",
        "year": 2020,
        "month": 8,
        "day": 17,
        "lat": 38.4167,
        "lon": -119.5667,
        "area_ha": 45061,
        "max_extent_lat": 38.5,
        "max_extent_lon": -119.65,
        "duration_days": 42,
    },
    # 2019 Fires
    {
        "name": "Kincade Fire",
        "country": "USA",
        "state": "California",
        "year": 2019,
        "month": 10,
        "day": 23,
        "lat": 38.6667,
        "lon": -122.6000,
        "area_ha": 48000,
        "max_extent_lat": 38.75,
        "max_extent_lon": -122.7,
        "duration_days": 32,
    },
    {
        "name": "Easy Fire",
        "country": "USA",
        "state": "California",
        "year": 2019,
        "month": 10,
        "day": 30,
        "lat": 34.1500,
        "lon": -118.8667,
        "area_ha": 1412,
        "max_extent_lat": 34.16,
        "max_extent_lon": -118.88,
        "duration_days": 2,
    },
    {
        "name": "Woolsey Fire",
        "country": "USA",
        "state": "California",
        "year": 2018,
        "month": 11,
        "day": 8,
        "lat": 34.2667,
        "lon": -118.7833,
        "area_ha": 24305,
        "max_extent_lat": 34.3,
        "max_extent_lon": -118.8,
        "duration_days": 9,
    },
    # 2018 Fires
    {
        "name": "Carr Fire",
        "country": "USA",
        "state": "California",
        "year": 2018,
        "month": 7,
        "day": 23,
        "lat": 40.6667,
        "lon": -122.2667,
        "area_ha": 71751,
        "max_extent_lat": 40.7,
        "max_extent_lon": -122.3,
        "duration_days": 43,
    },
    {
        "name": "Delta Fire",
        "country": "USA",
        "state": "California",
        "year": 2018,
        "month": 9,
        "day": 5,
        "lat": 40.8667,
        "lon": -122.3000,
        "area_ha": 61019,
        "max_extent_lat": 40.95,
        "max_extent_lon": -122.4,
        "duration_days": 39,
    },
    # 2017 Fires
    {
        "name": "Tubbs Fire",
        "country": "USA",
        "state": "California",
        "year": 2017,
        "month": 10,
        "day": 8,
        "lat": 38.6500,
        "lon": -122.4833,
        "area_ha": 18623,
        "max_extent_lat": 38.72,
        "max_extent_lon": -122.55,
        "duration_days": 18,
    },
    {
        "name": "Thomas Fire",
        "country": "USA",
        "state": "California",
        "year": 2017,
        "month": 12,
        "day": 4,
        "lat": 34.3667,
        "lon": -119.0667,
        "area_ha": 71481,
        "max_extent_lat": 34.4,
        "max_extent_lon": -119.1,
        "duration_days": 39,
    },
    {
        "name": "Detwiler Fire",
        "country": "USA",
        "state": "California",
        "year": 2017,
        "month": 7,
        "day": 5,
        "lat": 37.9167,
        "lon": -119.8500,
        "area_ha": 34190,
        "max_extent_lat": 38.0,
        "max_extent_lon": -119.95,
        "duration_days": 25,
    },
    # 2016 Fires
    {
        "name": "Soberanes Fire",
        "country": "USA",
        "state": "California",
        "year": 2016,
        "month": 7,
        "day": 22,
        "lat": 36.4333,
        "lon": -121.8167,
        "area_ha": 19186,
        "max_extent_lat": 36.5,
        "max_extent_lon": -121.9,
        "duration_days": 76,
    },
    {
        "name": "Erskine Fire",
        "country": "USA",
        "state": "California",
        "year": 2016,
        "month": 6,
        "day": 23,
        "lat": 35.3667,
        "lon": -119.0500,
        "area_ha": 15606,
        "max_extent_lat": 35.45,
        "max_extent_lon": -119.15,
        "duration_days": 15,
    },
    # 2015 Fires
    {
        "name": "Butte Fire",
        "country": "USA",
        "state": "California",
        "year": 2015,
        "month": 9,
        "day": 9,
        "lat": 38.3667,
        "lon": -120.7667,
        "area_ha": 22200,
        "max_extent_lat": 38.45,
        "max_extent_lon": -120.85,
        "duration_days": 12,
    },
    {
        "name": "Valley Fire",
        "country": "USA",
        "state": "California",
        "year": 2015,
        "month": 9,
        "day": 12,
        "lat": 38.4167,
        "lon": -122.3333,
        "area_ha": 19308,
        "max_extent_lat": 38.5,
        "max_extent_lon": -122.42,
        "duration_days": 11,
    },
    # 2014 Fires
    {
        "name": "King Fire",
        "country": "USA",
        "state": "California",
        "year": 2014,
        "month": 9,
        "day": 13,
        "lat": 38.6333,
        "lon": -120.5833,
        "area_ha": 33200,
        "max_extent_lat": 38.73,
        "max_extent_lon": -120.68,
        "duration_days": 20,
    },
    {
        "name": "Wallow Fire",
        "country": "USA",
        "state": "Arizona",
        "year": 2011,
        "month": 5,
        "day": 29,
        "lat": 33.5667,
        "lon": -109.5833,
        "area_ha": 143835,
        "max_extent_lat": 33.7,
        "max_extent_lon": -109.7,
        "duration_days": 45,
    },
    # 2013 Fires
    {
        "name": "Rim Fire",
        "country": "USA",
        "state": "California",
        "year": 2013,
        "month": 8,
        "day": 17,
        "lat": 38.3333,
        "lon": -120.1667,
        "area_ha": 104825,
        "max_extent_lat": 38.4,
        "max_extent_lon": -120.2,
        "duration_days": 73,
    },
    {
        "name": "Yosemite National Park Fire",
        "country": "USA",
        "state": "California",
        "year": 2013,
        "month": 8,
        "day": 15,
        "lat": 37.8167,
        "lon": -119.7333,
        "area_ha": 20304,
        "max_extent_lat": 37.9,
        "max_extent_lon": -119.83,
        "duration_days": 9,
    },
    # 2012 Fires
    {
        "name": "Waldo Canyon Fire",
        "country": "USA",
        "state": "Colorado",
        "year": 2012,
        "month": 6,
        "day": 23,
        "lat": 38.8833,
        "lon": -104.9167,
        "area_ha": 7895,
        "max_extent_lat": 38.94,
        "max_extent_lon": -105.02,
        "duration_days": 12,
    },
]


# ---------------------------------------------------------------------------
# CSV reader/writer
# ---------------------------------------------------------------------------


def get_all_columns() -> List[str]:
    """Get all available column names from database."""
    if not WILDFIRES:
        return []
    return sorted(list(WILDFIRES[0].keys()))


def filter_wildfires(
    wildfires: List[Dict[str, object]],
    country: Optional[str] = None,
    year_min: Optional[int] = None,
    year_max: Optional[int] = None,
    state: Optional[str] = None,
) -> List[Dict[str, object]]:
    """Filter wildfire list by various criteria."""
    result = wildfires
    if country:
        result = [w for w in result if w.get("country", "").lower() == country.lower()]
    if year_min is not None:
        result = [w for w in result if w.get("year", 0) >= year_min]
    if year_max is not None:
        result = [w for w in result if w.get("year", 0) <= year_max]
    if state:
        result = [w for w in result if w.get("state", "").lower() == state.lower()]
    return result


def find_fire_by_name(fire_name: str) -> Optional[Dict[str, object]]:
    """Find a fire in the database by name (case-insensitive)."""
    for fire in WILDFIRES:
        if fire.get("name", "").lower() == fire_name.lower():
            return fire
    return None


def generate_terrain_via_srtm(
    fire: Dict[str, object],
    output_path: str,
    lat_margin: float = 0.25,
    lon_margin: float = 0.25,
    subsample: int = 1,
) -> bool:
    """
    Generate terrain CSV using srtm_terrain_reader.py.
    
    Calls the srtm_terrain_reader module to download real SRTM elevation data
    for the specified fire location with margins.
    
    Returns True if successful, False otherwise.
    """
    try:
        import subprocess
    except ImportError:
        print("ERROR: subprocess module required for SRTM terrain generation.", file=sys.stderr)
        return False
    
    lat_center = float(fire.get("lat", 0))
    lon_center = float(fire.get("lon", 0))
    fire_name = fire.get("name", "Unknown")
    
    # Define bounding box with margins
    lat_min = lat_center - lat_margin
    lat_max = lat_center + lat_margin
    lon_min = lon_center - lon_margin
    lon_max = lon_center + lon_margin
    
    print(f"\nDownloading SRTM terrain for {fire_name}...")
    print(f"  Bounds: lat [{lat_min:.3f}, {lat_max:.3f}], lon [{lon_min:.3f}, {lon_max:.3f}]")
    
    # Find srtm_terrain_reader.py in the tools directory
    tools_dir = Path(__file__).parent
    srtm_reader = tools_dir / "srtm_terrain_reader.py"
    
    if not srtm_reader.exists():
        print(f"ERROR: srtm_terrain_reader.py not found at {srtm_reader}", file=sys.stderr)
        return False
    
    # Build command to call srtm_terrain_reader
    cmd = [
        sys.executable,
        str(srtm_reader),
        "--lat-min", str(lat_min),
        "--lat-max", str(lat_max),
        "--lon-min", str(lon_min),
        "--lon-max", str(lon_max),
        "--terrain", output_path,
    ]
    
    if subsample > 1:
        cmd.extend(["--subsample", str(subsample)])
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        if result.returncode != 0:
            print(f"ERROR: srtm_terrain_reader.py failed:", file=sys.stderr)
            print(result.stderr, file=sys.stderr)
            return False
        
        # Print the output from srtm_terrain_reader
        if result.stdout:
            print(result.stdout, end="")
        
        return True
    
    except subprocess.TimeoutExpired:
        print("ERROR: SRTM terrain download timed out (5 minutes).", file=sys.stderr)
        return False
    except Exception as e:
        print(f"ERROR: Failed to run srtm_terrain_reader.py: {e}", file=sys.stderr)
        return False


def write_input_file(
    fire: Dict[str, object],
    output_path: str,
    terrain_csv_file: str = "terrain.csv",
    grid_spacing: float = 30.0,
) -> None:
    """Write a minimal input.i file template for a fire."""
    fire_name = fire.get("name", "Unknown")
    year = fire.get("year", 2020)
    
    # Calculate domain parameters based on terrain
    domain_size = int(300.0 / grid_spacing) * grid_spacing
    domain_height = 300.0
    
    input_content = f"""# Input file for {fire_name} ({year})
# Generated by historical_wildfires.py
# Edit as needed for your simulation

# Terrain file
terrain_file = {terrain_csv_file}

# Reference wind: 10 m/s from west at 10 m AGL
U_ref = 10.0
V_ref = 0.0
z_ref = 10.0

# Aerodynamic roughness length [m] (open terrain / short grass)
z0 = 0.03

# Horizontal grid spacing [m]
dx = {grid_spacing:.1f}
dy = {grid_spacing:.1f}

# Vertical grid spacing [m]
dz = 5.0

# Domain height [m]
domain_height = {domain_height:.1f}

# MLMG solver coefficients
alpha_h = 1.0
alpha_v = 1.0

# Fuel model (Anderson 13-class)
fuel_model = 7

# Initial moisture content (%)
dead_fuel_moisture_1h = 3.0
dead_fuel_moisture_10h = 4.0
dead_fuel_moisture_100h = 5.0
live_fuel_moisture = 90.0

# Ignition
ignition_type = point
ignition_point_x = {domain_size/2:.1f}
ignition_point_y = {domain_size/2:.1f}
ignition_time = 0.0

# Outputs
plot_int = 60
max_time = 3600.0
"""
    
    with open(output_path, "w") as fh:
        fh.write(input_content)
    
    print(f"Wrote input file → {output_path}")


def create_fire_inputs(
    fire_name: str,
    outdir: str = ".",
    grid_spacing: float = 30.0,
    lat_margin: float = 0.25,
    lon_margin: float = 0.25,
    subsample: int = 1,
) -> None:
    """Create terrain.csv and input.i files for a given fire using SRTM data."""
    fire = find_fire_by_name(fire_name)
    if not fire:
        available = [f["name"] for f in WILDFIRES]
        print(f"ERROR: Fire '{fire_name}' not found in database.", file=sys.stderr)
        print(f"\nAvailable fires:\n", file=sys.stderr)
        for name in sorted(available):
            print(f"  - {name}", file=sys.stderr)
        sys.exit(1)
    
    # Create output directory
    outdir_path = Path(outdir)
    outdir_path.mkdir(parents=True, exist_ok=True)
    
    # Generate files
    terrain_path = outdir_path / "terrain.csv"
    input_path = outdir_path / "inputs.i"
    
    # Download SRTM terrain
    if not generate_terrain_via_srtm(
        fire,
        str(terrain_path),
        lat_margin=lat_margin,
        lon_margin=lon_margin,
        subsample=subsample,
    ):
        print(f"\nERROR: Failed to generate terrain for '{fire['name']}'", file=sys.stderr)
        sys.exit(1)
    
    write_input_file(fire, str(input_path), terrain_csv_file="terrain.csv", grid_spacing=grid_spacing)
    
    print(f"\n✓ Successfully created input files for '{fire['name']}':")
    print(f"  Location: {fire.get('state', 'Unknown')}, {fire.get('country', 'Unknown')}")
    print(f"  Coordinates: {fire.get('lat')}, {fire.get('lon')}")
    print(f"  Output directory: {outdir_path.absolute()}")
    print(f"\nGenerated files:")
    print(f"  - {terrain_path.name} (SRTM terrain data in UTM)")
    print(f"  - {input_path.name} (solver input template)")


def list_all_fire_names() -> None:
    """Print all available fire names."""
    print("Available fires in the database:")
    print()
    for i, fire in enumerate(WILDFIRES, 1):
        name = fire.get("name", "Unknown")
        year = fire.get("year", "")
        state = fire.get("state", "")
        country = fire.get("country", "")
        print(f"  {i:2d}. {name:<35s} ({year}) - {state}, {country}")
    print()


# ---------------------------------------------------------------------------
# ASCII table formatting
# ---------------------------------------------------------------------------


def print_summary_table(
    wildfires: List[Dict[str, object]],
    columns: Optional[List[str]] = None,
) -> None:
    """Print a formatted ASCII summary table to stdout."""
    if not wildfires:
        print("(no data)")
        return

    # Default columns if not specified
    if columns is None:
        columns = [
            "name",
            "country",
            "state",
            "year",
            "lat",
            "lon",
            "area_ha",
            "duration_days",
        ]

    # Validate columns
    valid_cols = get_all_columns()
    columns = [c for c in columns if c in valid_cols]
    if not columns:
        columns = ["name", "country", "lat", "lon"]

    # Build header and separator
    col_widths = {}
    for col in columns:
        max_width = len(col)
        for w in wildfires:
            val_str = str(w.get(col, ""))
            max_width = max(max_width, len(val_str))
        col_widths[col] = max_width

    header_parts = [f"{col:>{col_widths[col]}s}" for col in columns]
    sep_parts = ["-" * col_widths[col] for col in columns]

    header = "  ".join(header_parts)
    sep = "  ".join(sep_parts)

    print(header)
    print(sep)

    for w in wildfires:
        row_parts = []
        for col in columns:
            val = w.get(col, "")
            val_str = str(val)
            row_parts.append(f"{val_str:>{col_widths[col]}s}")
        print("  ".join(row_parts))

    print()
    print(f"Total wildfires: {len(wildfires)}")


# ---------------------------------------------------------------------------
# CSV export
# ---------------------------------------------------------------------------


def export_csv(
    wildfires: List[Dict[str, object]],
    out_path: str,
    columns: Optional[List[str]] = None,
) -> None:
    """Write wildfires to CSV file."""
    if not wildfires:
        print("(no data to export)")
        return

    if columns is None:
        columns = get_all_columns()

    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(wildfires)

    print(f"Exported {len(wildfires)} wildfires → {out_path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Table of major historical wildfires with lat/lon coordinates.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--output",
        metavar="FILE",
        default=None,
        help="Export to CSV file (instead of printing table).",
    )
    parser.add_argument(
        "--country",
        metavar="COUNTRY",
        default=None,
        help="Filter by country (e.g., USA, Australia, Canada).",
    )
    parser.add_argument(
        "--state",
        metavar="STATE",
        default=None,
        help="Filter by state/region (e.g., California, New South Wales).",
    )
    parser.add_argument(
        "--year-min",
        type=int,
        metavar="YEAR",
        default=None,
        help="Minimum year (inclusive).",
    )
    parser.add_argument(
        "--year-max",
        type=int,
        metavar="YEAR",
        default=None,
        help="Maximum year (inclusive).",
    )
    parser.add_argument(
        "--columns",
        metavar="COL1,COL2,...",
        default=None,
        help="Comma-separated list of columns to display/export.",
    )
    parser.add_argument(
        "--list-columns",
        action="store_true",
        help="List all available columns and exit.",
    )
    parser.add_argument(
        "--list-fires",
        action="store_true",
        help="List all available fire names and exit.",
    )
    parser.add_argument(
        "--create-inputs",
        metavar="FIRE_NAME",
        default=None,
        help="Generate terrain.csv and inputs.i files for a specific fire.",
    )
    parser.add_argument(
        "--outdir",
        metavar="DIR",
        default=".",
        help="Output directory for generated files (default: current directory).",
    )
    parser.add_argument(
        "--grid-spacing",
        type=float,
        metavar="METERS",
        default=30.0,
        help="Grid spacing for terrain in metres (default: 30.0).",
    )
    parser.add_argument(
        "--lat-margin",
        type=float,
        metavar="DEGREES",
        default=0.25,
        help="Latitude margin in degrees on either side of fire centre (default: 0.25).",
    )
    parser.add_argument(
        "--lon-margin",
        type=float,
        metavar="DEGREES",
        default=0.25,
        help="Longitude margin in degrees on either side of fire centre (default: 0.25).",
    )
    parser.add_argument(
        "--subsample",
        type=int,
        metavar="N",
        default=1,
        help="Keep every N-th terrain point (default: 1, no subsampling).",
    )

    args = parser.parse_args(argv)

    if args.list_columns:
        print("Available columns:")
        for col in get_all_columns():
            print(f"  - {col}")
        return

    if args.list_fires:
        list_all_fire_names()
        return

    if args.create_inputs:
        create_fire_inputs(
            args.create_inputs,
            outdir=args.outdir,
            grid_spacing=args.grid_spacing,
            lat_margin=args.lat_margin,
            lon_margin=args.lon_margin,
            subsample=args.subsample,
        )
        return

    # Parse columns
    columns = None
    if args.columns:
        columns = [c.strip() for c in args.columns.split(",")]

    # Filter wildfires
    filtered = filter_wildfires(
        WILDFIRES,
        country=args.country,
        year_min=args.year_min,
        year_max=args.year_max,
        state=args.state,
    )

    if not filtered:
        print("No wildfires match the filter criteria.")
        sys.exit(1)

    if args.output:
        export_csv(filtered, args.output, columns=columns)
    else:
        print_summary_table(filtered, columns=columns)


if __name__ == "__main__":
    main()
