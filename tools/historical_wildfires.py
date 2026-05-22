#!/usr/bin/env python3
"""
historical_wildfires.py – Table of major historical wildfires with lat/lon.

This tool provides a curated database of significant wildfires from the last
10+ years with location coordinates (latitude/longitude) that can be used for:

  1. Validating simulation results against real fire behavior
  2. Benchmarking fire spread models
  3. Comparison studies

Usage
-----
  # Print table to stdout
  python3 tools/historical_wildfires.py

  # Export to CSV
  python3 tools/historical_wildfires.py --output wildfires.csv

  # Filter by country/region
  python3 tools/historical_wildfires.py --country USA

  # Filter by year range
  python3 tools/historical_wildfires.py --year-min 2020 --year-max 2023

  # Display specific columns
  python3 tools/historical_wildfires.py --columns "name,country,year,lat,lon,area_ha"

References
----------
  Data compiled from NIFC (National Interagency Fire Center), InciWeb,
  Copernicus Emergency Management Service, and published case studies.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import Dict, List, Optional

# ---------------------------------------------------------------------------
# Wildfire Database
# ---------------------------------------------------------------------------

WILDFIRES: List[Dict[str, object]] = [
    # USA - Recent Major Fires
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
    # Australia
    {
        "name": "Black Summer Bushfires",
        "country": "Australia",
        "state": "New South Wales",
        "year": 2019,
        "month": 9,
        "day": 1,
        "lat": -33.8688,
        "lon": 151.2093,
        "area_ha": 1043000,
        "max_extent_lat": -33.9,
        "max_extent_lon": 151.2,
        "duration_days": 200,
    },
    {
        "name": "Tasmanian Bushfires",
        "country": "Australia",
        "state": "Tasmania",
        "year": 2016,
        "month": 1,
        "day": 3,
        "lat": -42.1667,
        "lon": 147.3333,
        "area_ha": 188000,
        "max_extent_lat": -42.17,
        "max_extent_lon": 147.33,
        "duration_days": 87,
    },
    # Canada
    {
        "name": "Park Fire (British Columbia)",
        "country": "Canada",
        "state": "British Columbia",
        "year": 2017,
        "month": 7,
        "day": 16,
        "lat": 54.1667,
        "lon": -120.3333,
        "area_ha": 239000,
        "max_extent_lat": 54.17,
        "max_extent_lon": -120.33,
        "duration_days": 178,
    },
    # Greece
    {
        "name": "Mati Wildfire",
        "country": "Greece",
        "state": "Attica",
        "year": 2018,
        "month": 7,
        "day": 23,
        "lat": 38.0333,
        "lon": 23.8667,
        "area_ha": 2400,
        "max_extent_lat": 38.04,
        "max_extent_lon": 23.87,
        "duration_days": 1,
    },
    # Portugal
    {
        "name": "Pedrógão Grande Fire",
        "country": "Portugal",
        "state": "Covilhã",
        "year": 2017,
        "month": 6,
        "day": 17,
        "lat": 40.0333,
        "lon": -7.3667,
        "area_ha": 50000,
        "max_extent_lat": 40.04,
        "max_extent_lon": -7.37,
        "duration_days": 4,
    },
    # Spain
    {
        "name": "Estepona Fire",
        "country": "Spain",
        "state": "Málaga",
        "year": 2012,
        "month": 9,
        "day": 2,
        "lat": 36.4333,
        "lon": -5.1667,
        "area_ha": 8150,
        "max_extent_lat": 36.44,
        "max_extent_lon": -5.17,
        "duration_days": 3,
    },
    # Russia
    {
        "name": "Russian Western Siberia Fires",
        "country": "Russia",
        "state": "Siberia",
        "year": 2021,
        "month": 7,
        "day": 1,
        "lat": 66.0,
        "lon": 93.0,
        "area_ha": 1000000,
        "max_extent_lat": 66.5,
        "max_extent_lon": 93.5,
        "duration_days": 180,
    },
    # Mediterranean region
    {
        "name": "Turkey Wildfires",
        "country": "Turkey",
        "state": "Antalya",
        "year": 2021,
        "month": 7,
        "day": 28,
        "lat": 36.5333,
        "lon": 30.7667,
        "area_ha": 400000,
        "max_extent_lat": 36.6,
        "max_extent_lon": 30.8,
        "duration_days": 150,
    },
    # Indonesia
    {
        "name": "Kalimantan Peatland Fires",
        "country": "Indonesia",
        "state": "Kalimantan",
        "year": 2015,
        "month": 8,
        "day": 1,
        "lat": -2.3333,
        "lon": 113.6667,
        "area_ha": 2600000,
        "max_extent_lat": -2.5,
        "max_extent_lon": 113.8,
        "duration_days": 180,
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

    args = parser.parse_args(argv)

    if args.list_columns:
        print("Available columns:")
        for col in get_all_columns():
            print(f"  - {col}")
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
