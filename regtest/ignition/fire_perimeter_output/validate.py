#!/usr/bin/env python3
"""
Validate fire_perimeter_output regression test.

Checks that:
1. fire_stats.csv was written and has the expected column headers.
2. At least one perimeter_NNNN.csv was written and contains X,Y columns.
3. At least one perimeter_NNNN.geojson was written and is valid JSON with
   a Polygon geometry.
"""
import json
import sys
from pathlib import Path
import glob


errors = []

# ── 1. fire_stats.csv ────────────────────────────────────────────────────────
stats_file = Path("fire_stats.csv")
if not stats_file.exists():
    errors.append("MISSING: fire_stats.csv")
else:
    with open(stats_file) as f:
        header = f.readline().strip()
    expected_cols = ["step", "time_s", "burned_area_ha", "perimeter_km"]
    for col in expected_cols:
        if col not in header:
            errors.append(f"fire_stats.csv missing column: {col}")
    # Check at least one data row
    with open(stats_file) as f:
        rows = f.readlines()
    if len(rows) < 2:
        errors.append("fire_stats.csv has no data rows")
    else:
        print(f"  fire_stats.csv: {len(rows)-1} data rows  ✓")

# ── 2. perimeter CSV ─────────────────────────────────────────────────────────
csv_files = sorted(glob.glob("perimeter_*.csv"))
if not csv_files:
    errors.append("MISSING: no perimeter_NNNN.csv files")
else:
    with open(csv_files[0]) as f:
        lines = [l.strip() for l in f if l.strip() and not l.startswith("#")]
    if not lines:
        errors.append(f"{csv_files[0]}: empty file")
    elif lines[0].lower() not in ("x,y", "x, y"):
        errors.append(f"{csv_files[0]}: expected header 'X,Y', got '{lines[0]}'")
    elif len(lines) < 2:
        errors.append(f"{csv_files[0]}: no coordinate rows")
    else:
        print(f"  {csv_files[0]}: {len(lines)-1} coordinate rows  ✓")

# ── 3. perimeter GeoJSON ──────────────────────────────────────────────────────
gjson_files = sorted(glob.glob("perimeter_*.geojson"))
if not gjson_files:
    errors.append("MISSING: no perimeter_NNNN.geojson files")
else:
    try:
        with open(gjson_files[0]) as f:
            gj = json.load(f)
        assert gj.get("type") == "FeatureCollection", "not a FeatureCollection"
        feats = gj.get("features", [])
        assert feats, "empty features list"
        geom = feats[0].get("geometry", {})
        assert geom.get("type") == "Polygon", (
            f"expected Polygon, got {geom.get('type')}")
        coords = geom.get("coordinates", [[]])[0]
        assert len(coords) >= 4, "polygon ring too short"
        print(f"  {gjson_files[0]}: {len(coords)} ring points  ✓")
    except Exception as e:
        errors.append(f"{gjson_files[0]}: {e}")

# ── Report ────────────────────────────────────────────────────────────────────
if errors:
    print("\nFAILED:")
    for e in errors:
        print(f"  ERROR: {e}")
    sys.exit(1)

print("\nAll fire_perimeter_output checks passed.")
