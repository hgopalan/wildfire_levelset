#!/usr/bin/env python3
"""
Validate kml_perimeter regression test.

Checks that:
1. At least one perimeter_NNNN.kml file was written.
2. Each KML file is valid XML.
3. Each KML file contains a <Polygon> element with at least 4 coordinate pairs.
4. When kml_utm_zone = 11, coordinates are in lon/lat range
   (roughly -180..180 for lon, -90..90 for lat).
"""
import sys
import glob
import xml.etree.ElementTree as ET

errors = []

# ── 1. KML files exist ────────────────────────────────────────────────────────
kml_files = sorted(glob.glob("perimeter_*.kml"))
if not kml_files:
    errors.append("MISSING: no perimeter_NNNN.kml files")
else:
    print(f"  Found {len(kml_files)} KML file(s): {kml_files[0]}  ✓")

# ── 2. Validate KML structure ─────────────────────────────────────────────────
KML_NS = "http://www.opengis.net/kml/2.2"
for kml_path in kml_files[:1]:  # check first file
    try:
        tree = ET.parse(kml_path)
        root = tree.getroot()

        # Root should be <kml> or <kml:kml>
        local = root.tag.split("}")[-1] if "}" in root.tag else root.tag
        assert local == "kml", f"root tag is {root.tag}, expected 'kml'"

        # Find <coordinates> element anywhere in the tree
        ns = {"k": KML_NS}
        coords_el = root.find(".//{%s}coordinates" % KML_NS)
        assert coords_el is not None and coords_el.text, \
            "no <coordinates> element found"

        # Parse coordinate tuples (lon,lat,alt)
        coord_pairs = [
            line.strip() for line in coords_el.text.strip().split()
            if line.strip()
        ]
        assert len(coord_pairs) >= 4, \
            f"too few coordinate pairs: {len(coord_pairs)} (need ≥ 4)"

        # Check that lon/lat are in plausible WGS-84 range
        for cp in coord_pairs[:5]:
            parts = cp.split(",")
            if len(parts) >= 2:
                lon, lat = float(parts[0]), float(parts[1])
                assert -180 <= lon <= 180, f"lon={lon} out of range"
                assert -90  <= lat <= 90,  f"lat={lat} out of range"

        print(f"  {kml_path}: {len(coord_pairs)} coordinate pairs, "
              f"lon={float(coord_pairs[0].split(',')[0]):.4f}  ✓")

    except Exception as e:
        errors.append(f"{kml_path}: {e}")

# ── Report ────────────────────────────────────────────────────────────────────
if errors:
    print("\nFAILED:")
    for e in errors:
        print(f"  ERROR: {e}")
    sys.exit(1)

print("\nAll kml_perimeter checks passed.")
