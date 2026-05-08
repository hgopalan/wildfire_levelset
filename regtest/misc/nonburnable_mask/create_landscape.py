#!/usr/bin/env python3
"""
create_nonburnable_landscape.py
Generate a synthetic LCP-format landscape for the nonburnable_mask regtest.

The landscape is 64x64 cells at ~6.25 m resolution covering a 400 m x 400 m
UTM Zone 11N domain.  The LEFT half (x < domain centre) is FM4 (chaparral,
code 4) and the RIGHT half is code 98 (Open Water, non-burnable).

The fire ignites on the left and should be blocked at the boundary.

Output: landscape_nb.lcp   (ASCII CSV: X Y ELEVATION SLOPE ASPECT FUEL_MODEL)
"""
import os
import sys

NX = 64
NY = 64
SPACING = 400.0 / NX  # ~6.25 m
X_ORIG  = 330000.0
Y_ORIG  = 3775000.0

OUTPUT = "landscape_nb.lcp"

rows = []
for j in range(NY):
    for i in range(NX):
        x = X_ORIG + (i + 0.5) * SPACING
        y = Y_ORIG + (j + 0.5) * SPACING
        # Flat terrain at 200 m elevation
        elev  = 200.0
        slope = 0.0
        aspect = 0.0
        # Left half = FM4 (code 4); right half = Open Water (code 98)
        fuel = 4 if (i < NX // 2) else 98
        rows.append((x, y, elev, slope, aspect, fuel))

with open(OUTPUT, "w") as fh:
    fh.write(
        "# Synthetic non-burnable mask landscape for regtest\n"
        "# Left half: FM4 (chaparral), right half: code 98 (Open Water)\n"
        "# Format: X Y ELEVATION SLOPE ASPECT FUEL_MODEL\n"
        f"# Grid: {NX}x{NY} at {SPACING:.2f} m spacing\n"
    )
    for (x, y, elev, slope, aspect, fuel) in rows:
        fh.write(f"{x:.2f} {y:.2f} {elev:.2f} {slope:.2f} {aspect:.2f} {fuel}\n")

n = len(rows)
print(f"Created landscape '{OUTPUT}' with {n} rows  "
      f"({NX//2}x{NY} burnable, {NX//2}x{NY} non-burnable)")
