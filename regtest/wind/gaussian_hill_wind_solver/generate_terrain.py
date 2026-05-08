#!/usr/bin/env python3
"""
generate_terrain.py
Generate a Gaussian hill terrain file for the wind_solver regression test.

Domain : 0–300 m × 0–300 m at 30 m spacing (11 × 11 = 121 data points)
Hill   : centre (150, 150) m, peak height 50 m, sigma 60 m

Output : terrain.csv  (whitespace-separated X Y Z)
"""
import math
import os

NX, NY = 11, 11
STEP = 30.0          # [m] grid spacing
CX, CY = 150.0, 150.0  # hill centre [m]
PEAK = 50.0          # peak elevation [m]
SIGMA = 60.0         # Gaussian half-width [m]
DENOM = 2.0 * SIGMA ** 2

out_path = os.path.join(os.path.dirname(__file__), "terrain.csv")

with open(out_path, "w") as fh:
    fh.write("# Gaussian hill terrain  X[m]  Y[m]  Z[m]\n")
    fh.write("# Domain: 0-300 x 0-300 m, peak=50 m at (150,150), sigma=60 m\n")
    for j in range(NY):
        y = j * STEP
        dy2 = (y - CY) ** 2
        for i in range(NX):
            x = i * STEP
            dx2 = (x - CX) ** 2
            z = PEAK * math.exp(-(dx2 + dy2) / DENOM)
            fh.write(f"{x:.1f} {y:.1f} {z:.6f}\n")

print(f"Wrote {NX * NY} terrain points to {out_path}")
