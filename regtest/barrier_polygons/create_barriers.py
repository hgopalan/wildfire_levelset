#!/usr/bin/env python3
"""
Setup step for barrier_polygons regression test.
Creates two small CSV barrier files:
  - barrier_north.csv  : east-west line across the northern part of the domain
  - barrier_east.csv   : short north-south stub blocking the downwind path
"""
from pathlib import Path

# North barrier: horizontal line at y=350 m from x=100 to x=400 m
# (domain is 500 x 500; fire starts at (100,250) so this is 100 m ahead)
north_pts = [(x, 350.0) for x in range(100, 410, 10)]
with open("barrier_north.csv", "w") as f:
    f.write("# X Y\n")
    for x, y in north_pts:
        f.write(f"{x:.1f} {y:.1f}\n")
print(f"Created barrier_north.csv  ({len(north_pts)} points)")

# East barrier: short vertical stub at x=350 m from y=200 to y=300 m
east_pts = [(350.0, y) for y in range(200, 310, 10)]
with open("barrier_east.csv", "w") as f:
    f.write("# X Y\n")
    for x, y in east_pts:
        f.write(f"{x:.1f} {y:.1f}\n")
print(f"Created barrier_east.csv  ({len(east_pts)} points)")
