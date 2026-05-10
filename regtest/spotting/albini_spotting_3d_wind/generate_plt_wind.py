#!/usr/bin/env python3
"""
generate_plt_wind.py
Generate a minimal single-level AMReX plotfile containing 3-D wind data (u, v, w)
for use in the albini_spotting_3d_wind regression test.

The file mimics the output of massconsistent_amr/src/wind_solver.cpp
(WriteSingleLevelPlotfile), producing a directory structure:

    plt_wind_3d/
        Header
        Level_0/
            Cell_H
            Cell_D_00000

Grid: 8 x 8 x 4  (coarse, fast to generate/read)
Domain: 329900 – 330500 E,  3774900 – 3775500 N,  0 – 40 m  (UTM Zone 11N)
Wind:   u = 5.0 m/s (westerly),  v = 0.5 m/s,  w = 0.0 m/s (uniform)
"""

import struct
import os

# ---- Grid & domain ----
NX, NY, NZ = 8, 8, 4
XMIN, XMAX = 329900.0, 330500.0
YMIN, YMAX = 3774900.0, 3775500.0
ZMIN, ZMAX = 0.0, 40.0

DX = (XMAX - XMIN) / NX
DY = (YMAX - YMIN) / NY
DZ = (ZMAX - ZMIN) / NZ

# Uniform 3-D wind (u, v, w)
U_WIND = 5.0
V_WIND = 0.5
W_WIND = 0.0

# Variable names (must include "u", "v", "w" for the reader)
VAR_NAMES = ["u", "v", "w", "vel_magnitude"]
NVAR = len(VAR_NAMES)

OUT_DIR = "plt_wind_3d"
LEVEL_DIR = os.path.join(OUT_DIR, "Level_0")

os.makedirs(LEVEL_DIR, exist_ok=True)

# ---- Build the full-domain FAB data ----
# Data layout: component-major, then z, y, x (Fortran column-major)
# i.e. flat index = comp * NX*NY*NZ + k*NX*NY + j*NX + i
N_CELLS = NX * NY * NZ
wind_vals = {
    "u": U_WIND,
    "v": V_WIND,
    "w": W_WIND,
    "vel_magnitude": (U_WIND**2 + V_WIND**2 + W_WIND**2)**0.5,
}

# Build binary data as doubles (8 bytes each)
data_bytes = bytearray()
for vname in VAR_NAMES:
    val = wind_vals[vname]
    for _ in range(N_CELLS):
        data_bytes += struct.pack("d", val)

# ---- FAB ASCII header (included at the start of the binary data file) ----
# Format used by AMReX WriteSingleLevelPlotfile:
#   FAB ((8,(64 11 52 0 1 12 0 1023)),(8,(64 11 52 0 1 12 0 1023))) ((lo0,lo1,lo2) (hi0,hi1,hi2) (0,0,0))
fab_header_str = (
    f"FAB ((8,(64 11 52 0 1 12 0 1023)),(8,(64 11 52 0 1 12 0 1023)))"
    f" (({0},{0},{0}) ({NX-1},{NY-1},{NZ-1}) ({0},{0},{0}))\n"
)
fab_header_bytes = fab_header_str.encode("ascii")

# ---- Write Level_0/Cell_D_00000 ----
data_file_path = os.path.join(LEVEL_DIR, "Cell_D_00000")
with open(data_file_path, "wb") as f:
    f.write(fab_header_bytes)
    f.write(bytes(data_bytes))

fab_offset = 0  # offset of the FAB header in Cell_D_00000

# ---- Write Level_0/Cell_H (VisMF header) ----
# AMReX VisMF-V1 format:
#   VisMF-V1
#   <nFAB>
#   <nComp>
#   <nGhost>
#   <box_descriptor> <nComp>  Cell_D_00000:<offset>
cell_h_path = os.path.join(LEVEL_DIR, "Cell_H")
box_str = f"(({0},{0},{0}) ({NX-1},{NY-1},{NZ-1}) ({0},{0},{0}))"
with open(cell_h_path, "w") as f:
    f.write("VisMF-V1\n")
    f.write(f"1\n")          # nFAB
    f.write(f"{NVAR}\n")     # nComp
    f.write("0\n")           # nGhost
    f.write(f"{box_str} {NVAR}  Cell_D_00000:{fab_offset}\n")

# ---- Write Header ----
# Standard AMReX plotfile header
header_path = os.path.join(OUT_DIR, "Header")
# bounding box for level 0
bbox_str = f"(({0},{0},{0}) ({NX-1},{NY-1},{NZ-1}) ({0},{0},{0}))"
with open(header_path, "w") as f:
    f.write("HyperCLaw-V1.1\n")
    f.write(f"{NVAR}\n")
    for vname in VAR_NAMES:
        f.write(f"{vname}\n")
    f.write("3\n")           # spacedim
    f.write("0.0\n")         # time
    f.write("0\n")           # finest_level
    f.write(f"{XMIN} {YMIN} {ZMIN}\n")  # prob_lo
    f.write(f"{XMAX} {YMAX} {ZMAX}\n")  # prob_hi
    # No refinement ratios (finest_level == 0)
    f.write(f"{bbox_str}\n") # level 0 bounding box
    f.write("0\n")           # step count for level 0
    f.write(f"{DX} {DY} {DZ}\n")  # cell sizes at level 0
    f.write("0\n")           # coordinate system (Cartesian)
    f.write("0\n")           # reserved
    f.write("Level_0/Cell\n")  # level 0 data path

print(f"Generated {OUT_DIR}/ with {NX}x{NY}x{NZ} grid")
print(f"  u={U_WIND} m/s, v={V_WIND} m/s, w={W_WIND} m/s (uniform)")
print(f"  Domain: ({XMIN},{YMIN},{ZMIN}) to ({XMAX},{YMAX},{ZMAX})")
