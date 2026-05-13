#!/usr/bin/env python3
"""generate_plt_wind.py — Synthetic AMReX plt wind file for ember_cascade_flux regression test.

Generates a minimal single-level AMReX plotfile containing 3-D wind data
(u, v, w) for use in the ``ember_cascade_flux`` regression test.

The file mimics the output of massconsistent_amr/src/wind_solver.cpp
(``WriteSingleLevelPlotfile``), producing the following directory layout::

    plt_wind_3d/
        Header
        Level_0/
            Cell_H
            Cell_D_00000

Grid: 8 × 8 × 4  (coarse, fast to generate/read)
Domain: 329 900 – 330 500 E,  3 774 900 – 3 775 500 N,  0 – 40 m  (UTM Zone 11N)
Wind:   u = 6.0 m/s (westerly),  v = 1.0 m/s (slight cross-wind),  w = 0.0 m/s

Wind speed is set slightly stronger than in the albini_spotting_3d_wind test so
that the ember cascade model produces a wider landing footprint to exercise the
Gaussian kernel over a meaningful fraction of the domain.

Usage::

    cd regtest/spotting/ember_cascade_flux
    python3 generate_plt_wind.py   # writes plt_wind_3d/ directory

The ``plt_wind_reader.H`` C++ reader will then open this directory at
simulation start when ``ember_cascade.use_3d_wind = 1`` and
``ember_cascade.plt_wind_file = plt_wind_3d`` are set in the inputs file.

Notes
-----
* Variable names *must* include ``"u"``, ``"v"``, and ``"w"`` exactly (case-
  sensitive) because that is what ``plt_wind_reader.H`` looks for when scanning
  the plotfile header.
* Data layout inside ``Cell_D_00000`` follows the AMReX / Fortran column-major
  convention: component-major with x as the fastest index,
  ``flat = comp*(nx*ny*nz) + k*(nx*ny) + j*nx + i``.
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

# Uniform 3-D wind
U_WIND = 6.0    # m/s eastward (dominant transport direction)
V_WIND = 1.0    # m/s northward (slight cross-wind)
W_WIND = 0.0    # m/s vertical

VAR_NAMES = ["u", "v", "w", "vel_magnitude"]
NVAR = len(VAR_NAMES)

OUT_DIR = "plt_wind_3d"
LEVEL_DIR = os.path.join(OUT_DIR, "Level_0")


def build_fab_data(nx, ny, nz, var_names, wind_vals):
    """Build the binary FAB payload for a uniform 3-D wind field.

    Parameters
    ----------
    nx, ny, nz : int
        Number of cells in each spatial direction.
    var_names : list[str]
        Ordered list of variable names.
    wind_vals : dict
        Mapping from variable name to its uniform scalar value.

    Returns
    -------
    bytes
        Raw binary payload (``nvar * nx * ny * nz`` IEEE-754 doubles).
    """
    n_cells = nx * ny * nz
    data = bytearray()
    for vname in var_names:
        val = wind_vals[vname]
        for _ in range(n_cells):
            data += struct.pack("d", val)
    return bytes(data)


def write_level0(level_dir, nx, ny, nz, var_names, wind_vals):
    """Write the ``Level_0/`` directory containing the FAB data and VisMF header."""
    os.makedirs(level_dir, exist_ok=True)
    nvar = len(var_names)

    fab_data = build_fab_data(nx, ny, nz, var_names, wind_vals)

    fab_header_str = (
        f"FAB ((8,(64 11 52 0 1 12 0 1023)),(8,(64 11 52 0 1 12 0 1023)))"
        f" (({0},{0},{0}) ({nx-1},{ny-1},{nz-1}) ({0},{0},{0}))\n"
    )
    fab_header_bytes = fab_header_str.encode("ascii")

    data_file_path = os.path.join(level_dir, "Cell_D_00000")
    with open(data_file_path, "wb") as f:
        f.write(fab_header_bytes)
        f.write(fab_data)

    box_str = f"(({0},{0},{0}) ({nx-1},{ny-1},{nz-1}) ({0},{0},{0}))"
    cell_h_path = os.path.join(level_dir, "Cell_H")
    with open(cell_h_path, "w") as f:
        f.write("VisMF-V1\n")
        f.write("1\n")
        f.write(f"{nvar}\n")
        f.write("0\n")
        f.write(f"{box_str} {nvar}  Cell_D_00000:0\n")


def write_plt_header(out_dir, nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax,
                     var_names):
    """Write the top-level AMReX plotfile ``Header`` file."""
    dx = (xmax - xmin) / nx
    dy = (ymax - ymin) / ny
    dz = (zmax - zmin) / nz
    nvar = len(var_names)

    bbox_str = f"(({0},{0},{0}) ({nx-1},{ny-1},{nz-1}) ({0},{0},{0}))"
    header_path = os.path.join(out_dir, "Header")
    with open(header_path, "w") as f:
        f.write("HyperCLaw-V1.1\n")
        f.write(f"{nvar}\n")
        for vname in var_names:
            f.write(f"{vname}\n")
        f.write("3\n")
        f.write("0.0\n")
        f.write("0\n")
        f.write(f"{xmin} {ymin} {zmin}\n")
        f.write(f"{xmax} {ymax} {zmax}\n")
        f.write(f"{bbox_str}\n")
        f.write("0\n")
        f.write(f"{dx} {dy} {dz}\n")
        f.write("0\n")
        f.write("0\n")
        f.write("Level_0/Cell\n")


def main():
    """Generate the synthetic plt wind directory for the ember_cascade_flux test.

    Writes ``plt_wind_3d/`` (relative to the current working directory)
    with a uniform u = 6 m/s, v = 1 m/s wind on an 8 × 8 × 4 grid.
    """
    os.makedirs(LEVEL_DIR, exist_ok=True)

    wind_vals = {
        "u": U_WIND,
        "v": V_WIND,
        "w": W_WIND,
        "vel_magnitude": (U_WIND**2 + V_WIND**2 + W_WIND**2) ** 0.5,
    }

    write_level0(LEVEL_DIR, NX, NY, NZ, VAR_NAMES, wind_vals)
    write_plt_header(OUT_DIR, NX, NY, NZ, XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX,
                     VAR_NAMES)

    print(f"Generated {OUT_DIR}/ with {NX}x{NY}x{NZ} grid")
    print(f"  u={U_WIND} m/s  v={V_WIND} m/s  w={W_WIND} m/s (uniform)")
    print(f"  Domain: ({XMIN},{YMIN},{ZMIN}) to ({XMAX},{YMAX},{ZMAX})")


if __name__ == "__main__":
    main()
