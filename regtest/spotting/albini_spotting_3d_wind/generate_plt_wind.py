#!/usr/bin/env python3
"""generate_plt_wind.py — Synthetic AMReX plt wind file for Albini 3-D spotting tests.

Generates a minimal single-level AMReX plotfile containing 3-D wind data
(u, v, w) for use in the ``albini_spotting_3d_wind`` regression test.

The file mimics the output of massconsistent_amr/src/wind_solver.cpp
(``WriteSingleLevelPlotfile``), producing the following directory layout::

    plt_wind_3d/
        Header
        Level_0/
            Cell_H
            Cell_D_00000

Grid: 8 × 8 × 4  (coarse, fast to generate/read)
Domain: 329 900 – 330 500 E,  3 774 900 – 3 775 500 N,  0 – 40 m  (UTM Zone 11N)
Wind:   u = 5.0 m/s (westerly),  v = 0.5 m/s,  w = 0.0 m/s (uniform)

Usage::

    cd regtest/spotting/albini_spotting_3d_wind
    python3 generate_plt_wind.py   # writes plt_wind_3d/ directory

The ``plt_wind_reader.H`` C++ reader will then open this directory at
simulation start when ``albini_spotting.use_3d_wind = 1`` and
``albini_spotting.plt_wind_file = plt_wind_3d`` are set in the inputs file.

GPU pathway note
----------------
In GPU builds (AMREX_USE_GPU) the reader uploads the wind data to the device
(``PltWindData::d_u2d`` / ``d_v2d``).  The Albini spotting routine then uses
``amrex::ParallelFor`` to interpolate the 3-D wind onto the fire-model 2-D
cell grid on the GPU before copying the result back to the host for the serial
trajectory integration step.  The CPU pathway (for non-GPU builds) reads the
same host-side arrays (``PltWindData::u2d`` / ``v2d``) in a plain double loop.
Both pathways produce bit-identical results for a given wind field.
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


def build_fab_data(nx, ny, nz, var_names, wind_vals):
    """Build the binary FAB payload for a uniform 3-D wind field.

    Data layout follows the AMReX / Fortran column-major convention used by
    ``WriteSingleLevelPlotfile``: component-major ordering with the fastest
    index in the x direction.

        flat_index = comp * (nx*ny*nz) + k*(nx*ny) + j*nx + i

    Parameters
    ----------
    nx, ny, nz : int
        Number of cells in each spatial direction.
    var_names : list of str
        Ordered list of variable names (e.g. ``["u", "v", "w"]``).
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
    """Write the ``Level_0/`` directory containing the FAB data and VisMF header.

    Creates two files inside *level_dir*:

    ``Cell_D_00000``
        ASCII FAB header line followed by the raw binary cell data.
    ``Cell_H``
        AMReX VisMF-V1 header pointing to ``Cell_D_00000``.

    Parameters
    ----------
    level_dir : str
        Path to the ``Level_0`` directory (created if absent).
    nx, ny, nz : int
        Grid dimensions.
    var_names : list of str
        Variable names in component order.
    wind_vals : dict
        Uniform wind values keyed by variable name.
    """
    os.makedirs(level_dir, exist_ok=True)
    nvar = len(var_names)

    # Build binary data
    fab_data = build_fab_data(nx, ny, nz, var_names, wind_vals)

    # ASCII FAB header (AMReX WriteSingleLevelPlotfile format)
    fab_header_str = (
        f"FAB ((8,(64 11 52 0 1 12 0 1023)),(8,(64 11 52 0 1 12 0 1023)))"
        f" (({0},{0},{0}) ({nx-1},{ny-1},{nz-1}) ({0},{0},{0}))\n"
    )
    fab_header_bytes = fab_header_str.encode("ascii")

    # Write Cell_D_00000
    data_file_path = os.path.join(level_dir, "Cell_D_00000")
    with open(data_file_path, "wb") as f:
        f.write(fab_header_bytes)
        f.write(fab_data)

    fab_offset = 0  # FAB header starts at byte 0 in Cell_D_00000

    # Write Cell_H (VisMF-V1 header)
    box_str = f"(({0},{0},{0}) ({nx-1},{ny-1},{nz-1}) ({0},{0},{0}))"
    cell_h_path = os.path.join(level_dir, "Cell_H")
    with open(cell_h_path, "w") as f:
        f.write("VisMF-V1\n")
        f.write(f"1\n")           # nFAB
        f.write(f"{nvar}\n")      # nComp
        f.write("0\n")            # nGhost
        f.write(f"{box_str} {nvar}  Cell_D_00000:{fab_offset}\n")


def write_plt_header(out_dir, nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax,
                     var_names):
    """Write the top-level AMReX plotfile ``Header`` file.

    Parameters
    ----------
    out_dir : str
        Root directory of the plotfile (e.g. ``plt_wind_3d/``).
    nx, ny, nz : int
        Grid cell counts.
    xmin, xmax : float
        Physical x-extent [m].
    ymin, ymax : float
        Physical y-extent [m].
    zmin, zmax : float
        Physical z-extent [m].
    var_names : list of str
        Variable names written in order to the header.
    """
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
        f.write("3\n")                          # spacedim
        f.write("0.0\n")                        # time
        f.write("0\n")                          # finest_level
        f.write(f"{xmin} {ymin} {zmin}\n")     # prob_lo
        f.write(f"{xmax} {ymax} {zmax}\n")     # prob_hi
        # No refinement ratios (finest_level == 0)
        f.write(f"{bbox_str}\n")                # level 0 bounding box
        f.write("0\n")                          # step count for level 0
        f.write(f"{dx} {dy} {dz}\n")           # cell sizes at level 0
        f.write("0\n")                          # coordinate system (Cartesian)
        f.write("0\n")                          # reserved
        f.write("Level_0/Cell\n")               # level 0 data path


def main():
    """Generate the synthetic plt wind directory used by the regression test.

    Writes ``plt_wind_3d/`` (relative to the current working directory)
    containing a uniform 3-D wind field on an 8 × 8 × 4 grid.
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
    print(f"  u={U_WIND} m/s, v={V_WIND} m/s, w={W_WIND} m/s (uniform)")
    print(f"  Domain: ({XMIN},{YMIN},{ZMIN}) to ({XMAX},{YMAX},{ZMAX})")


if __name__ == "__main__":
    main()
