# ===============================================================
# Usage:
#
# Example basic run:
#   python srtm_to_xyz_stl.py \
#       --lat-min 40 --lat-max 40.5 \
#       --lon-min -105 --lon-max -104.5
#
# Full custom filenames:
#   python srtm_to_xyz_stl.py \
#       --lat-min 39.5 --lat-max 40.2 \
#       --lon-min -106 --lon-max -105.2 \
#       --tif region.tif \
#       --xyz region.xyz \
#       --stl region.stl
#
# All outputs (TIF, XYZ, STL) are always written to the current directory.
# ===============================================================

import elevation
import rasterio
import numpy as np
import pyvista as pv
import argparse
import os
from math import floor


# -------------------------------
#  Force output to current folder
# -------------------------------
def force_current_directory(filename):
    filename = os.path.basename(filename)
    return os.path.abspath(filename)


# -------------------------------
#  Vectorized UTM Conversion
# -------------------------------
a = 6378137.0
f = 1 / 298.257223563
k0 = 0.9996
e2 = f * (2 - f)
e4 = e2 * e2
e6 = e4 * e2
e_prime_sq = e2 / (1 - e2)


def lonlat_to_utm_vectorized(lon, lat):
    lon = np.asarray(lon)
    lat = np.asarray(lat)

    zone_number = int(floor((np.mean(lon) + 180) / 6) + 1)
    northern = np.mean(lat) >= 0

    lon0 = (zone_number - 1) * 6 - 180 + 3
    lon_rad = np.radians(lon)
    lat_rad = np.radians(lat)
    lon0_rad = np.radians(lon0)

    N = a / np.sqrt(1 - e2 * np.sin(lat_rad) ** 2)
    T = np.tan(lat_rad) ** 2
    C = e_prime_sq * np.cos(lat_rad) ** 2
    A = np.cos(lat_rad) * (lon_rad - lon0_rad)

    M = (
        a
        * (
            (1 - e2 / 4 - 3 * e4 / 64 - 5 * e6 / 256) * lat_rad
            - (3 * e2 / 8 + 3 * e4 / 32 + 45 * e6 / 1024) * np.sin(2 * lat_rad)
            + (15 * e4 / 256 + 45 * e6 / 1024) * np.sin(4 * lat_rad)
            - (35 * e6 / 3072) * np.sin(6 * lat_rad)
        )
    )

    x = (
        k0
        * N
        * (
            A
            + (1 - T + C) * A**3 / 6
            + (5 - 18 * T + T**2 + 72 * C - 58 * e_prime_sq) * A**5 / 120
        )
        + 500000.0
    )

    y = k0 * (
        M
        + N
        * np.tan(lat_rad)
        * (
            A**2 / 2
            + (5 - T + 9 * C + 4 * C**2) * A**4 / 24
            + (61 - 58 * T + T**2 + 600 * C - 330 * e_prime_sq) * A**6 / 720
        )
    )

    if not northern:
        y += 10000000.0

    return x, y, zone_number


# -------------------------------
#  Download SRTM
# -------------------------------
def download_srtm(lat_min, lat_max, lon_min, lon_max, out_tif):
    elevation.clean()
    elevation.clip(
        bounds=(lon_min, lat_min, lon_max, lat_max),
        output=out_tif,
        product="SRTM1",
    )


# -------------------------------
#  TIFF → XYZ (UTM)
# -------------------------------
def tiff_to_xyz_utm(tif_path):
    with rasterio.open(tif_path) as src:
        z = src.read(1)
        z = np.where(z < 0, 0, z)

        h, w = z.shape
        xs, ys = np.meshgrid(np.arange(w), np.arange(h))

        lon, lat = rasterio.transform.xy(src.transform, ys, xs)
        lon = np.array(lon)
        lat = np.array(lat)

        utm_x, utm_y, _ = lonlat_to_utm_vectorized(lon, lat)

        return utm_x, utm_y, z


# -------------------------------
#  Save XYZ file
# -------------------------------
def save_xyz(utm_x, utm_y, z, xyz_output):
    data = np.column_stack([utm_x.ravel(), utm_y.ravel(), z.ravel()])
    np.savetxt(xyz_output, data, fmt="%.3f %.3f %.3f")
    print(f"Saved XYZ file: {xyz_output}")


# -------------------------------
#  XYZ → STL
# -------------------------------
def xyz_to_stl(utm_x, utm_y, z, stl_output):
    x1=utm_x.flatten(order='F')
    x2=utm_y.flatten(order='F')
    x3=z.flatten(order='F')
    data=np.column_stack([x1,x2,x3])
    mesh=pv.PolyData(data)
    mesh['elevation']=data[:,2]
    # Save STL
    mesh.save("terrain.vtk")



# -------------------------------
#  CLI Arguments
# -------------------------------
def parse_args():
    p = argparse.ArgumentParser(description="Download SRTM, convert to XYZ and STL")

    p.add_argument("--lat-min", type=float, required=True)
    p.add_argument("--lat-max", type=float, required=True)
    p.add_argument("--lon-min", type=float, required=True)
    p.add_argument("--lon-max", type=float, required=True)

    p.add_argument("--tif", type=str, default="srtm_clip.tif")
    p.add_argument("--xyz", type=str, default="terrain.xyz")
    p.add_argument("--stl", type=str, default="terrain.stl")

    return p.parse_args()


# -------------------------------
#  Main
# -------------------------------
def main():
    args = parse_args()

    out_tif = force_current_directory(args.tif)
    out_xyz = force_current_directory(args.xyz)
    out_stl = force_current_directory(args.stl)

    print("Downloading SRTM...")
    download_srtm(args.lat_min, args.lat_max, args.lon_min, args.lon_max, out_tif)

    print("Converting TIFF → UTM XYZ...")
    utm_x, utm_y, z = tiff_to_xyz_utm(out_tif)

    print("Saving XYZ...")
    save_xyz(utm_x, utm_y, z, out_xyz)

    print("Converting XYZ → STL...")
    xyz_to_stl(utm_x, utm_y, z, out_stl)

    print("All done.")


if __name__ == "__main__":
    main()