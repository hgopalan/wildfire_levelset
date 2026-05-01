#!/usr/bin/env python3
"""
run_regtest.py - Regression test for tools/terrain_wind_preprocess.py

Tests all major capabilities of the unified tool using only synthetic
data (no internet access or real WRF/SRTM/LANDFIRE files required):

  1. WRF wind extraction (basic)
  2. Bounding box read from WRF file (--lat-min/max/lon-min/max ignored)
  3. --no-terrain flag (terrain step skipped; wind still produced)
  4. --time-range T1:TN (multiple wind files)
  5. --time-index N (single time step)
  6. --interpolate-wind (wind on SRTM-resolution grid)
  7. SRTM terrain values not contaminated by WRF HGT_M
  8. Landscape from local raster files (--elev-file etc.)
  9. --subsample reduces output count

Run:
  python3 regtest/terrain_wind_preprocess/run_regtest.py
  python3 -m pytest regtest/terrain_wind_preprocess/run_regtest.py -v
"""

import os
import sys
import tempfile
import unittest

import numpy as np

# -----------------------------------------------------------------------
# Path setup – add tools/ to sys.path
# -----------------------------------------------------------------------
_THIS_DIR  = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.abspath(os.path.join(_THIS_DIR, "..", ".."))
_TOOLS_DIR = os.path.join(_REPO_ROOT, "tools")
if _TOOLS_DIR not in sys.path:
    sys.path.insert(0, _TOOLS_DIR)

import terrain_wind_preprocess as twp  # noqa: E402


# -----------------------------------------------------------------------
# Synthetic WRF netCDF helper
# -----------------------------------------------------------------------

def _make_wrf_nc(path, ny=10, nx=12, nz=2, n_times=3,
                 lat_sw=37.0, lon_sw=-120.0,
                 dlat=0.01, dlon=0.01):
    """Write a minimal WRF-style netCDF4 file."""
    try:
        import netCDF4 as nc
    except ImportError:
        raise unittest.SkipTest("netCDF4 not installed – skipping WRF regtest")

    ds = nc.Dataset(path, "w", format="NETCDF4")
    ds.createDimension("Time",             n_times)
    ds.createDimension("south_north",      ny)
    ds.createDimension("west_east",        nx)
    ds.createDimension("bottom_top",       nz)
    ds.createDimension("south_north_stag", ny + 1)
    ds.createDimension("west_east_stag",   nx + 1)

    lats = np.array([[lat_sw + j * dlat for _ in range(nx)]
                     for j in range(ny)], dtype=np.float32)
    lons = np.array([[lon_sw + i * dlon for i in range(nx)]
                     for _ in range(ny)], dtype=np.float32)

    xlat  = ds.createVariable("XLAT",  "f4", ("Time", "south_north", "west_east"))
    xlong = ds.createVariable("XLONG", "f4", ("Time", "south_north", "west_east"))
    for t in range(n_times):
        xlat [t] = lats
        xlong[t] = lons

    hgt_m = ds.createVariable("HGT_M", "f4", ("Time", "south_north", "west_east"))
    for t in range(n_times):
        for j in range(ny):
            hgt_m[t, j, :] = float(j * 10)  # 0, 10, 20, … m

    u_var = ds.createVariable(
        "U", "f4", ("Time", "bottom_top", "south_north", "west_east_stag"))
    for t in range(n_times):
        for k in range(nz):
            for j in range(ny):
                for i in range(nx + 1):
                    u_var[t, k, j, i] = float(k * 10 + t * 100 + i)

    v_var = ds.createVariable(
        "V", "f4", ("Time", "bottom_top", "south_north_stag", "west_east"))
    for t in range(n_times):
        for k in range(nz):
            for j in range(ny + 1):
                for i in range(nx):
                    v_var[t, k, j, i] = float(k * 10 + t * 100 + j)

    ds.close()
    return path


def _make_geotiff(path, data, west=-120.0, south=37.0, cellsize=0.001,
                  nodata=None, epsg=4326):
    """Write a small GeoTIFF from a 2-D numpy array (requires rasterio)."""
    try:
        import rasterio
        from rasterio.transform import from_bounds
        from rasterio.crs import CRS
    except ImportError:
        raise unittest.SkipTest("rasterio not installed – skipping GeoTIFF test")

    nrows, ncols = data.shape
    east  = west  + ncols * cellsize
    north = south + nrows * cellsize
    transform = from_bounds(west, south, east, north, ncols, nrows)
    crs = CRS.from_epsg(epsg)

    with rasterio.open(
        path, "w",
        driver="GTiff",
        height=nrows, width=ncols,
        count=1,
        dtype=data.dtype,
        crs=crs,
        transform=transform,
        nodata=nodata,
    ) as ds:
        ds.write(data, 1)
    return path


def _count(path):
    with open(path) as fh:
        return sum(1 for line in fh if not line.startswith("#") and line.strip())


# -----------------------------------------------------------------------
# Regression tests
# -----------------------------------------------------------------------

class RegressionTestTerrainWindPreprocess(unittest.TestCase):

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    # -------------------------------------------------------------------
    # Test 1: Basic WRF wind extraction
    # -------------------------------------------------------------------
    def test_01_wrf_wind_extraction_basic(self):
        """Wind file is produced with 4-column UTM-metre format."""
        nc_path = os.path.join(self.tmpdir, "wrf1.nc")
        _make_wrf_nc(nc_path, ny=6, nx=8, n_times=1)
        wout = os.path.join(self.tmpdir, "wind_basic.csv")

        twp.main([
            "--wrf-file", nc_path,
            "--wind",     wout,
            "--no-terrain",
            "--no-landscape",
        ])

        self.assertTrue(os.path.isfile(wout))
        with open(wout) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.split()
                self.assertEqual(len(parts), 4,
                                 f"Expected 4 columns, got: {line!r}")
                # UTM x should be metre-scale
                self.assertGreater(abs(float(parts[0])), 1e5)

    # -------------------------------------------------------------------
    # Test 2: Bounding box read from WRF file
    # -------------------------------------------------------------------
    def test_02_bbox_from_wrf_file(self):
        """read_wrf_bbox returns lat/lon range matching synthetic WRF grid."""
        ny, nx = 5, 7
        dlat = dlon = 0.01
        lat_sw, lon_sw = 35.0, -118.0
        nc_path = os.path.join(self.tmpdir, "wrf_bbox.nc")
        _make_wrf_nc(nc_path, ny=ny, nx=nx, n_times=1,
                     lat_sw=lat_sw, lon_sw=lon_sw,
                     dlat=dlat, dlon=dlon)

        lat_min, lat_max, lon_min, lon_max = twp.read_wrf_bbox(nc_path)

        self.assertAlmostEqual(lat_min, lat_sw, places=3)
        self.assertAlmostEqual(lat_max, lat_sw + (ny - 1) * dlat, places=3)
        self.assertAlmostEqual(lon_min, lon_sw, places=3)
        self.assertAlmostEqual(lon_max, lon_sw + (nx - 1) * dlon, places=3)

    # -------------------------------------------------------------------
    # Test 3: --no-terrain flag
    # -------------------------------------------------------------------
    def test_03_no_terrain_flag(self):
        """--no-terrain: terrain file absent; wind file present and non-empty."""
        nc_path = os.path.join(self.tmpdir, "wrf3.nc")
        _make_wrf_nc(nc_path, ny=5, nx=6, n_times=1)
        tout = os.path.join(self.tmpdir, "terrain3.xyz")
        wout = os.path.join(self.tmpdir, "wind3.csv")

        twp.main([
            "--wrf-file", nc_path,
            "--terrain",  tout,
            "--wind",     wout,
            "--no-terrain",
            "--no-landscape",
        ])

        self.assertFalse(os.path.isfile(tout),
                         "Terrain file must NOT be created with --no-terrain")
        self.assertTrue(os.path.isfile(wout))
        self.assertGreater(_count(wout), 0)

    # -------------------------------------------------------------------
    # Test 4: --time-range
    # -------------------------------------------------------------------
    def test_04_time_range_multiple_files(self):
        """--time-range 0:2 produces wind_t0.csv, wind_t1.csv, wind_t2.csv."""
        nc_path = os.path.join(self.tmpdir, "wrf4.nc")
        _make_wrf_nc(nc_path, ny=5, nx=6, n_times=3)
        base_wout = os.path.join(self.tmpdir, "wind4.csv")

        twp.main([
            "--wrf-file",   nc_path,
            "--wind",       base_wout,
            "--time-range", "0:2",
            "--no-terrain",
            "--no-landscape",
        ])

        for t in range(3):
            expected_path = twp._wind_output_path(base_wout, t)
            self.assertTrue(os.path.isfile(expected_path),
                            f"Expected time-step file not found: {expected_path}")
            self.assertGreater(_count(expected_path), 0)

    def test_04b_time_range_values_differ(self):
        """Wind values differ across time steps (different t offsets)."""
        ny, nx = 4, 5
        nc_path = os.path.join(self.tmpdir, "wrf4b.nc")
        _make_wrf_nc(nc_path, ny=ny, nx=nx, n_times=2, nz=1)
        base_wout = os.path.join(self.tmpdir, "wind4b.csv")

        twp.main([
            "--wrf-file",   nc_path,
            "--wind",       base_wout,
            "--time-range", "0:1",
            "--no-terrain",
            "--no-landscape",
        ])

        def _read_u(path):
            vals = []
            with open(path) as fh:
                for line in fh:
                    if line.startswith("#") or not line.strip():
                        continue
                    vals.append(float(line.split()[2]))
            return np.array(vals)

        u0 = _read_u(twp._wind_output_path(base_wout, 0))
        u1 = _read_u(twp._wind_output_path(base_wout, 1))
        # At t=1 the offset should be +100
        np.testing.assert_allclose(
            u1 - u0, np.full_like(u0, 100.0), atol=1e-3,
            err_msg="U difference between t=0 and t=1 should be 100",
        )

    # -------------------------------------------------------------------
    # Test 5: --time-index (single step)
    # -------------------------------------------------------------------
    def test_05_time_index_single_step(self):
        """--time-index 1 writes a single unsuffixed wind file."""
        nc_path = os.path.join(self.tmpdir, "wrf5.nc")
        _make_wrf_nc(nc_path, ny=5, nx=6, n_times=2)
        wout = os.path.join(self.tmpdir, "wind5.csv")

        twp.main([
            "--wrf-file",   nc_path,
            "--wind",       wout,
            "--time-index", "1",
            "--no-terrain",
            "--no-landscape",
        ])

        self.assertTrue(os.path.isfile(wout))
        # Suffixed version must NOT exist
        self.assertFalse(os.path.isfile(twp._wind_output_path(wout, 1)))

    def test_05b_time_index_u_values(self):
        """U values at t=1 include the t-offset of 100."""
        ny, nx = 3, 4
        nc_path = os.path.join(self.tmpdir, "wrf5b.nc")
        _make_wrf_nc(nc_path, ny=ny, nx=nx, n_times=2, nz=1)
        wout = os.path.join(self.tmpdir, "wind5b.csv")

        twp.main([
            "--wrf-file",   nc_path,
            "--wind",       wout,
            "--time-index", "1",
            "--no-terrain",
            "--no-landscape",
        ])

        u_vals = []
        with open(wout) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                u_vals.append(float(line.split()[2]))

        # At t=1, k=0: u_stag[j,i]=100+i → u_mass[j,i]=100+i+0.5
        expected = [100.0 + float(i) + 0.5
                    for _ in range(ny) for i in range(nx)]
        np.testing.assert_allclose(u_vals, expected, atol=1e-3)

    # -------------------------------------------------------------------
    # Test 6: --interpolate-wind (mocked SRTM grid)
    # -------------------------------------------------------------------
    def test_06_interpolate_wind(self):
        """--interpolate-wind: wind output has same number of rows as SRTM grid."""
        try:
            from scipy.interpolate import griddata  # noqa: F401
        except ImportError:
            self.skipTest("scipy not installed – skipping interpolate-wind test")

        ny_wrf, nx_wrf = 6, 8
        nc_path = os.path.join(self.tmpdir, "wrf6.nc")
        _make_wrf_nc(nc_path, ny=ny_wrf, nx=nx_wrf, n_times=1)
        tout = os.path.join(self.tmpdir, "terrain6.xyz")
        wout = os.path.join(self.tmpdir, "wind6.csv")

        # Build a finer synthetic SRTM grid matching the WRF domain
        lat_min, lat_max, lon_min, lon_max = twp.read_wrf_bbox(nc_path)
        srtm_ny, srtm_nx = 20, 24
        srtm_lats = np.linspace(lat_min, lat_max, srtm_ny)
        srtm_lons = np.linspace(lon_min, lon_max, srtm_nx)
        lon_2d, lat_2d = np.meshgrid(srtm_lons, srtm_lats)
        srtm_x, srtm_y = twp._latlon_to_utm(lat_2d, lon_2d)

        # Patch create_terrain_xyz_return_grid to return our synthetic grid
        original_fn = twp.create_terrain_xyz_return_grid

        def _fake_create(*args, **kwargs):
            twp.write_terrain_xyz(srtm_x, srtm_y,
                                  np.ones_like(srtm_x) * 150.0, tout)
            return srtm_x, srtm_y, np.ones_like(srtm_x) * 150.0

        twp.create_terrain_xyz_return_grid = _fake_create
        try:
            twp.main([
                "--wrf-file",       nc_path,
                "--terrain",        tout,
                "--wind",           wout,
                "--interpolate-wind",
                "--no-landscape",
            ])
        finally:
            twp.create_terrain_xyz_return_grid = original_fn

        self.assertTrue(os.path.isfile(wout))
        n_rows = _count(wout)
        self.assertEqual(n_rows, srtm_ny * srtm_nx,
                         f"Expected {srtm_ny * srtm_nx} rows, got {n_rows}")

    # -------------------------------------------------------------------
    # Test 7: SRTM terrain not contaminated by WRF HGT_M
    # -------------------------------------------------------------------
    def test_07_srtm_terrain_not_hgt_m(self):
        """When SRTM terrain is used, output Z differs from WRF HGT_M values.

        We mock create_terrain_xyz to avoid a real SRTM download and verify
        that the terrain file is populated from the mocked SRTM values (500 m),
        not from WRF HGT_M (0..40 m range for a 5×6 grid with j*10 pattern).
        """
        ny_wrf, nx_wrf = 5, 6
        nc_path = os.path.join(self.tmpdir, "wrf7.nc")
        _make_wrf_nc(nc_path, ny=ny_wrf, nx=nx_wrf, n_times=1)
        tout = os.path.join(self.tmpdir, "terrain7.xyz")

        # Build a synthetic SRTM grid with a fixed elevation clearly different
        # from WRF HGT_M values (which span 0..40 m for this grid)
        lat_min, lat_max, lon_min, lon_max = twp.read_wrf_bbox(nc_path)
        srtm_ny, srtm_nx = 8, 10
        lats = np.linspace(lat_min, lat_max, srtm_ny)
        lons = np.linspace(lon_min, lon_max, srtm_nx)
        lon_2d, lat_2d = np.meshgrid(lons, lats)
        srtm_x, srtm_y = twp._latlon_to_utm(lat_2d, lon_2d)
        srtm_z_const = 500.0  # clearly not WRF HGT_M values (0..40 m)

        # Patch create_terrain_xyz (called when --interpolate-wind is NOT set)
        # so no real SRTM download is triggered.
        original_fn = twp.create_terrain_xyz

        def _fake_create_xyz(output_path, **kwargs):
            z = np.full_like(srtm_x, srtm_z_const)
            twp.write_terrain_xyz(srtm_x, srtm_y, z, output_path)

        twp.create_terrain_xyz = _fake_create_xyz
        try:
            twp.main([
                "--wrf-file",  nc_path,
                "--terrain",   tout,
                "--no-landscape",
                "--no-wind",
            ])
        finally:
            twp.create_terrain_xyz = original_fn

        # All Z values in the terrain file should be 500 m (from SRTM mock),
        # not 0..40 m (WRF HGT_M range)
        self.assertTrue(os.path.isfile(tout))
        with open(tout) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                z = float(line.split()[2])
                self.assertAlmostEqual(z, srtm_z_const, places=1,
                                       msg=f"Z={z} is not the SRTM mock value")

    # -------------------------------------------------------------------
    # Test 8: Landscape from local raster files
    # -------------------------------------------------------------------
    def test_08_landscape_local_files(self):
        """create_landscape_from_files writes a valid 6-column LCP file."""
        try:
            import rasterio  # noqa: F401
        except ImportError:
            self.skipTest("rasterio not installed – skipping landscape test")

        nrows, ncols = 8, 10
        elev   = (100.0 + np.arange(nrows * ncols)
                  .reshape(nrows, ncols)).astype(np.float32)
        slope  = np.full((nrows, ncols), 15.0, dtype=np.float32)
        aspect = np.full((nrows, ncols), 225.0, dtype=np.float32)
        fuel   = np.full((nrows, ncols), 4, dtype=np.int16)

        e_path = _make_geotiff(os.path.join(self.tmpdir, "elev.tif"),   elev)
        s_path = _make_geotiff(os.path.join(self.tmpdir, "slope.tif"),  slope)
        a_path = _make_geotiff(os.path.join(self.tmpdir, "aspect.tif"), aspect)
        f_path = _make_geotiff(os.path.join(self.tmpdir, "fuel.tif"),   fuel)
        out    = os.path.join(self.tmpdir, "landscape8.lcp")

        twp.create_landscape_from_files(
            out, e_path, s_path, a_path, f_path,
            project_utm=False,
        )

        self.assertTrue(os.path.isfile(out))
        with open(out) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.split()
                self.assertEqual(len(parts), 6,
                                 f"Expected 6 columns, got: {line!r}")

    # -------------------------------------------------------------------
    # Test 9: --subsample reduces output count
    # -------------------------------------------------------------------
    def test_09_subsample(self):
        """--subsample 2 produces fewer wind rows than --subsample 1."""
        ny, nx = 12, 14
        nc_path = os.path.join(self.tmpdir, "wrf9.nc")
        _make_wrf_nc(nc_path, ny=ny, nx=nx, n_times=1)
        wout_full = os.path.join(self.tmpdir, "wind9_full.csv")
        wout_sub  = os.path.join(self.tmpdir, "wind9_sub.csv")

        twp.main([
            "--wrf-file", nc_path, "--wind", wout_full,
            "--no-terrain", "--no-landscape",
        ])
        twp.main([
            "--wrf-file", nc_path, "--wind", wout_sub,
            "--no-terrain", "--no-landscape",
            "--subsample", "2",
        ])

        n_full = _count(wout_full)
        n_sub  = _count(wout_sub)
        self.assertGreater(n_full, 0)
        self.assertLess(n_sub, n_full,
                        f"Expected fewer rows with --subsample 2: "
                        f"full={n_full}, sub={n_sub}")

    # -------------------------------------------------------------------
    # Test 10: _parse_time_range edge cases
    # -------------------------------------------------------------------
    def test_10_parse_time_range(self):
        self.assertEqual(twp._parse_time_range("0:0"), [0])
        self.assertEqual(twp._parse_time_range("2:5"), [2, 3, 4, 5])
        with self.assertRaises(ValueError):
            twp._parse_time_range("5:2")
        with self.assertRaises(ValueError):
            twp._parse_time_range("abc")


# -----------------------------------------------------------------------
# Entry point
# -----------------------------------------------------------------------

if __name__ == "__main__":
    unittest.main(verbosity=2)
