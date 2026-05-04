#!/usr/bin/env python3
"""
run_regtest.py - Regression test for the split preprocessing tools
  (supersede the deprecated tools/deprecated/terrain_wind_preprocess.py):
  tools/wrf_wind_reader.py, tools/srtm_terrain_reader.py, tools/landscape_writer.py

Tests all major capabilities using only synthetic data (no internet access or
real WRF/SRTM/LANDFIRE files required):

  1. WRF wind extraction (basic) — wrf_wind_reader
  2. Bounding box read from WRF file — wrf_wind_reader
  3. --no-inputs flag (terrain step skipped; wind still produced) — wrf_wind_reader
  4. --time-range T1:TN (multiple wind files) — wrf_wind_reader
  5. --time-index N (single time step) — wrf_wind_reader
  6. --terrain-file (wind interpolated onto existing terrain grid) — wrf_wind_reader
  7. SRTM terrain values independent of WRF HGT_M — srtm_terrain_reader
  8. Landscape from local raster files (--elev-file etc.) — landscape_writer
  9. --subsample reduces output count — wrf_wind_reader

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

import wrf_wind_reader as wrf         # noqa: E402
import srtm_terrain_reader as srtm    # noqa: E402
import landscape_writer as lw         # noqa: E402


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

class RegressionTestPreprocessTools(unittest.TestCase):

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

        wrf.main([
            "--wrf-file", nc_path,
            "--wind",     wout,
            "--no-inputs",
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
        """read_wrf_bbox returns center±0.45° bounding box from synthetic WRF grid."""
        ny, nx = 5, 7
        dlat = dlon = 0.01
        lat_sw, lon_sw = 35.0, -118.0
        nc_path = os.path.join(self.tmpdir, "wrf_bbox.nc")
        _make_wrf_nc(nc_path, ny=ny, nx=nx, n_times=1,
                     lat_sw=lat_sw, lon_sw=lon_sw,
                     dlat=dlat, dlon=dlon)

        lat_min, lat_max, lon_min, lon_max = wrf.read_wrf_bbox(nc_path)

        # Expected: center ± 0.45°
        center_lat = lat_sw + (ny - 1) * dlat / 2.0
        center_lon = lon_sw + (nx - 1) * dlon / 2.0
        self.assertAlmostEqual(lat_min, center_lat - wrf._WRF_BBOX_HALF_SPAN, places=3)
        self.assertAlmostEqual(lat_max, center_lat + wrf._WRF_BBOX_HALF_SPAN, places=3)
        self.assertAlmostEqual(lon_min, center_lon - wrf._WRF_BBOX_HALF_SPAN, places=3)
        self.assertAlmostEqual(lon_max, center_lon + wrf._WRF_BBOX_HALF_SPAN, places=3)

    # -------------------------------------------------------------------
    # Test 3: wind file present even when --no-inputs is used
    # -------------------------------------------------------------------
    def test_03_wind_without_inputs(self):
        """--no-inputs: inputs.i absent; wind file present and non-empty."""
        nc_path = os.path.join(self.tmpdir, "wrf3.nc")
        _make_wrf_nc(nc_path, ny=5, nx=6, n_times=1)
        inputs_out = os.path.join(self.tmpdir, "inputs3.i")
        wout = os.path.join(self.tmpdir, "wind3.csv")

        wrf.main([
            "--wrf-file", nc_path,
            "--inputs",   inputs_out,
            "--wind",     wout,
            "--no-inputs",
        ])

        self.assertFalse(os.path.isfile(inputs_out),
                         "inputs.i must NOT be created with --no-inputs")
        self.assertTrue(os.path.isfile(wout))
        self.assertGreater(_count(wout), 0)

    # -------------------------------------------------------------------
    # Test 4: --time-range
    # -------------------------------------------------------------------
    def test_04_time_range_multiple_files(self):
        """--time-range 0:2 produces wind4.csv, wind4_1.csv, wind4_2.csv."""
        nc_path = os.path.join(self.tmpdir, "wrf4.nc")
        _make_wrf_nc(nc_path, ny=5, nx=6, n_times=3)
        base_wout = os.path.join(self.tmpdir, "wind4.csv")

        wrf.main([
            "--wrf-file",   nc_path,
            "--wind",       base_wout,
            "--time-range", "0:2",
            "--no-inputs",
        ])

        for pos in range(3):
            expected_path = wrf._wind_output_path(base_wout, pos)
            self.assertTrue(os.path.isfile(expected_path),
                            f"Expected time-step file not found: {expected_path}")
            self.assertGreater(_count(expected_path), 0)

    def test_04b_time_range_values_differ(self):
        """Wind values differ across time steps (different t offsets)."""
        ny, nx = 4, 5
        nc_path = os.path.join(self.tmpdir, "wrf4b.nc")
        _make_wrf_nc(nc_path, ny=ny, nx=nx, n_times=2, nz=1)
        base_wout = os.path.join(self.tmpdir, "wind4b.csv")

        wrf.main([
            "--wrf-file",   nc_path,
            "--wind",       base_wout,
            "--time-range", "0:1",
            "--no-inputs",
        ])

        def _read_u(path):
            vals = []
            with open(path) as fh:
                for line in fh:
                    if line.startswith("#") or not line.strip():
                        continue
                    vals.append(float(line.split()[2]))
            return np.array(vals)

        u0 = _read_u(wrf._wind_output_path(base_wout, 0))
        u1 = _read_u(wrf._wind_output_path(base_wout, 1))
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

        wrf.main([
            "--wrf-file",   nc_path,
            "--wind",       wout,
            "--time-index", "1",
            "--no-inputs",
        ])

        self.assertTrue(os.path.isfile(wout))
        # Suffixed version must NOT exist (position=1 → wind5_1.csv)
        self.assertFalse(os.path.isfile(wrf._wind_output_path(wout, 1)))

    def test_05b_time_index_u_values(self):
        """U values at t=1 include the t-offset of 100."""
        ny, nx = 3, 4
        nc_path = os.path.join(self.tmpdir, "wrf5b.nc")
        _make_wrf_nc(nc_path, ny=ny, nx=nx, n_times=2, nz=1)
        wout = os.path.join(self.tmpdir, "wind5b.csv")

        wrf.main([
            "--wrf-file",   nc_path,
            "--wind",       wout,
            "--time-index", "1",
            "--no-inputs",
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
    # Test 6: --terrain-file (wind interpolated onto terrain grid)
    # -------------------------------------------------------------------
    def test_06_interpolate_wind_terrain_file(self):
        """--terrain-file: wind output has same number of rows as terrain grid."""
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
        lat_min, lat_max, lon_min, lon_max = wrf.read_wrf_bbox(nc_path)
        srtm_ny, srtm_nx = 20, 24
        srtm_lats = np.linspace(lat_min, lat_max, srtm_ny)
        srtm_lons = np.linspace(lon_min, lon_max, srtm_nx)
        lon_2d, lat_2d = np.meshgrid(srtm_lons, srtm_lats)
        srtm_x, srtm_y = wrf._latlon_to_utm(lat_2d, lon_2d)

        # Write a synthetic terrain XYZ file
        srtm.write_terrain_xyz(srtm_x, srtm_y,
                               np.ones_like(srtm_x) * 150.0, tout)

        wrf.main([
            "--wrf-file",    nc_path,
            "--wind",        wout,
            "--terrain-file", tout,
            "--no-inputs",
        ])

        self.assertTrue(os.path.isfile(wout))
        n_rows = _count(wout)
        self.assertEqual(n_rows, srtm_ny * srtm_nx,
                         f"Expected {srtm_ny * srtm_nx} rows, got {n_rows}")

    # -------------------------------------------------------------------
    # Test 7: SRTM terrain write_terrain_xyz works correctly
    # -------------------------------------------------------------------
    def test_07_write_terrain_xyz(self):
        """write_terrain_xyz produces a 3-column XYZ file with expected values."""
        ny, nx = 5, 6
        lat_min, lat_max = 37.0, 37.05
        lon_min, lon_max = -120.0, -119.95
        lats = np.linspace(lat_min, lat_max, ny)
        lons = np.linspace(lon_min, lon_max, nx)
        lon_2d, lat_2d = np.meshgrid(lons, lats)
        x_utm, y_utm = srtm._latlon_to_utm(lat_2d, lon_2d)
        z_const = 500.0
        z = np.full_like(x_utm, z_const)

        tout = os.path.join(self.tmpdir, "terrain7.xyz")
        srtm.write_terrain_xyz(x_utm, y_utm, z, tout)

        self.assertTrue(os.path.isfile(tout))
        with open(tout) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.split()
                self.assertEqual(len(parts), 3,
                                 f"Expected 3 columns, got: {line!r}")
                self.assertAlmostEqual(float(parts[2]), z_const, places=1,
                                       msg=f"Z={parts[2]} is not {z_const}")

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

        lw.create_landscape_from_files(
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

        wrf.main([
            "--wrf-file", nc_path, "--wind", wout_full,
            "--no-inputs",
        ])
        wrf.main([
            "--wrf-file", nc_path, "--wind", wout_sub,
            "--no-inputs",
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
        self.assertEqual(wrf._parse_time_range("0:0"), [0])
        self.assertEqual(wrf._parse_time_range("2:5"), [2, 3, 4, 5])
        with self.assertRaises(ValueError):
            wrf._parse_time_range("5:2")
        with self.assertRaises(ValueError):
            wrf._parse_time_range("abc")

    # -------------------------------------------------------------------
    # Test 11: write_inputs_file produces a valid FARSITE inputs.i
    # -------------------------------------------------------------------
    def test_11_write_inputs_file(self):
        """write_inputs_file writes a FARSITE inputs.i with expected keys."""
        out = os.path.join(self.tmpdir, "inputs_test.i")
        wrf.write_inputs_file(
            output_path=out,
            wind_base_file=os.path.join(self.tmpdir, "wind.csv"),
            multi_time=True,
            wind_time_spacing=3600.0,
            final_time=7200.0,
            domain_bounds=(330000.0, 3775000.0, 331000.0, 3776000.0),
        )
        self.assertTrue(os.path.isfile(out))
        content = open(out).read()
        self.assertIn("final_time = 7200.0", content)
        self.assertIn("use_time_dependent_wind = 1", content)
        self.assertIn("wind_time_spacing = 3600.0", content)
        self.assertIn("skip_levelset = 1", content)
        self.assertIn("farsite.enable = 1", content)
        self.assertIn("spotting.enable = 0", content)
        self.assertIn("crown.enable = 0", content)
        self.assertIn("albini_spotting.enable = 0", content)
        self.assertIn("prob_lo_x = 330000.00", content)
        self.assertIn("prob_hi_x = 331000.00", content)

    # -------------------------------------------------------------------
    # Test 12: read_wrf_time_spacing reads XTIME variable
    # -------------------------------------------------------------------
    def test_12_read_wrf_time_spacing_xtime(self):
        """read_wrf_time_spacing returns correct spacing from XTIME (minutes)."""
        try:
            import netCDF4 as nc
        except ImportError:
            self.skipTest("netCDF4 not installed")

        nc_path = os.path.join(self.tmpdir, "wrf_xtime.nc")
        ds = nc.Dataset(nc_path, "w", format="NETCDF4")
        ds.createDimension("Time", 3)
        xtime_var = ds.createVariable("XTIME", "f4", ("Time",))
        xtime_var[:] = [0.0, 60.0, 120.0]  # 60-minute spacing
        ds.close()

        spacing = wrf.read_wrf_time_spacing(nc_path)
        self.assertAlmostEqual(spacing, 3600.0, places=1)

    # -------------------------------------------------------------------
    # Test 13: inputs.i auto-generated with --time-range sets final_time
    # -------------------------------------------------------------------
    def test_13_auto_inputs_time_range(self):
        """--time-range 0:2 with XTIME spacing → final_time in inputs.i."""
        try:
            import netCDF4 as nc
        except ImportError:
            self.skipTest("netCDF4 not installed")

        nc_path = os.path.join(self.tmpdir, "wrf_tr.nc")
        # Create WRF file with XTIME (60-min spacing) and 3 time steps
        ds = nc.Dataset(nc_path, "w", format="NETCDF4")
        ds.createDimension("Time", 3)
        ds.createDimension("south_north", 4)
        ds.createDimension("west_east", 5)
        ds.createDimension("bottom_top", 1)
        ds.createDimension("south_north_stag", 5)
        ds.createDimension("west_east_stag", 6)

        import numpy as np
        xlat = ds.createVariable("XLAT", "f4", ("Time", "south_north", "west_east"))
        xlong = ds.createVariable("XLONG", "f4", ("Time", "south_north", "west_east"))
        for t in range(3):
            xlat[t] = np.array([[37.0 + j * 0.01 for _ in range(5)]
                                 for j in range(4)], dtype=np.float32)
            xlong[t] = np.array([[-120.0 + i * 0.01 for i in range(5)]
                                  for _ in range(4)], dtype=np.float32)
        u_var = ds.createVariable("U", "f4",
                                  ("Time", "bottom_top", "south_north", "west_east_stag"))
        v_var = ds.createVariable("V", "f4",
                                  ("Time", "bottom_top", "south_north_stag", "west_east"))
        u_var[:] = 1.0
        v_var[:] = 0.5
        xtime = ds.createVariable("XTIME", "f4", ("Time",))
        xtime[:] = [0.0, 60.0, 120.0]
        ds.close()

        wout = os.path.join(self.tmpdir, "wind_tr.csv")
        inputs_out = os.path.join(self.tmpdir, "inputs_tr.i")

        wrf.main([
            "--wrf-file",   nc_path,
            "--wind",       wout,
            "--time-range", "0:2",
            "--inputs",     inputs_out,
        ])

        self.assertTrue(os.path.isfile(inputs_out))
        content = open(inputs_out).read()
        # final_time = (2 - 0) * 3600 = 7200 s
        self.assertIn("final_time = 7200.0", content)
        self.assertIn("use_time_dependent_wind = 1", content)
        self.assertIn("wind_time_spacing = 3600.0", content)


# -----------------------------------------------------------------------
# Entry point
# -----------------------------------------------------------------------

if __name__ == "__main__":
    unittest.main(verbosity=2)
