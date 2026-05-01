#!/usr/bin/env python3
"""
test_terrain_wind_preprocess.py - Tests for tools/terrain_wind_preprocess.py

Tests:
  1.  _parse_time_range: valid and invalid inputs
  2.  _wind_output_path: time-stamped filename generation
  3.  read_wrf_bbox: bounding box extraction from a synthetic WRF netCDF
  4.  _destagger_u / _destagger_v: staggered-to-mass-point averaging
  5.  _latlon_to_utm: shape preservation and metre-scale values
  6.  extract_wrf_wind: full pipeline (synthetic WRF netCDF → UTM wind arrays)
  7.  interpolate_wind_to_grid: bilinear interpolation onto a finer grid
  8.  convert_wrf (legacy): terrain + wind files written from WRF HGT_M
  9.  --no-terrain flag: terrain step skipped, wind still extracted
  10. --interpolate-wind flag: wind output on SRTM grid (mocked SRTM download)
  11. --time-range flag: multiple wind files written
  12. CLI invocation: basic run, missing file, subsample, time-range

Run with:
  python3 tools/tests/test_terrain_wind_preprocess.py
  python3 -m pytest tools/tests/test_terrain_wind_preprocess.py -v
"""

import os
import subprocess
import sys
import tempfile
import unittest

import numpy as np

# ---------------------------------------------------------------------------
# Path setup
# ---------------------------------------------------------------------------

_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
_TOOLS_DIR = os.path.join(_REPO_ROOT, "tools")
sys.path.insert(0, _TOOLS_DIR)

import terrain_wind_preprocess as twp  # noqa: E402


# ---------------------------------------------------------------------------
# Synthetic WRF netCDF helper
# ---------------------------------------------------------------------------

def _make_wrf_nc(path, ny=10, nx=12, nz=3, n_times=2,
                 lat_sw=37.0, lon_sw=-120.0,
                 dlat=0.01, dlon=0.01):
    """Create a minimal WRF-style netCDF file for testing."""
    try:
        import netCDF4 as nc
    except ImportError:
        raise unittest.SkipTest("netCDF4 not installed – skipping WRF tests")

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
            hgt_m[t, j, :] = float(j * 10)

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


# ---------------------------------------------------------------------------
# Helper utilities
# ---------------------------------------------------------------------------

def _count_data_lines(path):
    with open(path) as fh:
        return sum(1 for l in fh if not l.startswith("#") and l.strip())


# ===========================================================================
# 1. _parse_time_range
# ===========================================================================

class TestParseTimeRange(unittest.TestCase):

    def test_single_step(self):
        self.assertEqual(twp._parse_time_range("3:3"), [3])

    def test_multi_step(self):
        self.assertEqual(twp._parse_time_range("0:4"), [0, 1, 2, 3, 4])

    def test_start_one(self):
        self.assertEqual(twp._parse_time_range("1:3"), [1, 2, 3])

    def test_zero_based(self):
        self.assertEqual(twp._parse_time_range("0:0"), [0])

    def test_invalid_format_no_colon(self):
        with self.assertRaises(ValueError):
            twp._parse_time_range("42")

    def test_invalid_format_two_colons(self):
        with self.assertRaises(ValueError):
            twp._parse_time_range("0:2:4")

    def test_invalid_non_integer(self):
        with self.assertRaises(ValueError):
            twp._parse_time_range("a:b")

    def test_t1_greater_than_tn(self):
        with self.assertRaises(ValueError):
            twp._parse_time_range("5:3")


# ===========================================================================
# 2. _wind_output_path
# ===========================================================================

class TestWindOutputPath(unittest.TestCase):

    def test_csv_extension(self):
        self.assertEqual(twp._wind_output_path("wind.csv", 0), "wind_t0.csv")

    def test_nonzero_index(self):
        self.assertEqual(twp._wind_output_path("wind.csv", 7), "wind_t7.csv")

    def test_path_with_directory(self):
        result = twp._wind_output_path("/tmp/output/wind.csv", 2)
        self.assertEqual(result, "/tmp/output/wind_t2.csv")

    def test_no_extension(self):
        self.assertEqual(twp._wind_output_path("wind", 1), "wind_t1")


# ===========================================================================
# 3. read_wrf_bbox
# ===========================================================================

class TestReadWrfBbox(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def test_bbox_extents(self):
        """Bbox should cover the expected lat/lon range."""
        ny, nx = 5, 6
        dlat = dlon = 0.01
        lat_sw, lon_sw = 37.0, -120.0
        path = os.path.join(self.tmpdir, "wrf_bbox.nc")
        _make_wrf_nc(path, ny=ny, nx=nx, n_times=1,
                     lat_sw=lat_sw, lon_sw=lon_sw, dlat=dlat, dlon=dlon)

        lat_min, lat_max, lon_min, lon_max = twp.read_wrf_bbox(path)

        self.assertAlmostEqual(lat_min, lat_sw, places=3)
        self.assertAlmostEqual(lat_max, lat_sw + (ny - 1) * dlat, places=3)
        self.assertAlmostEqual(lon_min, lon_sw, places=3)
        self.assertAlmostEqual(lon_max, lon_sw + (nx - 1) * dlon, places=3)

    def test_lat_min_less_than_max(self):
        path = os.path.join(self.tmpdir, "wrf_order.nc")
        _make_wrf_nc(path, n_times=1)
        lat_min, lat_max, lon_min, lon_max = twp.read_wrf_bbox(path)
        self.assertLess(lat_min, lat_max)
        self.assertLess(lon_min, lon_max)


# ===========================================================================
# 4. Destagger helpers
# ===========================================================================

class TestDestaggerU(unittest.TestCase):
    def test_shape(self):
        U_stag = np.ones((10, 13))
        u_mass = twp._destagger_u(U_stag)
        self.assertEqual(u_mass.shape, (10, 12))

    def test_average(self):
        U_stag = np.array([[0.0, 2.0]])
        np.testing.assert_allclose(twp._destagger_u(U_stag), np.array([[1.0]]))

    def test_known_values(self):
        U_stag = np.array([[1.0, 3.0, 5.0], [2.0, 4.0, 6.0]])
        expected = np.array([[2.0, 4.0], [3.0, 5.0]])
        np.testing.assert_allclose(twp._destagger_u(U_stag), expected)


class TestDestaggerV(unittest.TestCase):
    def test_shape(self):
        V_stag = np.ones((11, 12))
        v_mass = twp._destagger_v(V_stag)
        self.assertEqual(v_mass.shape, (10, 12))

    def test_average(self):
        V_stag = np.array([[0.0], [2.0]])
        np.testing.assert_allclose(twp._destagger_v(V_stag), np.array([[1.0]]))

    def test_known_values(self):
        V_stag = np.array([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]])
        expected = np.array([[2.0, 3.0], [4.0, 5.0]])
        np.testing.assert_allclose(twp._destagger_v(V_stag), expected)


# ===========================================================================
# 5. _latlon_to_utm
# ===========================================================================

class TestLatLonToUtm(unittest.TestCase):
    def test_shape_preserved(self):
        lats = np.array([[37.0, 37.1], [37.2, 37.3]])
        lons = np.array([[-120.0, -119.9], [-120.0, -119.9]])
        x, y = twp._latlon_to_utm(lats, lons)
        self.assertEqual(x.shape, lats.shape)
        self.assertEqual(y.shape, lats.shape)

    def test_utm_in_metres(self):
        lats = np.array([[37.0]])
        lons = np.array([[-120.0]])
        x, y = twp._latlon_to_utm(lats, lons)
        self.assertGreater(abs(float(x[0, 0])), 1e5)
        self.assertGreater(abs(float(y[0, 0])), 1e5)

    def test_east_increases_x(self):
        lats = np.array([[37.0, 37.0]])
        lons = np.array([[-120.5, -119.5]])
        x, _ = twp._latlon_to_utm(lats, lons)
        self.assertLess(float(x[0, 0]), float(x[0, 1]))

    def test_north_increases_y(self):
        lats = np.array([[36.0], [38.0]])
        lons = np.array([[-120.0], [-120.0]])
        _, y = twp._latlon_to_utm(lats, lons)
        self.assertLess(float(y[0, 0]), float(y[1, 0]))


# ===========================================================================
# 6. extract_wrf_wind
# ===========================================================================

class TestExtractWrfWind(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def _make_nc(self, **kwargs):
        path = os.path.join(self.tmpdir, "wrfout.nc")
        return _make_wrf_nc(path, **kwargs)

    def test_output_shapes(self):
        ny, nx = 8, 10
        nc_path = self._make_nc(ny=ny, nx=nx, n_times=1)
        x, y, u, v = twp.extract_wrf_wind(nc_path)
        self.assertEqual(x.shape, (ny, nx))
        self.assertEqual(u.shape, (ny, nx))
        self.assertEqual(v.shape, (ny, nx))

    def test_utm_coordinates_in_metres(self):
        nc_path = self._make_nc()
        x, y, _, _ = twp.extract_wrf_wind(nc_path)
        self.assertGreater(abs(float(x.ravel()[0])), 1e5)

    def test_destagger_u_values(self):
        """u_mass[j,i] = i + 0.5 for t=0, k=0."""
        ny, nx = 4, 6
        nc_path = self._make_nc(ny=ny, nx=nx, n_times=1, nz=1)
        _, _, u, _ = twp.extract_wrf_wind(nc_path, time_index=0, level=0)
        expected = np.array([[float(i) + 0.5 for i in range(nx)]
                              for _ in range(ny)])
        np.testing.assert_allclose(u, expected, rtol=1e-5)

    def test_destagger_v_values(self):
        """v_mass[j,i] = j + 0.5 for t=0, k=0."""
        ny, nx = 4, 6
        nc_path = self._make_nc(ny=ny, nx=nx, n_times=1, nz=1)
        _, _, _, v = twp.extract_wrf_wind(nc_path, time_index=0, level=0)
        expected = np.array([[float(j) + 0.5 for _ in range(nx)]
                              for j in range(ny)])
        np.testing.assert_allclose(v, expected, rtol=1e-5)

    def test_time_index_selection(self):
        """u_mass offset changes with time index."""
        ny, nx = 3, 4
        nc_path = self._make_nc(ny=ny, nx=nx, n_times=2, nz=1)
        _, _, u0, _ = twp.extract_wrf_wind(nc_path, time_index=0)
        _, _, u1, _ = twp.extract_wrf_wind(nc_path, time_index=1)
        # t=1 offset is +100 relative to t=0
        np.testing.assert_allclose(u1 - u0,
                                   np.full((ny, nx), 100.0), atol=1e-4)

    def test_level_selection(self):
        """u_mass offset changes with vertical level."""
        ny, nx, nz = 3, 4, 3
        nc_path = self._make_nc(ny=ny, nx=nx, n_times=1, nz=nz)
        _, _, u0, _ = twp.extract_wrf_wind(nc_path, level=0)
        _, _, u1, _ = twp.extract_wrf_wind(nc_path, level=1)
        np.testing.assert_allclose(u1 - u0,
                                   np.full((ny, nx), 10.0), atol=1e-4)

    def test_subsample_reduces_size(self):
        ny, nx = 10, 12
        nc_path = self._make_nc(ny=ny, nx=nx)
        x_full, _, _, _ = twp.extract_wrf_wind(nc_path, subsample=1)
        x_sub,  _, _, _ = twp.extract_wrf_wind(nc_path, subsample=2)
        self.assertLess(x_sub.size, x_full.size)


# ===========================================================================
# 7. interpolate_wind_to_grid
# ===========================================================================

class TestInterpolateWindToGrid(unittest.TestCase):

    def _make_uniform_grid(self, ny, nx, x0=0.0, y0=0.0, dx=1000.0, dy=1000.0):
        xs = np.array([[x0 + i * dx for i in range(nx)] for _ in range(ny)])
        ys = np.array([[y0 + j * dy for _ in range(nx)] for j in range(ny)])
        return xs, ys

    def test_shape_of_output(self):
        try:
            from scipy.interpolate import griddata  # noqa: F401
        except ImportError:
            self.skipTest("scipy not installed")

        wrf_x, wrf_y = self._make_uniform_grid(5, 6)
        u = np.ones((5, 6))
        v = np.zeros((5, 6))
        tgt_x, tgt_y = self._make_uniform_grid(10, 12, dx=500.0, dy=500.0)
        u_i, v_i = twp.interpolate_wind_to_grid(wrf_x, wrf_y, u, v, tgt_x, tgt_y)
        self.assertEqual(u_i.shape, tgt_x.shape)
        self.assertEqual(v_i.shape, tgt_y.shape)

    def test_constant_field_interpolated_exactly(self):
        """Interpolating a constant field should reproduce that constant."""
        try:
            from scipy.interpolate import griddata  # noqa: F401
        except ImportError:
            self.skipTest("scipy not installed")

        wrf_x, wrf_y = self._make_uniform_grid(8, 8, dx=1000.0, dy=1000.0)
        u = np.full((8, 8), 3.5)
        v = np.full((8, 8), -1.2)
        # Target grid that lies entirely inside the WRF domain
        tgt_x, tgt_y = self._make_uniform_grid(
            5, 5, x0=1500.0, y0=1500.0, dx=700.0, dy=700.0
        )
        u_i, v_i = twp.interpolate_wind_to_grid(wrf_x, wrf_y, u, v, tgt_x, tgt_y)
        np.testing.assert_allclose(u_i, 3.5, atol=1e-10)
        np.testing.assert_allclose(v_i, -1.2, atol=1e-10)


# ===========================================================================
# 8. convert_wrf (legacy)
# ===========================================================================

class TestConvertWrf(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def _make_nc(self, **kwargs):
        path = os.path.join(self.tmpdir, "wrfout.nc")
        return _make_wrf_nc(path, **kwargs)

    def test_output_files_created(self):
        nc = self._make_nc()
        tout = os.path.join(self.tmpdir, "t.csv")
        wout = os.path.join(self.tmpdir, "w.csv")
        twp.convert_wrf(nc, tout, wout)
        self.assertTrue(os.path.isfile(tout))
        self.assertTrue(os.path.isfile(wout))

    def test_terrain_row_count(self):
        ny, nx = 8, 10
        nc = self._make_nc(ny=ny, nx=nx)
        tout = os.path.join(self.tmpdir, "t_rows.csv")
        wout = os.path.join(self.tmpdir, "w_rows.csv")
        twp.convert_wrf(nc, tout, wout)
        self.assertEqual(_count_data_lines(tout), ny * nx)

    def test_terrain_three_columns(self):
        nc = self._make_nc()
        tout = os.path.join(self.tmpdir, "t_3col.csv")
        wout = os.path.join(self.tmpdir, "w_3col.csv")
        twp.convert_wrf(nc, tout, wout)
        with open(tout) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                self.assertEqual(len(line.split()), 3)

    def test_wind_four_columns(self):
        nc = self._make_nc()
        tout = os.path.join(self.tmpdir, "t_4col.csv")
        wout = os.path.join(self.tmpdir, "w_4col.csv")
        twp.convert_wrf(nc, tout, wout)
        with open(wout) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                self.assertEqual(len(line.split()), 4)

    def test_terrain_z_values(self):
        """Z values match synthetic HGT_M pattern (j * 10)."""
        ny, nx = 5, 6
        nc = self._make_nc(ny=ny, nx=nx, n_times=1)
        tout = os.path.join(self.tmpdir, "t_z.csv")
        wout = os.path.join(self.tmpdir, "w_z.csv")
        twp.convert_wrf(nc, tout, wout)
        z_vals = []
        with open(tout) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                z_vals.append(float(line.split()[2]))
        expected = [float(j * 10) for j in range(ny) for _ in range(nx)]
        np.testing.assert_allclose(z_vals, expected, atol=1e-3)


# ===========================================================================
# 9. main() - --no-terrain: terrain step skipped, wind still extracted
# ===========================================================================

class TestMainNoTerrain(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def _make_nc(self, **kwargs):
        path = os.path.join(self.tmpdir, "wrf_no_terrain.nc")
        return _make_wrf_nc(path, **kwargs)

    def test_no_terrain_skips_terrain_file(self):
        """--no-terrain: terrain file must NOT be created."""
        nc = self._make_nc()
        tout = os.path.join(self.tmpdir, "terrain.xyz")
        wout = os.path.join(self.tmpdir, "wind.csv")

        twp.main([
            "--wrf-file", nc,
            "--terrain", tout,
            "--wind",    wout,
            "--no-terrain",
            "--no-landscape",
        ])

        self.assertFalse(os.path.isfile(tout),
                         "terrain file should NOT exist with --no-terrain")

    def test_no_terrain_still_writes_wind(self):
        """--no-terrain: wind file should still be created."""
        nc = self._make_nc()
        wout = os.path.join(self.tmpdir, "wind_no_terrain.csv")

        twp.main([
            "--wrf-file", nc,
            "--wind",    wout,
            "--no-terrain",
            "--no-landscape",
        ])

        self.assertTrue(os.path.isfile(wout))
        self.assertGreater(_count_data_lines(wout), 0)

    def test_wind_file_has_four_columns(self):
        nc = self._make_nc()
        wout = os.path.join(self.tmpdir, "wind_4col.csv")
        twp.main([
            "--wrf-file", nc,
            "--wind",    wout,
            "--no-terrain",
            "--no-landscape",
        ])
        with open(wout) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                self.assertEqual(len(line.split()), 4)


# ===========================================================================
# 10. main() - --time-range: multiple wind files written
# ===========================================================================

class TestMainTimeRange(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def _make_nc(self, **kwargs):
        path = os.path.join(self.tmpdir, "wrf_trange.nc")
        return _make_wrf_nc(path, **kwargs)

    def test_multiple_wind_files_created(self):
        """--time-range 0:2 should produce wind_t0.csv, wind_t1.csv, wind_t2.csv."""
        nc = self._make_nc(n_times=3)
        base_wout = os.path.join(self.tmpdir, "wind.csv")

        twp.main([
            "--wrf-file",   nc,
            "--wind",       base_wout,
            "--time-range", "0:2",
            "--no-terrain",
            "--no-landscape",
        ])

        for t in range(3):
            expected = twp._wind_output_path(base_wout, t)
            self.assertTrue(os.path.isfile(expected),
                            f"Expected wind file not found: {expected}")

    def test_single_time_range_no_suffix(self):
        """A single-step --time-range should NOT produce a suffixed file."""
        nc = self._make_nc(n_times=2)
        base_wout = os.path.join(self.tmpdir, "wind_single.csv")

        twp.main([
            "--wrf-file",   nc,
            "--wind",       base_wout,
            "--time-range", "1:1",
            "--no-terrain",
            "--no-landscape",
        ])
        # single time-step: no _t1 suffix expected
        self.assertTrue(os.path.isfile(base_wout))

    def test_wind_values_differ_between_time_steps(self):
        """Different time steps should yield different wind values."""
        ny, nx = 4, 5
        nc = self._make_nc(ny=ny, nx=nx, n_times=2, nz=1)
        base_wout = os.path.join(self.tmpdir, "wind_diff.csv")

        twp.main([
            "--wrf-file",   nc,
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
        self.assertFalse(np.allclose(u0, u1),
                         "Wind values should differ across time steps")


# ===========================================================================
# 11. main() - --interpolate-wind (mocked SRTM)
# ===========================================================================

class TestMainInterpolateWind(unittest.TestCase):
    """Test --interpolate-wind using a monkey-patched SRTM download."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def _make_nc(self, ny=8, nx=10):
        path = os.path.join(self.tmpdir, "wrf_interp.nc")
        return _make_wrf_nc(path, ny=ny, nx=nx, n_times=1)

    def _mock_srtm_download_and_tiff(self, nc_path, tout, wout,
                                     srtm_ny=20, srtm_nx=25):
        """Run with a patched create_terrain_xyz_return_grid that returns a
        synthetic fine-resolution SRTM grid without hitting the network."""
        nc_lat_min, nc_lat_max, nc_lon_min, nc_lon_max = twp.read_wrf_bbox(nc_path)

        # Build synthetic SRTM grid
        srtm_lats = np.linspace(nc_lat_min, nc_lat_max, srtm_ny)
        srtm_lons = np.linspace(nc_lon_min, nc_lon_max, srtm_nx)
        lon_2d, lat_2d = np.meshgrid(srtm_lons, srtm_lats)
        srtm_x, srtm_y = twp._latlon_to_utm(lat_2d, lon_2d)
        srtm_z = np.ones_like(srtm_x) * 100.0

        # Write terrain file
        twp.write_terrain_xyz(srtm_x, srtm_y, srtm_z, tout)

        return srtm_x, srtm_y, srtm_z

    def test_interpolated_wind_on_srtm_grid(self):
        """With --interpolate-wind, wind output should match SRTM grid size."""
        try:
            from scipy.interpolate import griddata  # noqa: F401
        except ImportError:
            self.skipTest("scipy not installed")

        nc_path = self._make_nc(ny=8, nx=10)
        tout = os.path.join(self.tmpdir, "terrain.xyz")
        wout = os.path.join(self.tmpdir, "wind_interp.csv")

        srtm_ny, srtm_nx = 20, 25
        srtm_x, srtm_y, _ = self._mock_srtm_download_and_tiff(
            nc_path, tout, wout, srtm_ny=srtm_ny, srtm_nx=srtm_nx
        )

        # Patch the module-level function so main() skips the real download
        original_fn = twp.create_terrain_xyz_return_grid

        def fake_create(*args, **kwargs):
            return srtm_x, srtm_y, np.ones_like(srtm_x) * 100.0

        twp.create_terrain_xyz_return_grid = fake_create
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

        # Wind file should have srtm_ny * srtm_nx rows
        n_rows = _count_data_lines(wout)
        self.assertEqual(n_rows, srtm_ny * srtm_nx)


# ===========================================================================
# 12. CLI tests
# ===========================================================================

class TestCLI(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.script = os.path.join(_TOOLS_DIR, "terrain_wind_preprocess.py")

    def _make_nc(self, **kwargs):
        path = os.path.join(self.tmpdir, "wrf_cli.nc")
        return _make_wrf_nc(path, **kwargs)

    def test_cli_wind_only(self):
        """Basic CLI run: WRF wind extraction with --no-terrain."""
        nc = self._make_nc()
        wout = os.path.join(self.tmpdir, "w_cli.csv")
        result = subprocess.run(
            [sys.executable, self.script,
             "--wrf-file", nc,
             "--wind",     wout,
             "--no-terrain",
             "--no-landscape"],
            capture_output=True, text=True,
        )
        self.assertEqual(result.returncode, 0,
                         f"CLI failed:\n{result.stdout}\n{result.stderr}")
        self.assertTrue(os.path.isfile(wout))

    def test_cli_missing_wrf_file(self):
        """Non-existent WRF file should exit non-zero."""
        result = subprocess.run(
            [sys.executable, self.script,
             "--wrf-file", "no_such_file.nc",
             "--no-terrain", "--no-landscape"],
            capture_output=True, text=True,
        )
        self.assertNotEqual(result.returncode, 0)

    def test_cli_missing_bbox_without_wrf(self):
        """No --wrf-file and no bbox args should exit non-zero."""
        result = subprocess.run(
            [sys.executable, self.script, "--no-terrain", "--no-landscape"],
            capture_output=True, text=True,
        )
        self.assertNotEqual(result.returncode, 0)

    def test_cli_subsample_flag(self):
        """--subsample reduces output rows."""
        ny, nx = 8, 10
        nc = self._make_nc(ny=ny, nx=nx)
        wout_full = os.path.join(self.tmpdir, "w_full.csv")
        wout_sub  = os.path.join(self.tmpdir, "w_sub.csv")

        subprocess.run(
            [sys.executable, self.script,
             "--wrf-file", nc, "--wind", wout_full,
             "--no-terrain", "--no-landscape"],
            check=True, capture_output=True,
        )
        subprocess.run(
            [sys.executable, self.script,
             "--wrf-file", nc, "--wind", wout_sub,
             "--no-terrain", "--no-landscape",
             "--subsample", "2"],
            check=True, capture_output=True,
        )

        self.assertLess(_count_data_lines(wout_sub),
                        _count_data_lines(wout_full))

    def test_cli_time_range(self):
        """--time-range produces separate wind files for each time step."""
        nc = self._make_nc(n_times=3)
        base_wout = os.path.join(self.tmpdir, "w_range.csv")

        result = subprocess.run(
            [sys.executable, self.script,
             "--wrf-file",   nc,
             "--wind",       base_wout,
             "--time-range", "0:2",
             "--no-terrain", "--no-landscape"],
            capture_output=True, text=True,
        )
        self.assertEqual(result.returncode, 0,
                         f"CLI failed:\n{result.stdout}\n{result.stderr}")

        for t in range(3):
            expected = twp._wind_output_path(base_wout, t)
            self.assertTrue(os.path.isfile(expected),
                            f"Missing time-step file: {expected}")

    def test_cli_invalid_time_range(self):
        """Invalid --time-range should exit non-zero."""
        nc = self._make_nc()
        result = subprocess.run(
            [sys.executable, self.script,
             "--wrf-file",   nc,
             "--time-range", "bad",
             "--no-terrain", "--no-landscape"],
            capture_output=True, text=True,
        )
        self.assertNotEqual(result.returncode, 0)

    def test_cli_time_index(self):
        """--time-index 1 writes a single wind file (not time-stamped)."""
        nc = self._make_nc(n_times=2)
        wout = os.path.join(self.tmpdir, "w_tidx.csv")

        result = subprocess.run(
            [sys.executable, self.script,
             "--wrf-file",   nc,
             "--wind",       wout,
             "--time-index", "1",
             "--no-terrain", "--no-landscape"],
            capture_output=True, text=True,
        )
        self.assertEqual(result.returncode, 0,
                         f"CLI failed:\n{result.stdout}\n{result.stderr}")
        self.assertTrue(os.path.isfile(wout))

    def test_cli_no_wind_flag(self):
        """--no-wind should not produce a wind file."""
        nc = self._make_nc()
        wout = os.path.join(self.tmpdir, "w_nowind.csv")

        result = subprocess.run(
            [sys.executable, self.script,
             "--wrf-file", nc,
             "--wind",     wout,
             "--no-terrain", "--no-landscape",
             "--no-wind"],
            capture_output=True, text=True,
        )
        self.assertEqual(result.returncode, 0,
                         f"CLI failed:\n{result.stdout}\n{result.stderr}")
        self.assertFalse(os.path.isfile(wout),
                         "wind file should NOT exist with --no-wind")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    unittest.main(verbosity=2)
