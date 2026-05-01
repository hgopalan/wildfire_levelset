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

    def test_csv_extension_zero(self):
        """position=0 returns base path unchanged."""
        self.assertEqual(twp._wind_output_path("wind.csv", 0), "wind.csv")

    def test_csv_extension_nonzero(self):
        """position>0 inserts _N before extension."""
        self.assertEqual(twp._wind_output_path("wind.csv", 7), "wind_7.csv")

    def test_path_with_directory(self):
        result = twp._wind_output_path("/tmp/output/wind.csv", 2)
        self.assertEqual(result, "/tmp/output/wind_2.csv")

    def test_no_extension(self):
        self.assertEqual(twp._wind_output_path("wind", 1), "wind_1")


# ===========================================================================
# 3. read_wrf_bbox
# ===========================================================================

class TestReadWrfBbox(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def test_bbox_is_center_plus_minus_045(self):
        """Bbox should be centre ±0.45° regardless of domain extent."""
        ny, nx = 5, 6
        dlat = dlon = 0.01
        lat_sw, lon_sw = 37.0, -120.0
        path = os.path.join(self.tmpdir, "wrf_bbox.nc")
        _make_wrf_nc(path, ny=ny, nx=nx, n_times=1,
                     lat_sw=lat_sw, lon_sw=lon_sw, dlat=dlat, dlon=dlon)

        lat_min, lat_max, lon_min, lon_max = twp.read_wrf_bbox(path)

        center_lat = lat_sw + (ny - 1) * dlat / 2.0
        center_lon = lon_sw + (nx - 1) * dlon / 2.0
        self.assertAlmostEqual(lat_min, center_lat - 0.45, places=5)
        self.assertAlmostEqual(lat_max, center_lat + 0.45, places=5)
        self.assertAlmostEqual(lon_min, center_lon - 0.45, places=5)
        self.assertAlmostEqual(lon_max, center_lon + 0.45, places=5)

    def test_lat_min_less_than_max(self):
        path = os.path.join(self.tmpdir, "wrf_order.nc")
        _make_wrf_nc(path, n_times=1)
        lat_min, lat_max, lon_min, lon_max = twp.read_wrf_bbox(path)
        self.assertLess(lat_min, lat_max)
        self.assertLess(lon_min, lon_max)

    def test_span_is_0_9_degrees(self):
        """Each axis should span exactly 0.9° (2 × 0.45°)."""
        path = os.path.join(self.tmpdir, "wrf_span.nc")
        _make_wrf_nc(path, ny=10, nx=12, n_times=1)
        lat_min, lat_max, lon_min, lon_max = twp.read_wrf_bbox(path)
        self.assertAlmostEqual(lat_max - lat_min, 0.9, places=10)
        self.assertAlmostEqual(lon_max - lon_min, 0.9, places=10)


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

    def test_latlon_clip_reduces_size(self):
        """Providing lat/lon bounds should return fewer grid points."""
        ny, nx = 20, 20
        dlat = dlon = 0.05
        lat_sw, lon_sw = 37.0, -120.0
        nc_path = os.path.join(self.tmpdir, "wrf_clip.nc")
        _make_wrf_nc(nc_path, ny=ny, nx=nx, n_times=1,
                     lat_sw=lat_sw, lon_sw=lon_sw, dlat=dlat, dlon=dlon)

        # Clip to a sub-region that excludes some rows and columns
        center_lat = lat_sw + (ny - 1) * dlat / 2.0
        center_lon = lon_sw + (nx - 1) * dlon / 2.0
        x_full, _, _, _ = twp.extract_wrf_wind(nc_path)
        x_clip, _, _, _ = twp.extract_wrf_wind(
            nc_path,
            lat_min=center_lat - 0.1, lat_max=center_lat + 0.1,
            lon_min=center_lon - 0.1, lon_max=center_lon + 0.1,
        )
        self.assertLess(x_clip.size, x_full.size)

    def test_latlon_clip_values_within_bounds(self):
        """All returned UTM coordinates must correspond to points inside the clip box."""
        ny, nx = 20, 20
        dlat = dlon = 0.05
        lat_sw, lon_sw = 37.0, -120.0
        nc_path = os.path.join(self.tmpdir, "wrf_clip_vals.nc")
        _make_wrf_nc(nc_path, ny=ny, nx=nx, n_times=1,
                     lat_sw=lat_sw, lon_sw=lon_sw, dlat=dlat, dlon=dlon)

        center_lat = lat_sw + (ny - 1) * dlat / 2.0
        center_lon = lon_sw + (nx - 1) * dlon / 2.0
        clip_lat_min = center_lat - 0.45
        clip_lat_max = center_lat + 0.45
        clip_lon_min = center_lon - 0.45
        clip_lon_max = center_lon + 0.45

        x_clip, y_clip, u_clip, v_clip = twp.extract_wrf_wind(
            nc_path,
            lat_min=clip_lat_min, lat_max=clip_lat_max,
            lon_min=clip_lon_min, lon_max=clip_lon_max,
        )
        # Result must be non-empty 2-D arrays with matching shapes
        self.assertEqual(x_clip.ndim, 2)
        self.assertEqual(x_clip.shape, u_clip.shape)
        self.assertEqual(x_clip.shape, v_clip.shape)

        # Round-trip: convert UTM back to lat/lon and verify within clip bounds
        try:
            from pyproj import Transformer
            center_lon_approx = (clip_lon_min + clip_lon_max) / 2.0
            center_lat_approx = (clip_lat_min + clip_lat_max) / 2.0
            zone = int((center_lon_approx + 180.0) / 6.0) + 1
            epsg = 32600 + zone if center_lat_approx >= 0.0 else 32700 + zone
            t = Transformer.from_crs(f"EPSG:{epsg}", "EPSG:4326", always_xy=True)
            lons_back, lats_back = t.transform(x_clip.ravel(), y_clip.ravel())
            # All clipped points must lie within the requested bounds
            self.assertTrue(np.all(lats_back >= clip_lat_min - 1e-6))
            self.assertTrue(np.all(lats_back <= clip_lat_max + 1e-6))
            self.assertTrue(np.all(lons_back >= clip_lon_min - 1e-6))
            self.assertTrue(np.all(lons_back <= clip_lon_max + 1e-6))
        except ImportError:
            pass  # pyproj not installed; shape checks above are sufficient


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
# 9. main() - --no-terrain + --wrf-file: WRF HGT_M written as terrain file
# ===========================================================================

class TestMainNoTerrain(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def _make_nc(self, **kwargs):
        path = os.path.join(self.tmpdir, "wrf_no_terrain.nc")
        return _make_wrf_nc(path, **kwargs)

    def test_no_terrain_with_wrf_writes_terrain_file(self):
        """--no-terrain + --wrf-file: WRF HGT_M must be written as terrain file."""
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

        self.assertTrue(os.path.isfile(tout),
                        "terrain file SHOULD exist when --no-terrain + --wrf-file")

    def test_no_terrain_with_wrf_terrain_has_three_columns(self):
        """--no-terrain + --wrf-file: terrain file must have three columns (x y z)."""
        nc = self._make_nc()
        tout = os.path.join(self.tmpdir, "terrain_3col.xyz")

        twp.main([
            "--wrf-file", nc,
            "--terrain", tout,
            "--no-terrain",
            "--no-landscape",
            "--no-wind",
        ])

        self.assertTrue(os.path.isfile(tout))
        with open(tout) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                self.assertEqual(len(line.split()), 3,
                                 f"Expected 3 columns, got: {line!r}")

    def test_no_terrain_with_wrf_terrain_z_values(self):
        """--no-terrain + --wrf-file: terrain z values must match HGT_M (j*10)."""
        # _make_wrf_nc sets hgt_m[t, j, :] = j * 10 for all i
        nc = self._make_nc(ny=5, nx=4)
        tout = os.path.join(self.tmpdir, "terrain_z.xyz")

        twp.main([
            "--wrf-file", nc,
            "--terrain", tout,
            "--no-terrain",
            "--no-landscape",
            "--no-wind",
        ])

        z_values = []
        with open(tout) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.split()
                z_values.append(float(parts[2]))

        self.assertGreater(len(z_values), 0)
        # Heights must be non-negative multiples of 10 (0, 10, 20, 30, 40)
        for z in z_values:
            self.assertGreaterEqual(z, 0.0)
            self.assertAlmostEqual(z % 10, 0.0, places=3,
                                   msg=f"z={z} is not a multiple of 10")

    def test_no_terrain_still_writes_wind(self):
        """--no-terrain + --wrf-file: wind file should still be created."""
        nc = self._make_nc()
        tout = os.path.join(self.tmpdir, "terrain_wind_test.xyz")
        wout = os.path.join(self.tmpdir, "wind_no_terrain.csv")

        twp.main([
            "--wrf-file", nc,
            "--terrain", tout,
            "--wind",    wout,
            "--no-terrain",
            "--no-landscape",
        ])

        self.assertTrue(os.path.isfile(wout))
        self.assertGreater(_count_data_lines(wout), 0)

    def test_wind_file_has_four_columns(self):
        nc = self._make_nc()
        tout = os.path.join(self.tmpdir, "terrain_4col.xyz")
        wout = os.path.join(self.tmpdir, "wind_4col.csv")
        twp.main([
            "--wrf-file", nc,
            "--terrain", tout,
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
        """--time-range 0:2 should produce wind.csv, wind_1.csv, wind_2.csv."""
        nc = self._make_nc(n_times=3)
        base_wout = os.path.join(self.tmpdir, "wind.csv")
        tout = os.path.join(self.tmpdir, "terrain_trange.xyz")

        twp.main([
            "--wrf-file",   nc,
            "--wind",       base_wout,
            "--terrain",    tout,
            "--time-range", "0:2",
            "--no-terrain",
            "--no-landscape",
            "--no-inputs",
        ])

        for pos in range(3):
            expected = twp._wind_output_path(base_wout, pos)
            self.assertTrue(os.path.isfile(expected),
                            f"Expected wind file not found: {expected}")

    def test_single_time_range_no_suffix(self):
        """A single-step --time-range should NOT produce a suffixed file."""
        nc = self._make_nc(n_times=2)
        base_wout = os.path.join(self.tmpdir, "wind_single.csv")
        tout = os.path.join(self.tmpdir, "terrain_single.xyz")

        twp.main([
            "--wrf-file",   nc,
            "--wind",       base_wout,
            "--terrain",    tout,
            "--time-range", "1:1",
            "--no-terrain",
            "--no-landscape",
            "--no-inputs",
        ])
        # single time-step: no suffix expected
        self.assertTrue(os.path.isfile(base_wout))

    def test_wind_values_differ_between_time_steps(self):
        """Different time steps should yield different wind values."""
        ny, nx = 4, 5
        nc = self._make_nc(ny=ny, nx=nx, n_times=2, nz=1)
        base_wout = os.path.join(self.tmpdir, "wind_diff.csv")
        tout = os.path.join(self.tmpdir, "terrain_diff.xyz")

        twp.main([
            "--wrf-file",   nc,
            "--wind",       base_wout,
            "--terrain",    tout,
            "--time-range", "0:1",
            "--no-terrain",
            "--no-landscape",
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

        # position 0 → wind_diff.csv, position 1 → wind_diff_1.csv
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
             "--no-landscape",
             "--no-inputs"],
            capture_output=True, text=True,
            cwd=self.tmpdir,
        )
        self.assertEqual(result.returncode, 0,
                         f"CLI failed:\n{result.stdout}\n{result.stderr}")
        self.assertTrue(os.path.isfile(wout))

    def test_cli_missing_wrf_file(self):
        """Non-existent WRF file should exit non-zero."""
        result = subprocess.run(
            [sys.executable, self.script,
             "--wrf-file", "no_such_file.nc",
             "--no-terrain", "--no-landscape", "--no-inputs"],
            capture_output=True, text=True,
            cwd=self.tmpdir,
        )
        self.assertNotEqual(result.returncode, 0)

    def test_cli_missing_bbox_without_wrf(self):
        """No --wrf-file and no bbox args should exit non-zero."""
        result = subprocess.run(
            [sys.executable, self.script, "--no-terrain", "--no-landscape",
             "--no-inputs"],
            capture_output=True, text=True,
            cwd=self.tmpdir,
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
             "--no-terrain", "--no-landscape", "--no-inputs"],
            check=True, capture_output=True,
            cwd=self.tmpdir,
        )
        subprocess.run(
            [sys.executable, self.script,
             "--wrf-file", nc, "--wind", wout_sub,
             "--no-terrain", "--no-landscape", "--no-inputs",
             "--subsample", "2"],
            check=True, capture_output=True,
            cwd=self.tmpdir,
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
             "--no-terrain", "--no-landscape", "--no-inputs"],
            capture_output=True, text=True,
            cwd=self.tmpdir,
        )
        self.assertEqual(result.returncode, 0,
                         f"CLI failed:\n{result.stdout}\n{result.stderr}")

        for pos in range(3):
            expected = twp._wind_output_path(base_wout, pos)
            self.assertTrue(os.path.isfile(expected),
                            f"Missing time-step file: {expected}")

    def test_cli_invalid_time_range(self):
        """Invalid --time-range should exit non-zero."""
        nc = self._make_nc()
        result = subprocess.run(
            [sys.executable, self.script,
             "--wrf-file",   nc,
             "--time-range", "bad",
             "--no-terrain", "--no-landscape", "--no-inputs"],
            capture_output=True, text=True,
            cwd=self.tmpdir,
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
             "--no-terrain", "--no-landscape", "--no-inputs"],
            capture_output=True, text=True,
            cwd=self.tmpdir,
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
             "--no-wind", "--no-inputs"],
            capture_output=True, text=True,
            cwd=self.tmpdir,
        )
        self.assertEqual(result.returncode, 0,
                         f"CLI failed:\n{result.stdout}\n{result.stderr}")
        self.assertFalse(os.path.isfile(wout),
                         "wind file should NOT exist with --no-wind")


# ===========================================================================
# 13. write_inputs_file: fire box centered, --no-landscape FM4 default
# ===========================================================================

class TestWriteInputsFile(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def test_fire_box_centered(self):
        """Fire ignition box should be centred in the domain (45-55%)."""
        path = os.path.join(self.tmpdir, "inputs_center.i")
        bounds = (500000.0, 4000000.0, 503000.0, 4003000.0)  # 3 km × 3 km
        twp.write_inputs_file(path, domain_bounds=bounds)

        x_lo, y_lo, x_hi, y_hi = bounds
        w = x_hi - x_lo
        h = y_hi - y_lo
        expected_xmin = x_lo + 0.45 * w
        expected_xmax = x_lo + 0.55 * w
        expected_ymin = y_lo + 0.45 * h
        expected_ymax = y_lo + 0.55 * h

        with open(path) as fh:
            content = fh.read()
        self.assertIn(f"box_xmin = {expected_xmin:.2f}", content)
        self.assertIn(f"box_xmax = {expected_xmax:.2f}", content)
        self.assertIn(f"box_ymin = {expected_ymin:.2f}", content)
        self.assertIn(f"box_ymax = {expected_ymax:.2f}", content)

    def test_no_landscape_defaults_fm4(self):
        """When no landscape_file is given, inputs.i should use FM4 (S. California)."""
        path = os.path.join(self.tmpdir, "inputs_nolcp.i")
        twp.write_inputs_file(path, landscape_file=None)
        with open(path) as fh:
            content = fh.read()
        self.assertIn("rothermel.fuel_model = FM4", content)
        self.assertIn("Southern California", content)

    def test_landscape_file_no_default_comment(self):
        """When a landscape_file is given, no S.California default comment."""
        path = os.path.join(self.tmpdir, "inputs_withlcp.i")
        twp.write_inputs_file(path, landscape_file="/some/path/landscape.lcp")
        with open(path) as fh:
            content = fh.read()
        self.assertIn("rothermel.fuel_model = FM4", content)
        # The landscape file should suppress the 'defaulting to Southern California' note
        self.assertNotIn("defaulting to Southern California", content)


# ===========================================================================
# 14. main() --no-terrain: domain bounds in inputs.i from WRF UTM extents
# ===========================================================================

class TestMainNoTerrainDomainBounds(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def _make_nc(self, **kwargs):
        path = os.path.join(self.tmpdir, "wrf_domain_bounds.nc")
        return _make_wrf_nc(path, **kwargs)

    def test_inputs_has_domain_bounds_no_terrain(self):
        """--no-terrain + --wrf-file: inputs.i must contain prob_lo/hi and n_cell."""
        nc = self._make_nc(ny=10, nx=12, n_times=1)
        iout = os.path.join(self.tmpdir, "inputs_domain.i")
        tout = os.path.join(self.tmpdir, "terrain_domain.xyz")

        twp.main([
            "--wrf-file", nc,
            "--terrain",  tout,
            "--no-terrain",
            "--no-landscape",
            "--no-wind",
            "--inputs", iout,
        ])

        self.assertTrue(os.path.isfile(iout))
        with open(iout) as fh:
            content = fh.read()
        # Domain bounds must appear as uncommented lines
        self.assertRegex(content, r"(?m)^prob_lo_x\s*=")
        self.assertRegex(content, r"(?m)^prob_lo_y\s*=")
        self.assertRegex(content, r"(?m)^prob_hi_x\s*=")
        self.assertRegex(content, r"(?m)^prob_hi_y\s*=")
        self.assertRegex(content, r"(?m)^n_cell_x\s*=")
        self.assertRegex(content, r"(?m)^n_cell_y\s*=")

    def test_inputs_domain_bounds_are_metres(self):
        """Domain bounds from WRF UTM extents should be large (> 1e5 m)."""
        nc = self._make_nc(ny=10, nx=12, n_times=1)
        iout = os.path.join(self.tmpdir, "inputs_domain_metres.i")
        tout = os.path.join(self.tmpdir, "terrain_metres.xyz")

        twp.main([
            "--wrf-file", nc,
            "--terrain",  tout,
            "--no-terrain",
            "--no-landscape",
            "--no-wind",
            "--inputs", iout,
        ])

        import re
        with open(iout) as fh:
            content = fh.read()
        m = re.search(r"^prob_lo_x\s*=\s*([\d.]+)", content, re.MULTILINE)
        self.assertIsNotNone(m, "prob_lo_x not found")
        self.assertGreater(float(m.group(1)), 1e5,
                           "prob_lo_x should be in UTM metres (> 100 000 m)")

    def test_inputs_n_cell_at_30m_resolution(self):
        """n_cell_x / n_cell_y should equal round(span / 30) for WRF UTM extents."""
        nc = self._make_nc(ny=10, nx=12, n_times=1)
        iout = os.path.join(self.tmpdir, "inputs_ncell.i")
        tout = os.path.join(self.tmpdir, "terrain_ncell.xyz")

        twp.main([
            "--wrf-file", nc,
            "--terrain",  tout,
            "--no-terrain",
            "--no-landscape",
            "--no-wind",
            "--inputs", iout,
        ])

        import re
        with open(iout) as fh:
            content = fh.read()

        def _get(key):
            m = re.search(rf"^{key}\s*=\s*([\d.]+)", content, re.MULTILINE)
            self.assertIsNotNone(m, f"{key} not found")
            return float(m.group(1))

        x_lo = _get("prob_lo_x")
        x_hi = _get("prob_hi_x")
        y_lo = _get("prob_lo_y")
        y_hi = _get("prob_hi_y")
        n_cell_x = int(_get("n_cell_x"))
        n_cell_y = int(_get("n_cell_y"))

        expected_nx = max(1, round((x_hi - x_lo) / 30.0))
        expected_ny = max(1, round((y_hi - y_lo) / 30.0))
        self.assertEqual(n_cell_x, expected_nx)
        self.assertEqual(n_cell_y, expected_ny)

    def test_inputs_fire_box_centered(self):
        """Fire box in inputs.i should sit at 45-55% of domain in each axis."""
        nc = self._make_nc(ny=10, nx=12, n_times=1)
        iout = os.path.join(self.tmpdir, "inputs_firebox.i")
        tout = os.path.join(self.tmpdir, "terrain_firebox.xyz")

        twp.main([
            "--wrf-file", nc,
            "--terrain",  tout,
            "--no-terrain",
            "--no-landscape",
            "--no-wind",
            "--inputs", iout,
        ])

        import re
        with open(iout) as fh:
            content = fh.read()

        def _get(key):
            m = re.search(rf"^{key}\s*=\s*([\d.]+)", content, re.MULTILINE)
            self.assertIsNotNone(m, f"{key} not found in inputs.i")
            return float(m.group(1))

        x_lo = _get("prob_lo_x")
        x_hi = _get("prob_hi_x")
        y_lo = _get("prob_lo_y")
        y_hi = _get("prob_hi_y")
        box_xmin = _get("box_xmin")
        box_xmax = _get("box_xmax")
        box_ymin = _get("box_ymin")
        box_ymax = _get("box_ymax")

        w = x_hi - x_lo
        h = y_hi - y_lo
        self.assertAlmostEqual(box_xmin, x_lo + 0.45 * w, delta=1.0)
        self.assertAlmostEqual(box_xmax, x_lo + 0.55 * w, delta=1.0)
        self.assertAlmostEqual(box_ymin, y_lo + 0.45 * h, delta=1.0)
        self.assertAlmostEqual(box_ymax, y_lo + 0.55 * h, delta=1.0)


# ===========================================================================
# 15. _wrf_bbox_to_utm_domain_bounds
# ===========================================================================

class TestWrfBboxToUtmDomainBounds(unittest.TestCase):

    def test_returns_four_floats(self):
        result = twp._wrf_bbox_to_utm_domain_bounds(37.0, 37.9, -120.0, -119.1)
        self.assertEqual(len(result), 4)
        for v in result:
            self.assertIsInstance(v, float)

    def test_x_lo_less_than_x_hi(self):
        x_lo, y_lo, x_hi, y_hi = twp._wrf_bbox_to_utm_domain_bounds(
            37.0, 37.9, -120.0, -119.1
        )
        self.assertLess(x_lo, x_hi)
        self.assertLess(y_lo, y_hi)

    def test_utm_in_metres(self):
        x_lo, y_lo, x_hi, y_hi = twp._wrf_bbox_to_utm_domain_bounds(
            37.0, 37.9, -120.0, -119.1
        )
        for v in (x_lo, x_hi, y_lo, y_hi):
            self.assertGreater(abs(v), 1e5,
                               "Expected UTM coordinates to be in metres (> 1e5)")

    def test_wider_bbox_gives_wider_utm(self):
        """A wider lat/lon span should produce a wider UTM span."""
        _, _, x_hi_narrow, _ = twp._wrf_bbox_to_utm_domain_bounds(
            37.0, 37.1, -120.0, -119.9
        )
        x_lo_narrow, _, _, _ = twp._wrf_bbox_to_utm_domain_bounds(
            37.0, 37.1, -120.0, -119.9
        )
        _, _, x_hi_wide, _ = twp._wrf_bbox_to_utm_domain_bounds(
            37.0, 37.1, -120.0, -119.0
        )
        x_lo_wide, _, _, _ = twp._wrf_bbox_to_utm_domain_bounds(
            37.0, 37.1, -120.0, -119.0
        )
        self.assertGreater(x_hi_wide - x_lo_wide,
                           x_hi_narrow - x_lo_narrow)


# ===========================================================================
# 16. COG constants and download_landfire_cog (offline / unit tests)
# ===========================================================================

class TestCOGConstants(unittest.TestCase):
    """Verify the COG layer mapping constants are internally consistent."""

    def test_cog_layers_vintages(self):
        """_COG_LAYERS should define entries for 2014, 2016, and 2020."""
        for year in (2014, 2016, 2020):
            self.assertIn(year, twp._COG_LAYERS,
                          f"Missing _COG_LAYERS entry for vintage {year}")

    def test_cog_layers_keys(self):
        """Each vintage entry should have elev, slope, aspect, and fuel13."""
        for year, layers in twp._COG_LAYERS.items():
            for key in ("elev", "slope", "aspect", "fuel13"):
                self.assertIn(key, layers,
                              f"Missing key '{key}' for vintage {year}")

    def test_cog_paths_end_with_tif(self):
        """All COG paths should end with .tif."""
        for year, layers in twp._COG_LAYERS.items():
            for key, path in layers.items():
                self.assertTrue(
                    path.endswith(".tif"),
                    f"COG path for ({year}, {key}) does not end with .tif: {path}"
                )

    def test_cog_base_url_is_https(self):
        """_LANDFIRE_COG_BASE should be an HTTPS URL."""
        self.assertTrue(
            twp._LANDFIRE_COG_BASE.startswith("https://"),
            "_LANDFIRE_COG_BASE must use HTTPS"
        )

    def test_default_layers_match_cog_coverage(self):
        """_DEFAULT_LAYERS vintages with a COG mapping should align."""
        for year in twp._COG_LAYERS:
            self.assertIn(year, twp._DEFAULT_LAYERS,
                          f"_DEFAULT_LAYERS missing vintage {year}")


class TestDownloadLandFireCogOffline(unittest.TestCase):
    """Test download_landfire_cog without a real network connection.

    The function's rasterio call is monkey-patched so that no HTTP request
    is made.
    """

    def _make_fake_cog_bytes(self, shape=(4, 4), value=100.0):
        """Return minimal in-memory GeoTIFF bytes using rasterio MemoryFile."""
        try:
            import numpy as np
            import rasterio
            from rasterio.io import MemoryFile
            from rasterio.crs import CRS
            from rasterio.transform import from_bounds
        except ImportError:
            raise unittest.SkipTest("rasterio not installed")

        data = np.full(shape, value, dtype=np.float64)
        # Use a simple WGS-84 transform over a tiny area
        tf = from_bounds(-105.1, 40.0, -105.0, 40.1, shape[1], shape[0])
        with MemoryFile() as mf:
            with mf.open(driver="GTiff", height=shape[0], width=shape[1],
                         count=1, dtype=data.dtype,
                         crs=CRS.from_epsg(4326), transform=tf,
                         nodata=-9999.0) as ds:
                ds.write(data, 1)
            return mf.read()

    def test_download_landfire_cog_returns_four_layers(self):
        """download_landfire_cog should return exactly four layer entries."""
        try:
            import numpy as np
            import rasterio
            from rasterio.crs import CRS
            from rasterio.transform import from_bounds
        except ImportError:
            raise unittest.SkipTest("rasterio not installed")

        bbox = (-105.1, 40.0, -105.0, 40.1)
        shape = (4, 4)
        fake_data = np.ones(shape, dtype=np.float64)
        fake_tf = from_bounds(*bbox, shape[1], shape[0])
        fake_crs = CRS.from_epsg(4326)

        # Patch _read_landfire_cog to avoid any network call
        original = twp._read_landfire_cog
        twp._read_landfire_cog = lambda url, b: (fake_data, fake_tf, fake_crs, -9999.0)
        try:
            result = twp.download_landfire_cog(bbox, vintage=2020)
        finally:
            twp._read_landfire_cog = original

        self.assertEqual(len(result), 4)

    def test_download_landfire_cog_keys_match_default_layers(self):
        """The returned dict keys should match the 2020 default layer IDs."""
        try:
            import numpy as np
            import rasterio
            from rasterio.crs import CRS
            from rasterio.transform import from_bounds
        except ImportError:
            raise unittest.SkipTest("rasterio not installed")

        bbox = (-105.1, 40.0, -105.0, 40.1)
        shape = (4, 4)
        fake_data = np.ones(shape, dtype=np.float64)
        fake_tf = from_bounds(*bbox, shape[1], shape[0])
        fake_crs = CRS.from_epsg(4326)

        original = twp._read_landfire_cog
        twp._read_landfire_cog = lambda url, b: (fake_data, fake_tf, fake_crs, -9999.0)
        try:
            result = twp.download_landfire_cog(bbox, vintage=2020)
        finally:
            twp._read_landfire_cog = original

        expected_keys = set(twp._DEFAULT_LAYERS[2020].values())
        self.assertEqual(set(result.keys()), expected_keys)

    def test_download_landfire_cog_values_are_bytes(self):
        """Each value in the returned dict should be bytes (in-memory GeoTIFF)."""
        try:
            import numpy as np
            import rasterio
            from rasterio.crs import CRS
            from rasterio.transform import from_bounds
        except ImportError:
            raise unittest.SkipTest("rasterio not installed")

        bbox = (-105.1, 40.0, -105.0, 40.1)
        shape = (4, 4)
        fake_data = np.ones(shape, dtype=np.float64)
        fake_tf = from_bounds(*bbox, shape[1], shape[0])
        fake_crs = CRS.from_epsg(4326)

        original = twp._read_landfire_cog
        twp._read_landfire_cog = lambda url, b: (fake_data, fake_tf, fake_crs, -9999.0)
        try:
            result = twp.download_landfire_cog(bbox, vintage=2020)
        finally:
            twp._read_landfire_cog = original

        for lid, raw in result.items():
            self.assertIsInstance(raw, bytes,
                                  f"Value for layer '{lid}' is not bytes")
            # Minimal sanity: should be a valid GeoTIFF (starts with TIFF magic)
            self.assertIn(raw[:4], (b"II\x2a\x00", b"MM\x00\x2a"),
                          f"Layer '{lid}' bytes do not look like a GeoTIFF")

    def test_create_landscape_use_cog_false_uses_lfps(self):
        """create_landscape(use_cog=False) should call download_landfire, not COG."""
        calls = []

        def _fake_download_landfire(bbox, layer_ids, **kwargs):
            calls.append(("lfps", layer_ids))
            raise RuntimeError("stop here")  # abort after recording the call

        original = twp.download_landfire
        twp.download_landfire = _fake_download_landfire
        try:
            with self.assertRaises(RuntimeError):
                twp.create_landscape(
                    "/tmp/fake.lcp",
                    bbox=(-105.1, 40.0, -105.0, 40.1),
                    vintage=2020,
                    use_cog=False,
                )
        finally:
            twp.download_landfire = original

        self.assertEqual(len(calls), 1)
        self.assertEqual(calls[0][0], "lfps")

    def test_create_landscape_use_cog_true_uses_cog(self):
        """create_landscape(use_cog=True) should call download_landfire_cog first."""
        calls = []

        def _fake_download_landfire_cog(bbox, **kwargs):
            calls.append("cog")
            raise RuntimeError("stop here")

        original = twp.download_landfire_cog
        twp.download_landfire_cog = _fake_download_landfire_cog
        try:
            with self.assertRaises(Exception):
                twp.create_landscape(
                    "/tmp/fake.lcp",
                    bbox=(-105.1, 40.0, -105.0, 40.1),
                    vintage=2020,
                    use_cog=True,
                )
        finally:
            twp.download_landfire_cog = original

        self.assertIn("cog", calls)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    unittest.main(verbosity=2)
