#!/usr/bin/env python3
"""
Regression test for tools/wrf_to_terrain_wind.py

Tests:
  1. Destagger U (nx+1 columns) to mass-point centres
  2. Destagger V (ny+1 rows) to mass-point centres
  3. U/V already at mass points (no destagger needed)
  4. lat/lon → UTM projection preserves shape and produces metre-scale values
  5. Full pipeline: synthetic WRF-style netCDF → terrain.csv + wind.csv
  6. --subsample flag reduces output point count
  7. --time-index and --level flags select the correct slice
  8. CLI entry-point (subprocess invocation)

Run with:
  python3 tools/tests/test_wrf_to_terrain_wind.py
  python3 -m pytest tools/tests/test_wrf_to_terrain_wind.py -v
"""

import os
import subprocess
import sys
import tempfile
import unittest

import numpy as np

# Add repository root / tools/deprecated to path
_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
_TOOLS_DIR = os.path.join(_REPO_ROOT, "tools")
_DEPRECATED_DIR = os.path.join(_TOOLS_DIR, "deprecated")
sys.path.insert(0, _DEPRECATED_DIR)

import wrf_to_terrain_wind as w2tw  # noqa: E402


# ---------------------------------------------------------------------------
# Helper: create synthetic WRF-style netCDF file
# ---------------------------------------------------------------------------

def _make_wrf_nc(path, ny=10, nx=12, nz=3, n_times=2,
                 lat_sw=37.0, lon_sw=-120.0,
                 dlat=0.01, dlon=0.01):
    """Create a minimal WRF-style netCDF file.

    Variables written:
      XLAT   (Time, south_north, west_east)
      XLONG  (Time, south_north, west_east)
      HGT_M  (Time, south_north, west_east)
      U      (Time, bottom_top, south_north, west_east_stag)  — nx+1 columns
      V      (Time, bottom_top, south_north_stag, west_east)  — ny+1 rows
    """
    try:
        import netCDF4 as nc
    except ImportError:
        raise unittest.SkipTest("netCDF4 not installed – skipping WRF tests")

    ds = nc.Dataset(path, "w", format="NETCDF4")

    # Dimensions
    ds.createDimension("Time",           n_times)
    ds.createDimension("south_north",    ny)
    ds.createDimension("west_east",      nx)
    ds.createDimension("bottom_top",     nz)
    ds.createDimension("south_north_stag", ny + 1)
    ds.createDimension("west_east_stag",   nx + 1)

    # XLAT / XLONG
    lats = np.array([[lat_sw + j * dlat for _ in range(nx)] for j in range(ny)],
                    dtype=np.float32)
    lons = np.array([[lon_sw + i * dlon for i in range(nx)] for _ in range(ny)],
                    dtype=np.float32)

    xlat  = ds.createVariable("XLAT",  "f4", ("Time", "south_north", "west_east"))
    xlong = ds.createVariable("XLONG", "f4", ("Time", "south_north", "west_east"))
    for t in range(n_times):
        xlat[t, :, :] = lats
        xlong[t, :, :] = lons

    # HGT_M: simple terrain (increases northward)
    hgt_m = ds.createVariable("HGT_M", "f4", ("Time", "south_north", "west_east"))
    for t in range(n_times):
        for j in range(ny):
            hgt_m[t, j, :] = float(j * 10)  # 0, 10, 20, … m

    # U: staggered west-east, values = level * 10 + time * 100 + col
    u_var = ds.createVariable("U", "f4",
                               ("Time", "bottom_top", "south_north", "west_east_stag"))
    for t in range(n_times):
        for k in range(nz):
            for j in range(ny):
                for i in range(nx + 1):
                    u_var[t, k, j, i] = float(k * 10 + t * 100 + i)

    # V: staggered south-north, values = level * 10 + time * 100 + row
    v_var = ds.createVariable("V", "f4",
                               ("Time", "bottom_top", "south_north_stag", "west_east"))
    for t in range(n_times):
        for k in range(nz):
            for j in range(ny + 1):
                for i in range(nx):
                    v_var[t, k, j, i] = float(k * 10 + t * 100 + j)

    ds.close()
    return path


# ---------------------------------------------------------------------------
# Unit tests: individual functions
# ---------------------------------------------------------------------------

class TestDestaggerU(unittest.TestCase):
    def test_shape(self):
        U_stag = np.ones((10, 13))  # ny=10, nx+1=13
        u_mass = w2tw._destagger_u(U_stag)
        self.assertEqual(u_mass.shape, (10, 12))

    def test_average(self):
        """Destagger of [0, 2] → [1]."""
        U_stag = np.array([[0.0, 2.0]])
        u_mass = w2tw._destagger_u(U_stag)
        np.testing.assert_allclose(u_mass, np.array([[1.0]]))

    def test_known_values(self):
        """Check centre-cell averages for a small known matrix."""
        U_stag = np.array([
            [1.0, 3.0, 5.0],  # cols 0,1,2 (staggered nx+1=3)
            [2.0, 4.0, 6.0],
        ])
        expected = np.array([
            [2.0, 4.0],
            [3.0, 5.0],
        ])
        np.testing.assert_allclose(w2tw._destagger_u(U_stag), expected)


class TestDestaggerV(unittest.TestCase):
    def test_shape(self):
        V_stag = np.ones((11, 12))  # ny+1=11, nx=12
        v_mass = w2tw._destagger_v(V_stag)
        self.assertEqual(v_mass.shape, (10, 12))

    def test_average(self):
        V_stag = np.array([[0.0], [2.0]])
        v_mass = w2tw._destagger_v(V_stag)
        np.testing.assert_allclose(v_mass, np.array([[1.0]]))

    def test_known_values(self):
        V_stag = np.array([
            [1.0, 2.0],
            [3.0, 4.0],
            [5.0, 6.0],
        ])
        expected = np.array([
            [2.0, 3.0],
            [4.0, 5.0],
        ])
        np.testing.assert_allclose(w2tw._destagger_v(V_stag), expected)


class TestLatLonToUtm(unittest.TestCase):
    def test_shape_preserved(self):
        lats = np.array([[37.0, 37.1], [37.2, 37.3]])
        lons = np.array([[-120.0, -119.9], [-120.0, -119.9]])
        x, y = w2tw._latlon_to_utm(lats, lons)
        self.assertEqual(x.shape, lats.shape)
        self.assertEqual(y.shape, lats.shape)

    def test_utm_in_metres(self):
        """Projected coordinates should be in the hundreds-of-thousands range."""
        lats = np.array([[37.0]])
        lons = np.array([[-120.0]])
        x, y = w2tw._latlon_to_utm(lats, lons)
        self.assertGreater(abs(float(x[0, 0])), 1e5)
        self.assertGreater(abs(float(y[0, 0])), 1e5)

    def test_east_increases_x(self):
        """Moving east increases UTM x."""
        lats = np.array([[37.0, 37.0]])
        lons = np.array([[-120.5, -119.5]])
        x, _ = w2tw._latlon_to_utm(lats, lons)
        self.assertLess(float(x[0, 0]), float(x[0, 1]))

    def test_north_increases_y(self):
        """Moving north increases UTM y."""
        lats = np.array([[36.0], [38.0]])
        lons = np.array([[-120.0], [-120.0]])
        _, y = w2tw._latlon_to_utm(lats, lons)
        self.assertLess(float(y[0, 0]), float(y[1, 0]))


# ---------------------------------------------------------------------------
# Integration tests: full pipeline
# ---------------------------------------------------------------------------

class TestConvertWrf(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def _make_nc(self, **kwargs):
        path = os.path.join(self.tmpdir, "wrfout.nc")
        return _make_wrf_nc(path, **kwargs)

    def test_output_files_created(self):
        nc_path = self._make_nc()
        terrain_out = os.path.join(self.tmpdir, "terrain.csv")
        wind_out    = os.path.join(self.tmpdir, "wind.csv")

        w2tw.convert_wrf(nc_path, terrain_out, wind_out)

        self.assertTrue(os.path.isfile(terrain_out))
        self.assertTrue(os.path.isfile(wind_out))

    def test_terrain_row_count(self):
        """Terrain file should have ny*nx data rows."""
        ny, nx = 8, 10
        nc_path = self._make_nc(ny=ny, nx=nx)
        terrain_out = os.path.join(self.tmpdir, "terrain.csv")
        wind_out    = os.path.join(self.tmpdir, "wind.csv")

        w2tw.convert_wrf(nc_path, terrain_out, wind_out)

        with open(terrain_out) as fh:
            data_lines = [l for l in fh if not l.startswith("#") and l.strip()]
        self.assertEqual(len(data_lines), ny * nx)

    def test_terrain_column_count(self):
        """Each terrain row must have exactly 3 columns."""
        nc_path = self._make_nc()
        terrain_out = os.path.join(self.tmpdir, "t3col.csv")
        wind_out    = os.path.join(self.tmpdir, "w3col.csv")
        w2tw.convert_wrf(nc_path, terrain_out, wind_out)

        with open(terrain_out) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.split()
                self.assertEqual(len(parts), 3, f"Bad terrain line: {line!r}")

    def test_wind_column_count(self):
        """Each wind row must have exactly 4 columns."""
        nc_path = self._make_nc()
        terrain_out = os.path.join(self.tmpdir, "tw.csv")
        wind_out    = os.path.join(self.tmpdir, "ww.csv")
        w2tw.convert_wrf(nc_path, terrain_out, wind_out)

        with open(wind_out) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.split()
                self.assertEqual(len(parts), 4, f"Bad wind line: {line!r}")

    def test_utm_coordinates_in_metres(self):
        """UTM x/y values should exceed 1e5 (metres, not degrees)."""
        nc_path = self._make_nc()
        terrain_out = os.path.join(self.tmpdir, "utm_check.csv")
        wind_out    = os.path.join(self.tmpdir, "wutm.csv")
        w2tw.convert_wrf(nc_path, terrain_out, wind_out)

        with open(terrain_out) as fh:
            first_data = next(l for l in fh if not l.startswith("#") and l.strip())
        x_val = abs(float(first_data.split()[0]))
        self.assertGreater(x_val, 1e5)

    def test_destagger_u_values(self):
        """Destaggered U should equal average of adjacent staggered columns.

        In our synthetic file, U[t=0, k=0, j, i] = i (column index).
        After destaggering: u_mass[j, i] = 0.5*(i + i+1) = i + 0.5.
        """
        ny, nx = 4, 6
        nc_path = self._make_nc(ny=ny, nx=nx, n_times=1, nz=1)
        terrain_out = os.path.join(self.tmpdir, "tu.csv")
        wind_out    = os.path.join(self.tmpdir, "wu.csv")
        w2tw.convert_wrf(nc_path, terrain_out, wind_out, time_index=0, level=0)

        # Read U values from wind file
        u_vals = []
        with open(wind_out) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.split()
                u_vals.append(float(parts[2]))

        # Expected: for each row j, columns i=0..nx-1 → u = i + 0.5
        expected_u = [float(i) + 0.5 for _ in range(ny) for i in range(nx)]
        np.testing.assert_allclose(u_vals, expected_u, rtol=1e-5)

    def test_destagger_v_values(self):
        """Destaggered V should equal average of adjacent staggered rows.

        In our synthetic file, V[t=0, k=0, j, i] = j (row index).
        After destaggering: v_mass[j, i] = 0.5*(j + j+1) = j + 0.5.
        """
        ny, nx = 4, 6
        nc_path = self._make_nc(ny=ny, nx=nx, n_times=1, nz=1)
        terrain_out = os.path.join(self.tmpdir, "tv.csv")
        wind_out    = os.path.join(self.tmpdir, "wv.csv")
        w2tw.convert_wrf(nc_path, terrain_out, wind_out, time_index=0, level=0)

        v_vals = []
        with open(wind_out) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.split()
                v_vals.append(float(parts[3]))

        expected_v = [float(j) + 0.5 for j in range(ny) for _ in range(nx)]
        np.testing.assert_allclose(v_vals, expected_v, rtol=1e-5)

    def test_time_index_selection(self):
        """--time-index selects the correct time snapshot.

        Synthetic U[t, k=0, j, i] = t*100 + i.
        At t=0: u_mass = i + 0.5; at t=1: u_mass = 100 + i + 0.5.
        """
        ny, nx = 3, 4
        nc_path = self._make_nc(ny=ny, nx=nx, n_times=2, nz=1)

        for t_idx in [0, 1]:
            terrain_out = os.path.join(self.tmpdir, f"terrain_t{t_idx}.csv")
            wind_out    = os.path.join(self.tmpdir, f"wind_t{t_idx}.csv")
            w2tw.convert_wrf(nc_path, terrain_out, wind_out, time_index=t_idx)

            u_vals = []
            with open(wind_out) as fh:
                for line in fh:
                    if line.startswith("#") or not line.strip():
                        continue
                    u_vals.append(float(line.split()[2]))

            offset = t_idx * 100.0
            expected = [offset + float(i) + 0.5 for _ in range(ny) for i in range(nx)]
            np.testing.assert_allclose(u_vals, expected, rtol=1e-4,
                                       err_msg=f"time_index={t_idx}")

    def test_level_selection(self):
        """--level selects the correct vertical level.

        Synthetic U[t=0, k, j, i] = k*10 + i.
        Level k: u_mass = k*10 + i + 0.5.
        """
        ny, nx, nz = 3, 4, 3
        nc_path = self._make_nc(ny=ny, nx=nx, n_times=1, nz=nz)

        for lev in range(nz):
            terrain_out = os.path.join(self.tmpdir, f"terrain_lev{lev}.csv")
            wind_out    = os.path.join(self.tmpdir, f"wind_lev{lev}.csv")
            w2tw.convert_wrf(nc_path, terrain_out, wind_out, level=lev)

            u_vals = []
            with open(wind_out) as fh:
                for line in fh:
                    if line.startswith("#") or not line.strip():
                        continue
                    u_vals.append(float(line.split()[2]))

            offset = lev * 10.0
            expected = [offset + float(i) + 0.5 for _ in range(ny) for i in range(nx)]
            np.testing.assert_allclose(u_vals, expected, rtol=1e-4,
                                       err_msg=f"level={lev}")

    def test_subsample(self):
        """--subsample reduces the number of output rows."""
        ny, nx = 10, 12
        nc_path = self._make_nc(ny=ny, nx=nx)
        terrain_full = os.path.join(self.tmpdir, "terrain_full.csv")
        wind_full    = os.path.join(self.tmpdir, "wind_full.csv")
        terrain_sub  = os.path.join(self.tmpdir, "terrain_sub.csv")
        wind_sub     = os.path.join(self.tmpdir, "wind_sub.csv")

        w2tw.convert_wrf(nc_path, terrain_full, wind_full)
        w2tw.convert_wrf(nc_path, terrain_sub,  wind_sub,  subsample=2)

        def count(path):
            with open(path) as fh:
                return sum(1 for l in fh if not l.startswith("#") and l.strip())

        self.assertLess(count(terrain_sub), count(terrain_full))
        self.assertLess(count(wind_sub),    count(wind_full))

    def test_terrain_z_values(self):
        """Terrain Z should match the synthetic HGT_M pattern (j * 10)."""
        ny, nx = 5, 6
        nc_path = self._make_nc(ny=ny, nx=nx, n_times=1)
        terrain_out = os.path.join(self.tmpdir, "tz.csv")
        wind_out    = os.path.join(self.tmpdir, "wz.csv")
        w2tw.convert_wrf(nc_path, terrain_out, wind_out)

        z_vals = []
        with open(terrain_out) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                z_vals.append(float(line.split()[2]))

        # Expected: row j repeated nx times, value = j * 10
        expected = [float(j * 10) for j in range(ny) for _ in range(nx)]
        np.testing.assert_allclose(z_vals, expected, atol=1e-3)


# ---------------------------------------------------------------------------
# CLI tests
# ---------------------------------------------------------------------------

class TestCLI(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.script = os.path.join(_DEPRECATED_DIR, "wrf_to_terrain_wind.py")

    def _make_nc(self, **kwargs):
        path = os.path.join(self.tmpdir, "wrfout_cli.nc")
        return _make_wrf_nc(path, **kwargs)

    def test_cli_basic(self):
        """CLI invocation produces terrain and wind files."""
        nc = self._make_nc()
        tout = os.path.join(self.tmpdir, "t.csv")
        wout = os.path.join(self.tmpdir, "w.csv")

        result = subprocess.run(
            [sys.executable, self.script, nc, tout, wout],
            capture_output=True, text=True
        )
        self.assertEqual(result.returncode, 0,
                         f"CLI failed:\n{result.stdout}\n{result.stderr}")
        self.assertTrue(os.path.isfile(tout))
        self.assertTrue(os.path.isfile(wout))

    def test_cli_missing_input(self):
        """CLI exits non-zero when the WRF file does not exist."""
        result = subprocess.run(
            [sys.executable, self.script, "nofile.nc",
             "/tmp/t.csv", "/tmp/w.csv"],
            capture_output=True, text=True
        )
        self.assertNotEqual(result.returncode, 0)

    def test_cli_subsample_flag(self):
        """--subsample flag reduces output rows."""
        ny, nx = 8, 10
        nc = self._make_nc(ny=ny, nx=nx)
        tout_full = os.path.join(self.tmpdir, "tf.csv")
        wout_full = os.path.join(self.tmpdir, "wf.csv")
        tout_sub  = os.path.join(self.tmpdir, "ts.csv")
        wout_sub  = os.path.join(self.tmpdir, "ws.csv")

        subprocess.run([sys.executable, self.script, nc, tout_full, wout_full],
                       check=True, capture_output=True)
        subprocess.run([sys.executable, self.script, nc, tout_sub, wout_sub,
                        "--subsample", "2"],
                       check=True, capture_output=True)

        def count(path):
            with open(path) as fh:
                return sum(1 for l in fh if not l.startswith("#") and l.strip())

        self.assertLess(count(tout_sub), count(tout_full))


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    unittest.main(verbosity=2)
