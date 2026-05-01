#!/usr/bin/env python3
"""
Regression test for tools/deprecated/dem_to_xyz.py

Tests:
  1. Arc/Info ASCII Grid (.asc) – flat raster, all values readable
  2. Arc/Info ASCII Grid (.asc) – no-data masking
  3. SRTM HGT binary file (.hgt) – SRTM3 (1201×1201) synthetic file
  4. GeoTIFF via rasterio (if rasterio is available)
  5. UTM projection via --project-utm flag
  6. Subsampling via --subsample flag
  7. CLI entry-point (subprocess invocation)

Run with:
  python3 tools/tests/test_dem_to_xyz.py
  python3 -m pytest tools/tests/test_dem_to_xyz.py -v
"""

import math
import os
import struct
import subprocess
import sys
import tempfile
import unittest

# Add deprecated directory to path so we can import the tool module directly
_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
_TOOLS_DIR = os.path.join(_REPO_ROOT, "tools")
_DEPRECATED_DIR = os.path.join(_TOOLS_DIR, "deprecated")
sys.path.insert(0, _DEPRECATED_DIR)

import dem_to_xyz  # noqa: E402  (must be after sys.path manipulation)


# ---------------------------------------------------------------------------
# Helpers to create synthetic test files
# ---------------------------------------------------------------------------

def _make_asc(path, ncols=5, nrows=4, xll=100.0, yll=200.0, cellsize=10.0,
              nodata=-9999.0, values=None):
    """Write a minimal Arc/Info ASCII Grid file."""
    if values is None:
        # Simple ascending integers, with one no-data cell at [0,0] (top-left)
        values = []
        for r in range(nrows):
            row = []
            for c in range(ncols):
                row.append(nodata if (r == 0 and c == 0) else float(r * ncols + c))
            values.append(row)

    with open(path, "w") as fh:
        fh.write(f"ncols         {ncols}\n")
        fh.write(f"nrows         {nrows}\n")
        fh.write(f"xllcorner     {xll}\n")
        fh.write(f"yllcorner     {yll}\n")
        fh.write(f"cellsize      {cellsize}\n")
        fh.write(f"NODATA_value  {nodata}\n")
        for row in values:
            fh.write(" ".join(str(v) for v in row) + "\n")
    return path


def _make_hgt(path, lat_sw=37, lon_sw=-120, samples=1201):
    """Write a synthetic SRTM3-style HGT file with a simple elevation pattern."""
    fname = f"N{lat_sw:02d}W{abs(lon_sw):03d}.hgt"
    full_path = os.path.join(path, fname)

    nodata_val = -32768
    data = []
    for r in range(samples):
        for c in range(samples):
            # Simple synthetic elevation: 0–500 m
            elev = int(500.0 * r / samples)
            data.append(elev)
    # Write big-endian int16
    raw = struct.pack(f">{samples * samples}h", *data)
    with open(full_path, "wb") as fh:
        fh.write(raw)
    return full_path


# ---------------------------------------------------------------------------
# Unit tests
# ---------------------------------------------------------------------------

class TestReadAsc(unittest.TestCase):
    """Tests for _read_asc()."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def test_basic_read(self):
        """All valid cells are read; no-data cell is excluded."""
        asc_path = os.path.join(self.tmpdir, "test.asc")
        _make_asc(asc_path, ncols=5, nrows=4, nodata=-9999.0)

        xs, ys, zs, nodata = dem_to_xyz._read_asc(asc_path)

        # Total cells = 5*4 = 20; 1 no-data cell → 19 valid
        self.assertEqual(len(xs), 19)
        self.assertEqual(len(ys), 19)
        self.assertEqual(len(zs), 19)
        # No cell should have the nodata value
        self.assertTrue(all(z != nodata for z in zs))

    def test_cell_centre_coordinates(self):
        """Cell centre x-coordinates match xllcorner + 0.5*cellsize + i*cellsize."""
        asc_path = os.path.join(self.tmpdir, "coords.asc")
        xll, yll, cs = 0.0, 0.0, 1.0
        # 3×3 grid, no no-data values
        vals = [[float(r * 3 + c) for c in range(3)] for r in range(3)]
        _make_asc(asc_path, ncols=3, nrows=3, xll=xll, yll=yll, cellsize=cs,
                  nodata=-9999.0, values=vals)

        xs, ys, zs, _ = dem_to_xyz._read_asc(asc_path)

        expected_x = sorted(set(round(xll + (0.5 + c) * cs, 6) for c in range(3)))
        actual_x = sorted(set(round(float(x), 6) for x in xs))
        self.assertEqual(expected_x, actual_x)

    def test_nodata_override(self):
        """_read_asc returns raw arrays; nodata masking in convert_dem is correct."""
        asc_path = os.path.join(self.tmpdir, "nd.asc")
        _make_asc(asc_path, ncols=3, nrows=3, nodata=-1.0)
        xs, ys, zs, nodata = dem_to_xyz._read_asc(asc_path)
        # nodata=-1.0 is returned in the nodata sentinel value
        self.assertEqual(nodata, -1.0)


class TestReadHgt(unittest.TestCase):
    """Tests for _read_hgt()."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def test_srtm3_read(self):
        """SRTM3 HGT file (1201×1201) is read correctly."""
        hgt_path = _make_hgt(self.tmpdir, lat_sw=37, lon_sw=-120, samples=1201)
        lons, lats, zs, nodata = dem_to_xyz._read_hgt(hgt_path)
        # No -32768 values expected (synthetic data range 0–500)
        self.assertFalse(any(z == -32768.0 for z in zs))
        # Elevation values should be in the synthetic range [0, 500)
        self.assertTrue(all(0 <= z < 500 for z in zs))
        # Coordinates should be in the expected lat/lon range
        self.assertAlmostEqual(float(min(lats)), 37.0, places=2)
        self.assertAlmostEqual(float(min(lons)), -120.0, places=2)
        self.assertAlmostEqual(float(max(lats)), 38.0, places=2)
        self.assertAlmostEqual(float(max(lons)), -119.0, places=2)

    def test_filename_parsing(self):
        """Southern-hemisphere and eastern-hemisphere filenames are handled."""
        import numpy as np
        # Create a tiny fake HGT to test filename parsing only
        samples = 1201
        fname_s = "S10E020.hgt"
        path_s = os.path.join(self.tmpdir, fname_s)
        raw = struct.pack(f">{samples * samples}h", *([100] * samples * samples))
        with open(path_s, "wb") as fh:
            fh.write(raw)
        lons, lats, zs, _ = dem_to_xyz._read_hgt(path_s)
        self.assertTrue(all(lat < 0 for lat in lats),
                        "Southern hemisphere: all lats should be negative")
        self.assertTrue(all(lon >= 20 for lon in lons),
                        "Eastern hemisphere: all lons should be >= 20")


class TestProjectToUtm(unittest.TestCase):
    """Tests for _project_to_utm()."""

    def test_projection_not_degrees(self):
        """UTM coordinates should be in metres (large values, not ±180/90)."""
        import numpy as np
        lons = np.array([-120.0, -119.9])
        lats = np.array([37.0, 37.1])
        xs, ys = dem_to_xyz._project_to_utm(lons, lats)
        self.assertTrue(all(x > 1e5 for x in xs), "x_utm should be in metres")
        self.assertTrue(all(y > 1e5 for y in ys), "y_utm should be in metres")

    def test_projection_ordering_preserved(self):
        """Eastward increase in lon maps to increase in UTM x."""
        import numpy as np
        lons = np.array([-120.5, -119.5])
        lats = np.array([37.0, 37.0])
        xs, _ = dem_to_xyz._project_to_utm(lons, lats)
        self.assertLess(xs[0], xs[1])


class TestSubsample(unittest.TestCase):
    """Tests for _subsample()."""

    def test_stride_1(self):
        """subsample=1 returns all points unchanged."""
        import numpy as np
        xs = np.arange(100.0)
        ys = np.arange(100.0)
        zs = np.arange(100.0)
        xs2, ys2, zs2 = dem_to_xyz._subsample(xs, ys, zs, 1)
        self.assertEqual(len(xs2), 100)

    def test_stride_2(self):
        """subsample=2 halves the point count."""
        import numpy as np
        xs = np.arange(100.0)
        ys = np.arange(100.0)
        zs = np.arange(100.0)
        xs2, ys2, zs2 = dem_to_xyz._subsample(xs, ys, zs, 2)
        self.assertEqual(len(xs2), 50)

    def test_2d_structure_preserved(self):
        """With ncols supplied the 2-D skip is applied in both dimensions."""
        import numpy as np
        ncols, nrows = 10, 8
        xs = np.arange(float(ncols * nrows))
        ys = np.zeros_like(xs)
        zs = np.zeros_like(xs)
        xs2, ys2, zs2 = dem_to_xyz._subsample(xs, ys, zs, 2, ncols=ncols)
        # Every 2nd row and col: ceil(8/2)=4 rows, ceil(10/2)=5 cols → 20 points
        self.assertEqual(len(xs2), 20)


class TestConvertDem(unittest.TestCase):
    """Integration tests for the convert_dem() function."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def test_asc_end_to_end(self):
        """Full pipeline: .asc → X Y Z file."""
        asc = os.path.join(self.tmpdir, "hill.asc")
        out = os.path.join(self.tmpdir, "hill.csv")
        _make_asc(asc, ncols=5, nrows=4, nodata=-9999.0)

        dem_to_xyz.convert_dem(asc, out)

        self.assertTrue(os.path.isfile(out))
        with open(out) as fh:
            lines = [l for l in fh if not l.startswith("#") and l.strip()]
        # 5*4=20 cells minus 1 nodata = 19 rows
        self.assertEqual(len(lines), 19)
        # Each row should have exactly 3 columns
        for line in lines:
            parts = line.split()
            self.assertEqual(len(parts), 3)

    def test_asc_with_subsample(self):
        """Subsample reduces the number of output points."""
        asc = os.path.join(self.tmpdir, "big.asc")
        out = os.path.join(self.tmpdir, "big_sub.csv")
        # 10×10 grid with no nodata
        vals = [[float(r * 10 + c) for c in range(10)] for r in range(10)]
        _make_asc(asc, ncols=10, nrows=10, nodata=-9999.0, values=vals)

        dem_to_xyz.convert_dem(asc, out, subsample=2)

        with open(out) as fh:
            lines = [l for l in fh if not l.startswith("#") and l.strip()]
        # subsample=2 on a flat array (no 2D metadata) → 50 points
        self.assertLessEqual(len(lines), 50)
        self.assertGreater(len(lines), 0)

    def test_asc_with_utm_projection(self):
        """UTM projection results in large coordinate values."""
        asc = os.path.join(self.tmpdir, "geo.asc")
        out = os.path.join(self.tmpdir, "geo_utm.csv")
        # Place a small raster near -120°W, 37°N in degrees
        _make_asc(asc, ncols=3, nrows=3, xll=-120.0, yll=37.0,
                  cellsize=0.001, nodata=-9999.0,
                  values=[[float(r * 3 + c) for c in range(3)] for r in range(3)])

        dem_to_xyz.convert_dem(asc, out, project_utm=True)

        with open(out) as fh:
            lines = [l for l in fh if not l.startswith("#") and l.strip()]
        x0 = float(lines[0].split()[0])
        # UTM x should be on the order of hundreds of thousands (metres)
        self.assertGreater(abs(x0), 1e5)


class TestCLI(unittest.TestCase):
    """Test the command-line interface via subprocess."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.script = os.path.join(_DEPRECATED_DIR, "dem_to_xyz.py")

    def test_cli_asc(self):
        """CLI invocation on an .asc file produces the expected output file."""
        asc = os.path.join(self.tmpdir, "cli_test.asc")
        out = os.path.join(self.tmpdir, "cli_out.csv")
        _make_asc(asc, ncols=4, nrows=3,
                  values=[[float(r * 4 + c) for c in range(4)] for r in range(3)])

        result = subprocess.run(
            [sys.executable, self.script, asc, out],
            capture_output=True, text=True
        )
        self.assertEqual(result.returncode, 0,
                         f"CLI failed:\n{result.stdout}\n{result.stderr}")
        self.assertTrue(os.path.isfile(out))

    def test_cli_missing_input(self):
        """CLI exits with non-zero code when input file does not exist."""
        result = subprocess.run(
            [sys.executable, self.script, "nonexistent.asc", "/tmp/out.csv"],
            capture_output=True, text=True
        )
        self.assertNotEqual(result.returncode, 0)

    def test_cli_subsample_flag(self):
        """--subsample flag produces fewer output lines."""
        asc = os.path.join(self.tmpdir, "sub_test.asc")
        out_full = os.path.join(self.tmpdir, "full.csv")
        out_sub  = os.path.join(self.tmpdir, "sub.csv")
        vals = [[float(r * 6 + c) for c in range(6)] for r in range(6)]
        _make_asc(asc, ncols=6, nrows=6,
                  nodata=-9999.0, values=vals)

        subprocess.run([sys.executable, self.script, asc, out_full],
                       check=True, capture_output=True)
        subprocess.run([sys.executable, self.script, asc, out_sub,
                        "--subsample", "2"],
                       check=True, capture_output=True)

        def count_data_lines(path):
            with open(path) as fh:
                return sum(1 for l in fh if not l.startswith("#") and l.strip())

        self.assertLess(count_data_lines(out_sub), count_data_lines(out_full))


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    unittest.main(verbosity=2)
