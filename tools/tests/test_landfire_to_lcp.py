#!/usr/bin/env python3
"""
Regression/unit tests for tools/deprecated/landfire_to_lcp.py

Tests:
  1. assemble_landscape() with synthetic raster arrays
  2. _write_lcp() output format (correct columns, header, non-burnable filtering)
  3. create_landscape_from_files() end-to-end pipeline using GeoTIFF fixtures
  4. Non-burnable pixel filtering (codes 91-99 excluded by default)
  5. --keep-nonburnable flag includes non-burnable pixels with code 0
  6. subsample reduces the number of output points
  7. CLI --bbox validation (missing bbox, bad bounds)
  8. CLI local-files mode (--elev-file etc.) when rasterio is available
  9. CLI missing local-file argument produces non-zero exit

Run with:
  python3 tools/tests/test_landfire_to_lcp.py
  python3 -m pytest tools/tests/test_landfire_to_lcp.py -v
"""

import os
import subprocess
import sys
import tempfile
import unittest

import numpy as np

# Add deprecated directory to path so we can import the tool module directly
_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
_TOOLS_DIR = os.path.join(_REPO_ROOT, "tools")
_DEPRECATED_DIR = os.path.join(_TOOLS_DIR, "deprecated")
sys.path.insert(0, _DEPRECATED_DIR)

import landfire_to_lcp as lfp  # noqa: E402


# ---------------------------------------------------------------------------
# Helpers: synthetic raster fixtures
# ---------------------------------------------------------------------------

def _make_geotiff(path, data, west=0.0, south=0.0, cellsize=30.0,
                  nodata=None, epsg=4326):
    """Write a small GeoTIFF from a 2-D numpy array.

    Requires rasterio.  Raises ``unittest.SkipTest`` if rasterio is absent.
    """
    try:
        import rasterio
        from rasterio.transform import from_bounds
        from rasterio.crs import CRS
    except ImportError:
        raise unittest.SkipTest("rasterio not installed – skipping GeoTIFF tests")

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


def _make_synthetic_rasters(tmpdir, nrows=8, ncols=10, cellsize=30.0,
                             west=-120.0, south=37.0):
    """Create four minimal GeoTIFF rasters (elev, slope, aspect, fuel)."""
    import numpy as np

    elev   = (100.0 + np.arange(nrows * ncols).reshape(nrows, ncols)).astype(np.float32)
    slope  = (np.linspace(5.0, 30.0, nrows * ncols).reshape(nrows, ncols)).astype(np.float32)
    aspect = (np.linspace(0.0, 359.0, nrows * ncols).reshape(nrows, ncols)).astype(np.float32)
    # Mix of valid (1-13) and non-burnable codes
    fuel_vals = list(range(1, 14)) * (nrows * ncols // 13 + 1)
    fuel_vals = fuel_vals[:nrows * ncols]
    fuel = np.array(fuel_vals, dtype=np.int16).reshape(nrows, ncols)

    paths = {}
    paths["elev"]   = _make_geotiff(os.path.join(tmpdir, "elev.tif"),   elev,   west=west, south=south, cellsize=cellsize)
    paths["slope"]  = _make_geotiff(os.path.join(tmpdir, "slope.tif"),  slope,  west=west, south=south, cellsize=cellsize)
    paths["aspect"] = _make_geotiff(os.path.join(tmpdir, "aspect.tif"), aspect, west=west, south=south, cellsize=cellsize)
    paths["fuel"]   = _make_geotiff(os.path.join(tmpdir, "fuel.tif"),   fuel,   west=west, south=south, cellsize=cellsize)
    return paths, elev, slope, aspect, fuel


# ---------------------------------------------------------------------------
# Unit tests: _write_lcp()
# ---------------------------------------------------------------------------

class TestWriteLcp(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def test_output_file_created(self):
        xs = np.array([0.0, 1.0])
        ys = np.array([0.0, 1.0])
        elev  = np.array([100.0, 110.0])
        slope = np.array([5.0, 10.0])
        aspect = np.array([90.0, 180.0])
        fuel  = np.array([4.0, 7.0])
        out = os.path.join(self.tmpdir, "test.lcp")
        lfp._write_lcp(out, xs, ys, elev, slope, aspect, fuel)
        self.assertTrue(os.path.isfile(out))

    def test_header_present(self):
        out = os.path.join(self.tmpdir, "hdr.lcp")
        lfp._write_lcp(out, np.array([0.0]), np.array([0.0]),
                       np.array([100.0]), np.array([5.0]),
                       np.array([90.0]), np.array([4.0]))
        with open(out) as fh:
            content = fh.read()
        self.assertIn("#", content)

    def test_six_columns(self):
        xs     = np.array([10.0, 20.0, 30.0])
        ys     = np.array([10.0, 20.0, 30.0])
        elev   = np.array([150.0, 160.0, 170.0])
        slope  = np.array([10.0, 15.0, 20.0])
        aspect = np.array([90.0, 180.0, 270.0])
        fuel   = np.array([4.0, 5.0, 6.0])
        out = os.path.join(self.tmpdir, "cols.lcp")
        lfp._write_lcp(out, xs, ys, elev, slope, aspect, fuel)
        with open(out) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.split()
                self.assertEqual(len(parts), 6,
                                 f"Expected 6 columns, got: {line!r}")

    def test_fuel_model_integer(self):
        """Fuel model column should be an integer, not a float."""
        out = os.path.join(self.tmpdir, "fuel_int.lcp")
        lfp._write_lcp(out, np.array([0.0]), np.array([0.0]),
                       np.array([100.0]), np.array([5.0]),
                       np.array([90.0]), np.array([4.0]))
        with open(out) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                fuel_str = line.split()[-1]
                # Should parse as an int without raising
                int(fuel_str)


# ---------------------------------------------------------------------------
# Unit tests: assemble_landscape() with in-memory arrays
# ---------------------------------------------------------------------------

class TestAssembleLandscape(unittest.TestCase):

    def _make_transform(self, nrows, ncols, west=0.0, south=0.0, cellsize=30.0):
        """Return a simple Affine transform (no rasterio dependency)."""
        try:
            from rasterio.transform import from_bounds
            from rasterio.crs import CRS
        except ImportError:
            raise unittest.SkipTest("rasterio not installed")
        east  = west  + ncols * cellsize
        north = south + nrows * cellsize
        transform = from_bounds(west, south, east, north, ncols, nrows)
        crs = CRS.from_epsg(4326)
        return transform, crs

    def test_all_valid_pixels_returned(self):
        nrows, ncols = 4, 5
        elev   = np.full((nrows, ncols), 100.0)
        slope  = np.full((nrows, ncols), 10.0)
        aspect = np.full((nrows, ncols), 90.0)
        fuel   = np.ones((nrows, ncols), dtype=float) * 4  # FM4

        tf, crs = self._make_transform(nrows, ncols)
        xs, ys, e, s, a, f = lfp.assemble_landscape(
            elev, tf, crs, None,
            slope, tf, crs,
            aspect, tf, crs,
            fuel, tf, crs,
            project_utm=False,
            subsample=1,
        )
        self.assertEqual(len(xs), nrows * ncols)

    def test_nodata_excluded(self):
        nrows, ncols = 3, 4
        elev = np.full((nrows, ncols), 100.0)
        elev[0, 0] = -9999.0  # nodata
        slope  = np.full((nrows, ncols), 10.0)
        aspect = np.full((nrows, ncols), 90.0)
        fuel   = np.ones((nrows, ncols), dtype=float) * 4

        tf, crs = self._make_transform(nrows, ncols)
        xs, ys, e, s, a, f = lfp.assemble_landscape(
            elev, tf, crs, -9999.0,
            slope, tf, crs,
            aspect, tf, crs,
            fuel, tf, crs,
            project_utm=False,
        )
        self.assertEqual(len(xs), nrows * ncols - 1)

    def test_nonburnable_excluded_by_default(self):
        nrows, ncols = 3, 3
        elev   = np.full((nrows, ncols), 100.0)
        slope  = np.full((nrows, ncols), 10.0)
        aspect = np.full((nrows, ncols), 90.0)
        fuel   = np.array([[4, 91, 92],
                            [1, 93,  2],
                            [98, 99, 3]], dtype=float)

        tf, crs = self._make_transform(nrows, ncols)
        xs, ys, e, s, a, f = lfp.assemble_landscape(
            elev, tf, crs, None,
            slope, tf, crs,
            aspect, tf, crs,
            fuel, tf, crs,
            project_utm=False,
        )
        # Non-burnable codes: 91, 92, 93, 98, 99 → 5 excluded; valid: 4
        self.assertEqual(len(xs), 4)
        self.assertTrue(all(1 <= int(v) <= 13 for v in f))

    def test_nonburnable_kept_when_flag_set(self):
        nrows, ncols = 2, 2
        elev   = np.full((nrows, ncols), 100.0)
        slope  = np.full((nrows, ncols), 10.0)
        aspect = np.full((nrows, ncols), 90.0)
        fuel   = np.array([[4, 91], [93, 2]], dtype=float)

        tf, crs = self._make_transform(nrows, ncols)
        xs, ys, e, s, a, f = lfp.assemble_landscape(
            elev, tf, crs, None,
            slope, tf, crs,
            aspect, tf, crs,
            fuel, tf, crs,
            project_utm=False,
            keep_nonburnable=True,
        )
        # All 4 pixels returned (non-burnable not filtered)
        self.assertEqual(len(xs), 4)

    def test_subsample_reduces_points(self):
        nrows, ncols = 10, 12
        elev   = np.full((nrows, ncols), 100.0)
        slope  = np.full((nrows, ncols), 10.0)
        aspect = np.full((nrows, ncols), 90.0)
        fuel   = np.full((nrows, ncols), 4.0)

        tf, crs = self._make_transform(nrows, ncols)
        xs_full, *_ = lfp.assemble_landscape(
            elev, tf, crs, None,
            slope, tf, crs,
            aspect, tf, crs,
            fuel, tf, crs,
            project_utm=False, subsample=1,
        )
        xs_sub, *_ = lfp.assemble_landscape(
            elev, tf, crs, None,
            slope, tf, crs,
            aspect, tf, crs,
            fuel, tf, crs,
            project_utm=False, subsample=2,
        )
        self.assertLess(len(xs_sub), len(xs_full))
        self.assertGreater(len(xs_sub), 0)


# ---------------------------------------------------------------------------
# Integration tests: create_landscape_from_files()
# ---------------------------------------------------------------------------

class TestCreateLandscapeFromFiles(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def test_output_created(self):
        paths, *_ = _make_synthetic_rasters(self.tmpdir)
        out = os.path.join(self.tmpdir, "landscape.lcp")
        lfp.create_landscape_from_files(
            out,
            paths["elev"], paths["slope"], paths["aspect"], paths["fuel"],
            project_utm=False,
        )
        self.assertTrue(os.path.isfile(out))

    def test_data_rows_present(self):
        paths, *_ = _make_synthetic_rasters(self.tmpdir)
        out = os.path.join(self.tmpdir, "rows.lcp")
        lfp.create_landscape_from_files(
            out,
            paths["elev"], paths["slope"], paths["aspect"], paths["fuel"],
            project_utm=False,
        )
        with open(out) as fh:
            data_lines = [l for l in fh if not l.startswith("#") and l.strip()]
        self.assertGreater(len(data_lines), 0)

    def test_six_columns_in_output(self):
        paths, *_ = _make_synthetic_rasters(self.tmpdir)
        out = os.path.join(self.tmpdir, "6col.lcp")
        lfp.create_landscape_from_files(
            out,
            paths["elev"], paths["slope"], paths["aspect"], paths["fuel"],
            project_utm=False,
        )
        with open(out) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.split()
                self.assertEqual(len(parts), 6,
                                 f"Bad line: {line!r}")

    def test_subsample_reduces_output(self):
        paths, elev, *_ = _make_synthetic_rasters(self.tmpdir, nrows=10, ncols=10)
        out_full = os.path.join(self.tmpdir, "full.lcp")
        out_sub  = os.path.join(self.tmpdir, "sub.lcp")

        lfp.create_landscape_from_files(
            out_full,
            paths["elev"], paths["slope"], paths["aspect"], paths["fuel"],
            project_utm=False, subsample=1,
        )
        lfp.create_landscape_from_files(
            out_sub,
            paths["elev"], paths["slope"], paths["aspect"], paths["fuel"],
            project_utm=False, subsample=2,
        )

        def _count(p):
            with open(p) as fh:
                return sum(1 for l in fh if not l.startswith("#") and l.strip())

        self.assertLess(_count(out_sub), _count(out_full))

    def test_utm_projection(self):
        """When project_utm=True, coordinates should be in the metre range."""
        paths, *_ = _make_synthetic_rasters(
            self.tmpdir, west=-120.0, south=37.0, cellsize=0.0001
        )
        out = os.path.join(self.tmpdir, "utm.lcp")
        lfp.create_landscape_from_files(
            out,
            paths["elev"], paths["slope"], paths["aspect"], paths["fuel"],
            project_utm=True,
        )
        with open(out) as fh:
            first_data = next(
                (l for l in fh if not l.startswith("#") and l.strip()), None
            )
        self.assertIsNotNone(first_data)
        x_val = abs(float(first_data.split()[0]))
        # UTM x in metres should be on the order of hundreds of thousands
        self.assertGreater(x_val, 1e5)


# ---------------------------------------------------------------------------
# CLI tests
# ---------------------------------------------------------------------------

class TestCLI(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.script = os.path.join(_DEPRECATED_DIR, "landfire_to_lcp.py")

    def test_missing_bbox_exits_nonzero(self):
        """Omitting --bbox without local-file args should exit non-zero."""
        out = os.path.join(self.tmpdir, "no_bbox.lcp")
        result = subprocess.run(
            [sys.executable, self.script, out],
            capture_output=True, text=True,
        )
        self.assertNotEqual(result.returncode, 0)

    def test_bad_bbox_exits_nonzero(self):
        """MIN_LON >= MAX_LON should exit non-zero."""
        out = os.path.join(self.tmpdir, "bad_bbox.lcp")
        result = subprocess.run(
            [sys.executable, self.script,
             "--bbox", "-119.0", "37.0", "-120.0", "37.5",
             out],
            capture_output=True, text=True,
        )
        self.assertNotEqual(result.returncode, 0)

    def test_incomplete_local_files_exits_nonzero(self):
        """Specifying only --elev-file without the others should exit non-zero."""
        out = os.path.join(self.tmpdir, "partial.lcp")
        result = subprocess.run(
            [sys.executable, self.script,
             "--elev-file", "elev.tif",
             out],
            capture_output=True, text=True,
        )
        self.assertNotEqual(result.returncode, 0)

    def test_local_files_mode(self):
        """CLI local-files mode produces a valid landscape file."""
        paths, *_ = _make_synthetic_rasters(self.tmpdir)
        out = os.path.join(self.tmpdir, "cli_local.lcp")
        result = subprocess.run(
            [sys.executable, self.script,
             "--elev-file",   paths["elev"],
             "--slope-file",  paths["slope"],
             "--aspect-file", paths["aspect"],
             "--fuel-file",   paths["fuel"],
             "--no-utm",
             out],
            capture_output=True, text=True,
        )
        self.assertEqual(result.returncode, 0,
                         f"CLI failed:\n{result.stdout}\n{result.stderr}")
        self.assertTrue(os.path.isfile(out))
        with open(out) as fh:
            data_lines = [l for l in fh if not l.startswith("#") and l.strip()]
        self.assertGreater(len(data_lines), 0)

    def test_cli_subsample_flag(self):
        """--subsample flag produces fewer output rows."""
        paths, *_ = _make_synthetic_rasters(self.tmpdir, nrows=10, ncols=10)
        out_full = os.path.join(self.tmpdir, "cli_full.lcp")
        out_sub  = os.path.join(self.tmpdir, "cli_sub.lcp")

        base_args = [
            sys.executable, self.script,
            "--elev-file",   paths["elev"],
            "--slope-file",  paths["slope"],
            "--aspect-file", paths["aspect"],
            "--fuel-file",   paths["fuel"],
            "--no-utm",
        ]
        subprocess.run(base_args + [out_full], check=True, capture_output=True)
        subprocess.run(base_args + ["--subsample", "2", out_sub],
                       check=True, capture_output=True)

        def _count(p):
            with open(p) as fh:
                return sum(1 for l in fh if not l.startswith("#") and l.strip())

        self.assertLess(_count(out_sub), _count(out_full))


# ---------------------------------------------------------------------------
# Tests for internal helpers (no rasterio dependency)
# ---------------------------------------------------------------------------

class TestNonBurnableCodes(unittest.TestCase):
    """Ensure the expected set of non-burnable codes is defined."""

    def test_all_codes_present(self):
        for code in (91, 92, 93, 98, 99):
            self.assertIn(code, lfp._NONBURNABLE_CODES)

    def test_valid_fuel_not_in_set(self):
        for code in range(1, 14):
            self.assertNotIn(code, lfp._NONBURNABLE_CODES)


class TestDefaultLayerIds(unittest.TestCase):
    """Smoke-test that default layer IDs are defined for common vintages."""

    def test_2020_keys_exist(self):
        d = lfp._DEFAULT_LAYERS[2020]
        for key in ("elev", "slope", "aspect", "fuel13"):
            self.assertIn(key, d)
            self.assertIsInstance(d[key], str)
            self.assertGreater(len(d[key]), 0)

    def test_2016_keys_exist(self):
        d = lfp._DEFAULT_LAYERS[2016]
        for key in ("elev", "slope", "aspect", "fuel13"):
            self.assertIn(key, d)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    unittest.main(verbosity=2)
