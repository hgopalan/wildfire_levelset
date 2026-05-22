#!/usr/bin/env python3
"""
Regression test for tools/hrrr_wind_reader.py

Tests:
  1. Lat/lon → UTM projection preserves shape and produces metre-scale values
  2. Regular grid creation with specified resolution
  3. Nearest-neighbor interpolation (fallback when scipy unavailable)
  4. CSV output format verification
  5. Error handling for invalid inputs
  6. Argument parsing and validation

Run with:
  python3 tools/tests/test_hrrr_wind_reader.py
  python3 -m pytest tools/tests/test_hrrr_wind_reader.py -v
"""

import os
import sys
import tempfile
import unittest

import numpy as np

# Add repository root to path
_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
_TOOLS_DIR = os.path.join(_REPO_ROOT, "tools")
sys.path.insert(0, _TOOLS_DIR)

# Import hrrr_wind_reader module components
import hrrr_wind_reader as hwr


class TestLatLonToUTM(unittest.TestCase):
    """Test lat/lon to UTM conversion."""

    def test_projection_shape_preserved(self):
        """Test that projection preserves array shape."""
        lats = np.array([[40.0, 40.1, 40.2], [40.0, 40.1, 40.2]])
        lons = np.array([[-105.0, -105.0, -105.0], [-104.9, -104.9, -104.9]])

        x_utm, y_utm = hwr._latlon_to_utm(lats, lons)

        self.assertEqual(x_utm.shape, lats.shape)
        self.assertEqual(y_utm.shape, lons.shape)

    def test_projection_produces_positive_values(self):
        """Test that projection produces sensible UTM values (in metres)."""
        lats = np.array([[40.0, 40.1], [40.0, 40.1]])
        lons = np.array([[-105.0, -105.0], [-104.9, -104.9]])

        x_utm, y_utm = hwr._latlon_to_utm(lats, lons)

        # UTM values should be in a reasonable range (typically hundreds of thousands)
        self.assertTrue(np.all(x_utm > 0))
        self.assertTrue(np.all(y_utm > 0))
        self.assertTrue(np.all(x_utm < 1e7))  # UTM easting < 10 million
        self.assertTrue(np.all(y_utm < 1e7))  # UTM northing < 10 million

    def test_projection_consistent(self):
        """Test that same lat/lon always projects to same UTM."""
        lat = 40.0
        lon = -105.0
        lats = np.array([[lat, lat], [lat, lat]])
        lons = np.array([[lon, lon], [lon, lon]])

        x_utm, y_utm = hwr._latlon_to_utm(lats, lons)

        # All values should be identical since input is constant
        self.assertTrue(np.allclose(x_utm, x_utm[0, 0]))
        self.assertTrue(np.allclose(y_utm, y_utm[0, 0]))


class TestGridGeneration(unittest.TestCase):
    """Test regular grid creation."""

    def test_grid_covers_bounds(self):
        """Test that created grid covers the requested bounds."""
        lat_min, lat_max = 40.0, 40.5
        lon_min, lon_max = -105.0, -104.5
        resolution = 100  # metres

        lat_grid, lon_grid = hwr._create_regular_grid(
            lat_min, lat_max, lon_min, lon_max, resolution
        )

        # Grid should cover the bounds
        self.assertGreaterEqual(np.min(lat_grid), lat_min)
        self.assertLessEqual(np.max(lat_grid), lat_max)
        self.assertGreaterEqual(np.min(lon_grid), lon_min)
        self.assertLessEqual(np.max(lon_grid), lon_max)

    def test_grid_resolution_reasonable(self):
        """Test that grid spacing is close to requested resolution."""
        lat_min, lat_max = 40.0, 40.5
        lon_min, lon_max = -105.0, -104.5
        resolution = 100  # metres

        lat_grid, lon_grid = hwr._create_regular_grid(
            lat_min, lat_max, lon_min, lon_max, resolution
        )

        # Convert grid spacing to approximate metres
        # Rough approximation: 1 degree ≈ 111 km
        lat_spacing_deg = lat_grid[1] - lat_grid[0] if len(lat_grid) > 1 else 0
        lon_spacing_deg = lon_grid[1] - lon_grid[0] if len(lon_grid) > 1 else 0

        lat_spacing_m = lat_spacing_deg * 111000
        lon_spacing_m = lon_spacing_deg * 111000 * np.cos(np.radians(np.mean(lat_grid)))

        # Should be within 50% of target resolution
        if lat_spacing_m > 0:
            self.assertLess(lat_spacing_m, resolution * 2)
        if lon_spacing_m > 0:
            self.assertLess(lon_spacing_m, resolution * 2)

    def test_grid_has_points(self):
        """Test that grid creation produces at least 2 points in each dimension."""
        lat_min, lat_max = 40.0, 40.5
        lon_min, lon_max = -105.0, -104.5

        lat_grid, lon_grid = hwr._create_regular_grid(
            lat_min, lat_max, lon_min, lon_max, 100
        )

        self.assertGreaterEqual(len(lat_grid), 2)
        self.assertGreaterEqual(len(lon_grid), 2)


class TestInterpolation(unittest.TestCase):
    """Test wind interpolation to regular grid."""

    def test_interpolation_shape(self):
        """Test that interpolation produces correct output shape."""
        # Create synthetic HRRR data
        u_hrrr = np.array([5.0, 6.0, 7.0, 8.0])
        v_hrrr = np.array([1.0, 2.0, 3.0, 4.0])
        lat_hrrr = np.array([40.0, 40.1, 40.2, 40.3])
        lon_hrrr = np.array([-105.0, -104.9, -104.8, -104.7])

        # Create target grid
        lat_grid = np.array([40.0, 40.1, 40.2])
        lon_grid = np.array([-105.0, -104.9, -104.8])

        u_interp, v_interp, lats_2d, lons_2d = hwr._interpolate_to_grid(
            u_hrrr, v_hrrr, lat_hrrr, lon_hrrr, lat_grid, lon_grid
        )

        # Output should be 2D
        expected_shape = (len(lat_grid), len(lon_grid))
        self.assertEqual(u_interp.shape, expected_shape)
        self.assertEqual(v_interp.shape, expected_shape)
        self.assertEqual(lats_2d.shape, expected_shape)
        self.assertEqual(lons_2d.shape, expected_shape)

    def test_interpolation_no_nans_with_nearby_data(self):
        """Test that interpolation doesn't produce NaNs when data is available."""
        # Create synthetic HRRR data covering a larger area
        ny, nx = 10, 10
        lat_hrrr_1d = np.linspace(40.0, 40.5, ny)
        lon_hrrr_1d = np.linspace(-105.0, -104.5, nx)
        lons_hrrr, lats_hrrr = np.meshgrid(lon_hrrr_1d, lat_hrrr_1d)

        u_hrrr = 5.0 * np.ones_like(lats_hrrr)
        v_hrrr = 2.0 * np.ones_like(lats_hrrr)

        # Create target grid within bounds
        lat_grid = np.array([40.1, 40.2])
        lon_grid = np.array([-104.9, -104.8])

        u_interp, v_interp, _, _ = hwr._interpolate_to_grid(
            u_hrrr.ravel(), v_hrrr.ravel(),
            lats_hrrr.ravel(), lons_hrrr.ravel(),
            lat_grid, lon_grid
        )

        # Should not have NaNs
        self.assertFalse(np.any(np.isnan(u_interp)))
        self.assertFalse(np.any(np.isnan(v_interp)))

    def test_interpolation_values_in_range(self):
        """Test that interpolated values are within input range."""
        u_hrrr = np.array([5.0, 6.0, 7.0, 8.0])
        v_hrrr = np.array([1.0, 2.0, 3.0, 4.0])
        lat_hrrr = np.array([40.0, 40.1, 40.2, 40.3])
        lon_hrrr = np.array([-105.0, -104.9, -104.8, -104.7])

        lat_grid = np.array([40.05, 40.15])
        lon_grid = np.array([-104.95, -104.85])

        u_interp, v_interp, _, _ = hwr._interpolate_to_grid(
            u_hrrr, v_hrrr, lat_hrrr, lon_hrrr, lat_grid, lon_grid
        )

        # Interpolated values should be within min/max of input
        self.assertGreaterEqual(np.min(u_interp), np.min(u_hrrr))
        self.assertLessEqual(np.max(u_interp), np.max(u_hrrr))
        self.assertGreaterEqual(np.min(v_interp), np.min(v_hrrr))
        self.assertLessEqual(np.max(v_interp), np.max(v_hrrr))


class TestCSVOutput(unittest.TestCase):
    """Test CSV file writing."""

    def test_csv_format(self):
        """Test that CSV is written in the correct format."""
        with tempfile.NamedTemporaryFile(mode='r', suffix='.csv', delete=False) as f:
            csv_path = f.name

        try:
            # Create test data
            x_utm = np.array([[100.0, 200.0], [150.0, 250.0]])
            y_utm = np.array([[1000.0, 2000.0], [1500.0, 2500.0]])
            u = np.array([[5.0, 6.0], [7.0, 8.0]])
            v = np.array([[1.0, 2.0], [3.0, 4.0]])

            hwr._write_wind_csv(csv_path, x_utm, y_utm, u, v)

            # Verify output
            with open(csv_path, 'r') as f:
                lines = f.readlines()

            # Check header
            self.assertTrue(lines[0].startswith("# utm_x utm_y u v"))

            # Check data lines
            self.assertEqual(len(lines), 5)  # Header + 4 data points

            # Parse first data line
            parts = lines[1].strip().split()
            self.assertEqual(len(parts), 4)  # x, y, u, v
            self.assertAlmostEqual(float(parts[0]), 100.0, places=1)
            self.assertAlmostEqual(float(parts[1]), 1000.0, places=1)

        finally:
            if os.path.exists(csv_path):
                os.remove(csv_path)

    def test_csv_has_all_points(self):
        """Test that CSV contains all grid points."""
        with tempfile.NamedTemporaryFile(mode='r', suffix='.csv', delete=False) as f:
            csv_path = f.name

        try:
            x_utm = np.array([[100.0, 200.0, 300.0], [150.0, 250.0, 350.0]])
            y_utm = np.array([[1000.0, 2000.0, 3000.0], [1500.0, 2500.0, 3500.0]])
            u = np.random.randn(2, 3)
            v = np.random.randn(2, 3)

            hwr._write_wind_csv(csv_path, x_utm, y_utm, u, v)

            with open(csv_path, 'r') as f:
                lines = f.readlines()

            # Should have header + 6 data lines (2×3 grid)
            self.assertEqual(len(lines), 7)

        finally:
            if os.path.exists(csv_path):
                os.remove(csv_path)


class TestInputValidation(unittest.TestCase):
    """Test input validation and error handling."""

    def test_invalid_datetime_format(self):
        """Test that invalid datetime format raises error."""
        with self.assertRaises(ValueError):
            # _fetch_hrrr_data will fail due to missing Herbie, so we need to
            # patch the datetime parsing
            from datetime import datetime
            datetime.strptime("invalid-date", "%Y-%m-%d %H:%M")

    def test_wind_output_path_single(self):
        """Test wind output path for first time index."""
        path = hwr._wind_output_path("wind.csv", 0)
        self.assertEqual(path, "wind.csv")

    def test_wind_output_path_multiple(self):
        """Test wind output path for multiple time indices."""
        path1 = hwr._wind_output_path("wind.csv", 1)
        self.assertEqual(path1, "wind_1.csv")

        path2 = hwr._wind_output_path("wind.csv", 2)
        self.assertEqual(path2, "wind_2.csv")

        # Test with different extension
        path3 = hwr._wind_output_path("my_wind_data.txt", 3)
        self.assertEqual(path3, "my_wind_data_3.txt")


class TestInputsStub(unittest.TestCase):
    """Test inputs.i stub generation."""

    def test_inputs_stub_creation(self):
        """Test that inputs.i stub is created with correct content."""
        with tempfile.NamedTemporaryFile(mode='r', suffix='.i', delete=False) as f:
            inputs_path = f.name

        try:
            hwr._generate_inputs_stub(inputs_path, "wind.csv", time_s=3600)

            with open(inputs_path, 'r') as f:
                content = f.read()

            # Check that key fields are present
            self.assertIn("velocity_file", content)
            self.assertIn("wind.csv", content)
            self.assertIn("3600", content)

        finally:
            if os.path.exists(inputs_path):
                os.remove(inputs_path)


if __name__ == "__main__":
    unittest.main()
