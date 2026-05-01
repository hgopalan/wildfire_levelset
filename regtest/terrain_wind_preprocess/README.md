# terrain_wind_preprocess Regression Test

This regression test exercises all major capabilities of
`tools/terrain_wind_preprocess.py`:

| Test | Description |
|------|-------------|
| `test_wrf_wind_extraction` | Extract wind from a synthetic WRF netCDF → checks 4-column output, UTM coordinates |
| `test_bbox_from_wrf` | Bounding box is read from WRF file; `--lat-min/max/lon-min/max` are ignored |
| `test_no_terrain_flag` | `--no-terrain` skips terrain/landscape, wind still produced |
| `test_time_range` | `--time-range 0:2` writes `wind_t0.csv`, `wind_t1.csv`, `wind_t2.csv` |
| `test_time_index` | `--time-index 1` writes a single `wind.csv` from time step 1 |
| `test_interpolate_wind` | `--interpolate-wind` places wind values on the (mocked) SRTM terrain grid |
| `test_srtm_terrain_ignored_hgt_m` | Terrain Z values come from SRTM, not WRF HGT_M |
| `test_landscape_local_files` | Landscape created from local synthetic rasters (no LANDFIRE download) |
| `test_subsample` | `--subsample 2` reduces output point count |

## Running

```bash
# From the repository root
python3 regtest/terrain_wind_preprocess/run_regtest.py

# Or via pytest
python3 -m pytest regtest/terrain_wind_preprocess/run_regtest.py -v
```

All tests use synthetic data (synthetic WRF netCDF files and GeoTIFFs) so
they can run without an internet connection and without real WRF output or
SRTM/LANDFIRE downloads.

## Expected outputs

The `run_regtest.py` script generates all data in a temporary directory.  
Reference wind and terrain snapshots are embedded directly in the test script
and compared numerically with the tool output.
