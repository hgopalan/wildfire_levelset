# wildfire_levelset – Tools

Python utilities for converting geospatial data into the file formats
consumed by the wildfire level-set solver.

---

## Tools

### `terrain_wind_preprocess.py` — Unified terrain, landscape, and wind preprocessing

Downloads SRTM elevation data and LANDFIRE fuel/slope/aspect rasters for a
user-specified lat/lon bounding box and optionally extracts wind data from a
WRF-style netCDF file.  Also auto-generates an `inputs.i` file configured
for a FARSITE run with grid resolution set to 30 m.

**Outputs**

* **Terrain XYZ file** — `X Y Z` in UTM metres for `rothermel.terrain_file`
  (default: `terrain.xyz`)
* **Landscape LCP file** — `X Y ELEVATION SLOPE ASPECT FUEL_MODEL` for
  `rothermel.landscape_file` (default: `landscape.lcp`)
* **Wind CSV file** — `utm_x utm_y u v` in m/s (default: `wind.csv`;
  requires `--wrf-file`)
* **inputs.i file** — solver input file with `n_cell_x` and `n_cell_y` set
  to achieve 30 m grid resolution (default: `inputs.i`)

**Key behaviours**

* If `--wrf-file` is supplied, the bounding box is read from the WRF netCDF
  automatically; any `--lat-min/max/lon-min/max` values on the command line
  are **ignored**.
* Unless `--no-terrain` is given, SRTM elevation is used for the terrain and
  landscape files.
* If `--interpolate-wind` is given (requires `--wrf-file` and SRTM terrain),
  the WRF wind fields are interpolated to the SRTM terrain grid.
* `--time-range T1:TN` extracts WRF time steps T1 through TN **inclusive**
  (0-based indices).  When multiple time steps are requested, wind output
  filenames are derived from `--wind` by inserting `_tN` before the
  extension (e.g. `wind_t0.csv`, `wind_t1.csv`, …).
* A **local-files mode** is available: pass `--elev-file`, `--slope-file`,
  `--aspect-file`, and `--fuel-file` to use pre-existing rasters and skip
  the LFPS download entirely.

**Output formats**

```
# X Y Z (meters)                          ← terrain XYZ
x1  y1  z1
...

# X Y ELEVATION SLOPE ASPECT FUEL_MODEL   ← landscape LCP
x1  y1  elev1  slope1  aspect1  fuel1
...

# utm_x utm_y u v (m/s)                   ← wind CSV
utm_x1  utm_y1  u1  v1
...
```

**Usage**

```bash
# SRTM terrain + LANDFIRE landscape for a bounding box
python3 tools/terrain_wind_preprocess.py \
    --lat-min 40 --lat-max 40.5 \
    --lon-min -105 --lon-max -104.5

# Same but also extract wind from a WRF file (bbox comes from WRF)
python3 tools/terrain_wind_preprocess.py \
    --wrf-file wrfout_d01 \
    --wind wind.csv

# Skip SRTM and extract wind only
python3 tools/terrain_wind_preprocess.py \
    --wrf-file wrfout_d01 \
    --no-terrain \
    --wind wind.csv

# Interpolate WRF wind to SRTM resolution
python3 tools/terrain_wind_preprocess.py \
    --wrf-file wrfout_d01 \
    --wind wind.csv \
    --interpolate-wind

# Extract a range of time steps
python3 tools/terrain_wind_preprocess.py \
    --wrf-file wrfout_d01 \
    --wind wind.csv \
    --time-range 0:5

# Use local rasters instead of downloading from LANDFIRE
python3 tools/terrain_wind_preprocess.py \
    --lat-min 40 --lat-max 40.5 \
    --lon-min -105 --lon-max -104.5 \
    --elev-file elev.tif --slope-file slope.tif \
    --aspect-file aspect.tif --fuel-file fuel.tif

# Custom output filenames and skip inputs.i generation
python3 tools/terrain_wind_preprocess.py \
    --lat-min 39.5 --lat-max 40.2 \
    --lon-min -106 --lon-max -105.2 \
    --terrain region_terrain.xyz \
    --landscape region_landscape.lcp \
    --no-inputs
```

**Options**

| Flag | Description |
|------|-------------|
| `--wrf-file FILE` | WRF netCDF file; bbox is read from the file automatically |
| `--lat-min / --lat-max` | Latitude bounds (WGS-84 decimal degrees); required without `--wrf-file` |
| `--lon-min / --lon-max` | Longitude bounds (WGS-84 decimal degrees); required without `--wrf-file` |
| `--terrain FILE` | Output terrain XYZ file (default: `terrain.xyz`) |
| `--landscape FILE` | Output landscape LCP file (default: `landscape.lcp`) |
| `--wind FILE` | Output wind CSV file (default: `wind.csv`); requires `--wrf-file` |
| `--tif FILE` | Path for intermediate SRTM GeoTIFF (temporary if omitted) |
| `--no-terrain` | Skip the SRTM terrain XYZ and landscape steps |
| `--no-landscape` | Skip the LANDFIRE landscape LCP step |
| `--no-wind` | Skip the WRF wind extraction step |
| `--interpolate-wind` | Interpolate WRF wind fields to the SRTM terrain grid |
| `--time-range T1:TN` | Inclusive range of WRF time indices to extract (e.g. `0:4`) |
| `--time-index N` | Single WRF time index (default: 0) |
| `--level N` | WRF vertical level for U/V (default: 0 = lowest) |
| `--subsample N` | Keep every N-th point per dimension (default: 1) |
| `--vintage YEAR` | LANDFIRE data vintage year (default: 2020) |
| `--fuel-product ID` | Override LANDFIRE fuel model product ID |
| `--elev-product ID` | Override LANDFIRE elevation product ID |
| `--slope-product ID` | Override LANDFIRE slope product ID |
| `--aspect-product ID` | Override LANDFIRE aspect product ID |
| `--keep-nonburnable` | Include non-burnable LANDFIRE pixels in LCP output |
| `--cache-dir DIR` | Directory to cache downloaded LANDFIRE ZIP files |
| `--timeout N` | LANDFIRE API polling timeout in seconds (default: 300) |
| `--elev-file PATH` | Local elevation raster (skips LFPS download) |
| `--slope-file PATH` | Local slope raster |
| `--aspect-file PATH` | Local aspect raster |
| `--fuel-file PATH` | Local fuel model raster |
| `--inputs FILE` | Output `inputs.i` file (default: `inputs.i`) |
| `--no-inputs` | Skip automatic generation of the `inputs.i` file |

---

## Requirements

```bash
pip install numpy netCDF4 pyproj rasterio elevation requests
```

| Package | Used by |
|---------|---------|
| `numpy` | `terrain_wind_preprocess.py` |
| `netCDF4` | `terrain_wind_preprocess.py` (WRF wind extraction) |
| `pyproj` | `terrain_wind_preprocess.py` (UTM projection) |
| `rasterio` | `terrain_wind_preprocess.py` (raster I/O) |
| `elevation` | `terrain_wind_preprocess.py` (SRTM download) |
| `requests` | `terrain_wind_preprocess.py` (LANDFIRE API) |

---

## Running the Tests

```bash
# From the repository root
python3 -m pytest tools/tests/test_terrain_wind_preprocess.py -v
```

The test suite creates all necessary synthetic data files in temporary
directories; no external data downloads are required.
