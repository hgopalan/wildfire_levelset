# wildfire_levelset – Tools

Python utilities for converting geospatial data into the file formats
consumed by the wildfire level-set solver.

---

## Tools

### `dem_to_xyz.py` — Digital Elevation Map → X Y Z

Converts a Digital Elevation Map (DEM) raster into the space-separated
`X Y Z` terrain file format used by the solver (e.g. as
`rothermel.terrain_file`).

**Supported input formats**

| Format | Extension(s) | Notes |
|--------|-------------|-------|
| Arc/Info ASCII Grid | `.asc` | No extra library needed |
| SRTM HGT binary | `.hgt` | Filename must encode SW corner (e.g. `N37W120.hgt`) |
| GeoTIFF / generic raster | `.tif`, `.tiff`, `.img` | Requires `rasterio` |

**Output format** (matches solver convention)

```
# X Y Z (meters)
x1  y1  z1
x2  y2  z2
...
```

**Usage**

```bash
# Arc/Info ASCII Grid (projected, e.g. UTM)
python3 tools/dem_to_xyz.py dem/terrain.asc terrain.csv

# SRTM HGT — reproject lon/lat to UTM metres
python3 tools/dem_to_xyz.py dem/N37W120.hgt terrain.csv --project-utm

# GeoTIFF — keep every 4th point in each direction
python3 tools/dem_to_xyz.py dem/terrain.tif terrain.csv --subsample 4

# Override no-data value
python3 tools/dem_to_xyz.py dem/terrain.asc terrain.csv --nodata -9999
```

**Options**

| Flag | Description |
|------|-------------|
| `--nodata VALUE` | Override the no-data sentinel in the input |
| `--project-utm` | Reproject lon/lat → UTM metres (needs `pyproj`) |
| `--subsample N` | Keep every N-th point per dimension (default: 1) |

---

### `srtm_landfire_to_terrain.py` — Unified terrain + landscape preprocessing

Downloads SRTM elevation data and LANDFIRE fuel/slope/aspect rasters for a
user-specified lat/lon bounding box and writes:

* **Terrain XYZ file** — `X Y Z` in UTM metres for `rothermel.terrain_file`
* **Landscape LCP file** — `X Y ELEVATION SLOPE ASPECT FUEL_MODEL` for
  `rothermel.landscape_file`

Both outputs use UTM coordinates by default.  All four LANDFIRE layers
elevation, slope, aspect, fuel model) are fetched via the LANDFIRE Product
Service (LFPS) API and automatically resampled to a common grid.
Non-burnable pixels (FBFM13 codes 91–99) are excluded by default.

A **local-files mode** is also available: pass `--elev-file`, `--slope-file`,
`--aspect-file`, and `--fuel-file` to use pre-existing rasters and skip the
LFPS download entirely.

**Output formats**

```
# X Y Z (meters)                          ← terrain XYZ
x1  y1  z1
...

# X Y ELEVATION SLOPE ASPECT FUEL_MODEL   ← landscape LCP
x1  y1  elev1  slope1  aspect1  fuel1
...
```

**Usage**

```bash
# Download both terrain and landscape for a bounding box
python3 tools/srtm_landfire_to_terrain.py \
    --lat-min 40 --lat-max 40.5 \
    --lon-min -105 --lon-max -104.5

# Custom output filenames
python3 tools/srtm_landfire_to_terrain.py \
    --lat-min 39.5 --lat-max 40.2 \
    --lon-min -106 --lon-max -105.2 \
    --terrain region_terrain.xyz \
    --landscape region_landscape.lcp

# Skip the landscape step (terrain only)
python3 tools/srtm_landfire_to_terrain.py \
    --lat-min 40 --lat-max 40.5 \
    --lon-min -105 --lon-max -104.5 \
    --no-landscape

# Keep the intermediate SRTM GeoTIFF
python3 tools/srtm_landfire_to_terrain.py \
    --lat-min 40 --lat-max 40.5 \
    --lon-min -105 --lon-max -104.5 \
    --tif srtm_clip.tif

# Use local rasters instead of downloading from LANDFIRE
python3 tools/srtm_landfire_to_terrain.py \
    --lat-min 40 --lat-max 40.5 \
    --lon-min -105 --lon-max -104.5 \
    --elev-file elev.tif --slope-file slope.tif \
    --aspect-file aspect.tif --fuel-file fuel.tif

# Subsample to every 3rd point and select a specific LANDFIRE vintage
python3 tools/srtm_landfire_to_terrain.py \
    --lat-min 40 --lat-max 40.5 \
    --lon-min -105 --lon-max -104.5 \
    --subsample 3 --vintage 2016
```

**Options**

| Flag | Description |
|------|-------------|
| `--lat-min / --lat-max` | Latitude bounds (WGS-84 decimal degrees) |
| `--lon-min / --lon-max` | Longitude bounds (WGS-84 decimal degrees) |
| `--terrain FILE` | Output terrain XYZ file (default: `terrain.xyz`) |
| `--landscape FILE` | Output landscape LCP file (default: `landscape.lcp`) |
| `--tif FILE` | Path for intermediate SRTM GeoTIFF (temporary if omitted) |
| `--no-terrain` | Skip the SRTM terrain XYZ step |
| `--no-landscape` | Skip the LANDFIRE landscape LCP step |
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

---

### `wrf_to_terrain_wind.py` — WRF netCDF → terrain + wind files

Reads a WRF-style netCDF output file and produces two solver-compatible
ASCII files:

* **Terrain file** — `utm_x utm_y z`
* **Wind file** — `utm_x utm_y u v`

WRF uses a staggered Arakawa C-grid.  U and V are destaggered to mass-point
(cell-centre) locations before writing:

```
u_mass[j, i] = 0.5 * (U[j, i] + U[j, i+1])   # west-east destagger
v_mass[j, i] = 0.5 * (V[j, i] + V[j+1, i])   # south-north destagger
```

Latitude/longitude coordinates (`XLAT`, `XLONG`) are projected to UTM
metres using `pyproj`.

**WRF variables read**

| Variable | Shape | Description |
|----------|-------|-------------|
| `XLAT` | `(Time, ny, nx)` | Latitude at mass points |
| `XLONG` | `(Time, ny, nx)` | Longitude at mass points |
| `HGT_M` | `(Time, ny, nx)` | Terrain height [m] |
| `U` | `(Time, nz, ny, nx+1)` | U-wind on staggered west-east grid |
| `V` | `(Time, nz, ny+1, nx)` | V-wind on staggered south-north grid |

**Output formats**

```
# utm_x utm_y z (meters)        ← terrain file
utm_x1  utm_y1  z1
...

# utm_x utm_y u v (m/s)         ← wind file
utm_x1  utm_y1  u1  v1
...
```

**Usage**

```bash
# Basic (first time step, lowest model level)
python3 tools/wrf_to_terrain_wind.py wrfout_d01_2020-08-16_00:00:00 \
    terrain.csv wind.csv

# Select a specific time snapshot and vertical level
python3 tools/wrf_to_terrain_wind.py wrfout_d01 terrain.csv wind.csv \
    --time-index 2 --level 1

# Subsample (keep every 3rd point per dimension)
python3 tools/wrf_to_terrain_wind.py wrfout_d01 terrain.csv wind.csv \
    --subsample 3
```

**Options**

| Flag | Description |
|------|-------------|
| `--time-index N` | Time snapshot index (default: 0) |
| `--level N` | Vertical level for U/V (default: 0 = lowest) |
| `--subsample N` | Keep every N-th point per dimension (default: 1) |

---

### `utm_convert.py` — Bidirectional lat/lon ↔ UTM conversion

Pure-Python WGS-84 utility for converting between geographic (lat/lon) and
UTM (easting/northing) coordinates.  Uses `pyproj` as the backend when
available and falls back to a built-in Transverse Mercator formula otherwise
— no mandatory external dependencies.

**Module-level API**

```python
from tools.utm_convert import latlon_to_utm, utm_to_latlon

# lat/lon → UTM
easting, northing, zone_number, zone_letter = latlon_to_utm(34.10, -118.85)

# UTM → lat/lon
lat, lon = utm_to_latlon(330000, 3775000, zone_number=11, northern=True)
```

**Usage**

```bash
# Convert lat/lon → UTM (zone auto-detected)
python3 tools/utm_convert.py --to-utm 34.10 -118.85

# Convert lat/lon → UTM (force specific zone)
python3 tools/utm_convert.py --to-utm 34.10 -118.85 --zone 11

# Convert UTM → lat/lon
python3 tools/utm_convert.py --to-latlon 330000 3775000 --zone 11

# Southern hemisphere
python3 tools/utm_convert.py --to-latlon 330000 3775000 --zone 11 --south
```

**Options**

| Flag | Description |
|------|-------------|
| `--to-utm LAT LON` | Convert geographic coordinates to UTM |
| `--to-latlon EASTING NORTHING` | Convert UTM to geographic coordinates |
| `--zone N` | UTM zone number (1–60); auto-detected for `--to-utm` |
| `--south` | Treat coordinates as southern hemisphere (for `--to-latlon`) |

---

## Requirements

```bash
pip install numpy netCDF4 pyproj rasterio elevation requests
```

| Package | Used by |
|---------|---------|
| `numpy` | All tools |
| `netCDF4` | `wrf_to_terrain_wind.py` |
| `pyproj` | All tools (UTM projection; optional fallback in `utm_convert.py`) |
| `rasterio` | `dem_to_xyz.py`, `srtm_landfire_to_terrain.py` |
| `elevation` | `srtm_landfire_to_terrain.py` (SRTM download) |
| `requests` | `srtm_landfire_to_terrain.py` (LANDFIRE API) |

---

## Running the Tests

```bash
# From the repository root
python3 -m pytest tools/tests/ -v

# Or run individually
python3 tools/tests/test_dem_to_xyz.py
python3 tools/tests/test_wrf_to_terrain_wind.py
python3 tools/tests/test_landfire_to_lcp.py
```

The test suite creates all necessary synthetic data files in temporary
directories; no external data downloads are required.
