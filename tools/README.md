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

## Requirements

```bash
pip install numpy netCDF4 pyproj rasterio
```

| Package | Used by |
|---------|---------|
| `numpy` | Both tools |
| `netCDF4` | `wrf_to_terrain_wind.py` |
| `pyproj` | Both tools (`--project-utm` / UTM projection) |
| `rasterio` | `dem_to_xyz.py` (GeoTIFF reading) |

---

## Running the Tests

```bash
# From the repository root
python3 -m pytest tools/tests/ -v

# Or run individually
python3 tools/tests/test_dem_to_xyz.py
python3 tools/tests/test_wrf_to_terrain_wind.py
```

The test suite creates all necessary synthetic data files in temporary
directories; no external data downloads are required.
