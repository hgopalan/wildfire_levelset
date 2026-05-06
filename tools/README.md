# wildfire_levelset – Tools

Python utilities for converting geospatial data into the file formats
consumed by the wildfire level-set solver, and for post-processing simulation
results.

---

## Tools

### `wrf_wind_reader.py` — WRF Wind Extraction

Reads U and V 10-m (or lowest-level) wind components from a WRF-style netCDF
file and writes one or more CSV wind files for use as `velocity_file` in the
solver.

```bash
# Single time step
python3 tools/wrf_wind_reader.py --wrf-file wrfout_d01 --wind wind.csv

# Multiple time steps
python3 tools/wrf_wind_reader.py --wrf-file wrfout_d01 --wind wind.csv --time-range 0:4

# Interpolate wind onto an existing terrain grid
python3 tools/wrf_wind_reader.py --wrf-file wrfout_d01 --wind wind.csv \
    --terrain-file terrain.xyz
```

**Dependencies**: `netCDF4`, `numpy`, `pyproj`; `scipy` for `--terrain-file` or `--grid-resolution`.

---

### `srtm_terrain_reader.py` — SRTM Terrain Download

Downloads SRTM1 (30 m) elevation data for a lat/lon bounding box and writes a
UTM terrain XYZ file for `rothermel.terrain_file`.

```bash
python3 tools/srtm_terrain_reader.py \
    --lat-min 40 --lat-max 40.5 --lon-min -105 --lon-max -104.5
```

**Dependencies**: `elevation`, `rasterio`, `numpy`, `pyproj`, `scipy`.

---

### `landscape_writer.py` — LANDFIRE Landscape File Writer

Downloads LANDFIRE elevation, slope, aspect, and fuel model rasters for a
lat/lon bounding box and writes a FARSITE landscape file (`.lcp`) for
`rothermel.landscape_file`.

```bash
python3 tools/landscape_writer.py \
    --lat-min 40 --lat-max 40.5 --lon-min -105 --lon-max -104.5
```

**Dependencies**: `rasterio`, `numpy`, `pyproj`, `requests`.

---

### `farsite_weather_reader.py` — FARSITE Weather File Parser

Parses FARSITE `.wtr` RAWS weather-station files and writes time-stamped wind
CSV files for the solver.  Optionally generates `inputs.i` stubs and fuel
moisture files.

```bash
python3 tools/farsite_weather_reader.py --wtr fire.wtr --wind wind.csv
```

**Dependencies**: none (pure Python).

---

### `fuel_moisture_from_weather.py` — Equilibrium Fuel Moisture Calculator

Estimates 1-hr, 10-hr, and 100-hr dead fuel moisture from temperature and
relative humidity using the Nelson (2000) / Simard (1968) EMC model.

```bash
# Single observation
python3 tools/fuel_moisture_from_weather.py --temp 30 --rh 20

# Batch from a weather CSV
python3 tools/fuel_moisture_from_weather.py --csv weather.csv --solver-format
```

**Dependencies**: none (pure Python).

---

### `farsite_fms_reader.py` — FARSITE Fuel Moisture Scenario Parser

Reads a FARSITE `.fms` file specifying per-fuel-model dead and live fuel
moisture contents and writes solver-compatible moisture inputs.

```bash
# Print moisture summary and solver input blocks
python3 tools/farsite_fms_reader.py fire.fms

# Write a solver inputs.i moisture stub for fuel model 4
python3 tools/farsite_fms_reader.py fire.fms --fuel-model 4 \
    --write-inputs inputs_moisture.i

# Export all moisture values as CSV
python3 tools/farsite_fms_reader.py fire.fms --write-csv moisture.csv

# Query a single fuel model
python3 tools/farsite_fms_reader.py fire.fms --query 8
```

**Dependencies**: none (pure Python).

---

### `farsite_adj_reader.py` — FARSITE Fuel Adjustment File Parser

Reads per-fuel-model ROS multipliers from a FARSITE `.adj` file, generates
template files, and can apply adjustments to a Rothermel parameter CSV.

```bash
# Inspect an existing .adj file
python3 tools/farsite_adj_reader.py --adj fire.adj

# Generate a template for all FBFM13 models
python3 tools/farsite_adj_reader.py --generate --fuel-system 13 \
    --output template_13.adj

# Apply adjustments to a Rothermel CSV for calibration inspection
python3 tools/farsite_adj_reader.py --adj fire.adj \
    --apply-csv rothermel_params.csv \
    --output rothermel_adjusted.csv
```

**Dependencies**: none (pure Python).

---

### `farsite_fmd_reader.py` — FARSITE Fuel Moisture Schedule Parser

Reads, converts, queries, and generates FARSITE-compatible `.fmd` fuel
moisture schedule files.  The solver reads these via `fmd_file` and updates
moisture at every time step.

```bash
# Inspect a .fmd file
python3 tools/farsite_fmd_reader.py --fmd fire.fmd

# Convert to a flat CSV
python3 tools/farsite_fmd_reader.py --fmd fire.fmd --output moisture.csv

# Query interpolated moisture at t = 3600 s for FM4
python3 tools/farsite_fmd_reader.py --fmd fire.fmd \
    --query-time 3600 --fuel-model 4

# Generate a 24-hour constant-moisture template
python3 tools/farsite_fmd_reader.py --generate \
    --models 4 9 10 --start-month 7 --start-day 15 --hours 24 \
    --M-d1 8 --M-d10 10 --M-d100 15 --M-lh 90 --M-lw 120 \
    --output summer_moisture.fmd
```

**Dependencies**: none (pure Python).

---

### `plotfile_to_geotiff.py` — GIS Export

Converts AMReX 2-D plotfiles to GeoTIFF rasters and GeoJSON fire-perimeter
contours.  Multi-level AMR plotfiles are supported.

```bash
# Export all fire-behaviour variables
python3 tools/plotfile_to_geotiff.py plt0100 --outdir gis_out

# Batch-convert all plt#### directories
python3 tools/plotfile_to_geotiff.py --all --outdir gis_out \
    --utm-origin 450000 4200000 --epsg 32613
```

**Dependencies**: `rasterio`, `numpy`; `matplotlib` for GeoJSON perimeters.

---

### `perimeter_to_shapefile.py` — Fire Perimeter to Esri Shapefile

Converts wildfire_levelset fire perimeter output files (`.geojson`, `.csv`,
`.dat`) to Esri Shapefile format for import into GIS applications.

```bash
# Convert a single GeoJSON perimeter
python3 tools/perimeter_to_shapefile.py perimeter_0100.geojson

# Convert all GeoJSON perimeters with UTM zone 13N CRS
python3 tools/perimeter_to_shapefile.py perimeter_*.geojson --epsg 32613

# Batch convert everything in the current directory
python3 tools/perimeter_to_shapefile.py --all --outdir shapefiles/
```

**Dependencies**: `pyshp` (required); `shapely` for convex hull; `pyproj` for CRS WKT.

```bash
pip install pyshp shapely pyproj
```

---

### `utm_convert.py` — Lat/Lon ↔ UTM Coordinate Conversion

Bidirectional WGS-84 lat/lon ↔ UTM conversion.  Works as a standalone CLI
tool or as an importable Python module.

```bash
# Convert lat/lon → UTM (zone auto-detected)
python3 tools/utm_convert.py --to-utm 34.10 -118.85

# Convert UTM → lat/lon
python3 tools/utm_convert.py --to-latlon 330000 3775000 --zone 11
```

As a module:

```python
from tools.utm_convert import latlon_to_utm, utm_to_latlon

easting, northing, zone_number, zone_letter = latlon_to_utm(34.10, -118.85)
lat, lon = utm_to_latlon(330000, 3775000, zone_number=11, northern=True)
```

**Dependencies**: none (pure Python); `pyproj` optional for higher accuracy.

---

### `ensemble_burn_probability.py` — Ensemble Burn Probability Driver

FSPro-style ensemble driver.  Runs the solver *N* times with Latin hypercube
or random perturbations of wind speed, wind direction, and dead fuel moisture,
accumulates per-cell burn counts, and writes a probability map.

```bash
# Serial: 50 runs, ±20% wind speed, ±15° direction, ±2% moisture
python3 tools/ensemble_burn_probability.py \
    --exe ./wildfire_levelset \
    --inputs inputs.i \
    --n-runs 50 \
    --wind-speed-sigma 0.20 \
    --wind-dir-sigma 15.0 \
    --moisture-sigma 0.02

# MPI: each ensemble member uses 4 MPI ranks
python3 tools/ensemble_burn_probability.py \
    --exe ./wildfire_levelset \
    --inputs inputs.i \
    --n-runs 50 \
    --mpi-ranks 4
```

**Dependencies**: none required; `scipy` for Latin hypercube sampling (falls back to pure-Python LHS otherwise).

---

## Requirements

Install all optional Python dependencies with:

```bash
pip install numpy netCDF4 pyproj rasterio elevation requests scipy matplotlib shapely pyshp
```

| Package | Used by |
|---------|---------|
| `numpy` | `wrf_wind_reader.py`, `srtm_terrain_reader.py`, `landscape_writer.py`, `plotfile_to_geotiff.py` |
| `netCDF4` | `wrf_wind_reader.py` (WRF wind extraction) |
| `pyproj` | `wrf_wind_reader.py`, `srtm_terrain_reader.py`, `landscape_writer.py`, `utm_convert.py`, `perimeter_to_shapefile.py` |
| `rasterio` | `srtm_terrain_reader.py`, `landscape_writer.py`, `plotfile_to_geotiff.py` |
| `elevation` | `srtm_terrain_reader.py` (SRTM download) |
| `requests` | `landscape_writer.py` (LANDFIRE API) |
| `scipy` | `wrf_wind_reader.py` (IDW interpolation), `ensemble_burn_probability.py` (LHS) |
| `matplotlib` | `plotfile_to_geotiff.py` (GeoJSON perimeter contours) |
| `shapely` | `perimeter_to_shapefile.py` (convex hull) |
| `pyshp` | `perimeter_to_shapefile.py` (shapefile writing, **required**) |

---

## Deprecated

`tools/deprecated/` contains `terrain_wind_preprocess.py` — a unified legacy
script that previously combined terrain download, landscape generation, and WRF
wind extraction into a single tool.  It is superseded by the individual split
tools above and is retained only for reference.

---

## Running the Tests

```bash
# From the repository root
python3 -m pytest tools/tests/ -v
```

The test suite creates all necessary synthetic data files in temporary
directories; no external data downloads are required.
