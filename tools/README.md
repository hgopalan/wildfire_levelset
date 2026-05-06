# wildfire_levelset – Tools

Python utilities for converting geospatial data into the file formats
consumed by the wildfire level-set solver, and for post-processing simulation
results.

> **Deprecated tool**: `tools/deprecated/terrain_wind_preprocess.py` is a
> legacy unified script superseded by the individual tools below.  Use the
> specific tools listed here for all new workflows.

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
contours.  Multi-level AMR plotfiles are supported.  In addition to the stored
plotfile fields, **derived fields** are computed on-the-fly:

| Derived field    | Source fields    | Description                             |
|------------------|------------------|-----------------------------------------|
| `wind_speed`     | `velx`, `vely`   | Wind speed [m/s] (magnitude)            |
| `wind_direction` | `velx`, `vely`   | Met. wind direction [° from N, cw]      |

```bash
# Export all fire-behaviour variables (including wind_speed, wind_direction)
python3 tools/plotfile_to_geotiff.py plt0100 --outdir gis_out

# Export specific derived variables
python3 tools/plotfile_to_geotiff.py plt0100 --outdir gis_out \
    -v wind_speed wind_direction reaction_intensity

# Batch-convert all plt#### directories
python3 tools/plotfile_to_geotiff.py --all --outdir gis_out \
    --utm-origin 450000 4200000 --epsg 32613
```

**Dependencies**: `rasterio`, `numpy`; `matplotlib` for GeoJSON perimeters.

---

### `ensemble_burn_probability.py` — Ensemble Burn Probability Driver

FSPro-style ensemble driver.  Runs the solver *N* times with Latin hypercube
or random perturbations of wind speed, wind direction, and dead fuel moisture,
accumulates per-cell burn counts, and writes a probability map.

Optionally computes **conditional flame length exceedance** maps — the
probability that flame length exceeds a threshold M at each cell — using the
`flame_length` field from the solver plotfiles.

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

# Flame length exceedance at 0.5, 1.0, 2.0, 4.0 m thresholds
# (requires plot_int > 0 in inputs.i)
python3 tools/ensemble_burn_probability.py \
    --exe ./wildfire_levelset \
    --inputs inputs.i \
    --n-runs 50 \
    --fl-thresholds 0.5 1.0 2.0 4.0 \
    --fl-out-prefix fl_exceedance
```

**Outputs**:
- `burn_probability.csv` — P_burn per cell
- `fl_exceedance_0p5m.csv`, `fl_exceedance_1p0m.csv`, … — P(FL > M) per cell

**Dependencies**: none required; `scipy` for Latin hypercube sampling (falls back to pure-Python LHS otherwise).

---

### `plot_burn_probability.py` — Burn Probability Visualisation

Reads `burn_probability.csv` (or flame-length exceedance CSVs) from the
ensemble driver and produces heatmap plots.  Supports single-file, multi-file
panel layouts, and side-by-side difference maps.

```bash
# Basic heatmap
python3 tools/plot_burn_probability.py burn_probability.csv

# Custom colourmap and contour lines at 10%, 25%, 50%, 75%
python3 tools/plot_burn_probability.py burn_probability.csv \
    --cmap YlOrRd --vmin 0.1 --contours 0.10 0.25 0.50 0.75 \
    --out bp_map.png

# Compare two ensemble runs side by side with a difference panel
python3 tools/plot_burn_probability.py \
    --diff run_A/burn_probability.csv run_B/burn_probability.csv \
    --labels "Run A (dry)" "Run B (wet)" --out diff_map.png

# Plot multiple files in a panel row
python3 tools/plot_burn_probability.py \
    fl_exceedance_0p5m.csv fl_exceedance_2p0m.csv \
    --labels "P(FL>0.5 m)" "P(FL>2.0 m)" --out fl_panels.png

# Export as GeoTIFF (requires rasterio)
python3 tools/plot_burn_probability.py burn_probability.csv \
    --geotiff burn_prob.tif --epsg 32613
```

**Dependencies**: `matplotlib`, `numpy`; `rasterio` for `--geotiff` export.

---

### `isochrone_extractor.py` — Fire Arrival-Time Isochrones

Reads the `arrival_time` field from AMReX plotfiles and extracts fire
time-of-arrival isochrones as GeoJSON polygon features.  Each isochrone
marks the fire perimeter at a specific elapsed simulation time.

```bash
# 10-minute isochrones from a single plotfile
python3 tools/isochrone_extractor.py plt0100 --interval 600 --outdir iso_out

# Custom times [s]
python3 tools/isochrone_extractor.py plt0100 \
    --times 300 600 900 1200 1800 --outdir iso_out

# Batch all plt#### with UTM georeference
python3 tools/isochrone_extractor.py --all --interval 600 \
    --utm-origin 450000 4200000 --outdir iso_out
```

**Output**: `<outdir>/<plotfile>_isochrones.geojson` — GeoJSON FeatureCollection
of isochrone polygons with properties `time_s`, `time_min`, `label`.

**Dependencies**: `matplotlib`, `numpy`; imports the plotfile parser from
`plotfile_to_geotiff.py` when available.

---

### `fire_size_summary.py` — Fire Size Statistics

Reads the `fire_stats.csv` time series written by the solver (enabled with
`fire_stats_file = fire_stats.csv` in `inputs.i`) and prints a formatted ASCII
table of burned area, perimeter, active front cells, and emissions.  Optionally
saves matplotlib plots.

```bash
# Print ASCII table
python3 tools/fire_size_summary.py fire_stats.csv

# Print table and save PNG plots
python3 tools/fire_size_summary.py fire_stats.csv --plot --outdir fire_plots

# Export enriched CSV with time_min column
python3 tools/fire_size_summary.py fire_stats.csv --csv summary_min.csv
```

**Solver prerequisite**: set `fire_stats_file = fire_stats.csv` in `inputs.i`.

**Dependencies**: none required; `matplotlib` for `--plot`.

---

### `behavior_matrix.py` — Fuel Condition Fire Behavior Matrix

BehavePlus-style Rothermel (1972) fire behavior matrix tool.  Computes rate of
spread, fireline intensity, flame length, and reaction intensity across a grid
of wind speeds and dead fuel moisture values for any Anderson FBFM13 or Scott &
Burgan FBFM40 fuel model.  Pure Python — no solver required.

```bash
# Anderson FM 4 (chaparral), winds 0–10 m/s, moisture 4–20%
python3 tools/behavior_matrix.py --fuel-model 4 --fuel-system 13 \
    --wind-min 0 --wind-max 10 --wind-steps 11 \
    --moisture-min 0.04 --moisture-max 0.20 --moisture-steps 9 \
    --out fm4_matrix.csv

# Scott & Burgan SH5 (High load dry shrub) with slope correction
python3 tools/behavior_matrix.py --fuel-model 145 --fuel-system 40 \
    --slope 0.4 --out sh5_matrix.csv

# Save heatmap plots
python3 tools/behavior_matrix.py --fuel-model 4 --plot \
    --plot-out fm4_heatmap.png

# List available fuel model codes
python3 tools/behavior_matrix.py --list-fuels --fuel-system 40
```

**Outputs**: CSV matrix with columns `wind_m_s`, `moisture_pct`, `R_ros_m_min`,
`I_R_kW_m2`, `I_B_kW_m`, `L_f_m`, `phi_w`, `phi_s`.

**Dependencies**: none required (pure Python); `matplotlib`, `numpy` for `--plot`.

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

Optionally computes **conditional flame length exceedance** maps — the
probability that flame length exceeds a threshold M at each cell — using the
`flame_length` field from the solver plotfiles.

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

# Flame length exceedance at 0.5, 1.0, 2.0, 4.0 m thresholds
# (requires plot_int > 0 in inputs.i)
python3 tools/ensemble_burn_probability.py \
    --exe ./wildfire_levelset \
    --inputs inputs.i \
    --n-runs 50 \
    --fl-thresholds 0.5 1.0 2.0 4.0 \
    --fl-out-prefix fl_exceedance
```

**Outputs**:
- `burn_probability.csv` — P_burn per cell
- `fl_exceedance_0p5m.csv`, `fl_exceedance_1p0m.csv`, … — P(FL > M) per cell

**Dependencies**: none required; `scipy` for Latin hypercube sampling (falls back to pure-Python LHS otherwise).

---

## Requirements

Install all optional Python dependencies with:

```bash
pip install numpy netCDF4 pyproj rasterio elevation requests scipy matplotlib shapely pyshp
```

| Package | Used by |
|---------|---------|
| `numpy` | `wrf_wind_reader.py`, `srtm_terrain_reader.py`, `landscape_writer.py`, `plotfile_to_geotiff.py`, `plot_burn_probability.py`, `behavior_matrix.py` (optional) |
| `netCDF4` | `wrf_wind_reader.py` (WRF wind extraction) |
| `pyproj` | `wrf_wind_reader.py`, `srtm_terrain_reader.py`, `landscape_writer.py`, `utm_convert.py`, `perimeter_to_shapefile.py` |
| `rasterio` | `srtm_terrain_reader.py`, `landscape_writer.py`, `plotfile_to_geotiff.py`, `plot_burn_probability.py` (optional) |
| `elevation` | `srtm_terrain_reader.py` (SRTM download) |
| `requests` | `landscape_writer.py` (LANDFIRE API) |
| `scipy` | `wrf_wind_reader.py` (IDW interpolation), `ensemble_burn_probability.py` (LHS) |
| `matplotlib` | `plotfile_to_geotiff.py`, `isochrone_extractor.py`, `plot_burn_probability.py`, `fire_size_summary.py` (all optional) |
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
