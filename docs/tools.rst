Python Tools
============

The ``tools/`` directory contains Python utilities for converting geospatial data into the
file formats consumed by the wildfire level-set solver, as well as post-processing tools for
exporting simulation results to GIS formats.

Overview
--------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Script
     - Purpose
   * - ``wrf_wind_reader.py``
     - Extract U/V wind from a WRF netCDF file → solver wind CSV
   * - ``srtm_terrain_reader.py``
     - Download SRTM1 elevation data → UTM terrain XYZ file
   * - ``landscape_writer.py``
     - Download LANDFIRE rasters → FARSITE landscape file (``.lcp``)
   * - ``farsite_weather_reader.py``
     - Parse FARSITE ``.wtr`` RAWS weather files → time-stamped wind CSVs
   * - ``fuel_moisture_from_weather.py``
     - Estimate equilibrium dead-fuel moisture from temperature and RH
   * - ``farsite_fms_reader.py``
     - Parse FARSITE ``.fms`` fuel moisture scenario files → solver inputs
   * - ``farsite_adj_reader.py``
     - Inspect, generate, and apply FARSITE ``.adj`` fuel adjustment files
   * - ``farsite_fmd_reader.py``
     - Parse, query, convert, and generate FARSITE ``.fmd`` fuel moisture schedules
   * - ``plotfile_to_geotiff.py``
     - Convert AMReX plotfiles → GeoTIFF rasters and GeoJSON fire-perimeter contours
   * - ``perimeter_to_shapefile.py``
     - Convert fire perimeter files (``.geojson``, ``.csv``, ``.dat``) to Esri Shapefile
   * - ``utm_convert.py``
     - Bidirectional lat/lon ↔ UTM coordinate conversion (pure-Python or pyproj)
   * - ``ensemble_burn_probability.py``
     - FSPro-style ensemble burn probability driver (Latin hypercube sampling, MPI)
   * - ``plot_burn_probability.py``
     - Visualise ensemble burn probability maps (heatmaps, contours, difference maps)
   * - ``fire_size_summary.py``
     - Tabulate and plot fire size statistics over time from ``fire_stats.csv``;
       **NEW:** percentile statistics (10th/50th/90th) for FARSITE-style growth analysis
   * - ``isochrone_extractor.py``
     - Extract fire arrival-time isochrones from AMReX plotfiles → GeoJSON;
       **NEW:** matplotlib visualization with time labels overlaid on arrival time heatmap
   * - ``minimum_travel_path.py``
     - **NEW:** Extract minimum travel time (MTT) paths from ignition to destinations;
       follows steepest descent of arrival_time field (FARSITE path analysis)
   * - ``fire_period_analysis.py``
     - **NEW:** Classify burned cells as day/night based on burn period settings;
       spatial map and statistics (FARSITE/FSPro burn period concept)
   * - ``surface_fire_worksheet.py``
     - BehavePlus-style comprehensive single-point surface fire behavior worksheet
       (Rothermel 1972): reaction intensity, ROS, flame length, fire size, spotting,
       probability of ignition, crown fire assessment; batch CSV mode supported.
   * - ``values_at_risk.py``
     - FSPro-style values-at-risk (VAR) overlay: match asset inventory to burn
       probability map, compute expected loss per asset and category, produce
       prioritised at-risk list and optional Monte-Carlo exceedance curve.
   * - ``topographic_horizon_analysis.py``
     - Compute and visualise FARSITE topographic horizon angles from an XYZ
       terrain CSV or a built-in synthetic canyon profile.  Identifies which
       cells are in ridge shadow for a given solar elevation / azimuth.
       Also serves as executable documentation of the C++ horizon-scan algorithm.
   * - ``behavior_matrix.py``
     - Generate Rothermel fire behavior matrices across moisture / wind ranges
   * - ``crown_fire_worksheet.py``
     - Van Wagner (1977) crown fire initiation and active-crown ROS worksheet
   * - ``ignition_probability_table.py``
     - Anderson (1970) probability of ignition lookup tables
   * - ``satellite_goes_to_csv.py`` / Satellite Fire Detection Assimilation
     - Convert GOES NetCDF fire-detection granules to CSV; ingest GOES / VIIRS / CSV active-fire detections into the level-set initial condition or mid-simulation
   * - ``ember_cascade_analysis.py``
     - Post-process ``ember_cascade_flux`` / ``ember_cascade_ignition`` fields from
       AMReX plotfiles: summary statistics table, CSV export, and landing-flux
       density map visualisation (requires ``numpy``; ``matplotlib`` for plots)
   * - ``historical_wildfires.py``
     - Query a curated database of 29 major US historical wildfires (2009–2024);
       display tabular summary; export to CSV; generate ``terrain.csv`` and
       ``inputs.i`` files for a selected fire via automatic SRTM elevation
       download with configurable lat/lon bounds (±0.25° default).

A unified legacy tool (``terrain_wind_preprocess.py``) is retained in ``tools/deprecated/``
and superseded by the split tools above.

Dependencies
------------

Install with pip before using any tool:

.. code-block:: bash

   pip install numpy netCDF4 pyproj rasterio elevation requests scipy matplotlib shapely pyshp

Not all packages are required by every tool; see each section below.

``wrf_wind_reader.py`` — WRF Wind Extraction
----------------------------------------------

Reads U and V 10-m (or lowest-level) wind components from a WRF-style netCDF file and
writes one or more CSV wind files for use as ``velocity_file`` in the solver.

**Key workflow**

1. Opens the WRF netCDF and reads the simulation domain extent.
2. Reprojects each grid point from WRF geographic coordinates to UTM using pyproj.
3. Writes ``utm_x utm_y u v`` CSV rows (one per grid point per requested time step).
4. Optionally generates an ``inputs.i`` stub with ``velocity_file`` pre-configured.

**Typical usage**

.. code-block:: bash

   # Single time step, auto bounding box from WRF file
   python3 tools/wrf_wind_reader.py --wrf-file wrfout_d01 --wind wind.csv

   # Multiple time steps
   python3 tools/wrf_wind_reader.py --wrf-file wrfout_d01 --wind wind.csv --time-range 0:4

   # Interpolate wind onto an existing terrain XYZ grid
   python3 tools/wrf_wind_reader.py --wrf-file wrfout_d01 --wind wind.csv \
       --terrain-file terrain.xyz

**Key options**

* ``--wrf-file FILE`` — WRF netCDF input (required)
* ``--wind FILE`` — output wind CSV (default: ``wind.csv``)
* ``--time-range T1:TN`` — inclusive range of WRF time indices to extract
* ``--time-index N`` — single time index (default: 0)
* ``--terrain-file PATH`` — existing terrain XYZ; wind is IDW-interpolated onto its grid
* ``--grid-resolution M`` — target grid spacing in metres for regular-grid interpolation
* ``--subsample N`` — keep every N-th grid point (default: 1)

**Dependencies**: ``netCDF4``, ``numpy``, ``pyproj``; ``scipy`` for ``--terrain-file``
or ``--grid-resolution``.

``srtm_terrain_reader.py`` — SRTM Terrain Download
----------------------------------------------------

Downloads SRTM1 (30 m) elevation data for a lat/lon bounding box and writes a UTM terrain
XYZ file for ``rothermel.terrain_file``.

**Key workflow**

1. Downloads the SRTM tile(s) covering the bounding box via the ``elevation`` package.
2. Reprojects the raster to UTM using rasterio and pyproj.
3. Writes ``X Y Z`` rows in UTM metres.

**Typical usage**

.. code-block:: bash

   python3 tools/srtm_terrain_reader.py \
       --lat-min 40 --lat-max 40.5 --lon-min -105 --lon-max -104.5

   # With intermediate GeoTIFF and subsampling
   python3 tools/srtm_terrain_reader.py \
       --lat-min 40 --lat-max 40.5 --lon-min -105 --lon-max -104.5 \
       --tif srtm_raw.tif --subsample 2

**Key options**

* ``--lat-min/max``, ``--lon-min/max`` — bounding box (required)
* ``--terrain FILE`` — output terrain XYZ (default: ``terrain.xyz``)
* ``--tif FILE`` — keep intermediate GeoTIFF at this path
* ``--subsample N`` — keep every N-th point (default: 1)

**Dependencies**: ``elevation``, ``rasterio``, ``numpy``, ``pyproj``, ``scipy``.

``landscape_writer.py`` — LANDFIRE Landscape File Writer
---------------------------------------------------------

Downloads LANDFIRE elevation, slope, aspect, and fuel model rasters for a lat/lon bounding
box and writes a FARSITE landscape file (``.lcp``) for ``rothermel.landscape_file``.
Three remote sources are supported (AWS S3 COGs, Microsoft Planetary Computer, USGS REST API).
Local GeoTIFFs can be used to skip downloads entirely.

**Key workflow**

1. Downloads or reads LANDFIRE rasters for the specified bounding box.
2. Optionally derives terrain slope and aspect from SRTM data.
3. Writes an ASCII landscape file: ``X Y ELEVATION SLOPE ASPECT FUEL_MODEL``.

**Typical usage**

.. code-block:: bash

   # Download LANDFIRE via API
   python3 tools/landscape_writer.py \
       --lat-min 40 --lat-max 40.5 --lon-min -105 --lon-max -104.5

   # Use local rasters (no download)
   python3 tools/landscape_writer.py \
       --lat-min 40 --lat-max 40.5 --lon-min -105 --lon-max -104.5 \
       --elev-file elev.tif --slope-file slope.tif \
       --aspect-file aspect.tif --fuel-file fbfm13.tif

   # Scott & Burgan 40 fuel model system
   python3 tools/landscape_writer.py \
       --lat-min 40 --lat-max 40.5 --lon-min -105 --lon-max -104.5 \
       --fuel-model-type 40

**Key options**

* ``--lat-min/max``, ``--lon-min/max`` — bounding box (required)
* ``--landscape FILE`` — output ``.lcp`` file (default: ``landscape.lcp``)
* ``--fuel-model-type {13,40}`` — Anderson 13 or Scott & Burgan 40
* ``--vintage YEAR`` — LANDFIRE vintage year (default: 2020)
* ``--elev-file``, ``--slope-file``, ``--aspect-file``, ``--fuel-file`` — local rasters
* ``--srtm-slope-aspect`` — derive terrain from SRTM when only ``--fuel-file`` is given
* ``--subsample N`` — keep every N-th point (default: 1)

**Dependencies**: ``rasterio``, ``numpy``, ``pyproj``, ``requests``.

``farsite_weather_reader.py`` — FARSITE Weather File Parser
------------------------------------------------------------

Parses FARSITE ``.wtr`` RAWS weather-station files and writes time-stamped wind CSV files for
the solver.  Optionally generates ``inputs.i`` stubs and estimated fuel moisture files.

**Typical usage**

.. code-block:: bash

   # Single wind file from the first weather record
   python3 tools/farsite_weather_reader.py --wtr fire.wtr --wind wind.csv

   # Time-dependent wind snapshots for a 2 km × 2 km domain at 100 m resolution
   python3 tools/farsite_weather_reader.py \
       --wtr fire.wtr --wind wind.csv \
       --time-dependent --domain 2000 2000 --resolution 100

   # Generate solver input stub and estimated fuel moisture
   python3 tools/farsite_weather_reader.py \
       --wtr fire.wtr --wind wind.csv \
       --inputs-stub inputs.i --write-moisture

**Dependencies**: none (pure Python).

``fuel_moisture_from_weather.py`` — Equilibrium Fuel Moisture Calculator
--------------------------------------------------------------------------

Estimates 1-hr, 10-hr, and 100-hr dead fuel moisture content from temperature and relative
humidity using the Nelson (2000) / Simard (1968) Equilibrium Moisture Content (EMC) model,
as implemented in FARSITE, BehavePlus, and FLAMMAP.

**Typical usage**

.. code-block:: bash

   # Single observation
   python3 tools/fuel_moisture_from_weather.py --temp 30 --rh 20

   # Output in solver input format
   python3 tools/fuel_moisture_from_weather.py --temp 28 --rh 25 --solver-format

   # Batch from a weather CSV
   python3 tools/fuel_moisture_from_weather.py --csv weather.csv --solver-format

**Dependencies**: none (pure Python).

``plotfile_to_geotiff.py`` — GIS Export
-----------------------------------------

Converts AMReX 2-D plotfiles to GeoTIFF rasters and GeoJSON fire-perimeter contours for
import into QGIS, ArcGIS, or other GIS tools.  Multi-level AMR plotfiles (``finest_level > 0``)
are supported; finer-level data is composited onto the Level 0 base grid.

**Default fire-behaviour variables** exported when no ``-v`` filter is given:

``phi``, ``R``, ``fireline_intensity``, ``flame_length``, ``elevation``,
``slope``, ``aspect``, ``fuel_model``, ``fuel_consumption``,
``residual_fuel`` *(post-frontal burnout fraction)*, ``crown_fraction``,
``arrival_time``, ``reaction_intensity``, ``wind_speed``, ``wind_direction``.

**Typical usage**

.. code-block:: bash

   # Export all fire-behaviour variables
   python3 tools/plotfile_to_geotiff.py plt0100 --outdir gis_out

   # Export specific variables with UTM georeference
   python3 tools/plotfile_to_geotiff.py plt0100 \
       -v phi R fireline_intensity flame_length \
       --utm-origin 450000 4200000 --epsg 32613 --outdir gis_out

   # Post-frontal fuel state rasters
   python3 tools/plotfile_to_geotiff.py plt0100 \
       -v residual_fuel fuel_consumption arrival_time \
       --utm-origin 450000 4200000 --epsg 32613 --outdir gis_out

   # Batch-convert all plt#### directories
   python3 tools/plotfile_to_geotiff.py --all --outdir gis_out

**Key options**

* ``-v VAR [VAR …]`` — export only named variables (default: fire-behaviour subset)
* ``--all-vars`` — export every variable in the plotfile
* ``--all`` — process all ``plt####`` directories in the current working directory
* ``--utm-origin EASTING NORTHING`` — UTM origin to add to simulation coordinates
* ``--epsg CODE`` — EPSG CRS code (e.g. ``32613`` for UTM zone 13N)

**Output per plotfile**:

* ``plt####_<variable>.tif`` — single-band float32 GeoTIFF (deflate-compressed)
* ``plt####_fire_perimeter.geojson`` — :math:`\phi = 0` contour as GeoJSON LineString

**Dependencies**: ``rasterio``, ``numpy``; ``matplotlib`` for GeoJSON perimeters.

Combining the Tools: Typical Workflow
--------------------------------------

A typical field-scale simulation setup proceeds as follows:

1. **Download terrain and landscape data**:

   .. code-block:: bash

      python3 tools/srtm_terrain_reader.py \
          --lat-min 40 --lat-max 40.5 --lon-min -105 --lon-max -104.5

      python3 tools/landscape_writer.py \
          --lat-min 40 --lat-max 40.5 --lon-min -105 --lon-max -104.5

2. **Extract wind from a WRF forecast**:

   .. code-block:: bash

      python3 tools/wrf_wind_reader.py \
          --wrf-file wrfout_d01 --wind wind.csv --terrain-file terrain.xyz

3. **Run the solver** using the generated files:

   .. code-block:: text

      rothermel.terrain_file   = terrain.xyz
      rothermel.landscape_file = landscape.lcp
      velocity_file            = wind.csv

4. **Export results for GIS**:

   .. code-block:: bash

      python3 tools/plotfile_to_geotiff.py --all --outdir gis_out \
          --utm-origin 450000 4200000 --epsg 32613

``farsite_fms_reader.py`` — FARSITE Fuel Moisture Scenario Parser
------------------------------------------------------------------

Reads a FARSITE ``.fms`` file that specifies per-fuel-model dead and live fuel
moisture contents and writes solver-compatible moisture inputs.

**FMS file format**:

.. code-block:: text

   <num_fuel_models>
   <fuel_model_num>  <1hr%>  <10hr%>  <100hr%>  <live_herb%>  <live_wood%>

Moisture values are **integer percentages** (e.g. ``8`` → 8 %, which is 0.08 as
a fraction).  Non-burnable models (91–99) are silently skipped.

**Output modes**

1. **Print** (default): human-readable table and a ready-to-paste block of solver
   ``inputs.i`` parameters for every fuel model.
2. ``--write-inputs FILE``: write a solver ``inputs.i`` moisture stub for one model.
3. ``--write-csv FILE``: export all fuel models as a CSV (fractions and percentages).
4. ``--query NUM``: print the solver block for a single fuel model.

**Typical usage**

.. code-block:: bash

   # Print moisture summary
   python3 tools/farsite_fms_reader.py fire.fms

   # Write solver input stub using fuel model 4
   python3 tools/farsite_fms_reader.py fire.fms --fuel-model 4 \
       --write-inputs inputs_moisture.i

   # Export all moisture values as CSV
   python3 tools/farsite_fms_reader.py fire.fms --write-csv moisture.csv

   # Query a single fuel model
   python3 tools/farsite_fms_reader.py fire.fms --query 8

**Key options**

* ``FILE`` — input ``.fms`` file (positional, required)
* ``--fuel-model NUM`` — model to use for ``--write-inputs`` (default: first in file)
* ``--write-inputs FILE`` — write a solver ``inputs.i`` moisture stub
* ``--write-csv FILE`` — write all moisture values as a CSV table
* ``--query NUM`` — print only the entry for a specific fuel model

**Dependencies**: none (pure Python).

``farsite_adj_reader.py`` — FARSITE Fuel Adjustment File Parser
----------------------------------------------------------------

Reads per-fuel-model rate-of-spread (ROS) multipliers from a FARSITE ``.adj``
file, generates template files, and can apply adjustments to a Rothermel
parameter CSV.

**ADJ file format**:

.. code-block:: text

   # Comments start with # or !
   [optional: <num_fuel_models>]
   <fuel_model_number>  <adj_factor>

The ``adj_factor`` is a multiplicative scale applied to Rothermel R₀.

**Typical usage**

.. code-block:: bash

   # Inspect an existing .adj file
   python3 tools/farsite_adj_reader.py --adj fire.adj

   # Generate a template for all FBFM13 models with adj_factor = 1.0
   python3 tools/farsite_adj_reader.py --generate --fuel-system 13 \
       --output template_13.adj

   # Generate a template for specific models
   python3 tools/farsite_adj_reader.py --generate --models 4 9 10 \
       --output my_adj.adj

   # Apply adjustments to a Rothermel CSV for calibration inspection
   python3 tools/farsite_adj_reader.py --adj fire.adj \
       --apply-csv rothermel_params.csv \
       --output rothermel_adjusted.csv

**Key options**

* ``--adj FILE`` — path to ``.adj`` file to inspect or apply
* ``--generate`` — generate a template ``.adj`` file (requires ``--output``)
* ``--fuel-system {13|40}`` — fuel model system for template generation (default: 13)
* ``--models MODEL [MODEL …]`` — specific fuel model numbers to include in template
* ``--apply-csv FILE`` — Rothermel CSV to apply adjustments to
* ``--output / -o FILE`` — output file path

Integrate with the solver via ``inputs.i``:

.. code-block:: text

   fuel_adj_file  = fire.adj
   fuel_adj_model = 4

**Dependencies**: none (pure Python).

``farsite_fmd_reader.py`` — FARSITE Fuel Moisture Schedule Parser
------------------------------------------------------------------

Reads, converts, queries, and generates FARSITE-compatible ``.fmd`` fuel moisture
schedule files.  The solver reads these via ``fmd_file`` and updates moisture at
every time step.

**FMD file format**:

.. code-block:: text

   MONTH  DAY  HOUR  PRECIP  NUM_MODELS
   MODEL  1HR%  10HR%  100HR%  LHERB%  LWOOD%

**Output modes**

1. **Inspect** (default): print the full schedule as a table.
2. ``--output FILE`` (with ``--fmd``): convert the schedule to a flat CSV.
3. ``--query-time SECONDS --fuel-model N``: interpolated moisture at a simulation time.
4. ``--generate --output FILE``: write a constant-moisture ``.fmd`` template.

**Typical usage**

.. code-block:: bash

   # Inspect a .fmd file
   python3 tools/farsite_fmd_reader.py --fmd fire.fmd

   # Convert to a flat CSV for analysis
   python3 tools/farsite_fmd_reader.py --fmd fire.fmd --output moisture.csv

   # Query interpolated moisture at t = 3600 s for FM4
   python3 tools/farsite_fmd_reader.py --fmd fire.fmd \
       --query-time 3600 --fuel-model 4

   # Generate a 24-hour constant-moisture template
   python3 tools/farsite_fmd_reader.py --generate \
       --models 4 9 10 \
       --start-month 7 --start-day 15 --hours 24 \
       --M-d1 8 --M-d10 10 --M-d100 15 \
       --M-lh 90 --M-lw 120 \
       --output summer_moisture.fmd

**Key options**

* ``--fmd FILE`` — ``.fmd`` file to parse / inspect / convert
* ``--output / -o FILE`` — output CSV or ``.fmd`` template path
* ``--query-time SECONDS`` — query interpolated moisture at this simulation time
* ``--fuel-model N`` — fuel model code for ``--query-time`` (default: 0 = global)
* ``--generate`` — generate a constant-moisture ``.fmd`` template
* ``--models N [N …]`` — fuel model codes (default: FBFM13 1–13)
* ``--start-month M``, ``--start-day D`` — reference start date
* ``--hours N`` — number of hourly snapshots (default: 24)
* ``--M-d1``, ``--M-d10``, ``--M-d100``, ``--M-lh``, ``--M-lw`` — moisture [%]

Integrate with the solver via ``inputs.i``:

.. code-block:: text

   fmd_file        = summer_moisture.fmd
   fmd_start_month = 7
   fmd_start_day   = 15
   fmd_fuel_model  = 0

**Dependencies**: none (pure Python).

``perimeter_to_shapefile.py`` — Fire Perimeter to Shapefile Converter
----------------------------------------------------------------------

Converts wildfire_levelset fire perimeter output files to Esri Shapefile format
for import into QGIS, ArcGIS, and other GIS tools.

**Supported input formats**

* ``perimeter_NNNN.geojson`` — GeoJSON FeatureCollection with a Polygon
* ``perimeter_NNNN.csv`` — two-column CSV with X, Y coordinates
* ``phi_negative_NNNN.dat`` — whitespace-delimited X Y point cloud (φ < 0)
* ``phi_envelope_NNNN.dat`` — whitespace-delimited X Y convex hull vertices

**Typical usage**

.. code-block:: bash

   # Convert a single GeoJSON perimeter to shapefile
   python3 tools/perimeter_to_shapefile.py perimeter_0100.geojson

   # Convert all GeoJSON perimeters and assign UTM zone 13N (EPSG:32613)
   python3 tools/perimeter_to_shapefile.py perimeter_*.geojson --epsg 32613

   # Convert CSV perimeters, write to a different directory
   python3 tools/perimeter_to_shapefile.py perimeter_*.csv --outdir shapefiles/

   # Convert all phi_negative point clouds to Multipoint shapefiles
   python3 tools/perimeter_to_shapefile.py phi_negative_*.dat --point-cloud

   # Batch convert everything in the current directory
   python3 tools/perimeter_to_shapefile.py --all

**Key options**

* ``INPUT [INPUT …]`` — input perimeter file(s); glob patterns are supported
* ``--outdir DIR`` — output directory (default: same as input file)
* ``--epsg N`` — EPSG code for the ``.prj`` file (e.g. ``32613`` for UTM zone 13N)
* ``--point-cloud`` — write ``phi_negative`` ``.dat`` files as Multipoint shapefile
* ``--all`` — convert all supported perimeter files in the current directory
* ``--overwrite`` — overwrite existing ``.shp`` files

**Output per input file**: ``.shp``, ``.dbf``, ``.shx``, and (when ``--epsg`` is given) ``.prj``.

**Dependencies**: ``pyshp`` (required); ``shapely`` for convex hull; ``pyproj`` for rich CRS WKT.

.. code-block:: bash

   pip install pyshp shapely pyproj

``utm_convert.py`` — Lat/Lon ↔ UTM Coordinate Conversion
----------------------------------------------------------

Bidirectional WGS-84 lat/lon ↔ UTM conversion.  Works as a standalone CLI
tool or as an importable Python module.  Uses a pure-Python WGS-84 Transverse
Mercator implementation by default and automatically switches to the more
accurate pyproj backend when available.

**Module API**

.. code-block:: python

   from tools.utm_convert import latlon_to_utm, utm_to_latlon

   easting, northing, zone_number, zone_letter = latlon_to_utm(34.10, -118.85)
   lat, lon = utm_to_latlon(330000, 3775000, zone_number=11, northern=True)

**Typical CLI usage**

.. code-block:: bash

   # Convert lat/lon → UTM (zone auto-detected)
   python3 tools/utm_convert.py --to-utm 34.10 -118.85

   # Convert lat/lon → UTM (force specific zone)
   python3 tools/utm_convert.py --to-utm 34.10 -118.85 --zone 11

   # Convert UTM → lat/lon
   python3 tools/utm_convert.py --to-latlon 330000 3775000 --zone 11

   # Southern hemisphere
   python3 tools/utm_convert.py --to-latlon 330000 3775000 --zone 11 --south

**Key options**

* ``--to-utm LAT LON`` — convert lat/lon (decimal degrees) to UTM
* ``--to-latlon EASTING NORTHING`` — convert UTM metres to lat/lon
* ``--zone N`` — UTM zone number (1–60); auto-detected for ``--to-utm``, required for ``--to-latlon``
* ``--south`` — southern hemisphere (applies to ``--to-latlon`` only)

**Dependencies**: none (pure Python); ``pyproj`` optional for higher accuracy.

``ensemble_burn_probability.py`` — Ensemble Burn Probability Driver
--------------------------------------------------------------------

FSPro-style ensemble burn probability driver.  Runs the solver *N* times with
perturbed wind speed, wind direction, and dead fuel moisture (Latin hypercube
or random sampling), accumulates per-cell burn counts, and writes a probability
map.

**Perturbation model**

Three independent parameters are perturbed:

* wind speed factor (multiplicative, lognormal)
* wind direction offset (additive, normal, degrees)
* 1-hr dead fuel moisture offset (additive, normal, fraction)

**Typical usage**

.. code-block:: bash

   # Basic serial: 50 runs, ±20% wind speed, ±15° direction, ±2% moisture
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

   # Georeferenced probability map, parallel members
   python3 tools/ensemble_burn_probability.py \
       --exe ./wildfire_levelset \
       --inputs inputs.i \
       --n-runs 200 \
       --resolution 30 \
       --out burn_probability.csv \
       --geojson burn_probability.geojson \
       --jobs 2

   # Reproducible run with fixed seed
   python3 tools/ensemble_burn_probability.py \
       --exe ./wildfire_levelset \
       --inputs inputs.i \
       --n-runs 50 \
       --seed 42

**Key options**

* ``--exe FILE`` — solver executable path (default: ``./wildfire_levelset``)
* ``--inputs FILE`` — template ``inputs.i`` file (default: ``inputs.i``)
* ``--n-runs N`` — number of ensemble members (default: 50)
* ``--mpi-ranks N`` — MPI ranks per member; 0 = serial (default: 0)
* ``--mpirun CMD`` — MPI launcher command (default: ``mpirun``)
* ``--wind-speed-sigma S`` — std dev of multiplicative wind-speed factor (default: 0.20)
* ``--wind-dir-sigma D`` — std dev of wind-direction offset [deg] (default: 15.0)
* ``--moisture-sigma M`` — std dev of M_d1 offset [fraction] (default: 0.02)
* ``--resolution R`` — grid spacing for probability accumulation [m] (default: auto)
* ``--out FILE`` — output burn probability CSV (default: ``burn_probability.csv``)
* ``--geojson FILE`` — optional GeoJSON probability raster output
* ``--work-dir DIR`` — base working directory for scratch dirs (default: ``/tmp/ensemble``)
* ``--seed N`` — random seed (default: system clock)
* ``--sampler {lhs|random}`` — sampling method (default: ``lhs``)
* ``--jobs J`` — parallel solver runs (default: 1)
* ``--keep-runs`` — keep individual run directories after aggregation

**Output files**

* ``burn_probability.csv`` — ``X, Y, P_burn`` columns; one row per grid cell
* ``burn_probability.geojson`` — optional GeoJSON probability field (when ``--geojson`` is given)

**Dependencies**: none required; ``scipy`` for Latin hypercube sampling (falls back to pure-Python LHS otherwise).

``satellite_goes_to_csv.py`` / Satellite Fire Detection Assimilation
----------------------------------------------------------------------

**C++ module**: ``src/satellite_assimilation.H``
**Python helper**: ``tools/satellite_goes_to_csv.py`` (GOES NetCDF → CSV converter)

The satellite assimilation module ingests real-time active-fire detections from
GOES, VIIRS, or a pre-prepared CSV file, and merges them into the fire level-set
field ``phi`` as new ignition disks at the detected locations.  Detections can
be applied as the initial condition at ``t = 0``, re-applied during the
simulation at a configurable interval, or both.

Data Sources
~~~~~~~~~~~~

Three source modes are supported:

* **"file"** (recommended for HPC / offline) — reads a pre-downloaded
  lon/lat/confidence CSV.  Format::

      # longitude_deg, latitude_deg, confidence_pct
      -119.42, 37.85, 75
      -119.40, 37.86, 80

  Generate this file from NASA FIRMS with any of the workflows below.

* **"viirs"** — fetches the most-recent VIIRS-SNPP NRT active-fire CSV from
  the NASA FIRMS REST API.  Requires a free map key (register at
  https://firms.modaps.eosdis.nasa.gov/api/map_key/) and internet
  connectivity from the compute node.

* **"goes"** — downloads the latest NOAA GOES-16/17/18 ABI fire-detection
  granule from the public AWS S3 bucket and converts it to CSV using the
  optional ``satellite_goes_to_csv.py`` helper (requires ``netCDF4``).
  Falls back to the local cache if the helper is unavailable or the
  download fails.  Network errors are non-fatal.

Preparing a Detection CSV with Python (FIRMS)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A minimal Python snippet to download VIIRS NRT detections for a bounding
box and write the expected three-column CSV:

.. code-block:: python

   import csv, urllib.request

   MAP_KEY   = "YOUR_FIRMS_MAP_KEY"     # register at firms.modaps.eosdis.nasa.gov
   BBOX      = "-120,33,-114,42"        # lon_min,lat_min,lon_max,lat_max
   DAYS      = 1
   URL       = (f"https://firms.modaps.eosdis.nasa.gov/api/area/csv"
                f"/{MAP_KEY}/VIIRS_SNPP_NRT/{BBOX}/{DAYS}")

   with urllib.request.urlopen(URL) as resp:
       reader = csv.DictReader(resp.read().decode().splitlines())
       with open("fire_detections.csv", "w", newline="") as fh:
           fh.write("# longitude_deg,latitude_deg,confidence_pct\n")
           for row in reader:
               fh.write(f"{row['longitude']},{row['latitude']},{row['confidence']}\n")

Save the script and run it before the solver::

   python3 prepare_detections.py
   ./wildfire_levelset inputs.i

Coordinate Alignment
~~~~~~~~~~~~~~~~~~~~

Detections arrive in WGS-84 lon/lat degrees.  The module converts them to
UTM easting and northing [m] using a compact Karney (2011) Transverse
Mercator projection.  The simulation domain origin (``prob_lo_x``,
``prob_lo_y``) must be aligned with the UTM origin so that::

   sim_x = UTM_easting  - prob_lo_easting_m
   sim_y = UTM_northing - prob_lo_northing_m

Example for a domain spanning UTM zone 11N easting [500 000, 560 000] m and
northing [3 700 000, 3 760 000] m::

   geometry.prob_lo = 500000  3700000
   geometry.prob_hi = 560000  3760000

   satellite.utm_zone           = 11
   satellite.utm_northern       = 1
   satellite.prob_lo_easting_m  = 500000
   satellite.prob_lo_northing_m = 3700000

Input Parameters (prefix ``satellite.``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 38 10 52

   * - Parameter
     - Default
     - Description
   * - ``enable``
     - 0
     - 1 to activate the module
   * - ``source``
     - ``"file"``
     - Detection source: ``"file"``, ``"viirs"``, or ``"goes"``
   * - ``local_file``
     - ``""``
     - Path to pre-downloaded CSV (required for ``source="file"``)
   * - ``api_key``
     - ``""``
     - NASA FIRMS map key (required for ``source="viirs"``)
   * - ``viirs_url_base``
     - ``"https://firms.modaps.eosdis.nasa.gov/api/area/csv"``
     - Base URL for the NASA FIRMS REST API (override for proxies or mirrors)
   * - ``goes_product``
     - ``"ABI-L2-FDCF"``
     - GOES ABI fire product name
   * - ``goes_bucket``
     - ``"noaa-goes18"``
     - AWS S3 bucket (``noaa-goes16``, ``noaa-goes17``, or ``noaa-goes18``)
   * - ``bbox_lon_min``
     - \-180.0
     - Bounding box minimum longitude [deg]
   * - ``bbox_lon_max``
     - 180.0
     - Bounding box maximum longitude [deg]
   * - ``bbox_lat_min``
     - \-90.0
     - Bounding box minimum latitude [deg]
   * - ``bbox_lat_max``
     - 90.0
     - Bounding box maximum latitude [deg]
   * - ``utm_zone``
     - 10
     - UTM zone number 1–60
   * - ``utm_northern``
     - 1
     - 1 = northern hemisphere; 0 = southern
   * - ``prob_lo_easting_m``
     - 0.0
     - UTM easting of ``prob_lo_x`` [m]
   * - ``prob_lo_northing_m``
     - 0.0
     - UTM northing of ``prob_lo_y`` [m]
   * - ``fetch_interval_s``
     - 600.0
     - Re-fetch interval during simulation [s]
   * - ``use_as_ic``
     - 1
     - 1 to apply detections as initial condition at ``t = 0``
   * - ``use_mid_sim``
     - 1
     - 1 to re-apply new detections during the simulation run
   * - ``confidence_threshold``
     - 50
     - Minimum detection confidence [%]
   * - ``detection_radius_m``
     - 375.0
     - Radius of the ignition disk placed at each detection [m]
   * - ``local_cache_file``
     - ``""``
     - Path for caching fetched data; used as fallback on network failure
   * - ``suppress_if_burning``
     - 1
     - 1 = skip cells already burning (prevent extinguishing front)

Example — offline file-based assimilation::

   satellite.enable              = 1
   satellite.source              = file
   satellite.local_file          = fire_detections.csv
   satellite.utm_zone            = 11
   satellite.utm_northern        = 1
   satellite.prob_lo_easting_m   = 500000
   satellite.prob_lo_northing_m  = 3700000
   satellite.confidence_threshold = 60
   satellite.detection_radius_m  = 375.0
   satellite.use_as_ic           = 1
   satellite.use_mid_sim         = 0

Example — live VIIRS re-assimilation every 10 minutes::

   satellite.enable              = 1
   satellite.source              = viirs
   satellite.api_key             = ABCDEF123456
   satellite.bbox_lon_min        = -120.0
   satellite.bbox_lon_max        = -114.0
   satellite.bbox_lat_min        = 33.0
   satellite.bbox_lat_max        = 42.0
   satellite.utm_zone            = 11
   satellite.utm_northern        = 1
   satellite.prob_lo_easting_m   = 500000
   satellite.prob_lo_northing_m  = 3700000
   satellite.fetch_interval_s    = 600.0
   satellite.use_as_ic           = 1
   satellite.use_mid_sim         = 1
   satellite.local_cache_file    = satellite_cache.csv

**Test**: ``regtest/ignition/satellite_assimilation/`` (uses
``create_detections.py`` to generate a synthetic detection CSV).

**Dependencies**: ``curl`` system command (for VIIRS / GOES live modes);
``netCDF4`` Python package plus ``satellite_goes_to_csv.py`` helper script
(for GOES NetCDF conversion only; not required for ``"file"`` or ``"viirs"``
modes).


``ember_cascade_analysis.py`` — Ember Cascade Post-Processing
--------------------------------------------------------------

Post-processing tool for the flux-based ember cascade model
(``ember_cascade.enable = 1``).  Reads the ``ember_cascade_flux`` and
``ember_cascade_ignition`` fields from AMReX plotfiles and produces:

* A per-step summary table (max / mean landing flux, number of active cells,
  number of spot-fire ignitions).
* An optional CSV export for further analysis.
* Optional visualisation of the Gaussian landing-flux map alongside the fire
  perimeter contour.

The reader is dependency-free for plain AMReX plotfiles — it parses the ASCII
``Header`` and binary ``Cell_D_*`` FAB files directly without requiring the
AMReX Python bindings.

.. code-block:: bash

   # Summarise all steps
   python3 tools/ember_cascade_analysis.py --plt-dir .

   # Export statistics to CSV
   python3 tools/ember_cascade_analysis.py --plt-dir . --csv ember.csv

   # Visualise landing flux map at step 50
   python3 tools/ember_cascade_analysis.py --plt-dir . --plot 0050 \
       --output flux_step50.png

**Dependencies**: ``numpy`` (required); ``matplotlib`` (optional, for ``--plot``).

``historical_wildfires.py`` — Historical Wildfire Database & Simulation Setup
------------------------------------------------------------------------------

Provides a curated database of 19+ major historical wildfires (2012–2021) spanning
the USA, Australia, Canada, Europe, Russia, and Indonesia. Supports querying
the fire list, exporting to CSV, and automatic generation of simulation input files
(``terrain.csv`` and ``inputs.i``) with real SRTM elevation data.

**Key workflow**

1. Query the wildfire database and display as ASCII table.
2. Filter by country, state, year range, and other criteria.
3. Export selected fires to CSV for further analysis.
4. Generate ``terrain.csv`` (SRTM elevation) and ``inputs.i`` (solver template)
   for a selected fire at configurable lat/lon bounds (default ±0.25°).

**Typical usage**

.. code-block:: bash

   # List all available fires
   python3 tools/historical_wildfires.py --list-fires

   # Display all fires as a table
   python3 tools/historical_wildfires.py

   # Filter by country and year
   python3 tools/historical_wildfires.py --country USA --year-min 2020

   # Export to CSV
   python3 tools/historical_wildfires.py --country Australia --output aus_fires.csv

   # Generate terrain and input files for a fire
   python3 tools/historical_wildfires.py --create-inputs "Dixie Fire" \
       --outdir dixie_sim

   # With custom margins (±0.5° instead of default ±0.25°)
   python3 tools/historical_wildfires.py --create-inputs "Marshall Fire" \
       --outdir marshall_sim --lat-margin 0.5 --lon-margin 0.5

   # With SRTM subsampling (keep every 2nd point)
   python3 tools/historical_wildfires.py --create-inputs "Creek Fire" \
       --outdir creek_sim --subsample 2

**Key options**

* ``--list-fires`` — Print all available fire names with year, location.
* ``--output FILE`` — Export filtered fires to CSV (instead of printing table).
* ``--country COUNTRY`` — Filter by country name.
* ``--state STATE`` — Filter by state/region name.
* ``--year-min / --year-max YEAR`` — Filter by year range.
* ``--columns COL1,COL2,...`` — Select columns to display/export.
* ``--create-inputs FIRE_NAME`` — Generate ``terrain.csv`` and ``inputs.i``.
* ``--outdir DIR`` — Output directory for generated files (default: current).
* ``--lat-margin / --lon-margin DEGREES`` — Bounds margin (default: 0.25°).
* ``--subsample N`` — Keep every N-th SRTM grid point (default: 1).

**Generated files**

When ``--create-inputs`` is used, the tool calls ``srtm_terrain_reader.py`` to
download real SRTM1 elevation data and produces:

* ``terrain.csv`` — UTM-projected XYZ grid of elevation data for ``rothermel.terrain_file``.
* ``inputs.i`` — Solver input template with domain parameters, fuel model (Anderson class 7),
  fuel moisture, ignition point at domain centre, and basic output settings.

The ``inputs.i`` template includes comments and should be edited to customise
wind, weather, suppression, and other scenario parameters before running the solver.

**Database contents**

29 major US wildfires from 2009–2024:

* **California** (22 fires): Park Fire (2023), Coastal Fire (2023), Dixie Fire (2021),
  August Complex (2020), Creek Fire (2020), Apple Fire (2020), Slinkard Fire (2020),
  Kincade Fire (2019), Easy Fire (2019), Woolsey Fire (2018), Carr Fire (2018),
  Delta Fire (2018), Tubbs Fire (2017), Thomas Fire (2017), Detwiler Fire (2017),
  Soberanes Fire (2016), Erskine Fire (2016), Butte Fire (2015), Valley Fire (2015),
  King Fire (2014), Rim Fire (2013), Yosemite NP Fire (2013)

* **Oregon** (2 fires): Bootleg Fire (2021), Ditch Fire (2020)

* **Colorado** (2 fires): Marshall Fire (2021), Waldo Canyon Fire (2012)

* **New Mexico** (2 fires): Hermits Peak-Calf Canyon Fire (2022), Black Summer Fire (2022)

* **Arizona** (1 fire): Wallow Fire (2011)

Each entry includes fire name, state, ignition year/month/day, central lat/lon,
burned area (hectares), and duration (days).

**Dependencies**

``srtm_terrain_reader.py`` backend requires: ``elevation``, ``rasterio``,
``numpy``, ``pyproj``, ``scipy`` (install with ``pip install elevation rasterio ...``).

Plain table/filter operations have no dependencies.