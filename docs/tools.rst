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
   * - ``plotfile_to_geotiff.py``
     - Convert AMReX plotfiles → GeoTIFF rasters and GeoJSON fire-perimeter contours

A unified legacy tool (``terrain_wind_preprocess.py``) is retained in ``tools/deprecated/``
and superseded by the split tools above.

Dependencies
------------

Install with pip before using any tool:

.. code-block:: bash

   pip install numpy netCDF4 pyproj rasterio elevation requests scipy matplotlib

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

**Typical usage**

.. code-block:: bash

   # Export all fire-behaviour variables
   python3 tools/plotfile_to_geotiff.py plt0100 --outdir gis_out

   # Export specific variables with UTM georeference
   python3 tools/plotfile_to_geotiff.py plt0100 \
       -v phi R fireline_intensity flame_length \
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
