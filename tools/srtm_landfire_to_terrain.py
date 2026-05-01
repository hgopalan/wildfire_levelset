#!/usr/bin/env python3
"""
srtm_landfire_to_terrain.py - Unified terrain preprocessing tool for wildfire_levelset.

Downloads SRTM elevation data and LANDFIRE fuel/slope/aspect rasters for a
user-specified lat/lon bounding box and writes:
  1. A terrain XYZ file  (X Y Z in UTM metres) for ``rothermel.terrain_file``
  2. A landscape LCP file (X Y ELEVATION SLOPE ASPECT FUEL_MODEL) for
     ``rothermel.landscape_file``

Both outputs use UTM coordinates by default.

Usage:
  # Download both terrain and landscape
  python3 srtm_landfire_to_terrain.py \\
      --lat-min 40 --lat-max 40.5 \\
      --lon-min -105 --lon-max -104.5

  # Custom output filenames
  python3 srtm_landfire_to_terrain.py \\
      --lat-min 39.5 --lat-max 40.2 \\
      --lon-min -106 --lon-max -105.2 \\
      --terrain region_terrain.xyz \\
      --landscape region_landscape.lcp

  # Skip one of the outputs
  python3 srtm_landfire_to_terrain.py \\
      --lat-min 40 --lat-max 40.5 \\
      --lon-min -105 --lon-max -104.5 \\
      --no-landscape

  # Keep the intermediate SRTM GeoTIFF
  python3 srtm_landfire_to_terrain.py \\
      --lat-min 40 --lat-max 40.5 \\
      --lon-min -105 --lon-max -104.5 \\
      --tif srtm_clip.tif

Options:
  --lat-min / --lat-max     Latitude  bounds (WGS-84 decimal degrees)
  --lon-min / --lon-max     Longitude bounds (WGS-84 decimal degrees)
  --terrain FILE            Output terrain XYZ file  (default: terrain.xyz)
  --landscape FILE          Output landscape LCP file (default: landscape.lcp)
  --tif FILE                Path for intermediate SRTM GeoTIFF
                            (default: srtm_clip.tif; deleted after use unless
                            you specify an explicit path)
  --subsample N             Keep every N-th point in each dimension (default: 1)
  --no-terrain              Skip the SRTM terrain XYZ step
  --no-landscape            Skip the LANDFIRE landscape LCP step
  --keep-nonburnable        Include non-burnable LANDFIRE pixels in LCP output
  --vintage YEAR            LANDFIRE data vintage year (default: 2020)
  --fuel-product ID         Override LANDFIRE fuel model product ID
  --elev-product ID         Override LANDFIRE elevation product ID
  --cache-dir DIR           Directory to cache downloaded LANDFIRE ZIP files
  --timeout N               LANDFIRE API polling timeout in seconds (default: 300)
  --help                    Show this message and exit.

Terrain XYZ output format:
  # X Y Z (meters)
  x1  y1  z1
  x2  y2  z2
  ...

Landscape LCP output format:
  # X Y ELEVATION SLOPE ASPECT FUEL_MODEL
  x1  y1  elev1  slope1  aspect1  fuel1
  ...

Requires: elevation, rasterio, numpy, pyproj (and requests + pyproj for LANDFIRE)
Install : pip install elevation rasterio numpy pyproj requests
"""

import argparse
import io
import json
import math
import os
import sys
import tempfile
import time
import zipfile

# ---------------------------------------------------------------------------
# LANDFIRE LFPS API constants
# ---------------------------------------------------------------------------

_LFPS_BASE = (
    "https://lfps.usgs.gov/arcgis/rest/services/"
    "LandfireProductService/GPServer/LandfireProductService"
)
_SUBMIT_URL = f"{_LFPS_BASE}/submitJob"
_JOB_URL    = f"{_LFPS_BASE}/jobs/{{jobId}}"
_RESULT_URL = f"{_LFPS_BASE}/jobs/{{jobId}}/results/Output_File"

# Default vintage-to-product-ID mapping
_DEFAULT_LAYERS = {
    2020: {
        "elev":   "ELEV2020",
        "slope":  "SlpD2020",
        "aspect": "Asp2020",
        "fuel13": "F13_FBFM13",
    },
    2016: {
        "elev":   "ELEV2016",
        "slope":  "SlpD2016",
        "aspect": "Asp2016",
        "fuel13": "F13_FBFM13",
    },
    2014: {
        "elev":   "ELEV2014",
        "slope":  "SLP",
        "aspect": "ASP",
        "fuel13": "F13_FBFM13",
    },
}

# LANDFIRE FBFM13 non-burnable codes (outside the 1-13 Anderson range)
_NONBURNABLE_CODES = {91, 92, 93, 98, 99}


# ===========================================================================
# SRTM terrain helpers
# ===========================================================================

def download_srtm(lat_min, lat_max, lon_min, lon_max, out_tif):
    """Download SRTM1 elevation data clipped to the given bounding box.

    Parameters
    ----------
    lat_min, lat_max, lon_min, lon_max : float
        Bounding box in WGS-84 decimal degrees.
    out_tif : str
        Output path for the clipped GeoTIFF.
    """
    try:
        import elevation
    except ImportError:
        raise ImportError(
            "elevation is required to download SRTM data. "
            "Install with: pip install elevation"
        )
    elevation.clean()
    elevation.clip(
        bounds=(lon_min, lat_min, lon_max, lat_max),
        output=out_tif,
        product="SRTM1",
    )


def tiff_to_xyz_utm(tif_path, subsample=1):
    """Read a GeoTIFF and return (utm_x, utm_y, z) 2-D arrays in UTM metres.

    Parameters
    ----------
    tif_path : str
        Path to the input GeoTIFF.
    subsample : int
        Keep every *subsample*-th row and column (default: 1 = all).

    Returns
    -------
    utm_x, utm_y, z : numpy.ndarray (2-D)
    """
    import numpy as np
    try:
        import rasterio
        from rasterio.transform import xy as rasterio_xy
    except ImportError:
        raise ImportError(
            "rasterio is required to read GeoTIFF files. "
            "Install with: pip install rasterio"
        )
    try:
        from pyproj import Transformer
    except ImportError:
        raise ImportError(
            "pyproj is required for UTM projection. "
            "Install with: pip install pyproj"
        )

    with rasterio.open(tif_path) as src:
        z = src.read(1).astype(float)
        z = np.where(z < 0, 0.0, z)

        h, w = z.shape
        rows_idx, cols_idx = np.meshgrid(np.arange(h), np.arange(w), indexing="ij")
        lon, lat = rasterio_xy(src.transform, rows_idx.ravel(), cols_idx.ravel())
        lon = np.array(lon, dtype=np.float64).reshape(h, w)
        lat = np.array(lat, dtype=np.float64).reshape(h, w)

    # Determine UTM zone from centre
    center_lon = float(np.nanmean(lon))
    center_lat = float(np.nanmean(lat))
    zone = int((center_lon + 180.0) / 6.0) + 1
    epsg = (32600 + zone) if center_lat >= 0.0 else (32700 + zone)

    transformer = Transformer.from_crs("EPSG:4326", f"EPSG:{epsg}", always_xy=True)
    utm_x_flat, utm_y_flat = transformer.transform(lon.ravel(), lat.ravel())
    utm_x = np.array(utm_x_flat, dtype=np.float64).reshape(h, w)
    utm_y = np.array(utm_y_flat, dtype=np.float64).reshape(h, w)

    if subsample > 1:
        utm_x = utm_x[::subsample, ::subsample]
        utm_y = utm_y[::subsample, ::subsample]
        z     = z    [::subsample, ::subsample]

    return utm_x, utm_y, z


def write_terrain_xyz(utm_x, utm_y, z, output_path):
    """Write a terrain XYZ file compatible with ``rothermel.terrain_file``.

    Parameters
    ----------
    utm_x, utm_y : array-like
        Easting and northing in UTM metres.
    z : array-like
        Elevation in metres.
    output_path : str
        Destination file path.
    """
    import numpy as np

    data = np.column_stack([
        np.asarray(utm_x).ravel(),
        np.asarray(utm_y).ravel(),
        np.asarray(z).ravel(),
    ])
    with open(output_path, "w") as fh:
        fh.write("# X Y Z (meters)\n")
        for row in data:
            fh.write(f"{row[0]:.3f} {row[1]:.3f} {row[2]:.3f}\n")
    print(f"Wrote {len(data)} terrain points to '{output_path}'.")


def create_terrain_xyz(output_path, lat_min, lat_max, lon_min, lon_max,
                       tif_path=None, subsample=1):
    """Download SRTM data and write a UTM terrain XYZ file.

    Parameters
    ----------
    output_path : str
        Path for the output terrain XYZ file.
    lat_min, lat_max, lon_min, lon_max : float
        Bounding box in WGS-84 decimal degrees.
    tif_path : str or None
        Path for the intermediate GeoTIFF.  If *None* a temporary file is used
        and deleted after conversion.
    subsample : int
        Keep every N-th point in each dimension (default: 1).
    """
    _tmp_tif = tif_path is None
    if _tmp_tif:
        tmp = tempfile.NamedTemporaryFile(suffix=".tif", delete=False)
        tif_path = tmp.name
        tmp.close()

    try:
        print(f"Downloading SRTM data for bbox ({lat_min},{lon_min}) – ({lat_max},{lon_max}) …")
        download_srtm(lat_min, lat_max, lon_min, lon_max, tif_path)

        print("Converting SRTM GeoTIFF → UTM XYZ …")
        utm_x, utm_y, z = tiff_to_xyz_utm(tif_path, subsample=subsample)

        write_terrain_xyz(utm_x, utm_y, z, output_path)
    finally:
        if _tmp_tif and os.path.isfile(tif_path):
            os.remove(tif_path)


# ===========================================================================
# LANDFIRE landscape helpers
# ===========================================================================

def _submit_lfps_job(bbox, layer_ids, timeout_s=300):
    """Submit a LANDFIRE Product Service job and poll until completion.

    Returns the URL of the output ZIP file.
    """
    try:
        import requests
    except ImportError:
        raise ImportError(
            "requests is required to download LANDFIRE data. "
            "Install with: pip install requests"
        )

    min_lon, min_lat, max_lon, max_lat = bbox
    aoi = json.dumps({
        "xmin": min_lon, "ymin": min_lat,
        "xmax": max_lon, "ymax": max_lat,
        "spatialReference": {"wkid": 4326},
    })
    payload = {
        "Area_of_Interest": aoi,
        "Layer_list":       "|".join(layer_ids),
        "f":                "json",
    }

    print(f"Submitting LANDFIRE job for bbox {bbox} …")
    resp = _requests_post(_SUBMIT_URL, data=payload, timeout=60)
    resp.raise_for_status()
    job_info = resp.json()

    if "jobId" not in job_info:
        raise RuntimeError(f"LANDFIRE job submission failed: {job_info}")

    job_id = job_info["jobId"]
    print(f"  Job ID: {job_id}")

    status_url = _JOB_URL.format(jobId=job_id)
    deadline = time.time() + timeout_s
    poll_interval = 5

    while time.time() < deadline:
        time.sleep(poll_interval)
        sr = _requests_get(status_url, params={"f": "json"}, timeout=30)
        sr.raise_for_status()
        status = sr.json().get("jobStatus", "")
        print(f"  Status: {status}")

        if status == "esriJobSucceeded":
            break
        if status in ("esriJobFailed", "esriJobCancelled",
                      "esriJobCancelling", "esriJobDeleted"):
            raise RuntimeError(f"LANDFIRE job {job_id} ended with status: {status}")
        poll_interval = min(poll_interval * 1.5, 30)
    else:
        raise RuntimeError(f"LANDFIRE job {job_id} did not complete within {timeout_s}s")

    result_url = _RESULT_URL.format(jobId=job_id)
    rr = _requests_get(result_url, params={"f": "json"}, timeout=30)
    rr.raise_for_status()
    result = rr.json()

    download_url = (
        result.get("value", {}).get("url")
        or result.get("value")
    )
    if not download_url:
        raise RuntimeError(f"Could not parse download URL from LANDFIRE result: {result}")
    return download_url


def _requests_get(url, **kwargs):
    import requests
    return requests.get(url, **kwargs)


def _requests_post(url, **kwargs):
    import requests
    return requests.post(url, **kwargs)


def _download_zip(url, cache_dir=None):
    """Download a ZIP from *url* and return its bytes."""
    try:
        import requests
    except ImportError:
        raise ImportError(
            "requests is required to download LANDFIRE data. "
            "Install with: pip install requests"
        )

    print(f"Downloading output ZIP from {url} …")
    with _requests_get(url, stream=True, timeout=120) as r:
        r.raise_for_status()
        data = r.content

    if cache_dir:
        os.makedirs(cache_dir, exist_ok=True)
        fname = os.path.join(cache_dir, os.path.basename(url.split("?")[0]))
        with open(fname, "wb") as fh:
            fh.write(data)
        print(f"  Cached to {fname}")

    return data


def download_landfire(bbox, layer_ids, cache_dir=None, timeout_s=300):
    """Download LANDFIRE rasters and return a dict of {layer_id: bytes}.

    Parameters
    ----------
    bbox : tuple of float
        (min_lon, min_lat, max_lon, max_lat) in WGS-84 decimal degrees.
    layer_ids : list of str
        LANDFIRE layer product IDs to download.
    cache_dir : str or None
        Optional directory to cache the downloaded ZIP.
    timeout_s : int
        Maximum seconds to wait for the LFPS job to complete.

    Returns
    -------
    dict mapping layer_id (str) to raw GeoTIFF bytes.
    """
    download_url = _submit_lfps_job(bbox, layer_ids, timeout_s=timeout_s)
    zip_bytes = _download_zip(download_url, cache_dir=cache_dir)

    layer_map = {}
    with zipfile.ZipFile(io.BytesIO(zip_bytes)) as zf:
        names = zf.namelist()
        print(f"  ZIP contains: {names}")
        for name in names:
            lower = name.lower()
            if not lower.endswith(".tif") and not lower.endswith(".tiff"):
                continue
            raw = zf.read(name)
            matched = False
            for lid in layer_ids:
                if lid.lower() in lower:
                    layer_map[lid] = raw
                    matched = True
                    break
            if not matched:
                stem = os.path.splitext(os.path.basename(name))[0]
                layer_map[stem] = raw

    return layer_map


# ---------------------------------------------------------------------------
# Raster helpers
# ---------------------------------------------------------------------------

def _read_raster_bytes(raw_bytes):
    """Read a raster from raw bytes using rasterio."""
    import numpy as np
    try:
        import rasterio
        from rasterio.io import MemoryFile
    except ImportError:
        raise ImportError(
            "rasterio is required to process LANDFIRE rasters. "
            "Install with: pip install rasterio"
        )

    with MemoryFile(raw_bytes) as memfile:
        with memfile.open() as ds:
            data = ds.read(1).astype(np.float64)
            nodata = ds.nodata
            transform = ds.transform
            crs = ds.crs
    return data, transform, crs, nodata


def _read_raster_file(path, nodata_override=None):
    """Read a raster from a file path using rasterio."""
    import numpy as np
    try:
        import rasterio
    except ImportError:
        raise ImportError(
            "rasterio is required to process LANDFIRE rasters. "
            "Install with: pip install rasterio"
        )

    with rasterio.open(path) as ds:
        data = ds.read(1).astype(np.float64)
        nodata = nodata_override if nodata_override is not None else ds.nodata
        transform = ds.transform
        crs = ds.crs
    return data, transform, crs, nodata


def _raster_coords(data, transform):
    """Return (x_2d, y_2d) arrays of cell-centre coordinates."""
    import numpy as np
    from rasterio.transform import xy as rasterio_xy

    nrows, ncols = data.shape
    rows_idx, cols_idx = np.meshgrid(np.arange(nrows), np.arange(ncols),
                                     indexing="ij")
    xs, ys = rasterio_xy(transform, rows_idx.ravel(), cols_idx.ravel())
    xs = np.array(xs, dtype=np.float64).reshape(nrows, ncols)
    ys = np.array(ys, dtype=np.float64).reshape(nrows, ncols)
    return xs, ys


def _project_to_utm_landfire(xs, ys, src_crs):
    """Project *xs*, *ys* (any CRS) to UTM metres for LANDFIRE data."""
    try:
        from pyproj import Transformer
    except ImportError:
        raise ImportError(
            "pyproj is required for UTM reprojection. "
            "Install with: pip install pyproj"
        )
    import numpy as np

    # Resolve source authority string from rasterio CRS
    try:
        from rasterio.crs import CRS as RasterioCRS
        if isinstance(src_crs, RasterioCRS):
            src_epsg = src_crs.to_epsg()
            src_authority = f"EPSG:{src_epsg}" if src_epsg else src_crs.to_wkt()
        else:
            src_authority = str(src_crs)
    except Exception:
        src_authority = "EPSG:4326"

    # Determine UTM zone from centre (reproject to geographic first if needed)
    try:
        to_geo = Transformer.from_crs(src_authority, "EPSG:4326", always_xy=True)
        lons, lats = to_geo.transform(xs.ravel(), ys.ravel())
    except Exception:
        lons, lats = xs.ravel(), ys.ravel()

    center_lon = float(np.nanmean(lons))
    center_lat = float(np.nanmean(lats))
    zone = int((center_lon + 180.0) / 6.0) + 1
    epsg = (32600 + zone) if center_lat >= 0.0 else (32700 + zone)

    transformer = Transformer.from_crs(src_authority, f"EPSG:{epsg}", always_xy=True)
    x_utm, y_utm = transformer.transform(xs.ravel(), ys.ravel())
    return (np.array(x_utm, dtype=np.float64).reshape(xs.shape),
            np.array(y_utm, dtype=np.float64).reshape(xs.shape))


def _resample_to_grid(src_data, src_transform, src_crs,
                      ref_data, ref_transform, ref_crs,
                      resampling="nearest"):
    """Resample *src* raster to the grid of *ref* using rasterio.warp."""
    import numpy as np
    try:
        from rasterio.warp import reproject, Resampling
    except ImportError:
        raise ImportError("rasterio is required. Install with: pip install rasterio")

    resample_method = (Resampling.nearest if resampling == "nearest"
                       else Resampling.bilinear)

    nrows, ncols = ref_data.shape
    dst = np.zeros((nrows, ncols), dtype=np.float64)
    reproject(
        source=src_data,
        destination=dst,
        src_transform=src_transform,
        src_crs=src_crs,
        dst_transform=ref_transform,
        dst_crs=ref_crs,
        resampling=resample_method,
    )
    return dst


def assemble_landscape(elev_data, elev_transform, elev_crs, elev_nodata,
                       slope_data, slope_transform, slope_crs,
                       aspect_data, aspect_transform, aspect_crs,
                       fuel_data, fuel_transform, fuel_crs,
                       project_utm=True,
                       subsample=1,
                       keep_nonburnable=False):
    """Merge four rasters into flat arrays suitable for writing a .lcp file.

    All rasters are resampled to the elevation grid.  Non-burnable pixels
    (fuel codes 91-99) are excluded unless *keep_nonburnable* is True.

    Returns
    -------
    tuple of six 1-D numpy arrays: (xs, ys, elev, slope, aspect, fuel)
    """
    import numpy as np

    # Resample slope, aspect, and fuel to elevation grid if shapes differ
    if slope_data.shape != elev_data.shape:
        print("  Resampling slope to elevation grid …")
        slope_data = _resample_to_grid(
            slope_data, slope_transform, slope_crs,
            elev_data, elev_transform, elev_crs, resampling="bilinear")
    if aspect_data.shape != elev_data.shape:
        print("  Resampling aspect to elevation grid …")
        aspect_data = _resample_to_grid(
            aspect_data, aspect_transform, aspect_crs,
            elev_data, elev_transform, elev_crs, resampling="bilinear")
    if fuel_data.shape != elev_data.shape:
        print("  Resampling fuel model to elevation grid …")
        fuel_data = _resample_to_grid(
            fuel_data, fuel_transform, fuel_crs,
            elev_data, elev_transform, elev_crs, resampling="nearest")

    xs, ys = _raster_coords(elev_data, elev_transform)

    if project_utm:
        print("  Projecting coordinates to UTM metres …")
        xs, ys = _project_to_utm_landfire(xs, ys, elev_crs)

    fuel_int = np.round(fuel_data).astype(int)

    nodata_val = elev_nodata if elev_nodata is not None else -9999.0
    mask = elev_data != nodata_val
    mask &= np.isfinite(elev_data)

    if not keep_nonburnable:
        for code in _NONBURNABLE_CODES:
            mask &= fuel_int != code
        mask &= (fuel_int >= 1) & (fuel_int <= 13)

    if subsample > 1:
        skip = np.zeros(elev_data.shape, dtype=bool)
        skip[::subsample, ::subsample] = True
        mask &= skip

    return (xs[mask], ys[mask], elev_data[mask],
            slope_data[mask], aspect_data[mask], fuel_int[mask].astype(float))


def _write_lcp(path, xs, ys, elev, slope, aspect, fuel):
    """Write an ASCII landscape file (.lcp)."""
    with open(path, "w") as fh:
        fh.write(
            "# Landscape file generated by srtm_landfire_to_terrain.py\n"
            "# Format: X Y ELEVATION SLOPE ASPECT FUEL_MODEL\n"
            "# Units: X/Y metres, ELEVATION metres, SLOPE degrees,\n"
            "#        ASPECT degrees (0=North), FUEL_MODEL NFFL 1-13\n"
        )
        for x, y, z, s, a, f in zip(xs, ys, elev, slope, aspect, fuel):
            fh.write(f"{x:.2f} {y:.2f} {z:.2f} {s:.2f} {a:.2f} {int(f)}\n")
    print(f"Wrote {len(xs)} landscape points to '{path}'.")


# ---------------------------------------------------------------------------
# High-level landscape pipeline
# ---------------------------------------------------------------------------

def create_landscape(output_path,
                     bbox,
                     vintage=2020,
                     elev_product=None,
                     slope_product=None,
                     aspect_product=None,
                     fuel_product=None,
                     project_utm=True,
                     subsample=1,
                     keep_nonburnable=False,
                     cache_dir=None,
                     timeout_s=300):
    """Download LANDFIRE rasters and write an ASCII landscape file.

    Parameters
    ----------
    output_path : str
        Path for the output .lcp file.
    bbox : tuple of float
        (min_lon, min_lat, max_lon, max_lat) in WGS-84 decimal degrees.
    vintage : int
        LANDFIRE data vintage year (default: 2020).
    elev_product, slope_product, aspect_product, fuel_product : str or None
        Override the default LANDFIRE product IDs for each layer.
    project_utm : bool
        Reproject coordinates to UTM metres (default: True).
    subsample : int
        Keep every N-th pixel in each dimension (default: 1).
    keep_nonburnable : bool
        If True, include non-burnable pixels in output.
    cache_dir : str or None
        Directory to cache downloaded ZIP files.
    timeout_s : int
        Maximum seconds to wait for the LFPS API job to complete.
    """
    layer_defaults = _DEFAULT_LAYERS.get(vintage, _DEFAULT_LAYERS[2020])

    lid_elev   = elev_product   or layer_defaults["elev"]
    lid_slope  = slope_product  or layer_defaults["slope"]
    lid_aspect = aspect_product or layer_defaults["aspect"]
    lid_fuel   = fuel_product   or layer_defaults["fuel13"]

    layer_ids = [lid_elev, lid_slope, lid_aspect, lid_fuel]

    layer_map = download_landfire(bbox, layer_ids,
                                  cache_dir=cache_dir, timeout_s=timeout_s)

    def _get_bytes(lid):
        if lid in layer_map:
            return layer_map[lid]
        for key, val in layer_map.items():
            if lid.lower() in key.lower() or key.lower() in lid.lower():
                return val
        raise KeyError(
            f"Layer '{lid}' not found in downloaded data. "
            f"Available keys: {list(layer_map.keys())}"
        )

    print("Reading elevation raster …")
    elev_data, elev_tf, elev_crs, elev_nd = _read_raster_bytes(_get_bytes(lid_elev))
    print("Reading slope raster …")
    slope_data, slope_tf, slope_crs, _ = _read_raster_bytes(_get_bytes(lid_slope))
    print("Reading aspect raster …")
    asp_data, asp_tf, asp_crs, _        = _read_raster_bytes(_get_bytes(lid_aspect))
    print("Reading fuel model raster …")
    fuel_data, fuel_tf, fuel_crs, _     = _read_raster_bytes(_get_bytes(lid_fuel))

    print(f"Elevation grid: {elev_data.shape[0]}×{elev_data.shape[1]}")

    xs, ys, elev, slope, aspect, fuel = assemble_landscape(
        elev_data, elev_tf, elev_crs, elev_nd,
        slope_data, slope_tf, slope_crs,
        asp_data, asp_tf, asp_crs,
        fuel_data, fuel_tf, fuel_crs,
        project_utm=project_utm,
        subsample=subsample,
        keep_nonburnable=keep_nonburnable,
    )

    _write_lcp(output_path, xs, ys, elev, slope, aspect, fuel)


def create_landscape_from_files(output_path,
                                 elev_path,
                                 slope_path,
                                 aspect_path,
                                 fuel_path,
                                 project_utm=True,
                                 subsample=1,
                                 keep_nonburnable=False,
                                 elev_nodata=None):
    """Create a landscape file from local raster files (no download needed).

    Parameters
    ----------
    output_path : str
        Path for the output .lcp file.
    elev_path, slope_path, aspect_path, fuel_path : str
        Paths to the elevation, slope, aspect, and fuel model rasters.
    project_utm : bool
        Reproject coordinates to UTM metres (default: True).
    subsample : int
        Keep every N-th pixel in each dimension (default: 1).
    keep_nonburnable : bool
        If True, include non-burnable pixels in output.
    elev_nodata : float or None
        Override the no-data value for the elevation raster.
    """
    print("Reading elevation raster …")
    elev_data, elev_tf, elev_crs, elev_nd = _read_raster_file(
        elev_path, nodata_override=elev_nodata)
    print("Reading slope raster …")
    slope_data, slope_tf, slope_crs, _ = _read_raster_file(slope_path)
    print("Reading aspect raster …")
    asp_data, asp_tf, asp_crs, _        = _read_raster_file(aspect_path)
    print("Reading fuel model raster …")
    fuel_data, fuel_tf, fuel_crs, _     = _read_raster_file(fuel_path)

    print(f"Elevation grid: {elev_data.shape[0]}×{elev_data.shape[1]}")

    xs, ys, elev, slope, aspect, fuel = assemble_landscape(
        elev_data, elev_tf, elev_crs, elev_nd,
        slope_data, slope_tf, slope_crs,
        asp_data, asp_tf, asp_crs,
        fuel_data, fuel_tf, fuel_crs,
        project_utm=project_utm,
        subsample=subsample,
        keep_nonburnable=keep_nonburnable,
    )

    _write_lcp(output_path, xs, ys, elev, slope, aspect, fuel)


# ===========================================================================
# CLI
# ===========================================================================

def _build_parser():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # Bounding box (required)
    parser.add_argument("--lat-min", type=float, required=True,
                        help="Minimum latitude (WGS-84 decimal degrees)")
    parser.add_argument("--lat-max", type=float, required=True,
                        help="Maximum latitude (WGS-84 decimal degrees)")
    parser.add_argument("--lon-min", type=float, required=True,
                        help="Minimum longitude (WGS-84 decimal degrees)")
    parser.add_argument("--lon-max", type=float, required=True,
                        help="Maximum longitude (WGS-84 decimal degrees)")

    # Output files
    parser.add_argument("--terrain", default="terrain.xyz",
                        help="Output terrain XYZ file (default: terrain.xyz)")
    parser.add_argument("--landscape", default="landscape.lcp",
                        help="Output landscape LCP file (default: landscape.lcp)")
    parser.add_argument("--tif", default=None,
                        help="Path for the intermediate SRTM GeoTIFF "
                             "(temporary file used and removed if not specified)")

    # Control which outputs to produce
    parser.add_argument("--no-terrain", action="store_true",
                        help="Skip the SRTM terrain XYZ step")
    parser.add_argument("--no-landscape", action="store_true",
                        help="Skip the LANDFIRE landscape LCP step")

    # Common options
    parser.add_argument("--subsample", type=int, default=1, metavar="N",
                        help="Keep every N-th point in each dimension (default: 1)")

    # LANDFIRE-specific options
    parser.add_argument("--vintage", type=int, default=2020, metavar="YEAR",
                        help="LANDFIRE data vintage year (default: 2020)")
    parser.add_argument("--elev-product", default=None, metavar="ID",
                        help="Override LANDFIRE elevation product ID")
    parser.add_argument("--slope-product", default=None, metavar="ID",
                        help="Override LANDFIRE slope product ID")
    parser.add_argument("--aspect-product", default=None, metavar="ID",
                        help="Override LANDFIRE aspect product ID")
    parser.add_argument("--fuel-product", default=None, metavar="ID",
                        help="Override LANDFIRE fuel model product ID")
    parser.add_argument("--keep-nonburnable", action="store_true",
                        help="Include non-burnable LANDFIRE pixels in LCP output")
    parser.add_argument("--cache-dir", default=None, metavar="DIR",
                        help="Directory to cache downloaded LANDFIRE ZIP files")
    parser.add_argument("--timeout", type=int, default=300, metavar="N",
                        help="LANDFIRE API polling timeout in seconds (default: 300)")

    # Local-file landscape mode (skip LANDFIRE download)
    local = parser.add_argument_group(
        "local files mode (use pre-existing rasters instead of downloading)"
    )
    local.add_argument("--elev-file", default=None, metavar="PATH",
                       help="Local elevation raster (skips LFPS download)")
    local.add_argument("--slope-file", default=None, metavar="PATH",
                       help="Local slope raster")
    local.add_argument("--aspect-file", default=None, metavar="PATH",
                       help="Local aspect raster")
    local.add_argument("--fuel-file", default=None, metavar="PATH",
                       help="Local fuel model raster")

    return parser


def main(argv=None):
    parser = _build_parser()
    args = parser.parse_args(argv)

    if args.subsample < 1:
        print("ERROR: --subsample must be >= 1", file=sys.stderr)
        sys.exit(1)

    bbox = (args.lon_min, args.lat_min, args.lon_max, args.lat_max)

    # Validate bounding box
    if args.lon_min >= args.lon_max or args.lat_min >= args.lat_max:
        print(
            "ERROR: bounding box must satisfy lon-min < lon-max and lat-min < lat-max",
            file=sys.stderr,
        )
        sys.exit(1)

    # --- Terrain XYZ step ---
    if not args.no_terrain:
        tif_path = (os.path.abspath(args.tif)
                    if args.tif is not None else None)
        out_terrain = os.path.abspath(args.terrain)
        create_terrain_xyz(
            out_terrain,
            lat_min=args.lat_min,
            lat_max=args.lat_max,
            lon_min=args.lon_min,
            lon_max=args.lon_max,
            tif_path=tif_path,
            subsample=args.subsample,
        )

    # --- Landscape LCP step ---
    if not args.no_landscape:
        out_lcp = os.path.abspath(args.landscape)
        local_files = [args.elev_file, args.slope_file,
                       args.aspect_file, args.fuel_file]
        using_local = any(f is not None for f in local_files)

        if using_local:
            missing = []
            for name, val in [("--elev-file",   args.elev_file),
                               ("--slope-file",  args.slope_file),
                               ("--aspect-file", args.aspect_file),
                               ("--fuel-file",   args.fuel_file)]:
                if val is None:
                    missing.append(name)
            if missing:
                print(
                    f"ERROR: local-files mode requires all four rasters. "
                    f"Missing: {', '.join(missing)}",
                    file=sys.stderr,
                )
                sys.exit(1)
            create_landscape_from_files(
                out_lcp,
                args.elev_file,
                args.slope_file,
                args.aspect_file,
                args.fuel_file,
                subsample=args.subsample,
                keep_nonburnable=args.keep_nonburnable,
            )
        else:
            create_landscape(
                out_lcp,
                bbox=bbox,
                vintage=args.vintage,
                elev_product=args.elev_product,
                slope_product=args.slope_product,
                aspect_product=args.aspect_product,
                fuel_product=args.fuel_product,
                subsample=args.subsample,
                keep_nonburnable=args.keep_nonburnable,
                cache_dir=args.cache_dir,
                timeout_s=args.timeout,
            )

        # Report UTM extents
        if os.path.isfile(out_lcp):
            xs, ys = [], []
            with open(out_lcp) as fh:
                for line in fh:
                    if line.startswith("#") or not line.strip():
                        continue
                    parts = line.split()
                    if len(parts) >= 2:
                        try:
                            xs.append(float(parts[0]))
                            ys.append(float(parts[1]))
                        except ValueError:
                            pass
            if xs and ys:
                print(
                    f"\nUTM extents of landscape data (use for prob_lo / prob_hi):\n"
                    f"  prob_lo_x = {min(xs):.0f}\n"
                    f"  prob_lo_y = {min(ys):.0f}\n"
                    f"  prob_hi_x = {max(xs):.0f}\n"
                    f"  prob_hi_y = {max(ys):.0f}"
                )

    print("\nDone.")


if __name__ == "__main__":
    main()
