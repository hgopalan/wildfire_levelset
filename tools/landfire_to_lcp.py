#!/usr/bin/env python3
"""
landfire_to_lcp.py - Download LANDFIRE rasters and create a landscape file.

Downloads LANDFIRE elevation, slope, aspect, and fuel model rasters for a
user-specified bounding box via the LANDFIRE Product Service (LFPS) API and
assembles them into an ASCII landscape file (.lcp) compatible with
wildfire_levelset.

Output landscape file format (X Y ELEVATION SLOPE ASPECT FUEL_MODEL):
  # X Y ELEVATION SLOPE ASPECT FUEL_MODEL
  x1  y1  elev1  slope1  aspect1  fuel1
  x2  y2  elev2  slope2  aspect2  fuel2
  ...

Where:
  - X, Y          : coordinates in metres (UTM by default)
  - ELEVATION     : elevation above sea level in metres
  - SLOPE         : slope angle in degrees (0-90)
  - ASPECT        : slope aspect in degrees (0-360, 0=North, 90=East)
  - FUEL_MODEL    : Anderson 13-model NFFL fuel model number (1-13)

Non-burnable LANDFIRE classes (urban, water, agriculture, snow/ice, barren)
are mapped to 0 and excluded from the output by default (--keep-nonburnable
keeps them as 0).

Usage:
  python3 landfire_to_lcp.py --bbox MIN_LON MIN_LAT MAX_LON MAX_LAT output.lcp
  python3 landfire_to_lcp.py --bbox -120.5 37.0 -120.0 37.5 landscape.lcp
  python3 landfire_to_lcp.py --bbox -120.5 37.0 -120.0 37.5 landscape.lcp \\
      --subsample 4 --project-utm

Options:
  --bbox MIN_LON MIN_LAT MAX_LON MAX_LAT  Bounding box in decimal degrees (WGS-84)
  --subsample N      Keep every N-th point in each dimension (default: 1)
  --project-utm      Project coordinates to UTM metres (default: True)
  --no-utm           Keep native raster coordinates (no UTM reprojection)
  --fuel-product P   LANDFIRE fuel model product ID (default: F13_FBFM13)
  --elev-product P   LANDFIRE elevation product ID (default: ELEV2020)
  --vintage YEAR     LANDFIRE data vintage year (default: 2020)
  --cache-dir DIR    Directory to save downloaded ZIP files (default: system temp)
  --timeout N        API polling timeout in seconds (default: 300)
  --keep-nonburnable Include non-burnable pixels (fuel code 0) in output
  --help             Show this message and exit.

LANDFIRE product IDs vary by vintage.  Common values for 2020 data:
  Elevation : ELEV2020
  Slope     : SlpD2020
  Aspect    : Asp2020
  FBFM13    : F13_FBFM13
  FBFM40    : F40_FBFM40

Requires: rasterio, numpy (and pyproj for UTM reprojection)
Install : pip install rasterio numpy pyproj
"""

import argparse
import io
import json
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
        "elev":    "ELEV2020",
        "slope":   "SlpD2020",
        "aspect":  "Asp2020",
        "fuel13":  "F13_FBFM13",
    },
    2016: {
        "elev":    "ELEV2016",
        "slope":   "SlpD2016",
        "aspect":  "Asp2016",
        "fuel13":  "F13_FBFM13",
    },
    2014: {
        "elev":    "ELEV2014",
        "slope":   "SLP",
        "aspect":  "ASP",
        "fuel13":  "F13_FBFM13",
    },
}

# LANDFIRE FBFM13 non-burnable codes (outside the 1-13 Anderson range)
_NONBURNABLE_CODES = {91, 92, 93, 98, 99}

# ---------------------------------------------------------------------------
# LFPS API helpers
# ---------------------------------------------------------------------------

def _submit_lfps_job(bbox, layer_ids, timeout_s=300):
    """Submit a LANDFIRE Product Service job and poll until completion.

    Parameters
    ----------
    bbox : tuple of float
        (min_lon, min_lat, max_lon, max_lat) in WGS-84 decimal degrees.
    layer_ids : list of str
        LANDFIRE layer product IDs to request (e.g. ["ELEV2020", "SlpD2020"]).
    timeout_s : int
        Maximum seconds to wait for the job to complete.

    Returns
    -------
    str
        URL pointing to the output ZIP file.

    Raises
    ------
    RuntimeError
        If the job fails or times out.
    ImportError
        If the ``requests`` package is not available.
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
    layer_list = "|".join(layer_ids)

    payload = {
        "Area_of_Interest": aoi,
        "Layer_list":       layer_list,
        "f":                "json",
    }

    print(f"Submitting LANDFIRE job for bbox {bbox} …")
    resp = requests.post(_SUBMIT_URL, data=payload, timeout=60)
    resp.raise_for_status()
    job_info = resp.json()

    if "jobId" not in job_info:
        raise RuntimeError(
            f"LANDFIRE job submission failed: {job_info}"
        )

    job_id = job_info["jobId"]
    print(f"  Job ID: {job_id}")

    # Poll until done
    status_url = _JOB_URL.format(jobId=job_id)
    deadline = time.time() + timeout_s
    # poll_interval grows from 5 s up to 30 s (exponential back-off) across
    # iterations; it is initialised once here, outside the loop, so each poll
    # builds on the previous interval rather than resetting.
    poll_interval = 5  # seconds between polls

    while time.time() < deadline:
        time.sleep(poll_interval)
        sr = requests.get(status_url, params={"f": "json"}, timeout=30)
        sr.raise_for_status()
        status = sr.json().get("jobStatus", "")
        print(f"  Status: {status}")

        if status == "esriJobSucceeded":
            break
        if status in ("esriJobFailed", "esriJobCancelled",
                      "esriJobCancelling", "esriJobDeleted"):
            raise RuntimeError(
                f"LANDFIRE job {job_id} ended with status: {status}"
            )
        # Increase poll interval gradually to reduce server load
        poll_interval = min(poll_interval * 1.5, 30)
    else:
        raise RuntimeError(
            f"LANDFIRE job {job_id} did not complete within {timeout_s}s"
        )

    # Retrieve output file URL
    result_url = _RESULT_URL.format(jobId=job_id)
    rr = requests.get(result_url, params={"f": "json"}, timeout=30)
    rr.raise_for_status()
    result = rr.json()

    download_url = (
        result.get("value", {}).get("url")
        or result.get("value")
    )
    if not download_url:
        raise RuntimeError(
            f"Could not parse download URL from LANDFIRE result: {result}"
        )

    return download_url


def _download_zip(url, cache_dir=None):
    """Download a ZIP file from *url* and return its bytes.

    If *cache_dir* is provided the ZIP is also saved there for future reuse.
    """
    try:
        import requests
    except ImportError:
        raise ImportError(
            "requests is required to download LANDFIRE data. "
            "Install with: pip install requests"
        )

    print(f"Downloading output ZIP from {url} …")
    with requests.get(url, stream=True, timeout=120) as r:
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
    dict mapping layer_id (str) to file bytes (bytes).
        Keys are the layer IDs supplied; values are the raw bytes of each
        GeoTIFF extracted from the ZIP.
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
            # Heuristically match zip entry to requested layer ID
            matched = False
            for lid in layer_ids:
                if lid.lower() in lower:
                    layer_map[lid] = raw
                    matched = True
                    break
            if not matched:
                # Use the filename stem as key
                stem = os.path.splitext(os.path.basename(name))[0]
                layer_map[stem] = raw

    return layer_map


# ---------------------------------------------------------------------------
# Raster reading helpers
# ---------------------------------------------------------------------------

def _read_raster_bytes(raw_bytes):
    """Read a raster from raw bytes using rasterio.

    Returns (data_2d, transform, crs, nodata) where *data_2d* is a 2-D
    numpy float64 array, *transform* is an Affine transform, *crs* is a
    CRS object, and *nodata* is the no-data sentinel (or None).
    """
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
    """Return (x_2d, y_2d) arrays of cell centre coordinates."""
    import numpy as np
    from rasterio.transform import xy as rasterio_xy

    nrows, ncols = data.shape
    rows_idx, cols_idx = np.meshgrid(np.arange(nrows), np.arange(ncols),
                                     indexing="ij")
    xs, ys = rasterio_xy(transform,
                         rows_idx.ravel(),
                         cols_idx.ravel())
    xs = np.array(xs, dtype=np.float64).reshape(nrows, ncols)
    ys = np.array(ys, dtype=np.float64).reshape(nrows, ncols)
    return xs, ys


def _project_to_utm(xs, ys, src_crs):
    """Project *xs*, *ys* (any CRS) to UTM metres.

    The target UTM zone is chosen from the centre of the data.
    Returns (x_utm, y_utm) in metres.
    """
    try:
        from pyproj import Transformer, CRS
    except ImportError:
        raise ImportError(
            "pyproj is required for UTM reprojection. "
            "Install with: pip install pyproj"
        )
    import numpy as np

    # Ensure we have a CRS object
    try:
        from rasterio.crs import CRS as RasterioCRS
        if isinstance(src_crs, RasterioCRS):
            src_epsg = src_crs.to_epsg()
            if src_epsg:
                src_authority = f"EPSG:{src_epsg}"
            else:
                src_authority = src_crs.to_wkt()
        else:
            src_authority = str(src_crs)
    except Exception:
        src_authority = "EPSG:4326"

    # Determine UTM zone from centre
    # First project to geographic if not already
    try:
        to_geo = Transformer.from_crs(src_authority, "EPSG:4326", always_xy=True)
        lons, lats = to_geo.transform(xs.ravel(), ys.ravel())
    except Exception:
        lons, lats = xs.ravel(), ys.ravel()

    center_lon = float(np.nanmean(lons))
    center_lat = float(np.nanmean(lats))
    zone = int((center_lon + 180.0) / 6.0) + 1
    hemisphere = "north" if center_lat >= 0.0 else "south"
    epsg = 32600 + zone if hemisphere == "north" else 32700 + zone

    transformer = Transformer.from_crs(src_authority, f"EPSG:{epsg}",
                                       always_xy=True)
    x_utm, y_utm = transformer.transform(xs.ravel(), ys.ravel())
    return (np.array(x_utm, dtype=np.float64).reshape(xs.shape),
            np.array(y_utm, dtype=np.float64).reshape(xs.shape))


# ---------------------------------------------------------------------------
# LCP assembly
# ---------------------------------------------------------------------------

def _resample_to_grid(src_data, src_transform, src_crs,
                      ref_data, ref_transform, ref_crs,
                      resampling="nearest"):
    """Resample *src* raster to the grid of *ref* raster using rasterio.warp."""
    import numpy as np
    try:
        import rasterio
        from rasterio.warp import reproject, Resampling
        from rasterio.io import MemoryFile
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

    # Resample slope, aspect, and fuel to elevation grid
    if slope_data.shape != elev_data.shape:
        print("  Resampling slope to elevation grid …")
        slope_data = _resample_to_grid(
            slope_data, slope_transform, slope_crs,
            elev_data, elev_transform, elev_crs,
            resampling="bilinear",
        )
    if aspect_data.shape != elev_data.shape:
        print("  Resampling aspect to elevation grid …")
        aspect_data = _resample_to_grid(
            aspect_data, aspect_transform, aspect_crs,
            elev_data, elev_transform, elev_crs,
            resampling="bilinear",
        )
    if fuel_data.shape != elev_data.shape:
        print("  Resampling fuel model to elevation grid …")
        fuel_data = _resample_to_grid(
            fuel_data, fuel_transform, fuel_crs,
            elev_data, elev_transform, elev_crs,
            resampling="nearest",
        )

    # Build x/y coordinate arrays
    xs, ys = _raster_coords(elev_data, elev_transform)

    if project_utm:
        print("  Projecting coordinates to UTM metres …")
        xs, ys = _project_to_utm(xs, ys, elev_crs)

    # Clamp fuel codes and mark non-burnable
    fuel_int = np.round(fuel_data).astype(int)

    # Build valid-pixel mask
    nodata_val = elev_nodata if elev_nodata is not None else -9999.0
    mask = elev_data != nodata_val
    mask &= np.isfinite(elev_data)

    if not keep_nonburnable:
        for code in _NONBURNABLE_CODES:
            mask &= fuel_int != code
        # Also exclude fuel codes outside 1-13
        mask &= (fuel_int >= 1) & (fuel_int <= 13)

    # Apply subsample
    if subsample > 1:
        skip = np.zeros(elev_data.shape, dtype=bool)
        skip[::subsample, ::subsample] = True
        mask &= skip

    xs_flat     = xs[mask]
    ys_flat     = ys[mask]
    elev_flat   = elev_data[mask]
    slope_flat  = slope_data[mask]
    aspect_flat = aspect_data[mask]
    fuel_flat   = fuel_int[mask].astype(float)

    return xs_flat, ys_flat, elev_flat, slope_flat, aspect_flat, fuel_flat


# ---------------------------------------------------------------------------
# Writer
# ---------------------------------------------------------------------------

def _write_lcp(path, xs, ys, elev, slope, aspect, fuel):
    """Write an ASCII landscape file (.lcp)."""
    with open(path, "w") as fh:
        fh.write(
            "# Landscape file generated by landfire_to_lcp.py\n"
            "# Format: X Y ELEVATION SLOPE ASPECT FUEL_MODEL\n"
            "# Units: X/Y metres, ELEVATION metres, SLOPE degrees,\n"
            "#        ASPECT degrees (0=North), FUEL_MODEL NFFL 1-13\n"
        )
        for x, y, z, s, a, f in zip(xs, ys, elev, slope, aspect, fuel):
            fh.write(
                f"{x:.2f} {y:.2f} {z:.2f} {s:.2f} {a:.2f} {int(f)}\n"
            )
    print(f"Wrote {len(xs)} landscape points to '{path}'.")


# ---------------------------------------------------------------------------
# High-level pipeline
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
        If None, the defaults for *vintage* are used.
    project_utm : bool
        Reproject coordinates to UTM metres (default: True).
    subsample : int
        Keep every N-th pixel in each dimension (default: 1 = all pixels).
    keep_nonburnable : bool
        If True, include non-burnable pixels (fuel code 0) in output.
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

    # Download from LANDFIRE
    layer_map = download_landfire(bbox, layer_ids,
                                  cache_dir=cache_dir, timeout_s=timeout_s)

    def _get_bytes(lid):
        if lid in layer_map:
            return layer_map[lid]
        # Try case-insensitive partial match
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
        Keep every N-th pixel in each dimension (default: 1 = all pixels).
    keep_nonburnable : bool
        If True, include non-burnable pixels (fuel code 0) in output.
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


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _build_parser():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # Required arguments
    parser.add_argument(
        "--bbox", nargs=4, type=float,
        metavar=("MIN_LON", "MIN_LAT", "MAX_LON", "MAX_LAT"),
        help="Bounding box in WGS-84 decimal degrees",
    )
    parser.add_argument(
        "output",
        help="Output landscape file (.lcp)",
    )

    # Optional arguments
    parser.add_argument(
        "--vintage", type=int, default=2020,
        help="LANDFIRE data vintage year (default: 2020)",
    )
    parser.add_argument(
        "--elev-product", default=None, metavar="ID",
        help="LANDFIRE elevation product ID (overrides vintage default)",
    )
    parser.add_argument(
        "--slope-product", default=None, metavar="ID",
        help="LANDFIRE slope product ID (overrides vintage default)",
    )
    parser.add_argument(
        "--aspect-product", default=None, metavar="ID",
        help="LANDFIRE aspect product ID (overrides vintage default)",
    )
    parser.add_argument(
        "--fuel-product", default=None, metavar="ID",
        help="LANDFIRE fuel model product ID (overrides vintage default)",
    )
    parser.add_argument(
        "--subsample", type=int, default=1, metavar="N",
        help="Keep every N-th point in each dimension (default: 1)",
    )
    parser.add_argument(
        "--no-utm", action="store_true",
        help="Do not reproject coordinates to UTM metres",
    )
    parser.add_argument(
        "--keep-nonburnable", action="store_true",
        help="Include non-burnable pixels (fuel code 0) in output",
    )
    parser.add_argument(
        "--cache-dir", default=None, metavar="DIR",
        help="Directory to cache downloaded ZIP files",
    )
    parser.add_argument(
        "--timeout", type=int, default=300, metavar="N",
        help="API polling timeout in seconds (default: 300)",
    )

    # Local-files mode (no download)
    local_group = parser.add_argument_group(
        "local files mode (skip download, use pre-existing rasters)"
    )
    local_group.add_argument(
        "--elev-file", default=None, metavar="PATH",
        help="Path to local elevation raster (skips LFPS download)",
    )
    local_group.add_argument(
        "--slope-file", default=None, metavar="PATH",
        help="Path to local slope raster",
    )
    local_group.add_argument(
        "--aspect-file", default=None, metavar="PATH",
        help="Path to local aspect raster",
    )
    local_group.add_argument(
        "--fuel-file", default=None, metavar="PATH",
        help="Path to local fuel model raster",
    )

    return parser


def main(argv=None):
    parser = _build_parser()
    args = parser.parse_args(argv)

    if args.subsample < 1:
        print("ERROR: --subsample must be >= 1", file=sys.stderr)
        sys.exit(1)

    local_files = [args.elev_file, args.slope_file,
                   args.aspect_file, args.fuel_file]
    using_local = any(f is not None for f in local_files)

    if using_local:
        # Local-files mode: all four must be supplied
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
            args.output,
            args.elev_file,
            args.slope_file,
            args.aspect_file,
            args.fuel_file,
            project_utm=not args.no_utm,
            subsample=args.subsample,
            keep_nonburnable=args.keep_nonburnable,
        )
    else:
        # Download mode: bbox is required
        if args.bbox is None:
            print(
                "ERROR: --bbox is required unless using --elev-file / "
                "--slope-file / --aspect-file / --fuel-file.",
                file=sys.stderr,
            )
            parser.print_help(sys.stderr)
            sys.exit(1)

        bbox = tuple(args.bbox)
        min_lon, min_lat, max_lon, max_lat = bbox
        if min_lon >= max_lon or min_lat >= max_lat:
            print(
                "ERROR: --bbox must satisfy MIN_LON < MAX_LON and "
                "MIN_LAT < MAX_LAT",
                file=sys.stderr,
            )
            sys.exit(1)

        create_landscape(
            args.output,
            bbox=bbox,
            vintage=args.vintage,
            elev_product=args.elev_product,
            slope_product=args.slope_product,
            aspect_product=args.aspect_product,
            fuel_product=args.fuel_product,
            project_utm=not args.no_utm,
            subsample=args.subsample,
            keep_nonburnable=args.keep_nonburnable,
            cache_dir=args.cache_dir,
            timeout_s=args.timeout,
        )


if __name__ == "__main__":
    main()
