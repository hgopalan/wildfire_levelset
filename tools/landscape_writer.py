#!/usr/bin/env python3
"""
landscape_writer.py - Download LANDFIRE rasters and write a FARSITE landscape (.lcp) file.

Fetches LANDFIRE elevation, slope, aspect, and fuel model rasters for a
lat/lon bounding box and writes an ASCII landscape file (``.lcp``) suitable
for use as ``rothermel.landscape_file`` in the wildfire_levelset solver.

Three remote sources are supported (tried in order):
  1. **cog**  – Cloud Optimized GeoTIFFs from the public LANDFIRE AWS S3 bucket.
  2. **mspc** – Microsoft Planetary Computer STAC catalog (Azure Blob COGs).
  3. **lfps** – Legacy LANDFIRE Product Service REST API hosted by USGS.

Alternatively, local GeoTIFF files may be supplied for all four layers
(elevation, slope, aspect, fuel model) to skip any network download.

Usage examples
--------------
  # Download LANDFIRE (API) for a bounding box
  python3 landscape_writer.py \\
      --lat-min 40 --lat-max 40.5 \\
      --lon-min -105 --lon-max -104.5

  # Use local raster files (no download)
  python3 landscape_writer.py \\
      --lat-min 40 --lat-max 40.5 \\
      --lon-min -105 --lon-max -104.5 \\
      --elev-file elev.tif \\
      --slope-file slope.tif \\
      --aspect-file aspect.tif \\
      --fuel-file fbfm13.tif

  # Use local fuel raster with SRTM-derived elevation/slope/aspect
  python3 landscape_writer.py \\
      --lat-min 40 --lat-max 40.5 \\
      --lon-min -105 --lon-max -104.5 \\
      --fuel-file fuel.tif \\
      --srtm-slope-aspect

  # Scott & Burgan 40 fuel model system
  python3 landscape_writer.py \\
      --lat-min 40 --lat-max 40.5 \\
      --lon-min -105 --lon-max -104.5 \\
      --fuel-model-type 40

Options
-------
  --lat-min / --lat-max   Latitude  bounds (WGS-84 degrees; required).
  --lon-min / --lon-max   Longitude bounds (WGS-84 degrees; required).
  --landscape FILE        Output landscape LCP file (default: landscape.lcp).
  --subsample N           Keep every N-th point (default: 1).
  --vintage YEAR          LANDFIRE vintage year (default: 2020).
  --fuel-model-type {13,40}
                          Fuel model system (default: 13).
  --keep-nonburnable      Include non-burnable LANDFIRE pixels.
  --cache-dir DIR         Cache directory for LANDFIRE ZIP files (lfps only).
  --timeout N             LANDFIRE API polling timeout in seconds (default: 300).
  --use-lfps              Force the legacy LFPS REST API.
  --sources LIST          Comma-separated source priority list (cog,mspc,lfps).
  --no-vintage-fallback   Disable automatic vintage fallback.
  --elev-product ID       Override LANDFIRE elevation product ID.
  --slope-product ID      Override LANDFIRE slope product ID.
  --aspect-product ID     Override LANDFIRE aspect product ID.
  --fuel-product ID       Override LANDFIRE fuel model product ID.
  --elev-file PATH        Local elevation raster (skips download).
  --slope-file PATH       Local slope raster.
  --aspect-file PATH      Local aspect raster.
  --fuel-file PATH        Local fuel model raster.
  --srtm-slope-aspect     Derive elevation/slope/aspect from SRTM when
                          --fuel-file is given without --elev-file etc.
  --tif FILE              Intermediate SRTM GeoTIFF path (for --srtm-slope-aspect).
  --help                  Show this message and exit.

Requires: pip install rasterio numpy pyproj requests
Optional: pip install elevation (for --srtm-slope-aspect)
          pip install pystac-client planetary-computer (for mspc source)
"""

import argparse
import io
import json
import math
import os
import re
import sys
import tempfile
import time
import zipfile

# ---------------------------------------------------------------------------
# Supported LANDFIRE source names
# ---------------------------------------------------------------------------

_LANDFIRE_SOURCES = ("cog", "mspc", "lfps")
_VINTAGE_FALLBACK_ORDER = (2020, 2016, 2014)

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

_LFPS_SSL_WARNED = False

# ---------------------------------------------------------------------------
# LANDFIRE Cloud Optimized GeoTIFF (COG) constants
# ---------------------------------------------------------------------------

_LANDFIRE_COG_BASE = "https://s3.amazonaws.com/landfire/"

_COG_LAYERS = {
    2020: {
        "elev":   "US_200/US_200ELEV.tif",
        "slope":  "US_200/US_200SlpD.tif",
        "aspect": "US_200/US_200Asp.tif",
        "fuel13": "US_200/US_200FBFM13.tif",
        "fuel40": "US_200/US_200FBFM40.tif",
    },
    2016: {
        "elev":   "US_140/US_140ELEV.tif",
        "slope":  "US_140/US_140SlpD.tif",
        "aspect": "US_140/US_140Asp.tif",
        "fuel13": "US_140/US_140FBFM13.tif",
        "fuel40": "US_140/US_140FBFM40.tif",
    },
    2014: {
        "elev":   "US_130/US_130ELEV.tif",
        "slope":  "US_130/US_130SLP.tif",
        "aspect": "US_130/US_130ASP.tif",
        "fuel13": "US_130/US_130FBFM13.tif",
        "fuel40": "US_130/US_130FBFM40.tif",
    },
}

_DEFAULT_LAYERS = {
    2020: {
        "elev":   "ELEV2020",
        "slope":  "SlpD2020",
        "aspect": "Asp2020",
        "fuel13": "F13_FBFM13",
        "fuel40": "F40_FBFM40",
    },
    2016: {
        "elev":   "ELEV2016",
        "slope":  "SlpD2016",
        "aspect": "Asp2016",
        "fuel13": "F13_FBFM13",
        "fuel40": "F40_FBFM40",
    },
    2014: {
        "elev":   "ELEV2014",
        "slope":  "SLP",
        "aspect": "ASP",
        "fuel13": "F13_FBFM13",
        "fuel40": "F40_FBFM40",
    },
}

_NONBURNABLE_CODES = {91, 92, 93, 98, 99}

# ---------------------------------------------------------------------------
# Microsoft Planetary Computer (MSPC) LANDFIRE constants
# ---------------------------------------------------------------------------

_MSPC_STAC_ENDPOINT = "https://planetarycomputer.microsoft.com/api/stac/v1"
_MSPC_LANDFIRE_COLLECTION = "landfire"

_MSPC_ASSET_KEYS = {
    2020: {
        "elev":   ["US_200ELEV",   "200ELEV",   "elev"],
        "slope":  ["US_200SlpD",   "200SlpD",   "slope"],
        "aspect": ["US_200Asp",    "200Asp",    "aspect"],
        "fuel13": ["US_200FBFM13", "200FBFM13", "fuel13"],
        "fuel40": ["US_200FBFM40", "200FBFM40", "fuel40"],
    },
    2016: {
        "elev":   ["US_140ELEV",   "140ELEV",   "elev"],
        "slope":  ["US_140SlpD",   "140SlpD",   "slope"],
        "aspect": ["US_140Asp",    "140Asp",    "aspect"],
        "fuel13": ["US_140FBFM13", "140FBFM13", "fuel13"],
        "fuel40": ["US_140FBFM40", "140FBFM40", "fuel40"],
    },
    2014: {
        "elev":   ["US_130ELEV",   "130ELEV",   "elev"],
        "slope":  ["US_130SLP",    "130SLP",    "slope"],
        "aspect": ["US_130ASP",    "130ASP",    "aspect"],
        "fuel13": ["US_130FBFM13", "130FBFM13", "fuel13"],
        "fuel40": ["US_130FBFM40", "130FBFM40", "fuel40"],
    },
}


# ===========================================================================
# Shared projection helper
# ===========================================================================

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

    try:
        from rasterio.crs import CRS as RasterioCRS
        if isinstance(src_crs, RasterioCRS):
            src_epsg = src_crs.to_epsg()
            src_authority = f"EPSG:{src_epsg}" if src_epsg else src_crs.to_wkt()
        else:
            src_authority = str(src_crs)
    except Exception:
        src_authority = "EPSG:4326"

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


# ===========================================================================
# LANDFIRE LFPS API helpers
# ===========================================================================

def _submit_lfps_job(bbox, layer_ids, timeout_s=300):
    """Submit a LANDFIRE Product Service job and poll until completion."""
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
    resp = _requests_post(_SUBMIT_URL, data=payload, params={"f": "json"}, timeout=60)
    resp.raise_for_status()
    job_info = _parse_json_response(resp, "LANDFIRE job submission")

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
        status = _parse_json_response(sr, "LANDFIRE job status").get("jobStatus", "")
        print(f"  Status: {status}")

        if status == "esriJobSucceeded":
            break
        if status in ("esriJobFailed", "esriJobCancelled",
                      "esriJobCancelling", "esriJobDeleted"):
            raise RuntimeError(
                f"LANDFIRE job {job_id} ended with status: {status}")
        poll_interval = min(poll_interval * 1.5, 30)
    else:
        raise RuntimeError(
            f"LANDFIRE job {job_id} did not complete within {timeout_s}s")

    result_url = _RESULT_URL.format(jobId=job_id)
    rr = _requests_get(result_url, params={"f": "json"}, timeout=30)
    rr.raise_for_status()
    result = _parse_json_response(rr, "LANDFIRE job result")

    download_url = (
        result.get("value", {}).get("url")
        or result.get("value")
    )
    if not download_url:
        raise RuntimeError(
            f"Could not parse download URL from LANDFIRE result: {result}")
    return download_url


def _parse_json_response(resp, context="LANDFIRE API"):
    """Parse JSON from *resp*, raising RuntimeError on empty or invalid body."""
    text = resp.text.strip() if resp.text else ""
    if not text:
        raise RuntimeError(
            f"{context} returned an empty response "
            f"(HTTP {resp.status_code}). "
            "The LFPS service may be temporarily unavailable; "
            "retry later or check https://lfps.usgs.gov for service status."
        )
    content_type = resp.headers.get("Content-Type", "")
    if "text/html" in content_type or text.lstrip().startswith("<"):
        m = re.search(r"<title[^>]*>([^<]+)</title>", text, re.IGNORECASE)
        title_hint = f" (page title: {m.group(1).strip()!r})" if m else ""
        raise RuntimeError(
            f"{context} returned an HTML page instead of JSON "
            f"(HTTP {resp.status_code}){title_hint}. "
            "The LFPS service may be temporarily unavailable or is returning "
            "a maintenance/error page. Retry later or check "
            "https://lfps.usgs.gov for service status."
        )
    try:
        return json.loads(text)
    except json.JSONDecodeError as exc:
        raise RuntimeError(
            f"{context} returned non-JSON response "
            f"(HTTP {resp.status_code}): {text[:500]!r}"
        ) from exc


def _requests_get(url, **kwargs):
    import requests
    if url.startswith(_LFPS_BASE):
        kwargs.setdefault("verify", False)
        _suppress_lfps_ssl_warning()
    return requests.get(url, **kwargs)


def _requests_post(url, **kwargs):
    import requests
    if url.startswith(_LFPS_BASE):
        kwargs.setdefault("verify", False)
        _suppress_lfps_ssl_warning()
    return requests.post(url, **kwargs)


def _suppress_lfps_ssl_warning():
    global _LFPS_SSL_WARNED
    if _LFPS_SSL_WARNED:
        return
    _LFPS_SSL_WARNED = True
    print(
        "WARNING: SSL certificate verification disabled for lfps.usgs.gov "
        "(self-signed certificate in chain). "
        "Ensure you trust this host before proceeding."
    )
    try:
        import urllib3
        urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
    except (ImportError, AttributeError):
        pass


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
    """Download LANDFIRE rasters and return a dict of {layer_id: bytes}."""
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


def _read_landfire_cog(cog_url, bbox):
    """Read a LANDFIRE COG clipped to *bbox* (min_lon, min_lat, max_lon, max_lat)."""
    import numpy as np
    try:
        import rasterio
        from rasterio.windows import from_bounds
    except ImportError:
        raise ImportError(
            "rasterio is required to read LANDFIRE COGs. "
            "Install with: pip install rasterio"
        )

    min_lon, min_lat, max_lon, max_lat = bbox

    env_opts = {
        "GDAL_HTTP_MAX_RETRY": "3",
        "GDAL_HTTP_RETRY_DELAY": "2",
        "GDAL_DISABLE_READDIR_ON_OPEN": "EMPTY_DIR",
        "CPL_VSIL_CURL_ALLOWED_EXTENSIONS": ".tif,.tiff",
    }

    print(f"  Reading COG window from {cog_url} …")
    with rasterio.Env(**env_opts):
        with rasterio.open(cog_url) as ds:
            src_crs = ds.crs
            nodata = ds.nodata

            if src_crs.to_epsg() == 4326:
                win_min_x, win_min_y = min_lon, min_lat
                win_max_x, win_max_y = max_lon, max_lat
            else:
                try:
                    from pyproj import Transformer
                    t = Transformer.from_crs(
                        "EPSG:4326", src_crs, always_xy=True
                    )
                    corners_x, corners_y = t.transform(
                        [min_lon, max_lon, min_lon, max_lon],
                        [min_lat, min_lat, max_lat, max_lat],
                    )
                    win_min_x, win_max_x = min(corners_x), max(corners_x)
                    win_min_y, win_max_y = min(corners_y), max(corners_y)
                except Exception as exc:
                    raise RuntimeError(
                        f"Failed to transform bbox to COG CRS ({src_crs}): {exc}"
                    ) from exc

            window = from_bounds(
                win_min_x, win_min_y, win_max_x, win_max_y, ds.transform
            )
            win_transform = ds.window_transform(window)
            data = ds.read(1, window=window).astype(np.float64)

    return data, win_transform, src_crs, nodata


def download_landfire_cog(bbox, vintage=2020,
                          lid_elev=None, lid_slope=None,
                          lid_aspect=None, lid_fuel=None,
                          fuel_type="13"):
    """Download LANDFIRE rasters from Cloud Optimized GeoTIFFs on AWS S3."""
    try:
        import rasterio
        from rasterio.io import MemoryFile
    except ImportError:
        raise ImportError(
            "rasterio is required to read LANDFIRE COGs. "
            "Install with: pip install rasterio"
        )

    cog_vintage = vintage if vintage in _COG_LAYERS else 2020
    if vintage not in _COG_LAYERS:
        print(
            f"WARNING: No COG layer mapping for vintage {vintage}; "
            f"using {cog_vintage}.",
            file=sys.stderr,
        )

    layer_defaults = _DEFAULT_LAYERS.get(vintage, _DEFAULT_LAYERS[2020])
    cog_paths = _COG_LAYERS[cog_vintage]

    fuel_key = "fuel40" if fuel_type == "40" else "fuel13"
    layers = {
        lid_elev   or layer_defaults["elev"]:      cog_paths["elev"],
        lid_slope  or layer_defaults["slope"]:     cog_paths["slope"],
        lid_aspect or layer_defaults["aspect"]:    cog_paths["aspect"],
        lid_fuel   or layer_defaults[fuel_key]:    cog_paths[fuel_key],
    }

    result = {}
    for lid, rel_path in layers.items():
        cog_url = _LANDFIRE_COG_BASE + rel_path
        data, transform, crs, nodata = _read_landfire_cog(cog_url, bbox)

        with MemoryFile() as mf:
            with mf.open(
                driver="GTiff",
                height=data.shape[0],
                width=data.shape[1],
                count=1,
                dtype=data.dtype,
                crs=crs,
                transform=transform,
                nodata=nodata,
            ) as ds:
                ds.write(data, 1)
            result[lid] = mf.read()

    return result


def download_landfire_mspc(bbox, vintage=2020,
                           lid_elev=None, lid_slope=None,
                           lid_aspect=None, lid_fuel=None,
                           fuel_type="13"):
    """Download LANDFIRE rasters from Microsoft Planetary Computer (MSPC)."""
    try:
        import pystac_client
    except ImportError:
        raise ImportError(
            "pystac-client is required to access Microsoft Planetary Computer. "
            "Install with: pip install pystac-client planetary-computer"
        )
    try:
        import planetary_computer
    except ImportError:
        raise ImportError(
            "planetary-computer is required to sign MSPC asset URLs. "
            "Install with: pip install planetary-computer"
        )
    try:
        import rasterio
        from rasterio.io import MemoryFile
    except ImportError:
        raise ImportError(
            "rasterio is required to read LANDFIRE COGs. "
            "Install with: pip install rasterio"
        )

    mspc_vintage = vintage if vintage in _MSPC_ASSET_KEYS else 2020
    if vintage not in _MSPC_ASSET_KEYS:
        print(
            f"WARNING: No MSPC asset key mapping for vintage {vintage}; "
            f"using {mspc_vintage}.",
            file=sys.stderr,
        )

    layer_defaults = _DEFAULT_LAYERS.get(vintage, _DEFAULT_LAYERS[2020])
    asset_key_candidates = _MSPC_ASSET_KEYS[mspc_vintage]

    fuel_key = "fuel40" if fuel_type == "40" else "fuel13"
    lid_map = {
        lid_elev   or layer_defaults["elev"]:      "elev",
        lid_slope  or layer_defaults["slope"]:     "slope",
        lid_aspect or layer_defaults["aspect"]:    "aspect",
        lid_fuel   or layer_defaults[fuel_key]:    fuel_key,
    }

    print("Querying Microsoft Planetary Computer STAC for LANDFIRE items …")
    catalog = pystac_client.Client.open(
        _MSPC_STAC_ENDPOINT,
        modifier=planetary_computer.sign_inplace,
    )

    min_lon, min_lat, max_lon, max_lat = bbox
    search = catalog.search(
        collections=[_MSPC_LANDFIRE_COLLECTION],
        bbox=[min_lon, min_lat, max_lon, max_lat],
        max_items=10,
    )
    items = list(search.items())
    if not items:
        raise RuntimeError(
            f"No LANDFIRE items found in the MSPC STAC catalog for "
            f"bbox {bbox}.  The collection '{_MSPC_LANDFIRE_COLLECTION}' may "
            "not cover this region or the catalog may be temporarily "
            "unavailable.  Try a different --sources option."
        )

    print(f"  Found {len(items)} MSPC LANDFIRE item(s); using first matching.")

    def _find_asset_url(item, candidates):
        for key in candidates:
            if key in item.assets:
                return item.assets[key].href
            for asset_key in item.assets:
                if asset_key.lower() == key.lower():
                    return item.assets[asset_key].href
        available = list(item.assets.keys())
        raise RuntimeError(
            f"Could not find any of the expected asset keys "
            f"{candidates!r} in MSPC LANDFIRE item '{item.id}'. "
            f"Available assets: {available}"
        )

    result = {}
    for lid, layer_key in lid_map.items():
        candidates = asset_key_candidates[layer_key]
        asset_url = None
        last_exc = None
        for item in items:
            try:
                asset_url = _find_asset_url(item, candidates)
                break
            except RuntimeError as exc:
                last_exc = exc
                continue

        if asset_url is None:
            raise RuntimeError(
                f"None of the {len(items)} MSPC item(s) contained assets "
                f"matching {candidates!r} for layer '{lid}'."
            ) from last_exc

        print(f"  Reading MSPC COG for layer '{lid}' from {asset_url[:80]} …")
        data, transform, crs, nodata = _read_landfire_cog(asset_url, bbox)

        with MemoryFile() as mf:
            with mf.open(
                driver="GTiff",
                height=data.shape[0],
                width=data.shape[1],
                count=1,
                dtype=data.dtype,
                crs=crs,
                transform=transform,
                nodata=nodata,
            ) as ds:
                ds.write(data, 1)
            result[lid] = mf.read()

    return result


# ===========================================================================
# Raster I/O helpers
# ===========================================================================

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


def _read_raster_file_clipped(path, lat_min, lat_max, lon_min, lon_max,
                               nodata_override=None):
    """Read only the portion of a raster that overlaps a WGS-84 bounding box."""
    import math as _math
    import numpy as np
    try:
        import rasterio
        from rasterio import windows as rio_windows
    except ImportError:
        raise ImportError(
            "rasterio is required to process LANDFIRE rasters. "
            "Install with: pip install rasterio"
        )

    with rasterio.open(path) as ds:
        src_crs = ds.crs

        corners_lon = [lon_min, lon_min, lon_max, lon_max]
        corners_lat = [lat_min, lat_max, lat_min, lat_max]

        if src_crs is not None and src_crs.to_epsg() != 4326:
            try:
                from pyproj import Transformer
                tf = Transformer.from_crs("EPSG:4326", src_crs, always_xy=True)
                xs_proj, ys_proj = tf.transform(corners_lon, corners_lat)
            except Exception:
                xs_proj, ys_proj = corners_lon, corners_lat
        else:
            xs_proj, ys_proj = corners_lon, corners_lat

        x_min_proj = min(xs_proj)
        x_max_proj = max(xs_proj)
        y_min_proj = min(ys_proj)
        y_max_proj = max(ys_proj)

        win = rio_windows.from_bounds(
            x_min_proj, y_min_proj, x_max_proj, y_max_proj,
            transform=ds.transform,
        )

        full_win = rio_windows.Window(0, 0, ds.width, ds.height)
        try:
            win = win.intersection(full_win)
        except rio_windows.WindowError:
            empty = np.empty((0, 0), dtype=np.float64)
            return empty, ds.transform, src_crs, nodata_override or ds.nodata

        col_off = max(0, int(_math.floor(win.col_off)))
        row_off = max(0, int(_math.floor(win.row_off)))
        col_end = min(ds.width,  int(_math.ceil(win.col_off + win.width)))
        row_end = min(ds.height, int(_math.ceil(win.row_off + win.height)))

        win_int = rio_windows.Window(
            col_off=col_off,
            row_off=row_off,
            width=col_end - col_off,
            height=row_end - row_off,
        )

        data = ds.read(1, window=win_int).astype(np.float64)
        nodata = nodata_override if nodata_override is not None else ds.nodata
        transform = ds.window_transform(win_int)
        crs = src_crs

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


# ===========================================================================
# Landscape assembly
# ===========================================================================

def assemble_landscape(elev_data, elev_transform, elev_crs, elev_nodata,
                       slope_data, slope_transform, slope_crs,
                       aspect_data, aspect_transform, aspect_crs,
                       fuel_data, fuel_transform, fuel_crs,
                       project_utm=True,
                       subsample=1,
                       keep_nonburnable=False,
                       fuel_type="13"):
    """Merge four rasters into flat arrays suitable for writing a .lcp file.

    Parameters
    ----------
    fuel_type : str
        Fuel model system used in *fuel_data*: ``"13"`` for Anderson 13
        (FBFM13, codes 1–13) or ``"40"`` for Scott & Burgan 40 (FBFM40,
        codes 101–204).  Determines which pixel values are kept as burnable
        when *keep_nonburnable* is ``False``.
    """
    import numpy as np

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
        if fuel_type == "40":
            mask &= (
                ((fuel_int >= 101) & (fuel_int <= 109)) |
                ((fuel_int >= 121) & (fuel_int <= 124)) |
                ((fuel_int >= 141) & (fuel_int <= 149)) |
                ((fuel_int >= 161) & (fuel_int <= 165)) |
                ((fuel_int >= 181) & (fuel_int <= 189)) |
                ((fuel_int >= 201) & (fuel_int <= 204))
            )
        else:
            mask &= (fuel_int >= 1) & (fuel_int <= 13)

    if subsample > 1:
        skip = np.zeros(elev_data.shape, dtype=bool)
        skip[::subsample, ::subsample] = True
        mask &= skip

    return (xs[mask], ys[mask], elev_data[mask],
            slope_data[mask], aspect_data[mask], fuel_int[mask].astype(float))


def _write_lcp(path, xs, ys, elev, slope, aspect, fuel, fuel_type="13"):
    """Write an ASCII landscape file (.lcp)."""
    fuel_system = ("FBFM40 (Scott & Burgan 40)" if fuel_type == "40"
                   else "FBFM13 (Anderson 13, NFFL 1-13)")
    with open(path, "w") as fh:
        fh.write(
            "# Landscape file generated by landscape_writer.py\n"
            "# Format: X Y ELEVATION SLOPE ASPECT FUEL_MODEL\n"
            "# Units: X/Y metres, ELEVATION metres, SLOPE degrees,\n"
            f"#        ASPECT degrees (0=North), FUEL_MODEL {fuel_system}\n"
        )
        for x, y, z, s, a, f in zip(xs, ys, elev, slope, aspect, fuel):
            fh.write(f"{x:.2f} {y:.2f} {z:.2f} {s:.2f} {a:.2f} {int(f)}\n")
    print(f"Wrote {len(xs)} landscape points to '{path}'.")


# ===========================================================================
# High-level landscape creation
# ===========================================================================

def create_landscape(output_path, bbox, vintage=2020,
                     elev_product=None, slope_product=None,
                     aspect_product=None, fuel_product=None,
                     project_utm=True, subsample=1,
                     keep_nonburnable=False, cache_dir=None,
                     timeout_s=300, use_cog=True,
                     sources=None, vintage_fallback=True,
                     fuel_type="13"):
    """Download LANDFIRE rasters and write an ASCII landscape file.

    Sources are tried in the order given by *sources* (default
    ``['cog', 'mspc', 'lfps']``).  Each source is attempted in turn; if a
    source fails, the next one is tried automatically.

    When all sources fail for the requested vintage and *vintage_fallback* is
    ``True`` (the default), the function retries with the next older vintage
    in the sequence ``2020 → 2016 → 2014``.
    """
    if fuel_type not in ("13", "40"):
        raise ValueError(
            f"fuel_type must be '13' or '40', got: {fuel_type!r}"
        )
    print(f"Fuel model system: {'FBFM13 (Anderson 13)' if fuel_type == '13' else 'FBFM40 (Scott & Burgan 40)'}")
    _have_overrides = any((elev_product, slope_product,
                           aspect_product, fuel_product))

    if vintage_fallback and not _have_overrides and use_cog:
        if vintage in _VINTAGE_FALLBACK_ORDER:
            _fb_start = _VINTAGE_FALLBACK_ORDER.index(vintage)
            vintages_to_try = list(_VINTAGE_FALLBACK_ORDER[_fb_start:])
        else:
            vintages_to_try = [vintage]
    else:
        vintages_to_try = [vintage]

    if not use_cog:
        _fixed_sources = ["lfps"]
    elif _have_overrides:
        print(
            "WARNING: Custom --*-product IDs are not supported by the COG or "
            "MSPC downloaders; using LFPS API only.",
            file=sys.stderr,
        )
        _fixed_sources = ["lfps"]
    else:
        _fixed_sources = list(sources) if sources else ["cog", "mspc", "lfps"]

    layer_map = None
    last_exc = None

    for v_idx, v_attempt in enumerate(vintages_to_try):
        if v_idx > 0:
            print(
                f"WARNING: All LANDFIRE sources failed for vintage "
                f"{vintages_to_try[v_idx - 1]}; retrying with older "
                f"vintage {v_attempt} …",
                file=sys.stderr,
            )

        layer_defaults = _DEFAULT_LAYERS.get(v_attempt, _DEFAULT_LAYERS[2020])
        lid_elev   = elev_product   or layer_defaults["elev"]
        lid_slope  = slope_product  or layer_defaults["slope"]
        lid_aspect = aspect_product or layer_defaults["aspect"]
        fuel_key   = "fuel40" if fuel_type == "40" else "fuel13"
        lid_fuel   = fuel_product   or layer_defaults[fuel_key]
        layer_ids  = [lid_elev, lid_slope, lid_aspect, lid_fuel]

        effective_sources = list(_fixed_sources)
        if v_attempt not in _COG_LAYERS and "cog" in effective_sources:
            print(
                f"WARNING: No COG layer mapping for vintage {v_attempt}; "
                "skipping 'cog' source.",
                file=sys.stderr,
            )
            effective_sources = [s for s in effective_sources if s != "cog"]
        if v_attempt not in _MSPC_ASSET_KEYS and "mspc" in effective_sources:
            print(
                f"WARNING: No MSPC asset key mapping for vintage {v_attempt}; "
                "skipping 'mspc' source.",
                file=sys.stderr,
            )
            effective_sources = [s for s in effective_sources if s != "mspc"]

        if not effective_sources:
            raise RuntimeError(
                "No LANDFIRE download sources remain after filtering.  "
                "Pass --sources lfps to use the USGS LFPS API."
            )

        vintage_layer_map = None

        for source in effective_sources:
            try:
                if source == "cog":
                    print("Fetching LANDFIRE layers from COG (AWS S3 / HTTPS) …")
                    vintage_layer_map = download_landfire_cog(
                        bbox, vintage=v_attempt,
                        lid_elev=lid_elev, lid_slope=lid_slope,
                        lid_aspect=lid_aspect, lid_fuel=lid_fuel,
                        fuel_type=fuel_type,
                    )
                elif source == "mspc":
                    print(
                        "Fetching LANDFIRE layers from Microsoft Planetary "
                        "Computer (MSPC) …"
                    )
                    vintage_layer_map = download_landfire_mspc(
                        bbox, vintage=v_attempt,
                        lid_elev=lid_elev, lid_slope=lid_slope,
                        lid_aspect=lid_aspect, lid_fuel=lid_fuel,
                        fuel_type=fuel_type,
                    )
                elif source == "lfps":
                    print("Fetching LANDFIRE layers from USGS LFPS API …")
                    vintage_layer_map = download_landfire(
                        bbox, layer_ids,
                        cache_dir=cache_dir,
                        timeout_s=timeout_s,
                    )
                else:
                    print(
                        f"WARNING: Unknown LANDFIRE source '{source}'; "
                        "skipping.",
                        file=sys.stderr,
                    )
                    continue

                break

            except Exception as exc:
                remaining = effective_sources[
                    effective_sources.index(source) + 1:]
                if remaining:
                    print(
                        f"WARNING: LANDFIRE source '{source}' failed "
                        f"({type(exc).__name__}: {exc}); "
                        f"trying next source: {remaining[0]} …",
                        file=sys.stderr,
                    )
                else:
                    print(
                        f"WARNING: LANDFIRE source '{source}' failed "
                        f"({type(exc).__name__}: {exc}); "
                        "no more sources to try.",
                        file=sys.stderr,
                    )
                last_exc = exc
                vintage_layer_map = None

        if vintage_layer_map is not None:
            layer_map = vintage_layer_map
            break

    if layer_map is None:
        vintages_tried = ", ".join(str(v) for v in vintages_to_try)
        raise RuntimeError(
            f"All LANDFIRE download sources failed for vintage(s) "
            f"{vintages_tried}. "
            "Check your internet connection or try a different --sources "
            "combination.  Last error: " + str(last_exc)
        ) from last_exc

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
    asp_data, asp_tf, asp_crs, _ = _read_raster_bytes(_get_bytes(lid_aspect))
    print("Reading fuel model raster …")
    fuel_data, fuel_tf, fuel_crs, _ = _read_raster_bytes(_get_bytes(lid_fuel))

    print(f"Elevation grid: {elev_data.shape[0]}×{elev_data.shape[1]}")

    xs, ys, elev, slope, aspect, fuel = assemble_landscape(
        elev_data, elev_tf, elev_crs, elev_nd,
        slope_data, slope_tf, slope_crs,
        asp_data, asp_tf, asp_crs,
        fuel_data, fuel_tf, fuel_crs,
        project_utm=project_utm, subsample=subsample,
        keep_nonburnable=keep_nonburnable,
        fuel_type=fuel_type,
    )
    _write_lcp(output_path, xs, ys, elev, slope, aspect, fuel, fuel_type=fuel_type)


def create_landscape_from_files(output_path, elev_path, slope_path,
                                 aspect_path, fuel_path,
                                 project_utm=True, subsample=1,
                                 keep_nonburnable=False, elev_nodata=None,
                                 fuel_type="13"):
    """Create a landscape file from local raster files (no download needed)."""
    print("Reading elevation raster …")
    elev_data, elev_tf, elev_crs, elev_nd = _read_raster_file(
        elev_path, nodata_override=elev_nodata)
    print("Reading slope raster …")
    slope_data, slope_tf, slope_crs, _ = _read_raster_file(slope_path)
    print("Reading aspect raster …")
    asp_data, asp_tf, asp_crs, _ = _read_raster_file(aspect_path)
    print("Reading fuel model raster …")
    fuel_data, fuel_tf, fuel_crs, _ = _read_raster_file(fuel_path)

    print(f"Elevation grid: {elev_data.shape[0]}×{elev_data.shape[1]}")

    xs, ys, elev, slope, aspect, fuel = assemble_landscape(
        elev_data, elev_tf, elev_crs, elev_nd,
        slope_data, slope_tf, slope_crs,
        asp_data, asp_tf, asp_crs,
        fuel_data, fuel_tf, fuel_crs,
        project_utm=project_utm, subsample=subsample,
        keep_nonburnable=keep_nonburnable,
        fuel_type=fuel_type,
    )
    _write_lcp(output_path, xs, ys, elev, slope, aspect, fuel, fuel_type=fuel_type)


def create_landscape_srtm_with_fuel(output_path, lat_min, lat_max,
                                     lon_min, lon_max, fuel_path,
                                     tif_path=None, project_utm=True,
                                     subsample=1, keep_nonburnable=False,
                                     fuel_type="13"):
    """Create a landscape file using SRTM-derived terrain and a local fuel raster.

    Downloads SRTM elevation data for the given bounding box, derives slope
    and aspect via finite differences, then interpolates the SRTM data onto
    the fuel raster's grid.  This is the backend for ``--srtm-slope-aspect``.
    """
    from srtm_terrain_reader import (
        download_srtm, _compute_slope_aspect_from_srtm_tif,
    )

    _tmp_tif = tif_path is None
    if _tmp_tif:
        tmp = tempfile.NamedTemporaryFile(suffix=".tif", delete=False)
        tif_path = tmp.name
        tmp.close()

    try:
        print(
            f"Downloading SRTM data for bbox "
            f"({lat_min},{lon_min}) – ({lat_max},{lon_max}) …"
        )
        download_srtm(lat_min, lat_max, lon_min, lon_max, tif_path)

        print("Computing slope and aspect from SRTM elevation …")
        elev_data, slope_data, aspect_data, elev_tf, elev_crs, elev_nd = \
            _compute_slope_aspect_from_srtm_tif(tif_path)
    finally:
        if _tmp_tif and os.path.isfile(tif_path):
            os.remove(tif_path)

    print(f"SRTM elevation grid: {elev_data.shape[0]}×{elev_data.shape[1]}")
    print("Reading fuel model raster …")
    fuel_data, fuel_tf, fuel_crs, _ = _read_raster_file_clipped(
        fuel_path, lat_min, lat_max, lon_min, lon_max)

    print(f"Fuel grid (clipped to bbox): {fuel_data.shape[0]}×{fuel_data.shape[1]}")
    print("Interpolating SRTM elevation, slope, and aspect to fuel grid …")
    elev_data = _resample_to_grid(
        elev_data, elev_tf, elev_crs,
        fuel_data, fuel_tf, fuel_crs, resampling="bilinear")
    slope_data = _resample_to_grid(
        slope_data, elev_tf, elev_crs,
        fuel_data, fuel_tf, fuel_crs, resampling="bilinear")
    aspect_data = _resample_to_grid(
        aspect_data, elev_tf, elev_crs,
        fuel_data, fuel_tf, fuel_crs, resampling="bilinear")

    xs, ys, elev, slope, aspect, fuel = assemble_landscape(
        elev_data, fuel_tf, fuel_crs, elev_nd,
        slope_data, fuel_tf, fuel_crs,
        aspect_data, fuel_tf, fuel_crs,
        fuel_data, fuel_tf, fuel_crs,
        project_utm=project_utm, subsample=subsample,
        keep_nonburnable=keep_nonburnable,
        fuel_type=fuel_type,
    )
    _write_lcp(output_path, xs, ys, elev, slope, aspect, fuel, fuel_type=fuel_type)


# ===========================================================================
# CLI
# ===========================================================================

def _build_parser():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument("--lat-min", type=float, required=True,
                        help="Minimum latitude (WGS-84 decimal degrees)")
    parser.add_argument("--lat-max", type=float, required=True,
                        help="Maximum latitude (WGS-84 decimal degrees)")
    parser.add_argument("--lon-min", type=float, required=True,
                        help="Minimum longitude (WGS-84 decimal degrees)")
    parser.add_argument("--lon-max", type=float, required=True,
                        help="Maximum longitude (WGS-84 decimal degrees)")

    parser.add_argument("--landscape", default="landscape.lcp",
                        help="Output landscape LCP file (default: landscape.lcp)")
    parser.add_argument("--subsample", type=int, default=1, metavar="N",
                        help="Keep every N-th point in each dimension (default: 1)")

    # LANDFIRE options
    parser.add_argument("--vintage", type=int, default=2020, metavar="YEAR",
                        help="LANDFIRE data vintage year (default: 2020)")
    parser.add_argument(
        "--fuel-model-type", default="13", choices=["13", "40"],
        metavar="{13,40}",
        help=(
            "Fuel model system: '13' for Anderson 13 (FBFM13, default) or "
            "'40' for Scott & Burgan 40 (FBFM40)."
        ),
    )
    parser.add_argument("--keep-nonburnable", action="store_true",
                        help="Include non-burnable LANDFIRE pixels in LCP output")
    parser.add_argument("--elev-product",   default=None, metavar="ID",
                        help="Override LANDFIRE elevation product ID")
    parser.add_argument("--slope-product",  default=None, metavar="ID",
                        help="Override LANDFIRE slope product ID")
    parser.add_argument("--aspect-product", default=None, metavar="ID",
                        help="Override LANDFIRE aspect product ID")
    parser.add_argument("--fuel-product",   default=None, metavar="ID",
                        help="Override LANDFIRE fuel model product ID")
    parser.add_argument("--cache-dir", default=None, metavar="DIR",
                        help="Directory to cache downloaded LANDFIRE ZIP files "
                             "(only used with --use-lfps)")
    parser.add_argument("--timeout", type=int, default=300, metavar="N",
                        help="LANDFIRE API polling timeout in seconds (default: 300)")
    parser.add_argument(
        "--use-lfps", action="store_true",
        help="Force the legacy LFPS REST API; equivalent to --sources lfps.",
    )
    parser.add_argument(
        "--sources", default=None, metavar="LIST",
        help=(
            "Comma-separated priority list of LANDFIRE sources to try in "
            "order.  Valid tokens: 'cog', 'mspc', 'lfps'.  "
            "Default: 'cog,mspc,lfps'."
        ),
    )
    parser.add_argument(
        "--no-vintage-fallback", action="store_true",
        help="Disable automatic vintage fallback (2020 → 2016 → 2014).",
    )

    # Local-file mode
    local = parser.add_argument_group(
        "local files mode (use pre-existing rasters instead of downloading)"
    )
    local.add_argument("--elev-file",   default=None, metavar="PATH",
                       help="Local elevation raster (skips download)")
    local.add_argument("--slope-file",  default=None, metavar="PATH",
                       help="Local slope raster")
    local.add_argument("--aspect-file", default=None, metavar="PATH",
                       help="Local aspect raster")
    local.add_argument("--fuel-file",   default=None, metavar="PATH",
                       help="Local fuel model raster")
    local.add_argument(
        "--srtm-slope-aspect", action="store_true",
        help=(
            "When --fuel-file is given without --elev-file/--slope-file/"
            "--aspect-file, derive elevation, slope, and aspect from SRTM "
            "data.  Requires the 'elevation' package (pip install elevation)."
        ),
    )
    local.add_argument("--tif", default=None, metavar="FILE",
                       help="Intermediate SRTM GeoTIFF path (for --srtm-slope-aspect)")

    return parser


def main(argv=None):
    parser = _build_parser()
    args = parser.parse_args(argv)

    if args.lat_min >= args.lat_max or args.lon_min >= args.lon_max:
        print(
            "ERROR: bounding box must satisfy "
            "lon-min < lon-max and lat-min < lat-max",
            file=sys.stderr,
        )
        sys.exit(1)

    if args.subsample < 1:
        print("ERROR: --subsample must be >= 1", file=sys.stderr)
        sys.exit(1)

    if args.srtm_slope_aspect and args.fuel_file is None:
        print(
            "ERROR: --srtm-slope-aspect requires --fuel-file to be specified.",
            file=sys.stderr,
        )
        sys.exit(1)

    # Resolve LANDFIRE source priority list
    if args.use_lfps:
        landfire_sources = ["lfps"]
    elif args.sources is not None:
        raw_sources = [s.strip().lower() for s in args.sources.split(",")
                       if s.strip()]
        invalid = [s for s in raw_sources if s not in _LANDFIRE_SOURCES]
        if invalid:
            print(
                f"ERROR: Unknown --sources token(s): {invalid}. "
                f"Valid tokens are: {list(_LANDFIRE_SOURCES)}",
                file=sys.stderr,
            )
            sys.exit(1)
        landfire_sources = raw_sources
    else:
        landfire_sources = list(_LANDFIRE_SOURCES)

    bbox = (args.lon_min, args.lat_min, args.lon_max, args.lat_max)
    out_lcp = os.path.abspath(args.landscape)

    # Determine which mode to use
    using_srtm_terrain = (
        args.srtm_slope_aspect
        and args.fuel_file is not None
        and args.elev_file is None
        and args.slope_file is None
        and args.aspect_file is None
    )
    local_files = [args.elev_file, args.slope_file, args.aspect_file, args.fuel_file]
    using_local = any(f is not None for f in local_files) and not using_srtm_terrain

    if using_srtm_terrain:
        tif_path = os.path.abspath(args.tif) if args.tif is not None else None
        create_landscape_srtm_with_fuel(
            out_lcp,
            lat_min=args.lat_min, lat_max=args.lat_max,
            lon_min=args.lon_min, lon_max=args.lon_max,
            fuel_path=args.fuel_file,
            tif_path=tif_path,
            subsample=args.subsample,
            keep_nonburnable=args.keep_nonburnable,
            fuel_type=args.fuel_model_type,
        )
    elif using_local:
        missing = [
            name for name, val in [
                ("--elev-file",   args.elev_file),
                ("--slope-file",  args.slope_file),
                ("--aspect-file", args.aspect_file),
                ("--fuel-file",   args.fuel_file),
            ] if val is None
        ]
        if missing:
            print(
                f"ERROR: local-files mode requires all four rasters. "
                f"Missing: {', '.join(missing)}",
                file=sys.stderr,
            )
            sys.exit(1)
        create_landscape_from_files(
            out_lcp,
            args.elev_file, args.slope_file,
            args.aspect_file, args.fuel_file,
            subsample=args.subsample,
            keep_nonburnable=args.keep_nonburnable,
            fuel_type=args.fuel_model_type,
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
            use_cog=not args.use_lfps,
            sources=landfire_sources,
            vintage_fallback=not args.no_vintage_fallback,
            fuel_type=args.fuel_model_type,
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
                f"\nUTM extents of landscape data "
                f"(use for prob_lo / prob_hi):\n"
                f"  prob_lo_x = {min(xs):.0f}\n"
                f"  prob_lo_y = {min(ys):.0f}\n"
                f"  prob_hi_x = {max(xs):.0f}\n"
                f"  prob_hi_y = {max(ys):.0f}\n"
            )


if __name__ == "__main__":
    main()
