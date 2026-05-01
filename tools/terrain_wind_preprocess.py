#!/usr/bin/env python3
"""
terrain_wind_preprocess.py - Unified terrain, landscape, and wind preprocessing tool.

Merges srtm_landfire_to_terrain.py and wrf_to_terrain_wind.py into a single
script that can download SRTM elevation data, LANDFIRE fuel/slope/aspect
rasters, and extract wind data from a WRF-style netCDF file.

Key behaviours
--------------
* If ``--wrf-file`` is supplied, the bounding box (lat/lon min/max) is derived
  from the WRF netCDF automatically: the domain centre is computed and the
  bounding box is clipped to **±0.45 degrees** in each direction (total span
  0.9° × 0.9°, within the SRTM 1-degree-tile download limit).  Any
  ``--lat-min/max/lon-min/max`` values on the command line are **ignored**.
* Unless ``--no-terrain`` is given, SRTM elevation is used for the terrain
  and landscape files.  When ``--no-terrain`` *and* ``--wrf-file`` are both
  given, the ``HGT`` (or ``HGT_M``) variable from the WRF file is written as
  the terrain output instead of downloading SRTM data.
* If ``--interpolate-wind`` is given (requires ``--wrf-file`` and SRTM terrain
  to be active), the WRF wind fields are interpolated to the SRTM terrain grid
  and the resulting wind file shares the same (x, y) points as the terrain
  file.
* ``--time-range T1:TN`` extracts WRF time steps T1 through TN **inclusive**
  (0-based indices).  When multiple time steps are requested the wind output
  filenames are derived from ``--wind`` by inserting ``_tN`` before the
  extension, e.g. ``wind_t0.csv``, ``wind_t1.csv``, …
* ``--time-index N`` (single time step, default 0) is the legacy flag and
  takes precedence when ``--time-range`` is not supplied.

Usage examples
--------------
  # SRTM terrain + LANDFIRE landscape for a bounding box
  python3 terrain_wind_preprocess.py \\
      --lat-min 40 --lat-max 40.5 \\
      --lon-min -105 --lon-max -104.5

  # Same but also extract wind from a WRF file (bbox comes from WRF)
  python3 terrain_wind_preprocess.py \\
      --wrf-file wrfout_d01 \\
      --wind wind.csv

  # Skip SRTM (use WRF HGT_M for terrain) and extract wind only
  python3 terrain_wind_preprocess.py \\
      --wrf-file wrfout_d01 \\
      --no-terrain \\
      --wind wind.csv

  # Interpolate WRF wind to SRTM resolution
  python3 terrain_wind_preprocess.py \\
      --wrf-file wrfout_d01 \\
      --wind wind.csv \\
      --interpolate-wind

  # Extract a range of time steps
  python3 terrain_wind_preprocess.py \\
      --wrf-file wrfout_d01 \\
      --wind wind.csv \\
      --time-range 0:5

Options
-------
  --wrf-file FILE       WRF netCDF file; if given the lat/lon bbox is read
                        from the file automatically.
  --lat-min / --lat-max Latitude  bounds (WGS-84 degrees); required when
                        --wrf-file is not given.
  --lon-min / --lon-max Longitude bounds (WGS-84 degrees); required when
                        --wrf-file is not given.
  --terrain FILE        Output terrain XYZ file   (default: terrain.xyz)
  --landscape FILE      Output landscape LCP file  (default: landscape.lcp)
  --wind FILE           Output wind CSV file       (default: wind.csv)
                        Requires --wrf-file.
  --tif FILE            Intermediate SRTM GeoTIFF (temp file if not given).
  --subsample N         Keep every N-th point (default: 1).
  --no-terrain          Skip SRTM terrain and landscape steps.
  --no-landscape        Skip the LANDFIRE landscape step.
  --no-wind             Skip the wind extraction step.
  --interpolate-wind    Interpolate WRF wind to SRTM terrain grid.
  --time-range T1:TN    Time index range (inclusive) to extract from WRF.
  --time-index N        Single WRF time index (default: 0).
  --level N             WRF vertical level for wind (default: 0 = lowest).
  --keep-nonburnable    Include non-burnable LANDFIRE pixels.
  --vintage YEAR        LANDFIRE vintage year (default: 2020).
  --fuel-product ID     Override LANDFIRE fuel model product ID.
  --elev-product ID     Override LANDFIRE elevation product ID.
  --slope-product ID    Override LANDFIRE slope product ID.
  --aspect-product ID   Override LANDFIRE aspect product ID.
  --cache-dir DIR       Cache directory for LANDFIRE ZIP files.
  --timeout N           LANDFIRE API polling timeout in seconds (default: 300).
  --elev-file PATH      Local elevation raster (skips download).
  --slope-file PATH     Local slope raster.
  --aspect-file PATH    Local aspect raster.
  --fuel-file PATH      Local fuel model raster.
  --use-lfps            Force the legacy LFPS REST API for LANDFIRE download
                        instead of the default COG (S3/HTTPS) approach.  Use
                        this only if the S3 COG endpoint is unreachable.
  --help                Show this message and exit.

LANDFIRE data source
--------------------
By default, LANDFIRE elevation, slope, aspect and fuel-model rasters are
fetched as Cloud Optimized GeoTIFFs (COGs) from the public LANDFIRE AWS S3
bucket (``https://s3.amazonaws.com/landfire/``).  Only the pixels within the
requested bounding box are downloaded, so the transfer is fast and avoids the
LFPS REST API which can fail with SSL certificate errors.  Pass ``--use-lfps``
to revert to the legacy API if needed.
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

# ---------------------------------------------------------------------------
# LANDFIRE Cloud Optimized GeoTIFF (COG) constants
# ---------------------------------------------------------------------------

#: Base HTTPS URL for the public LANDFIRE AWS S3 bucket.
_LANDFIRE_COG_BASE = "https://s3.amazonaws.com/landfire/"

#: Mapping from vintage year → layer key → relative COG path on S3.
#: These correspond to LANDFIRE national products hosted as COGs; only the
#: pixels inside the requested bounding box are transferred.
_COG_LAYERS = {
    2020: {
        "elev":   "US_200/US_200ELEV.tif",
        "slope":  "US_200/US_200SlpD.tif",
        "aspect": "US_200/US_200Asp.tif",
        "fuel13": "US_200/US_200FBFM13.tif",
    },
    2016: {
        "elev":   "US_140/US_140ELEV.tif",
        "slope":  "US_140/US_140SlpD.tif",
        "aspect": "US_140/US_140Asp.tif",
        "fuel13": "US_140/US_140FBFM13.tif",
    },
    2014: {
        "elev":   "US_130/US_130ELEV.tif",
        "slope":  "US_130/US_130SLP.tif",
        "aspect": "US_130/US_130ASP.tif",
        "fuel13": "US_130/US_130FBFM13.tif",
    },
}

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

# LANDFIRE FBFM13 non-burnable codes
_NONBURNABLE_CODES = {91, 92, 93, 98, 99}


# ===========================================================================
# Shared projection helpers
# ===========================================================================

def _latlon_to_utm(lats, lons):
    """Project (lat, lon) arrays (WGS-84) to UTM metres.

    The UTM zone is determined from the centre of the domain.
    Returns (x_utm, y_utm) both in metres.
    """
    try:
        from pyproj import Transformer
    except ImportError:
        raise ImportError(
            "pyproj is required for UTM projection. "
            "Install with: pip install pyproj"
        )
    import numpy as np

    lats = np.asarray(lats, dtype=np.float64)
    lons = np.asarray(lons, dtype=np.float64)

    center_lat = float(np.mean(lats))
    center_lon = float(np.mean(lons))
    zone = int((center_lon + 180.0) / 6.0) + 1
    hemisphere = "north" if center_lat >= 0.0 else "south"
    epsg = 32600 + zone if hemisphere == "north" else 32700 + zone

    transformer = Transformer.from_crs("EPSG:4326", f"EPSG:{epsg}", always_xy=True)
    x_utm, y_utm = transformer.transform(lons.ravel(), lats.ravel())
    return (
        np.array(x_utm, dtype=np.float64).reshape(lats.shape),
        np.array(y_utm, dtype=np.float64).reshape(lats.shape),
    )


# ===========================================================================
# SRTM terrain helpers
# ===========================================================================

def download_srtm(lat_min, lat_max, lon_min, lon_max, out_tif):
    """Download SRTM1 elevation data clipped to the given bounding box."""
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


def _fix_srtm_zeros(z):
    """Replace spurious zero-elevation pixels with interpolated values.

    SRTM tiles sometimes contain random zero values at tile borders that are
    not genuine sea-level points.  A zero is treated as spurious when it is
    surrounded by cells whose elevation is > 0.  Such cells are filled by
    linear interpolation from the nearest valid (non-zero) neighbours using
    ``scipy.interpolate.griddata``.

    Parameters
    ----------
    z : numpy.ndarray (2-D, float)
        Elevation grid, potentially containing spurious zeros.

    Returns
    -------
    numpy.ndarray
        Copy of *z* with spurious zeros replaced.
    """
    import numpy as np

    zero_mask = z == 0.0
    if not np.any(zero_mask):
        return z

    # A zero is spurious only when at least one immediate neighbour is > 0
    from scipy.ndimage import generic_filter

    spurious = generic_filter(
        z,
        lambda p: float((p[len(p) // 2] == 0.0) and np.any(p > 0.0)),
        size=3,
        mode="nearest",
    ).astype(bool)

    if not np.any(spurious):
        return z

    print(f"  Fixing {int(np.sum(spurious))} spurious zero elevation pixels …")
    z_out = z.copy()
    h, w = z.shape
    rows, cols = np.mgrid[0:h, 0:w]

    valid = ~spurious
    try:
        from scipy.interpolate import griddata

        pts_valid = np.column_stack([rows[valid].ravel(), cols[valid].ravel()])
        pts_spurious = np.column_stack(
            [rows[spurious].ravel(), cols[spurious].ravel()]
        )
        z_out[spurious] = griddata(
            pts_valid, z[valid].ravel(), pts_spurious, method="linear", fill_value=0.0
        )
    except ImportError:
        # Fallback: replace each spurious cell with the mean of its valid neighbours
        for r, c in zip(rows[spurious].ravel(), cols[spurious].ravel()):
            neighbours = []
            for dr in (-1, 0, 1):
                for dc in (-1, 0, 1):
                    if dr == 0 and dc == 0:
                        continue
                    rr, cc = r + dr, c + dc
                    if 0 <= rr < h and 0 <= cc < w and z[rr, cc] > 0.0:
                        neighbours.append(z[rr, cc])
            if neighbours:
                z_out[r, c] = float(np.mean(neighbours))

    return z_out


def _smooth_terrain_border(z, fraction=0.2):
    """Smooth the outer *fraction* of the terrain grid to reduce tile-seam artefacts.

    A Gaussian-smoothed version of the full grid is blended into the border
    region.  The blend weight transitions linearly from 1 (fully smoothed) at
    the grid edges to 0 (original values) at the interior boundary of the
    border band, so there is no hard discontinuity.

    Parameters
    ----------
    z : numpy.ndarray (2-D, float)
        Elevation grid.
    fraction : float
        Width of the border band expressed as a fraction of the grid dimension
        (default 0.2 → outer 20 % on each side).

    Returns
    -------
    numpy.ndarray
        Smoothed copy of *z*.
    """
    import numpy as np

    try:
        from scipy.ndimage import gaussian_filter
    except ImportError:
        print(
            "WARNING: scipy not available; skipping terrain border smoothing.",
            file=sys.stderr,
        )
        return z

    h, w = z.shape
    border_h = max(1, int(h * fraction))
    border_w = max(1, int(w * fraction))

    # Gaussian sigma proportional to border width for effective smoothing
    sigma = max(1.0, min(border_h, border_w) / 4.0)
    z_smooth = gaussian_filter(z, sigma=sigma)

    # Build a 2-D blending weight: 1.0 at edges, 0.0 at interior boundary
    weight = np.zeros((h, w), dtype=np.float64)

    # Vertical ramp (top and bottom)
    ramp_v = np.linspace(1.0, 0.0, border_h, endpoint=False)
    weight[:border_h, :] = np.maximum(weight[:border_h, :],
                                       ramp_v[:, np.newaxis])
    weight[-border_h:, :] = np.maximum(weight[-border_h:, :],
                                        ramp_v[::-1, np.newaxis])

    # Horizontal ramp (left and right)
    ramp_h = np.linspace(1.0, 0.0, border_w, endpoint=False)
    weight[:, :border_w] = np.maximum(weight[:, :border_w],
                                       ramp_h[np.newaxis, :])
    weight[:, -border_w:] = np.maximum(weight[:, -border_w:],
                                        ramp_h[np.newaxis, ::-1])

    return z * (1.0 - weight) + z_smooth * weight


def tiff_to_xyz_utm(tif_path, subsample=1):
    """Read a GeoTIFF and return (utm_x, utm_y, z) 2-D arrays in UTM metres."""
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

    # Fix spurious zero-elevation pixels caused by SRTM tile-seam artefacts
    z = _fix_srtm_zeros(z)

    # Smooth the outer 20 % border to reduce tile-seam elevation discontinuities
    z = _smooth_terrain_border(z, fraction=0.2)

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
    """Write a terrain XYZ file compatible with ``rothermel.terrain_file``."""
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
    """Download SRTM data and write a UTM terrain XYZ file."""
    _tmp_tif = tif_path is None
    if _tmp_tif:
        tmp = tempfile.NamedTemporaryFile(suffix=".tif", delete=False)
        tif_path = tmp.name
        tmp.close()

    try:
        print(f"Downloading SRTM data for bbox "
              f"({lat_min},{lon_min}) – ({lat_max},{lon_max}) …")
        download_srtm(lat_min, lat_max, lon_min, lon_max, tif_path)

        print("Converting SRTM GeoTIFF → UTM XYZ …")
        utm_x, utm_y, z = tiff_to_xyz_utm(tif_path, subsample=subsample)

        write_terrain_xyz(utm_x, utm_y, z, output_path)
    finally:
        if _tmp_tif and os.path.isfile(tif_path):
            os.remove(tif_path)


def create_terrain_xyz_return_grid(lat_min, lat_max, lon_min, lon_max,
                                   output_path, tif_path=None, subsample=1):
    """Download SRTM, write terrain XYZ, and return (utm_x, utm_y, z) arrays.

    Like ``create_terrain_xyz`` but also returns the 2-D grid arrays so they
    can be used for wind interpolation.
    """
    _tmp_tif = tif_path is None
    if _tmp_tif:
        tmp = tempfile.NamedTemporaryFile(suffix=".tif", delete=False)
        tif_path = tmp.name
        tmp.close()

    try:
        print(f"Downloading SRTM data for bbox "
              f"({lat_min},{lon_min}) – ({lat_max},{lon_max}) …")
        download_srtm(lat_min, lat_max, lon_min, lon_max, tif_path)

        print("Converting SRTM GeoTIFF → UTM XYZ …")
        utm_x, utm_y, z = tiff_to_xyz_utm(tif_path, subsample=subsample)

        write_terrain_xyz(utm_x, utm_y, z, output_path)
    finally:
        if _tmp_tif and os.path.isfile(tif_path):
            os.remove(tif_path)

    return utm_x, utm_y, z


# ===========================================================================
# LANDFIRE landscape helpers
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
            raise RuntimeError(
                f"LANDFIRE job {job_id} ended with status: {status}")
        poll_interval = min(poll_interval * 1.5, 30)
    else:
        raise RuntimeError(
            f"LANDFIRE job {job_id} did not complete within {timeout_s}s")

    result_url = _RESULT_URL.format(jobId=job_id)
    rr = _requests_get(result_url, params={"f": "json"}, timeout=30)
    rr.raise_for_status()
    result = rr.json()

    download_url = (
        result.get("value", {}).get("url")
        or result.get("value")
    )
    if not download_url:
        raise RuntimeError(
            f"Could not parse download URL from LANDFIRE result: {result}")
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
    """Read a LANDFIRE COG clipped to *bbox* (min_lon, min_lat, max_lon, max_lat).

    Opens the remote Cloud Optimized GeoTIFF via HTTPS and downloads only the
    window of pixels that falls within *bbox*.  No API polling or ZIP
    extraction is required.

    Parameters
    ----------
    cog_url : str
        Full HTTPS URL of the COG file.
    bbox : tuple of float
        ``(min_lon, min_lat, max_lon, max_lat)`` in WGS-84 degrees.

    Returns
    -------
    tuple
        ``(data, transform, crs, nodata)`` where *data* is a 2-D float64
        array, *transform* is a rasterio ``Affine`` for the windowed region,
        *crs* is the dataset CRS, and *nodata* is the nodata sentinel value.
    """
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

            # Transform the WGS-84 bbox to the dataset's native CRS for windowing.
            # LANDFIRE COGs use Albers Equal Area (EPSG:5070), not WGS-84.
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
                          lid_aspect=None, lid_fuel=None):
    """Download LANDFIRE rasters from Cloud Optimized GeoTIFFs on AWS S3.

    Reads the four LANDFIRE layers (elevation, slope, aspect, fuel model)
    from the public LANDFIRE AWS S3 bucket via HTTPS windowed reads.  Only
    pixels within *bbox* are transferred, so no job submission or ZIP
    extraction is needed and the download avoids the LFPS SSL endpoint.

    Parameters
    ----------
    bbox : tuple of float
        ``(min_lon, min_lat, max_lon, max_lat)`` in WGS-84 degrees.
    vintage : int
        LANDFIRE data vintage year (2014, 2016, or 2020; default 2020).
        If the vintage is not in ``_COG_LAYERS``, 2020 is used as a fallback.
    lid_elev, lid_slope, lid_aspect, lid_fuel : str or None
        Layer IDs to use as keys in the returned dict.  Defaults to the
        standard LFPS product IDs for *vintage* so that the return value is
        a drop-in replacement for :func:`download_landfire`.

    Returns
    -------
    dict
        ``{layer_id: bytes}`` where each value is an in-memory GeoTIFF,
        compatible with :func:`_read_raster_bytes`.

    Raises
    ------
    ValueError
        If no COG layer mapping exists for the requested vintage.
    RuntimeError
        If the bbox-to-CRS transformation or the rasterio read fails.
    """
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

    layers = {
        lid_elev   or layer_defaults["elev"]:   cog_paths["elev"],
        lid_slope  or layer_defaults["slope"]:  cog_paths["slope"],
        lid_aspect or layer_defaults["aspect"]: cog_paths["aspect"],
        lid_fuel   or layer_defaults["fuel13"]: cog_paths["fuel13"],
    }

    result = {}
    for lid, rel_path in layers.items():
        cog_url = _LANDFIRE_COG_BASE + rel_path
        data, transform, crs, nodata = _read_landfire_cog(cog_url, bbox)

        # Encode as an in-memory GeoTIFF so _read_raster_bytes() can consume
        # the result without modification.
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
    """Merge four rasters into flat arrays suitable for writing a .lcp file."""
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
            "# Landscape file generated by terrain_wind_preprocess.py\n"
            "# Format: X Y ELEVATION SLOPE ASPECT FUEL_MODEL\n"
            "# Units: X/Y metres, ELEVATION metres, SLOPE degrees,\n"
            "#        ASPECT degrees (0=North), FUEL_MODEL NFFL 1-13\n"
        )
        for x, y, z, s, a, f in zip(xs, ys, elev, slope, aspect, fuel):
            fh.write(f"{x:.2f} {y:.2f} {z:.2f} {s:.2f} {a:.2f} {int(f)}\n")
    print(f"Wrote {len(xs)} landscape points to '{path}'.")


def create_landscape(output_path, bbox, vintage=2020,
                     elev_product=None, slope_product=None,
                     aspect_product=None, fuel_product=None,
                     project_utm=True, subsample=1,
                     keep_nonburnable=False, cache_dir=None,
                     timeout_s=300, use_cog=True):
    """Download LANDFIRE rasters and write an ASCII landscape file.

    By default (*use_cog=True*) layers are fetched from the public LANDFIRE
    AWS S3 COG endpoint, which avoids the LFPS REST API and its occasional
    SSL certificate verification failures.  Set *use_cog=False* (or pass
    ``--use-lfps`` on the command line) to use the legacy LFPS API instead.

    If any ``*_product`` override is provided, those custom layer IDs cannot
    be resolved to a COG path; the function automatically falls back to the
    LFPS API in that case.
    """
    layer_defaults = _DEFAULT_LAYERS.get(vintage, _DEFAULT_LAYERS[2020])

    lid_elev   = elev_product   or layer_defaults["elev"]
    lid_slope  = slope_product  or layer_defaults["slope"]
    lid_aspect = aspect_product or layer_defaults["aspect"]
    lid_fuel   = fuel_product   or layer_defaults["fuel13"]

    # Custom product IDs cannot be mapped to a COG path; fall back to LFPS.
    _have_overrides = any([elev_product, slope_product,
                           aspect_product, fuel_product])
    _use_cog = use_cog and not _have_overrides and vintage in _COG_LAYERS
    if use_cog and _have_overrides:
        print(
            "WARNING: Custom --*-product IDs are not supported by the COG "
            "downloader; falling back to LFPS API.",
            file=sys.stderr,
        )
    if use_cog and vintage not in _COG_LAYERS:
        print(
            f"WARNING: No COG layer mapping for vintage {vintage}; "
            "falling back to LFPS API.",
            file=sys.stderr,
        )

    layer_ids = [lid_elev, lid_slope, lid_aspect, lid_fuel]

    if _use_cog:
        print("Fetching LANDFIRE layers from COG (HTTPS/S3) …")
        try:
            layer_map = download_landfire_cog(
                bbox, vintage=vintage,
                lid_elev=lid_elev, lid_slope=lid_slope,
                lid_aspect=lid_aspect, lid_fuel=lid_fuel,
            )
        except Exception as exc:
            print(
                f"WARNING: COG download failed ({exc}); "
                "retrying via LFPS API …",
                file=sys.stderr,
            )
            layer_map = download_landfire(bbox, layer_ids,
                                          cache_dir=cache_dir,
                                          timeout_s=timeout_s)
    else:
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
    )
    _write_lcp(output_path, xs, ys, elev, slope, aspect, fuel)


def create_landscape_from_files(output_path, elev_path, slope_path,
                                 aspect_path, fuel_path,
                                 project_utm=True, subsample=1,
                                 keep_nonburnable=False, elev_nodata=None):
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
    )
    _write_lcp(output_path, xs, ys, elev, slope, aspect, fuel)


# ===========================================================================
# WRF helpers
# ===========================================================================

def _open_nc(path):
    """Open a netCDF4 file and return the dataset."""
    try:
        import netCDF4 as nc
    except ImportError:
        raise ImportError(
            "netCDF4 is required to read WRF output files. "
            "Install with: pip install netCDF4"
        )
    return nc.Dataset(path, "r")


def _get_var(ds, *candidates):
    """Return the first variable found in *candidates*; raise if none found."""
    for name in candidates:
        if name in ds.variables:
            return ds.variables[name]
    raise KeyError(
        f"None of the expected variables {candidates} found in file. "
        f"Available variables: {list(ds.variables.keys())}"
    )


def _destagger_u(U_stag):
    """Average staggered U (shape ny, nx+1) to mass points (shape ny, nx)."""
    return 0.5 * (U_stag[:, :-1] + U_stag[:, 1:])


def _destagger_v(V_stag):
    """Average staggered V (shape ny+1, nx) to mass points (shape ny, nx)."""
    return 0.5 * (V_stag[:-1, :] + V_stag[1:, :])


def _write_wind(path, x_utm, y_utm, u, v):
    """Write wind file: utm_x  utm_y  u  v (m/s, one row per point)."""
    with open(path, "w") as fh:
        fh.write("# utm_x utm_y u v (m/s)\n")
        for x, y, uv, vv in zip(
            x_utm.ravel(), y_utm.ravel(), u.ravel(), v.ravel()
        ):
            fh.write(f"{x:.2f} {y:.2f} {float(uv):.6f} {float(vv):.6f}\n")


_WRF_BBOX_HALF_SPAN = 0.45  # degrees; SRTM download limit is 1°×1°


def read_wrf_bbox(wrf_path):
    """Read a clipped lat/lon bounding box from a WRF netCDF file.

    The bounding box is centred on the domain centre and extends
    ``_WRF_BBOX_HALF_SPAN`` degrees (0.45°) in each direction, keeping the
    total span below the SRTM 1-degree-tile download limit.

    Parameters
    ----------
    wrf_path : str
        Path to the WRF output netCDF file.

    Returns
    -------
    tuple of float
        (lat_min, lat_max, lon_min, lon_max) in WGS-84 decimal degrees,
        each being ``center ± 0.45``.
    """
    import numpy as np

    ds = _open_nc(wrf_path)
    try:
        lat_var = _get_var(ds, "XLAT", "lat", "latitude", "LAT")
        lon_var = _get_var(ds, "XLONG", "XLON", "lon", "longitude", "LON")

        lat_arr = np.array(lat_var[:], dtype=np.float64)
        lon_arr = np.array(lon_var[:], dtype=np.float64)
    finally:
        ds.close()

    # Collapse any time dimension
    if lat_arr.ndim == 3:
        lat_arr = lat_arr[0]
        lon_arr = lon_arr[0]

    center_lat = 0.5 * (float(np.min(lat_arr)) + float(np.max(lat_arr)))
    center_lon = 0.5 * (float(np.min(lon_arr)) + float(np.max(lon_arr)))

    return (
        center_lat - _WRF_BBOX_HALF_SPAN,
        center_lat + _WRF_BBOX_HALF_SPAN,
        center_lon - _WRF_BBOX_HALF_SPAN,
        center_lon + _WRF_BBOX_HALF_SPAN,
    )


def _parse_time_range(time_range_str):
    """Parse a 'T1:TN' string and return a list of integer indices.

    The range is **inclusive** on both ends: '0:3' → [0, 1, 2, 3].

    Parameters
    ----------
    time_range_str : str
        A string of the form 'T1:TN' (e.g. '0:4').

    Returns
    -------
    list of int

    Raises
    ------
    ValueError
        If the string cannot be parsed or T1 > TN.
    """
    parts = time_range_str.split(":")
    if len(parts) != 2:
        raise ValueError(
            f"--time-range must be in the form 'T1:TN', got: {time_range_str!r}"
        )
    try:
        t1, tn = int(parts[0]), int(parts[1])
    except ValueError:
        raise ValueError(
            f"--time-range indices must be integers, got: {time_range_str!r}"
        )
    if t1 > tn:
        raise ValueError(
            f"--time-range T1 must be <= TN, got {t1}:{tn}"
        )
    return list(range(t1, tn + 1))


def _wind_output_path(base_path, position):
    """Return the wind output path for *position* (0-based) within a time range.

    position=0 returns *base_path* unchanged (e.g. ``wind.csv``).
    position>0 inserts ``_<position>`` before the file extension
    (e.g. position=1 → ``wind_1.csv``, position=2 → ``wind_2.csv``).

    This matches the naming convention expected by ``velocity_field.H`` for
    time-dependent wind-field loading (base file at t=0, ``base_N.csv`` at
    subsequent time steps).
    """
    if position == 0:
        return base_path
    root, ext = os.path.splitext(base_path)
    return f"{root}_{position}{ext}"


def read_wrf_time_spacing(wrf_path):
    """Return the time spacing (seconds) between consecutive WRF snapshots.

    Tries ``XTIME`` (float, minutes from start) first, then ``Times``
    (character array in ``YYYY-MM-DD_HH:MM:SS`` format).  Falls back to
    3600 s (1 hour) if neither variable is present or only one snapshot
    exists.

    Parameters
    ----------
    wrf_path : str
        Path to the WRF output netCDF file.

    Returns
    -------
    float
        Time spacing in seconds.
    """
    import numpy as np

    ds = _open_nc(wrf_path)
    try:
        if "XTIME" in ds.variables:
            xtime = np.array(ds.variables["XTIME"][:], dtype=np.float64)
            if xtime.size > 1:
                return float((xtime[1] - xtime[0]) * 60.0)

        if "Times" in ds.variables:
            times_var = ds.variables["Times"]
            if times_var.shape[0] > 1:
                try:
                    from datetime import datetime

                    def _parse_wrf_time(chars):
                        s = "".join(c.decode("utf-8") if isinstance(c, bytes)
                                    else c for c in chars).strip().replace("_", "T")
                        for fmt in ("%Y-%m-%dT%H:%M:%S", "%Y-%m-%d %H:%M:%S"):
                            try:
                                return datetime.strptime(s, fmt)
                            except ValueError:
                                continue
                        return None

                    t0 = _parse_wrf_time(times_var[0])
                    t1 = _parse_wrf_time(times_var[1])
                    if t0 is not None and t1 is not None:
                        return float((t1 - t0).total_seconds())
                except Exception:
                    pass
    finally:
        ds.close()

    print("WARNING: Could not determine WRF time spacing; defaulting to 3600 s",
          file=sys.stderr)
    return 3600.0


def _wrf_bbox_to_utm_domain_bounds(lat_min, lat_max, lon_min, lon_max):
    """Convert a WRF lat/lon bounding box to UTM domain bounds.

    Projects the four corners of the lat/lon bounding box to UTM metres and
    returns the enclosing axis-aligned rectangle.

    Parameters
    ----------
    lat_min, lat_max : float
        Latitude bounds in WGS-84 decimal degrees.
    lon_min, lon_max : float
        Longitude bounds in WGS-84 decimal degrees.

    Returns
    -------
    tuple of float
        ``(x_lo, y_lo, x_hi, y_hi)`` in UTM metres.
    """
    import numpy as np

    corner_lats = np.array([lat_min, lat_min, lat_max, lat_max])
    corner_lons = np.array([lon_min, lon_max, lon_min, lon_max])
    x_utm, y_utm = _latlon_to_utm(corner_lats, corner_lons)
    return (float(np.min(x_utm)), float(np.min(y_utm)),
            float(np.max(x_utm)), float(np.max(y_utm)))


def _read_domain_bounds(file_path):
    """Read x/y min/max from a terrain XYZ or landscape LCP file.

    Returns
    -------
    tuple of float or None
        (x_min, y_min, x_max, y_max), or ``None`` if the file cannot be read.
    """
    xs, ys = [], []
    try:
        with open(file_path) as fh:
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
    except OSError:
        return None
    if not xs:
        return None
    return min(xs), min(ys), max(xs), max(ys)


def write_inputs_file(output_path,
                      terrain_file=None,
                      landscape_file=None,
                      wind_base_file=None,
                      multi_time=False,
                      wind_time_spacing=None,
                      final_time=None,
                      domain_bounds=None):
    """Write an ``inputs.i`` file configured for a pure FARSITE run.

    The generated file enables FARSITE elliptical fire spread and disables
    firebrand spotting, crown-fire initiation, and level-set advection.

    Parameters
    ----------
    output_path : str
        Destination path for the ``inputs.i`` file.
    terrain_file : str or None
        Path to the terrain XYZ file (``rothermel.terrain_file``).
    landscape_file : str or None
        Path to the landscape LCP file (``rothermel.landscape_file``).
    wind_base_file : str or None
        Base path for the wind CSV file (``velocity_file``).
    multi_time : bool
        Whether multiple wind time steps were extracted
        (enables ``use_time_dependent_wind``).
    wind_time_spacing : float or None
        Wind time-step spacing in seconds for ``wind_time_spacing``.
    final_time : float or None
        Simulation stop time in seconds for ``final_time``.
    domain_bounds : tuple of float or None
        ``(x_min, y_min, x_max, y_max)`` in UTM metres.  Used to set
        ``prob_lo`` / ``prob_hi`` and a default box ignition.  When
        ``None``, placeholder comments are written instead.
    """
    bounds = domain_bounds  # may be None

    with open(output_path, "w") as fh:
        fh.write(
            "# inputs.i – FARSITE run (no spotting, no crown fire)\n"
            "# Auto-generated by terrain_wind_preprocess.py\n"
            "#\n"
            "# Edit the ignition box extents and fuel model before running.\n\n"
        )

        # --- Grid & domain -------------------------------------------
        fh.write("# Grid & domain\n")
        if bounds is not None:
            x_lo, y_lo, x_hi, y_hi = bounds
            n_cell_x = max(1, round((x_hi - x_lo) / 30.0))
            n_cell_y = max(1, round((y_hi - y_lo) / 30.0))
            fh.write(f"n_cell_x = {n_cell_x}\n")
            fh.write(f"n_cell_y = {n_cell_y}\n")
            fh.write(f"max_grid_size = {max(1, max(n_cell_x, n_cell_y) // 2)}\n")
            fh.write(f"prob_lo_x = {x_lo:.2f}\n")
            fh.write(f"prob_lo_y = {y_lo:.2f}\n")
            fh.write(f"prob_hi_x = {x_hi:.2f}\n")
            fh.write(f"prob_hi_y = {y_hi:.2f}\n")
        else:
            fh.write("# n_cell_x = <round((x_max - x_min) / 30)>  # 30 m resolution\n")
            fh.write("# n_cell_y = <round((y_max - y_min) / 30)>  # 30 m resolution\n")
            fh.write("# max_grid_size = <n_cell_x or n_cell_y / 2>\n")
            fh.write("# prob_lo_x = <x_min from terrain/landscape>\n")
            fh.write("# prob_lo_y = <y_min from terrain/landscape>\n")
            fh.write("# prob_hi_x = <x_max from terrain/landscape>\n")
            fh.write("# prob_hi_y = <y_max from terrain/landscape>\n")
        fh.write("\n")

        # --- Time & output -------------------------------------------
        fh.write("# Time & output\n")
        if final_time is not None and final_time > 0.0:
            fh.write(f"final_time = {final_time:.1f}\n")
        else:
            fh.write("# final_time = <simulation_duration_seconds>\n")
        fh.write("cfl = 0.5\n")
        fh.write("plot_int = 10\n")
        fh.write("reinit_int = -1\n\n")

        # --- Wind / velocity -----------------------------------------
        fh.write("# Wind\n")
        if wind_base_file is not None:
            fh.write(f"velocity_file = {os.path.basename(wind_base_file)}\n")
        else:
            fh.write("# velocity_file = wind.csv\n")
        if multi_time:
            fh.write("use_time_dependent_wind = 1\n")
            if wind_time_spacing is not None:
                fh.write(f"wind_time_spacing = {wind_time_spacing:.1f}\n")
        else:
            fh.write("use_time_dependent_wind = 0\n")
        fh.write("\n")

        # --- Ignition source (box) -----------------------------------
        fh.write("# Initial ignition source (adjust extents to your fire location)\n")
        fh.write("source_type = box\n")
        if bounds is not None:
            x_lo, y_lo, x_hi, y_hi = bounds
            w = x_hi - x_lo
            h = y_hi - y_lo
            fh.write(f"box_xmin = {x_lo + 0.45 * w:.2f}\n")
            fh.write(f"box_xmax = {x_lo + 0.55 * w:.2f}\n")
            fh.write(f"box_ymin = {y_lo + 0.45 * h:.2f}\n")
            fh.write(f"box_ymax = {y_lo + 0.55 * h:.2f}\n")
        else:
            fh.write("# box_xmin / box_xmax / box_ymin / box_ymax = <UTM metres>\n")
        fh.write("\n")

        # --- Rothermel / fuel ----------------------------------------
        fh.write("# Rothermel fuel model (edit as needed)\n")
        if landscape_file is None:
            fh.write(
                "# No landscape file: defaulting to Southern California chaparral "
                "(NFFL FM4)\n"
            )
        fh.write("rothermel.fuel_model = FM4\n")
        fh.write("rothermel.M_f = 0.08\n")
        if terrain_file is not None:
            fh.write(f"rothermel.terrain_file = {os.path.basename(terrain_file)}\n")
        if landscape_file is not None:
            fh.write(
                f"rothermel.landscape_file = {os.path.basename(landscape_file)}\n"
            )
        fh.write("\n")

        # --- FARSITE model (no spotting, no crown) -------------------
        fh.write("# FARSITE ellipse model (Richards 1990)\n")
        fh.write("skip_levelset = 1\n")
        fh.write("farsite.enable = 1\n")
        fh.write("farsite.use_anderson_LW = 1\n")
        fh.write("farsite.phi_threshold = 0.1\n\n")

        fh.write("# Disabled sub-models\n")
        fh.write("spotting.enable = 0\n")
        fh.write("crown.enable = 0\n")
        fh.write("albini_spotting.enable = 0\n")

    print(f"Wrote inputs file: '{output_path}'.")


def interpolate_wind_to_grid(wrf_x, wrf_y, u, v, target_x, target_y):
    """Interpolate WRF wind components (u, v) onto a target (x, y) grid.

    Uses scipy's linear ``griddata`` interpolation.  Points outside the WRF
    domain are filled with 0.0.

    Parameters
    ----------
    wrf_x, wrf_y : array-like, shape (ny_wrf, nx_wrf)
        UTM coordinates of the WRF mass-point grid (metres).
    u, v : array-like, shape (ny_wrf, nx_wrf)
        Wind components at WRF mass-point locations (m/s).
    target_x, target_y : array-like, shape (ny_srtm, nx_srtm)
        UTM coordinates of the target (e.g. SRTM) grid (metres).

    Returns
    -------
    u_interp, v_interp : numpy.ndarray, shape matching *target_x*
    """
    import numpy as np
    try:
        from scipy.interpolate import griddata
    except ImportError:
        raise ImportError(
            "scipy is required for wind interpolation. "
            "Install with: pip install scipy"
        )

    wrf_pts = np.column_stack([
        np.asarray(wrf_x).ravel(),
        np.asarray(wrf_y).ravel(),
    ])
    target_pts = np.column_stack([
        np.asarray(target_x).ravel(),
        np.asarray(target_y).ravel(),
    ])
    target_shape = np.asarray(target_x).shape

    u_interp = griddata(wrf_pts, np.asarray(u).ravel(),
                        target_pts, method="linear", fill_value=0.0)
    v_interp = griddata(wrf_pts, np.asarray(v).ravel(),
                        target_pts, method="linear", fill_value=0.0)

    return u_interp.reshape(target_shape), v_interp.reshape(target_shape)


def extract_wrf_wind(wrf_path, time_index=0, level=0, subsample=1,
                     lat_min=None, lat_max=None,
                     lon_min=None, lon_max=None):
    """Extract wind and coordinate data from a WRF netCDF file.

    Parameters
    ----------
    wrf_path : str
        Path to the WRF output netCDF file.
    time_index : int
        Time snapshot index to extract (default 0).
    level : int
        Vertical model level for wind (default 0 = lowest level).
    subsample : int
        Keep every *subsample*-th point in each dimension (default 1).
    lat_min, lat_max, lon_min, lon_max : float or None
        Optional lat/lon bounding box (WGS-84 degrees).  When all four are
        provided the WRF grid is clipped to the tightest rectangular sub-grid
        whose grid points all fall within the specified bounds before UTM
        projection and subsampling.

    Returns
    -------
    x_utm, y_utm, u_mass, v_mass : numpy.ndarray (2-D)
        UTM coordinates and destaggered wind components.
    """
    import numpy as np

    ds = _open_nc(wrf_path)

    lat_var = _get_var(ds, "XLAT", "lat", "latitude", "LAT")
    lon_var = _get_var(ds, "XLONG", "XLON", "lon", "longitude", "LON")

    lat_arr = np.array(lat_var[:], dtype=np.float64)
    lon_arr = np.array(lon_var[:], dtype=np.float64)

    if lat_arr.ndim == 3:
        lat_2d = lat_arr[0]
        lon_2d = lon_arr[0]
    else:
        lat_2d = lat_arr
        lon_2d = lon_arr

    # U wind (staggered west-east)
    u_var = _get_var(ds, "U", "u", "U10", "u10")
    u_arr = np.array(u_var[:], dtype=np.float64)
    if u_arr.ndim == 4:
        U_2d_stag = u_arr[time_index, level]
    elif u_arr.ndim == 3:
        U_2d_stag = u_arr[time_index]
    else:
        U_2d_stag = u_arr

    # V wind (staggered south-north)
    v_var = _get_var(ds, "V", "v", "V10", "v10")
    v_arr = np.array(v_var[:], dtype=np.float64)
    if v_arr.ndim == 4:
        V_2d_stag = v_arr[time_index, level]
    elif v_arr.ndim == 3:
        V_2d_stag = v_arr[time_index]
    else:
        V_2d_stag = v_arr

    ds.close()

    ny, nx = lat_2d.shape

    # Destagger U and V to mass-point locations
    if U_2d_stag.shape == (ny, nx + 1):
        u_mass = _destagger_u(U_2d_stag)
    elif U_2d_stag.shape == (ny, nx):
        u_mass = U_2d_stag
    else:
        raise ValueError(
            f"Unexpected U shape {U_2d_stag.shape} for domain ({ny}, {nx})"
        )

    if V_2d_stag.shape == (ny + 1, nx):
        v_mass = _destagger_v(V_2d_stag)
    elif V_2d_stag.shape == (ny, nx):
        v_mass = V_2d_stag
    else:
        raise ValueError(
            f"Unexpected V shape {V_2d_stag.shape} for domain ({ny}, {nx})"
        )

    # Clip to lat/lon bounds if provided
    clip_bounds = (lat_min is not None and lat_max is not None
                   and lon_min is not None and lon_max is not None)
    if clip_bounds:
        mask = ((lat_2d >= lat_min) & (lat_2d <= lat_max) &
                (lon_2d >= lon_min) & (lon_2d <= lon_max))
        rows, cols = np.where(mask)
        if rows.size == 0:
            raise ValueError(
                f"No WRF grid points fall within the lat/lon bounds "
                f"lat=[{lat_min}, {lat_max}] lon=[{lon_min}, {lon_max}]."
            )
        r0, r1 = int(rows.min()), int(rows.max()) + 1
        c0, c1 = int(cols.min()), int(cols.max()) + 1
        lat_2d = lat_2d[r0:r1, c0:c1]
        lon_2d = lon_2d[r0:r1, c0:c1]
        u_mass = u_mass[r0:r1, c0:c1]
        v_mass = v_mass[r0:r1, c0:c1]
        ny, nx = lat_2d.shape
        print(f"Clipped WRF wind grid to lat=[{lat_min:.4f}, {lat_max:.4f}] "
              f"lon=[{lon_min:.4f}, {lon_max:.4f}]: {ny}×{nx} points.")

    # Project to UTM
    print(f"Projecting {ny}×{nx} WRF grid to UTM …")
    x_utm, y_utm = _latlon_to_utm(lat_2d, lon_2d)

    # Subsample
    if subsample > 1:
        sl = slice(None, None, subsample)
        x_utm  = x_utm [sl, sl]
        y_utm  = y_utm [sl, sl]
        u_mass = u_mass[sl, sl]
        v_mass = v_mass[sl, sl]
        ny_out, nx_out = x_utm.shape
        print(f"Subsampled to {ny_out}×{nx_out} = {ny_out * nx_out} WRF points.")

    return x_utm, y_utm, u_mass, v_mass


def convert_wrf(wrf_path, terrain_out, wind_out,
                time_index=0, level=0, subsample=1):
    """Extract terrain and wind data from a WRF netCDF file.

    .. deprecated::
        Use the unified ``main()`` entry point instead.  Kept for backward
        compatibility.  Writes WRF HGT_M as the terrain file.
    """
    import numpy as np

    ds = _open_nc(wrf_path)

    lat_var = _get_var(ds, "XLAT", "lat", "latitude", "LAT")
    lon_var = _get_var(ds, "XLONG", "XLON", "lon", "longitude", "LON")

    lat_arr = np.array(lat_var[:], dtype=np.float64)
    lon_arr = np.array(lon_var[:], dtype=np.float64)

    if lat_arr.ndim == 3:
        lat_2d = lat_arr[0]
        lon_2d = lon_arr[0]
    else:
        lat_2d = lat_arr
        lon_2d = lon_arr

    hgt_var = _get_var(ds, "HGT_M", "HGT", "hgt", "TERRAIN", "terrain", "elevation")
    hgt_arr = np.array(hgt_var[:], dtype=np.float64)
    hgt_2d = hgt_arr[time_index] if hgt_arr.ndim == 3 else hgt_arr

    u_var = _get_var(ds, "U", "u", "U10", "u10")
    u_arr = np.array(u_var[:], dtype=np.float64)
    if u_arr.ndim == 4:
        U_2d_stag = u_arr[time_index, level]
    elif u_arr.ndim == 3:
        U_2d_stag = u_arr[time_index]
    else:
        U_2d_stag = u_arr

    v_var = _get_var(ds, "V", "v", "V10", "v10")
    v_arr = np.array(v_var[:], dtype=np.float64)
    if v_arr.ndim == 4:
        V_2d_stag = v_arr[time_index, level]
    elif v_arr.ndim == 3:
        V_2d_stag = v_arr[time_index]
    else:
        V_2d_stag = v_arr

    ds.close()

    ny, nx = lat_2d.shape

    if U_2d_stag.shape == (ny, nx + 1):
        u_mass = _destagger_u(U_2d_stag)
    elif U_2d_stag.shape == (ny, nx):
        u_mass = U_2d_stag
    else:
        raise ValueError(f"Unexpected U shape {U_2d_stag.shape}")

    if V_2d_stag.shape == (ny + 1, nx):
        v_mass = _destagger_v(V_2d_stag)
    elif V_2d_stag.shape == (ny, nx):
        v_mass = V_2d_stag
    else:
        raise ValueError(f"Unexpected V shape {V_2d_stag.shape}")

    print(f"Projecting {ny}×{nx} grid to UTM …")
    x_utm, y_utm = _latlon_to_utm(lat_2d, lon_2d)

    if subsample > 1:
        sl = slice(None, None, subsample)
        hgt_2d = hgt_2d[sl, sl]
        u_mass = u_mass[sl, sl]
        v_mass = v_mass[sl, sl]
        x_utm  = x_utm [sl, sl]
        y_utm  = y_utm [sl, sl]
        ny_out, nx_out = x_utm.shape
        print(f"Subsampled to {ny_out}×{nx_out} = {ny_out * nx_out} points.")

    # Write terrain (WRF HGT_M) – legacy behaviour
    with open(terrain_out, "w") as fh:
        fh.write("# utm_x utm_y z (meters)\n")
        for x, y, zv in zip(x_utm.ravel(), y_utm.ravel(), hgt_2d.ravel()):
            fh.write(f"{x:.2f} {y:.2f} {float(zv):.6f}\n")
    print(f"Wrote terrain file: '{terrain_out}' ({x_utm.size} points)")

    _write_wind(wind_out, x_utm, y_utm, u_mass, v_mass)
    print(f"Wrote wind file:    '{wind_out}' ({x_utm.size} points)")


def extract_wrf_terrain(wrf_path, output_path, subsample=1,
                        lat_min=None, lat_max=None,
                        lon_min=None, lon_max=None):
    """Extract terrain height from a WRF netCDF file and write an XYZ file.

    Parameters
    ----------
    wrf_path : str
        Path to the WRF output netCDF file.
    output_path : str
        Path for the output terrain XYZ file (utm_x utm_y z).
    subsample : int
        Keep every *subsample*-th point in each dimension (default 1).
    lat_min, lat_max, lon_min, lon_max : float or None
        Optional lat/lon bounding box for clipping the WRF grid.

    Returns
    -------
    x_utm, y_utm : numpy.ndarray (2-D)
        UTM coordinates of the terrain grid points.
    """
    import numpy as np

    ds = _open_nc(wrf_path)
    try:
        lat_var = _get_var(ds, "XLAT", "lat", "latitude", "LAT")
        lon_var = _get_var(ds, "XLONG", "XLON", "lon", "longitude", "LON")
        hgt_var = _get_var(ds, "HGT_M", "HGT", "hgt", "TERRAIN",
                           "terrain", "elevation")
        lat_arr = np.array(lat_var[:], dtype=np.float64)
        lon_arr = np.array(lon_var[:], dtype=np.float64)
        hgt_arr = np.array(hgt_var[:], dtype=np.float64)
    finally:
        ds.close()

    lat_2d = lat_arr[0] if lat_arr.ndim == 3 else lat_arr
    lon_2d = lon_arr[0] if lon_arr.ndim == 3 else lon_arr
    hgt_2d = hgt_arr[0] if hgt_arr.ndim == 3 else hgt_arr

    # Clip to lat/lon bounds if provided
    clip_bounds = (lat_min is not None and lat_max is not None
                   and lon_min is not None and lon_max is not None)
    if clip_bounds:
        mask = ((lat_2d >= lat_min) & (lat_2d <= lat_max) &
                (lon_2d >= lon_min) & (lon_2d <= lon_max))
        rows, cols = np.where(mask)
        if rows.size == 0:
            raise ValueError(
                f"No WRF grid points fall within the lat/lon bounds "
                f"lat=[{lat_min}, {lat_max}] lon=[{lon_min}, {lon_max}]."
            )
        r0, r1 = int(rows.min()), int(rows.max()) + 1
        c0, c1 = int(cols.min()), int(cols.max()) + 1
        lat_2d = lat_2d[r0:r1, c0:c1]
        lon_2d = lon_2d[r0:r1, c0:c1]
        hgt_2d = hgt_2d[r0:r1, c0:c1]
        ny, nx = lat_2d.shape
        print(f"Clipped WRF terrain grid to lat=[{lat_min:.4f}, {lat_max:.4f}] "
              f"lon=[{lon_min:.4f}, {lon_max:.4f}]: {ny}×{nx} points.")

    ny, nx = lat_2d.shape
    print(f"Projecting {ny}×{nx} WRF terrain grid to UTM …")
    x_utm, y_utm = _latlon_to_utm(lat_2d, lon_2d)

    if subsample > 1:
        sl = slice(None, None, subsample)
        x_utm = x_utm[sl, sl]
        y_utm = y_utm[sl, sl]
        hgt_2d = hgt_2d[sl, sl]
        ny_out, nx_out = x_utm.shape
        print(f"Subsampled to {ny_out}×{nx_out} = {ny_out * nx_out} points.")

    with open(output_path, "w") as fh:
        fh.write("# utm_x utm_y z (meters)\n")
        for xv, yv, zv in zip(x_utm.ravel(), y_utm.ravel(), hgt_2d.ravel()):
            fh.write(f"{xv:.2f} {yv:.2f} {float(zv):.6f}\n")
    print(f"Wrote WRF terrain file: '{output_path}' ({x_utm.size} points)")

    return x_utm, y_utm


# ===========================================================================
# CLI
# ===========================================================================

def _build_parser():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # WRF input (optional; if given, bbox is read from the file)
    parser.add_argument(
        "--wrf-file", default=None, metavar="FILE",
        help=(
            "WRF netCDF file.  If given, the lat/lon bounding box is centred "
            "on the domain centre and clipped to ±0.45° in each direction; "
            "--lat-min/max/lon-min/max are ignored."
        ),
    )

    # Bounding box (required when --wrf-file is not given)
    parser.add_argument("--lat-min", type=float, default=None,
                        help="Minimum latitude  (WGS-84 decimal degrees)")
    parser.add_argument("--lat-max", type=float, default=None,
                        help="Maximum latitude  (WGS-84 decimal degrees)")
    parser.add_argument("--lon-min", type=float, default=None,
                        help="Minimum longitude (WGS-84 decimal degrees)")
    parser.add_argument("--lon-max", type=float, default=None,
                        help="Maximum longitude (WGS-84 decimal degrees)")

    # Output files
    parser.add_argument("--terrain",   default="terrain.xyz",
                        help="Output terrain XYZ file   (default: terrain.xyz)")
    parser.add_argument("--landscape", default="landscape.lcp",
                        help="Output landscape LCP file  (default: landscape.lcp)")
    parser.add_argument("--wind",      default="wind.csv",
                        help="Output wind CSV file       (default: wind.csv); "
                             "requires --wrf-file")
    parser.add_argument("--tif", default=None, metavar="FILE",
                        help="Path for the intermediate SRTM GeoTIFF "
                             "(temporary file used and removed if not specified)")

    # Control which outputs to produce
    parser.add_argument("--no-terrain",   action="store_true",
                        help="Skip the SRTM terrain XYZ and landscape steps")
    parser.add_argument("--no-landscape", action="store_true",
                        help="Skip the LANDFIRE landscape LCP step")
    parser.add_argument("--no-wind",      action="store_true",
                        help="Skip the WRF wind extraction step")

    # Wind-specific options
    parser.add_argument(
        "--interpolate-wind", action="store_true",
        help=(
            "Interpolate WRF wind fields to the SRTM terrain grid.  "
            "Requires --wrf-file and that --no-terrain is NOT set."
        ),
    )
    parser.add_argument(
        "--time-range", default=None, metavar="T1:TN",
        help=(
            "Inclusive range of WRF time indices to extract, e.g. '0:4' "
            "extracts indices 0,1,2,3,4.  Overrides --time-index when given."
        ),
    )
    parser.add_argument(
        "--time-index", type=int, default=0, metavar="N",
        help="Single WRF time index to extract (default: 0)",
    )
    parser.add_argument(
        "--level", type=int, default=0, metavar="N",
        help="WRF vertical level for U/V (default: 0 = lowest)",
    )

    # Common options
    parser.add_argument("--subsample", type=int, default=1, metavar="N",
                        help="Keep every N-th point in each dimension (default: 1)")

    # LANDFIRE-specific options
    parser.add_argument("--vintage", type=int, default=2020, metavar="YEAR",
                        help="LANDFIRE data vintage year (default: 2020)")
    parser.add_argument("--elev-product",   default=None, metavar="ID",
                        help="Override LANDFIRE elevation product ID")
    parser.add_argument("--slope-product",  default=None, metavar="ID",
                        help="Override LANDFIRE slope product ID")
    parser.add_argument("--aspect-product", default=None, metavar="ID",
                        help="Override LANDFIRE aspect product ID")
    parser.add_argument("--fuel-product",   default=None, metavar="ID",
                        help="Override LANDFIRE fuel model product ID")
    parser.add_argument("--keep-nonburnable", action="store_true",
                        help="Include non-burnable LANDFIRE pixels in LCP output")
    parser.add_argument("--cache-dir", default=None, metavar="DIR",
                        help="Directory to cache downloaded LANDFIRE ZIP files "
                             "(only used with --use-lfps)")
    parser.add_argument("--timeout", type=int, default=300, metavar="N",
                        help="LANDFIRE API polling timeout in seconds "
                             "(default: 300; only used with --use-lfps)")
    parser.add_argument(
        "--use-lfps", action="store_true",
        help=(
            "Use the legacy LFPS REST API to download LANDFIRE data instead "
            "of the default COG (S3/HTTPS) approach.  This may be needed if "
            "the S3 endpoint is unreachable, but is susceptible to LFPS SSL "
            "certificate verification failures."
        ),
    )

    # Local-file landscape mode
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

    # inputs.i generation
    parser.add_argument(
        "--inputs", default="inputs.i", metavar="FILE",
        help="Output inputs.i file for a FARSITE run (default: inputs.i)",
    )
    parser.add_argument(
        "--no-inputs", action="store_true",
        help="Skip automatic generation of the inputs.i file",
    )

    return parser


def main(argv=None):
    parser = _build_parser()
    args = parser.parse_args(argv)

    # -----------------------------------------------------------------------
    # Validate / resolve bounding box
    # -----------------------------------------------------------------------
    if args.wrf_file is not None:
        if not os.path.isfile(args.wrf_file):
            print(f"ERROR: WRF file not found: {args.wrf_file}", file=sys.stderr)
            sys.exit(1)
        lat_min, lat_max, lon_min, lon_max = read_wrf_bbox(args.wrf_file)
        print(f"Bounding box from WRF file (centre ±0.45°): "
              f"lat=[{lat_min:.4f}, {lat_max:.4f}] "
              f"lon=[{lon_min:.4f}, {lon_max:.4f}]")
    else:
        # Bounding box required from CLI
        missing = [name for name, val in [
            ("--lat-min", args.lat_min), ("--lat-max", args.lat_max),
            ("--lon-min", args.lon_min), ("--lon-max", args.lon_max),
        ] if val is None]
        if missing:
            print(
                f"ERROR: {', '.join(missing)} are required when --wrf-file is "
                f"not given.",
                file=sys.stderr,
            )
            sys.exit(1)
        lat_min = args.lat_min
        lat_max = args.lat_max
        lon_min = args.lon_min
        lon_max = args.lon_max

        if lon_min >= lon_max or lat_min >= lat_max:
            print(
                "ERROR: bounding box must satisfy "
                "lon-min < lon-max and lat-min < lat-max",
                file=sys.stderr,
            )
            sys.exit(1)

    if args.subsample < 1:
        print("ERROR: --subsample must be >= 1", file=sys.stderr)
        sys.exit(1)

    bbox = (lon_min, lat_min, lon_max, lat_max)

    # -----------------------------------------------------------------------
    # Resolve time indices for WRF
    # -----------------------------------------------------------------------
    if args.time_range is not None:
        try:
            time_indices = _parse_time_range(args.time_range)
        except ValueError as exc:
            print(f"ERROR: {exc}", file=sys.stderr)
            sys.exit(1)
    else:
        time_indices = [args.time_index]

    # -----------------------------------------------------------------------
    # SRTM terrain XYZ step
    # -----------------------------------------------------------------------
    srtm_x = srtm_y = None  # kept for wind interpolation if needed

    if not args.no_terrain:
        tif_path = (os.path.abspath(args.tif)
                    if args.tif is not None else None)
        out_terrain = os.path.abspath(args.terrain)

        if args.interpolate_wind and args.wrf_file is not None:
            # We need the SRTM grid arrays for later interpolation
            srtm_x, srtm_y, _ = create_terrain_xyz_return_grid(
                lat_min=lat_min, lat_max=lat_max,
                lon_min=lon_min, lon_max=lon_max,
                output_path=out_terrain,
                tif_path=tif_path,
                subsample=args.subsample,
            )
        else:
            create_terrain_xyz(
                out_terrain,
                lat_min=lat_min, lat_max=lat_max,
                lon_min=lon_min, lon_max=lon_max,
                tif_path=tif_path,
                subsample=args.subsample,
            )
    elif args.no_terrain and args.wrf_file is not None:
        # When --no-terrain is set together with --wrf-file, write the WRF
        # HGT / HGT_M field as the terrain output instead of downloading SRTM.
        out_terrain = os.path.abspath(args.terrain)
        extract_wrf_terrain(
            args.wrf_file,
            out_terrain,
            subsample=args.subsample,
            lat_min=lat_min,
            lat_max=lat_max,
            lon_min=lon_min,
            lon_max=lon_max,
        )

    # -----------------------------------------------------------------------
    # LANDFIRE landscape LCP step
    # -----------------------------------------------------------------------
    if not args.no_terrain and not args.no_landscape:
        out_lcp = os.path.abspath(args.landscape)
        local_files = [args.elev_file, args.slope_file,
                       args.aspect_file, args.fuel_file]
        using_local = any(f is not None for f in local_files)

        if using_local:
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

    # -----------------------------------------------------------------------
    # WRF wind extraction step
    # -----------------------------------------------------------------------
    wind_base_out = None   # set below if wind was extracted
    multi_time = False

    if args.wrf_file is not None and not args.no_wind:
        multi_time = len(time_indices) > 1
        wind_base_out = os.path.abspath(args.wind)

        for pos, t_idx in enumerate(time_indices):
            print(f"\nExtracting WRF wind for time index {t_idx} …")
            wrf_x, wrf_y, u_mass, v_mass = extract_wrf_wind(
                args.wrf_file,
                time_index=t_idx,
                level=args.level,
                subsample=args.subsample,
                lat_min=lat_min,
                lat_max=lat_max,
                lon_min=lon_min,
                lon_max=lon_max,
            )

            wind_out = (_wind_output_path(wind_base_out, pos)
                        if multi_time else wind_base_out)

            if args.interpolate_wind:
                if args.no_terrain:
                    print(
                        "WARNING: --interpolate-wind has no effect when "
                        "--no-terrain is set; writing WRF-resolution wind.",
                        file=sys.stderr,
                    )
                    _write_wind(wind_out, wrf_x, wrf_y, u_mass, v_mass)
                elif srtm_x is None:
                    # SRTM grid not available (no --interpolate-wind was
                    # triggered with grid capture above); fall back gracefully
                    print(
                        "WARNING: SRTM grid not available for interpolation; "
                        "writing WRF-resolution wind.",
                        file=sys.stderr,
                    )
                    _write_wind(wind_out, wrf_x, wrf_y, u_mass, v_mass)
                else:
                    print(
                        f"Interpolating wind from WRF grid "
                        f"({wrf_x.size} pts) to SRTM grid "
                        f"({srtm_x.size} pts) …"
                    )
                    u_interp, v_interp = interpolate_wind_to_grid(
                        wrf_x, wrf_y, u_mass, v_mass, srtm_x, srtm_y
                    )
                    _write_wind(wind_out, srtm_x, srtm_y, u_interp, v_interp)
                    print(f"Wrote wind file: '{wind_out}' "
                          f"({srtm_x.size} points on SRTM grid, "
                          f"utm_x utm_y u v, time index {t_idx})")
                    continue
            else:
                _write_wind(wind_out, wrf_x, wrf_y, u_mass, v_mass)

            print(f"Wrote wind file: '{wind_out}' ({wrf_x.size} points, "
                  f"utm_x utm_y u v, time index {t_idx})")

    # -----------------------------------------------------------------------
    # Auto-generate inputs.i for a FARSITE run
    # -----------------------------------------------------------------------
    if not getattr(args, "no_inputs", False):
        # Determine domain bounds:
        # 1. Prefer reading from landscape or terrain file when terrain is active.
        # 2. When --no-terrain + --wrf-file, compute from WRF lat/lon bbox → UTM.
        domain_bounds = None
        for candidate in [
            os.path.abspath(args.landscape) if not args.no_terrain and not args.no_landscape else None,
            os.path.abspath(args.terrain) if not args.no_terrain else None,
        ]:
            if candidate is not None and os.path.isfile(candidate):
                domain_bounds = _read_domain_bounds(candidate)
                if domain_bounds is not None:
                    break

        if domain_bounds is None and args.no_terrain and args.wrf_file is not None:
            # No terrain file: derive domain extents directly from the WRF
            # lat/lon bounding box converted to UTM at 30 m resolution.
            print("Computing domain bounds from WRF bbox UTM extents "
                  "(--no-terrain mode) …")
            domain_bounds = _wrf_bbox_to_utm_domain_bounds(
                lat_min, lat_max, lon_min, lon_max
            )
            x_lo, y_lo, x_hi, y_hi = domain_bounds
            n_x = max(1, round((x_hi - x_lo) / 30.0))
            n_y = max(1, round((y_hi - y_lo) / 30.0))
            print(f"  UTM extents: x=[{x_lo:.0f}, {x_hi:.0f}]  "
                  f"y=[{y_lo:.0f}, {y_hi:.0f}]")
            print(f"  n_cell_x = {n_x}  n_cell_y = {n_y}  (30 m resolution)")

        # Determine WRF time spacing and final_time
        wind_time_spacing = None
        final_time = None

        if args.wrf_file is not None:
            wind_time_spacing = read_wrf_time_spacing(args.wrf_file)

        if args.time_range is not None and wind_time_spacing is not None:
            t_indices = _parse_time_range(args.time_range)
            if len(t_indices) > 1:
                t1 = t_indices[0]
                tn = t_indices[-1]
                final_time = (tn - t1) * wind_time_spacing

        terrain_out = (
            os.path.abspath(args.terrain)
            if (not args.no_terrain or args.wrf_file is not None)
            else None
        )
        landscape_out = (os.path.abspath(args.landscape)
                         if not args.no_terrain and not args.no_landscape else None)

        write_inputs_file(
            output_path=os.path.abspath(args.inputs),
            terrain_file=terrain_out if (terrain_out and os.path.isfile(terrain_out)) else None,
            landscape_file=landscape_out if (landscape_out and os.path.isfile(landscape_out)) else None,
            wind_base_file=wind_base_out,
            multi_time=multi_time,
            wind_time_spacing=wind_time_spacing,
            final_time=final_time,
            domain_bounds=domain_bounds,
        )


if __name__ == "__main__":
    main()
