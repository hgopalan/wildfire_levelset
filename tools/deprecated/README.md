# Deprecated Tools

The scripts in this directory have been superseded by the unified
`tools/srtm_landfire_to_terrain.py` tool, which integrates SRTM terrain
download, DEM conversion, and LANDFIRE landscape creation into a single
command.

| Deprecated script | Functionality now in |
|---|---|
| `srtm_to_xyz_stl.py` | `tools/srtm_landfire_to_terrain.py` |
| `dem_to_xyz.py`       | `tools/srtm_landfire_to_terrain.py` |
| `landfire_to_lcp.py`  | `tools/srtm_landfire_to_terrain.py` |

These scripts are kept here for reference and backward compatibility.
They are **not** actively maintained.  Please use
`tools/srtm_landfire_to_terrain.py` for new workflows.
