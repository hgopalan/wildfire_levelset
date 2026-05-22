# HRRR Wind Reader Example

This document provides a practical example of how to use the `hrrr_wind_reader.py` tool
to extract wind data from NOAA's High Resolution Rapid Refresh (HRRR) model and prepare
it for use with the wildfire_levelset solver.

## Installation

First, install the required dependencies:

```bash
pip install herbie-data numpy pyproj scipy
```

- **herbie-data**: For fetching HRRR data from NOAA
- **numpy**: Numerical computations
- **pyproj**: Coordinate system transformations (UTM projection)
- **scipy**: Optional, for better interpolation quality

## Basic Usage

### Single Time Snapshot

Extract wind data for a specific date and time:

```bash
python3 tools/hrrr_wind_reader.py \
    --lat-min 40.0 \
    --lat-max 40.5 \
    --lon-min -105.0 \
    --lon-max -104.5 \
    --date-time "2023-08-15 14:00" \
    --wind wind.csv
```

This creates:
- `wind.csv`: Wind data on a regular 100 m grid (default)
- `inputs.i`: Minimal solver configuration stub

### Custom Grid Resolution

To use a finer grid (e.g., 30 m resolution):

```bash
python3 tools/hrrr_wind_reader.py \
    --lat-min 40.0 \
    --lat-max 40.5 \
    --lon-min -105.0 \
    --lon-max -104.5 \
    --date-time "2023-08-15 14:00" \
    --wind wind.csv \
    --grid-resolution 30
```

### Multiple Time Steps

Extract wind at multiple forecast times (e.g., every hour for 4 hours):

```bash
python3 tools/hrrr_wind_reader.py \
    --lat-min 40.0 \
    --lat-max 40.5 \
    --lon-min -105.0 \
    --lon-max -104.5 \
    --date-time "2023-08-15 14:00" \
    --wind wind.csv \
    --time-indices 0 1 2 3
```

This creates:
- `wind.csv`: Analysis time (f00)
- `wind_1.csv`: +1 hour forecast (f01)
- `wind_2.csv`: +2 hour forecast (f02)
- `wind_3.csv`: +3 hour forecast (f03)

## Understanding the Output

The `wind.csv` file has the format:

```
# utm_x utm_y u v (m/s)
<easting_1> <northing_1> <u_1> <v_1>
<easting_2> <northing_2> <u_2> <v_2>
...
```

Where:
- **utm_x**: Easting coordinate in metres (UTM projection)
- **utm_y**: Northing coordinate in metres (UTM projection)
- **u**: Zonal wind component (m/s) — positive is eastward
- **v**: Meridional wind component (m/s) — positive is northward

The solver computes wind speed magnitude and direction from the U and V components.

## Workflow Example: Fire Simulation Setup

Here's a complete workflow to set up a fire simulation with HRRR wind data:

```bash
# 1. Get terrain data (see srtm_terrain_reader.py)
python3 tools/srtm_terrain_reader.py \
    --lat-min 40.0 --lat-max 40.5 \
    --lon-min -105.0 --lon-max -104.5 \
    --output terrain.xyz

# 2. Get landscape data (see landscape_writer.py)
python3 tools/landscape_writer.py \
    --lat-min 40.0 --lat-max 40.5 \
    --lon-min -105.0 --lon-max -104.5 \
    --output landscape.lcp

# 3. Extract HRRR wind data
python3 tools/hrrr_wind_reader.py \
    --lat-min 40.0 --lat-max 40.5 \
    --lon-min -105.0 --lon-max -104.5 \
    --date-time "2023-08-15 14:00" \
    --wind wind.csv

# 4. Create fuel moisture file (example with constant moisture)
cat > moisture.csv << EOF
fuel_model,moisture_pct
1,8
2,9
3,10
EOF

# 5. Update inputs.i stub with additional parameters
# (See inputs.i documentation for full configuration)

# 6. Run the solver
./wildfire_levelset inputs.i
```

## Common Issues and Troubleshooting

### Issue: `herbie-data` module not found

**Solution**: Install herbie-data explicitly:

```bash
pip install --upgrade herbie-data
```

### Issue: No HRRR data available for the requested date/time

**Solution**: 
- HRRR data is typically available for recent dates (usually up to ~25 days in the past)
- Use a more recent date, e.g., yesterday:
  ```bash
  # Get yesterday's date in YYYY-MM-DD format
  python3 -c "from datetime import datetime, timedelta; \
    print((datetime.now() - timedelta(days=1)).strftime('%Y-%m-%d 12:00'))"
  ```

### Issue: NaN values in wind.csv

**Solution**:
- This typically indicates the grid is outside the HRRR data bounds
- Verify your lat/lon bounds are valid (latitude: -90 to 90, longitude: -180 to 180)
- Ensure the bounding box is contiguous (min < max)

### Issue: Slow data fetching

**Solution**:
- HRRR data must be downloaded from NOAA servers, which may be slow
- For operational use, consider downloading data in advance and caching it locally
- Use herbie's local caching features (see herbie documentation)

## Advanced Options

### Skip inputs.i Generation

```bash
python3 tools/hrrr_wind_reader.py \
    --lat-min 40.0 --lat-max 40.5 \
    --lon-min -105.0 --lon-max -104.5 \
    --date-time "2023-08-15 14:00" \
    --wind wind.csv \
    --no-inputs
```

### Custom inputs.i Output Name

```bash
python3 tools/hrrr_wind_reader.py \
    --lat-min 40.0 --lat-max 40.5 \
    --lon-min -105.0 --lon-max -104.5 \
    --date-time "2023-08-15 14:00" \
    --wind wind.csv \
    --inputs solver_config.i
```

## Data Quality Notes

- HRRR provides 10-meter wind components from the surface-level analysis
- Wind data are interpolated onto a regular grid; values between HRRR grid points use bilinear interpolation (or nearest-neighbor if scipy is not installed)
- HRRR data are typically available with ~1 km horizontal resolution
- Time-dependent forecasts (f01, f02, etc.) allow the solver to simulate evolving wind patterns

## References

- NOAA HRRR: https://www.ncei.noaa.gov/products/high-resolution-rapid-refresh
- Herbie documentation: https://herbie.readthedocs.io/
- wildfire_levelset solver: https://github.com/hgopalan/wildfire_levelset
