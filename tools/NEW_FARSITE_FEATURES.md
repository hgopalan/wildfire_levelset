# New FARSITE-Compatible Features

This document describes four new "easy-to-implement" FARSITE features that have been added to the wildfire_levelset toolkit.

## Feature 3: Fire Size Percentile Statistics

**Enhancement to:** `tools/fire_size_summary.py`

### What it does
Adds percentile statistics (10th, 50th/median, 90th percentile) to the fire size summary output, similar to FARSITE's statistical fire growth analysis.

### Usage
```bash
# Standard usage - percentiles are automatically displayed
python3 tools/fire_size_summary.py fire_stats.csv
```

### Output Example
```
==================================================================================
PERCENTILE STATISTICS (Fire Growth Distribution)
==================================================================================
Metric                         10th %   50th % (Median)         90th %
--------------------------------------------------------------------------------
Burned Area [ha]               12.3456           45.6789        123.4567
Perimeter [km]                  1.2345            3.4567          8.9012
Active Front Cells            123.0000          456.0000        789.0000
Head ROS [m/s]                  0.0012            0.0045          0.0123
```

### Dependencies
- Requires `numpy` for percentile calculations
- Gracefully degrades if numpy is not available

---

## Feature 4: Isochrone Time Labels with Visualization

**Enhancement to:** `tools/isochrone_extractor.py`

### What it does
Adds optional visualization of fire arrival-time isochrones with time labels, similar to FARSITE's isochrone display. Isochrones are already labeled in the GeoJSON output; this adds matplotlib-based plotting.

### Usage
```bash
# Extract isochrones with visualization
python3 tools/isochrone_extractor.py plt0100 --interval 600 \
    --plot isochrones.png --outdir iso_out

# Process all plotfiles with plots
python3 tools/isochrone_extractor.py --all --interval 600 --plot auto
```

### Features
- Overlays isochrone contours on arrival time heatmap
- Labels each isochrone with its time (e.g., "t=10.0 min")
- Color-coded isochrones with legend
- Saves as PNG with 150 DPI resolution

---

## Feature 5: Minimum Travel Path Tool

**New tool:** `tools/minimum_travel_path.py`

### What it does
Extracts minimum travel time (MTT) paths from ignition points to specified destinations. Follows the steepest descent gradient of the arrival_time field, similar to FARSITE's path analysis capabilities.

### Usage
```bash
# Single destination point
python3 tools/minimum_travel_path.py plt0100 --dest 5000 3000 \
    --output path.csv

# Multiple destinations with visualization
python3 tools/minimum_travel_path.py plt0100 \
    --dest 5000 3000 6000 3500 \
    --output paths.csv --plot paths.png

# Customize path tracing
python3 tools/minimum_travel_path.py plt0100 --dest 5000 3000 \
    --output path.csv --max-steps 20000
```

### Output
**CSV file** with columns:
- `path_id` - Path identifier (0, 1, 2, ...)
- `x`, `y` - Physical coordinates [m]
- `arrival_time_s` - Arrival time in seconds
- `arrival_time_min` - Arrival time in minutes

**Plot** (optional):
- Paths overlaid on arrival time contours
- Start points (destinations) marked with circles
- End points (ignition) marked with stars
- Color-coded by path ID

### Use Cases
- Understanding fire spread paths
- Identifying critical control points
- Analyzing firebreak placement effectiveness
- Visualizing dominant fire progression routes

---

## Feature 6: Fire Period Analysis (Day/Night Burning)

**New tool:** `tools/fire_period_analysis.py`

### What it does
Classifies burned cells as "day" or "night" ignition based on burn period settings, similar to FARSITE's burn period concept. Analyzes when fire spread occurred relative to the diurnal burning window.

### Usage
```bash
# Read burn period from inputs file
python3 tools/fire_period_analysis.py plt0100 --inputs inputs.i \
    --output burn_period.csv

# Explicit burn period parameters
python3 tools/fire_period_analysis.py plt0100 \
    --start-hour 10.0 --end-hour 20.0 --sim-start-hour 12.0 \
    --output burn_period.csv --plot burn_period.png

# With visualization
python3 tools/fire_period_analysis.py plt0100 --inputs inputs.i \
    --output burn_period.csv --plot burn_period_map.png
```

### Burn Period Parameters
- `--start-hour`: Burn period start (e.g., 10.0 = 10:00 AM)
- `--end-hour`: Burn period end (e.g., 20.0 = 8:00 PM)
- `--sim-start-hour`: Simulation start time in local clock

These are read from the inputs file under `burn_period.*` parameters if available.

### Output
**CSV file** with statistics:
```csv
total_burned_area_ha,day_burned_area_ha,night_burned_area_ha,day_fraction,night_fraction,...
1234.56,987.65,246.91,0.80,0.20,...
```

**Plot** (optional):
- Spatial map showing day (orange) vs night (blue) burning
- Summary statistics overlay
- Grid with physical coordinates

### Interpretation
- **Day burning**: Fire spread occurred during the burn period window (typically 10:00-20:00)
- **Night burning**: Fire spread occurred outside the burn period window
- Higher day fraction indicates fire behaved according to typical diurnal patterns
- Night burning may indicate extreme fire conditions or long-range spotting

---

## Installation

All tools require Python 3.7+ with numpy:
```bash
pip install numpy matplotlib
```

Matplotlib is optional for visualization features but recommended.

---

## References

These features align with FARSITE capabilities:
- **Percentile statistics**: Similar to FARSITE fire growth percentile reporting
- **Isochrone visualization**: Matches FARSITE isochrone display
- **MTT paths**: Based on Finney (2002) minimum travel time methods
- **Burn period**: Implements FARSITE/FSPro burn period concept

### Citations
- Finney, M.A. (1998). FARSITE: Fire Area Simulator. USDA Forest Service RMRS-RP-4.
- Finney, M.A. (2002). Fire growth using minimum travel time methods. Canadian Journal of Forest Research, 32(8), 1420-1424.
