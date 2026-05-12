# Satellite Fire Detection Assimilation Regression Test

## Purpose

Tests the **satellite fire detection assimilation** feature (`satellite.*`):
reading a pre-generated CSV of lon/lat active-fire detections, converting
them to simulation-domain coordinates via UTM projection, and applying them
as both an **initial condition** (t = 0) and a **mid-simulation re-marking
event** (at t ≈ 300 s).

This mirrors the operational workflow of ingesting GOES or VIIRS active-fire
detections into a running simulation to correct the fire perimeter based on
the latest satellite observations.

## Configuration

| Parameter | Value |
|-----------|-------|
| Domain | 20 km × 20 km (UTM Zone 10 N) |
| Grid | 64 × 64 cells |
| Wind | 5 m/s eastward |
| Fuel | FM4 (Southern California chaparral) |
| Propagation | `levelset` |
| Simulation time | 700 s |
| Satellite source | `file` (pre-generated CSV) |
| Detection radius | 500 m per point |
| Fetch interval | 300 s |
| Use as IC | yes (`satellite.use_as_ic = 1`) |
| Mid-sim re-marking | yes (`satellite.use_mid_sim = 1`) |
| Suppress if burning | yes (`satellite.suppress_if_burning = 1`) |

## Setup

Run the Python setup script to create the detection CSV:

```bash
python3 create_detections.py
./build/levelset regtest/ignition/satellite_assimilation/inputs.i
```

## Expected Behaviour

1. At t = 0, `fire_detections.csv` is read and three fire disks of radius
   500 m are merged into `phi` at simulated positions (3 km, 5 km),
   (10 km, 10 km), and (16 km, 14 km).
2. Fire spreads eastward driven by 5 m/s wind.
3. At t ≈ 300 s, the CSV is re-read and the same detections are re-applied.
   Since `suppress_if_burning = 1`, only unburned cells are affected.
4. The simulation runs to t = 700 s without error.

## Detection CSV format

```
# lon_deg,lat_deg,confidence_pct
-119.970000,37.045000,85
-119.900000,37.090000,75
-119.840000,37.126000,90
```

- Lines starting with `#` are comments.
- `confidence_pct` must be ≥ `satellite.confidence_threshold` (default 50 %)
  to be included.
- The coordinate projection uses the UTM zone and origin offsets specified
  in `inputs.i` (`satellite.utm_zone`, `satellite.prob_lo_easting_m`, etc.).

## File-mode vs real-time modes

In this test `satellite.source = file` is used so no network access is
required.  For real-time operation, change `source` to `viirs` (with a
valid NASA FIRMS `api_key`) or `goes` (public, no key required).

## Coordinate alignment

The inputs file sets:

```
satellite.prob_lo_easting_m  = 570000.0
satellite.prob_lo_northing_m = 4100000.0
```

These must match the UTM easting/northing corresponding to the simulation's
(`prob_lo_x`, `prob_lo_y`) = (0, 0) corner, so that

```
sim_x = UTM_easting  - prob_lo_easting_m
sim_y = UTM_northing - prob_lo_northing_m
```

maps detections from geographic space to the simulation domain.
