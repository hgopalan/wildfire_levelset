# FARSITE Barrier Polygon / Firebreak Regression Test

## Purpose

Tests the **barrier polygon** feature (`barrier_files`): reading multiple CSV
files of barrier vertex coordinates, finding the nearest grid cells, and
extinguishing fire (setting phi > 0) in those cells every time step.

This models firebreaks, fuel breaks, roads, or any physical barrier without
using the AMReX Embedded Boundary (EB) framework.

## Configuration

| Parameter | Value |
|-----------|-------|
| Domain | 500 m × 500 m × 500 m |
| Grid | 64³ cells |
| Wind | 3 m/s eastward |
| Fuel | FM4 (Southern California chaparral) |
| Propagation | `levelset` |
| Simulation time | 600 s |
| Barrier 1 | `barrier_north.csv` – east–west line at y = 350 m |
| Barrier 2 | `barrier_east.csv`  – north–south stub at x = 350 m |

## Setup

Run the Python setup script to create the barrier CSV files:

```bash
python3 create_barriers.py
./build/levelset inputs.i
```

## Expected Behaviour

- Fire starts as a sphere at (100, 250) m and spreads predominantly eastward.
- The northern barrier line (y = 350 m) prevents fire from spreading into the
  northern half of the domain above y = 350 m.
- The eastern stub (x = 350 m, 200 ≤ y ≤ 300 m) creates a gap in the
  eastward path, visible as a notch in the fire perimeter.
- Both barriers appear as thin unburned strips in the `phi` plotfile.

## Input Format

Each barrier CSV file contains `X Y [Z]` coordinates (one per line).
Comment lines (starting with `#`) and empty lines are ignored.
Multiple files are specified as a space-separated list:

    barrier_files = barrier_north.csv barrier_east.csv

## Run Command

```bash
cd regtest/barrier_polygons
python3 create_barriers.py
../../build/levelset inputs.i
```
