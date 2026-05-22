#!/usr/bin/env python3
"""
Setup step for multi_wtr_spatial regtest.

Creates:
  1. A station list file (stations.csv) with 3 weather stations at different
     locations.
  2. Three separate .wtr files with spatially-varying temperature and RH,
     creating a gradient across the domain.

This test validates that:
  * apply_multi_wtr_TRH_to_spatial() correctly interpolates T and RH using IDW
  * Spatially-varying EMC is computed from the interpolated weather fields
  * The fire spread varies based on local moisture content

Station layout (UTM Zone 11N):
  Station 1: SW corner (330000, 3775000) – Hot & Dry  (T=40°C, RH=15%)
  Station 2: SE corner (331000, 3775000) – Moderate   (T=30°C, RH=30%)
  Station 3: North     (330500, 3776000) – Cool & Wet (T=20°C, RH=60%)

This creates a T/RH gradient from hot/dry southwest to cool/wet north.
"""
from pathlib import Path

# Station list: id, x_m, y_m, wtr_file
stations = [
    (1, 330000.0, 3775000.0, "station1.wtr"),  # SW: Hot & Dry
    (2, 331000.0, 3775000.0, "station2.wtr"),  # SE: Moderate
    (3, 330500.0, 3776000.0, "station3.wtr"),  # N:  Cool & Wet
]

# Write station list
with open("stations.csv", "w") as f:
    f.write("# station_id, x_m, y_m, wtr_file\n")
    for sid, x, y, wtr in stations:
        f.write(f"{sid}, {x:.1f}, {y:.1f}, {wtr}\n")

print(f"Created stations.csv with {len(stations)} weather stations")

# Create WTR files with different T/RH conditions
# All stations have same wind (10 mph from 270°) but different T/RH

wtr_template = """# Synthetic FARSITE .wtr file – Station {sid}
# MONTH DAY HOUR TEMP_F RH_PCT PRECIP_IN WSPD_MPH WDIR_DEG
7  1  700  {T1}  {RH1}  0.00  10  270
7  1  800  {T2}  {RH2}  0.00  10  270
7  1  900  {T3}  {RH3}  0.00  10  270
7  1 1000  {T4}  {RH4}  0.00  10  270
7  1 1100  {T5}  {RH5}  0.00  10  270
7  1 1200  {T6}  {RH6}  0.00  10  270
"""

# Station 1: Hot & Dry (SW corner)
# T: 40°C = 104°F, RH: 15%
s1_data = {
    "sid": 1,
    "T1": 104, "T2": 104, "T3": 104, "T4": 104, "T5": 104, "T6": 104,
    "RH1": 15, "RH2": 15, "RH3": 15, "RH4": 15, "RH5": 15, "RH6": 15,
}

# Station 2: Moderate (SE corner)
# T: 30°C = 86°F, RH: 30%
s2_data = {
    "sid": 2,
    "T1": 86, "T2": 86, "T3": 86, "T4": 86, "T5": 86, "T6": 86,
    "RH1": 30, "RH2": 30, "RH3": 30, "RH4": 30, "RH5": 30, "RH6": 30,
}

# Station 3: Cool & Wet (North)
# T: 20°C = 68°F, RH: 60%
s3_data = {
    "sid": 3,
    "T1": 68, "T2": 68, "T3": 68, "T4": 68, "T5": 68, "T6": 68,
    "RH1": 60, "RH2": 60, "RH3": 60, "RH4": 60, "RH5": 60, "RH6": 60,
}

for fname, data in [("station1.wtr", s1_data),
                     ("station2.wtr", s2_data),
                     ("station3.wtr", s3_data)]:
    with open(fname, "w") as f:
        f.write(wtr_template.format(**data))
    print(f"Created {fname}")

print("\nTest setup complete. Run wildfire_levelset with inputs.i")
