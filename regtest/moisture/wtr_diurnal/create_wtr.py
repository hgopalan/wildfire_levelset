#!/usr/bin/env python3
"""
Setup step for wtr_diurnal regtest.

Creates a synthetic FARSITE .wtr file (fire.wtr) covering a single summer day
(July 1) with hourly records from 07:00 to 19:00.  The file exercises:

  * Time-varying wind speed and direction fed through the WTR parser
    (wtr_weather.H → schedule_setup.H).
  * Diurnal temperature / relative humidity cycle used to update dead-fuel
    moisture via the Nelson (2000) model (fuel_moisture_scheduler.H).
  * Zero precipitation in every record (no rain-driven wetting in this test;
    see wtr_rain_wetting for the precipitation case).

FARSITE .wtr format (whitespace-delimited):
  MONTH  DAY  HOUR  TEMP_F  RH_PCT  PRECIP_IN  WIND_SPEED_MPH  WIND_DIR_DEG

The diurnal pattern used:
  * Temperature: peaks at 14:00 (100 °F → 37.8 °C), cools toward morning/evening.
  * Relative humidity: minimum at 14:00 (15 %), higher at night/morning (50–65 %).
  * Wind: gentle morning (5 mph / 270°), strengthening afternoon (12 mph / 250°).
  * Precipitation: 0 in/hr throughout.
"""
from pathlib import Path

# MONTH  DAY  HOUR  TEMP_F  RH_%  PRECIP_IN  WSPD_MPH  WDIR_DEG
records = [
    (7, 1,  700,  72, 60, 0.00,  5, 270),
    (7, 1,  800,  76, 55, 0.00,  6, 265),
    (7, 1,  900,  82, 48, 0.00,  7, 260),
    (7, 1, 1000,  88, 38, 0.00,  8, 255),
    (7, 1, 1100,  94, 28, 0.00,  9, 255),
    (7, 1, 1200,  98, 20, 0.00, 10, 250),
    (7, 1, 1300, 100, 17, 0.00, 11, 250),
    (7, 1, 1400, 100, 15, 0.00, 12, 250),
    (7, 1, 1500,  98, 17, 0.00, 11, 255),
    (7, 1, 1600,  95, 22, 0.00, 10, 260),
    (7, 1, 1700,  90, 28, 0.00,  8, 265),
    (7, 1, 1800,  84, 35, 0.00,  6, 270),
    (7, 1, 1900,  78, 45, 0.00,  5, 275),
]

lines = ["# Synthetic FARSITE .wtr file – wtr_diurnal regtest\n",
         "# MONTH DAY HOUR TEMP_F RH_PCT PRECIP_IN WSPD_MPH WDIR_DEG\n"]
for m, d, h, T, rh, p, ws, wd in records:
    lines.append(f"{m:2d} {d:2d} {h:4d}  {T:3d}  {rh:3d}  {p:.2f}  {ws:2d}  {wd:3d}\n")

out = Path("fire.wtr")
out.write_text("".join(lines))
print(f"Created {out}  ({len(records)} hourly records)")
