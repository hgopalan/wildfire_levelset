#!/usr/bin/env python3
"""
Setup step for wtr_rain_wetting regtest.

Creates a synthetic FARSITE .wtr file (fire.wtr) that begins dry, then
transitions to a period of steady precipitation (rain rate > threshold),
exercising the rain-driven fuel-moisture wetting model (precipitation_moisture.H
coupled with fuel_moisture_scheduler.H).

FARSITE .wtr format (whitespace-delimited):
  MONTH  DAY  HOUR  TEMP_F  RH_PCT  PRECIP_IN  WIND_SPEED_MPH  WIND_DIR_DEG

Weather scenario used:
  * Hours 0–2  (07:00–09:00):  clear and dry – wind 6 mph, 0 precip.
  * Hours 3–6  (10:00–13:00):  steady rain 0.15 in/hr (≈ 3.8 mm/hr)
    which exceeds the default 0.25 mm/hr wetting threshold; moisture
    climbs toward the saturation ceiling (M_sat = 1.20).
  * Hour  7    (14:00):        rain ends, clear again; moisture starts
    drying back toward equilibrium under the diurnal T/RH schedule.

The test verifies that:
  1. load_wtr_weather() correctly parses non-zero PRECIP_IN records.
  2. get_precip_at_time() returns positive rates during the rain window.
  3. The solver does not abort when wtr_file is combined with
     diurnal_moisture.enable and precipitation wetting parameters.
"""
from pathlib import Path

# MONTH  DAY  HOUR  TEMP_F  RH_%  PRECIP_IN  WSPD_MPH  WDIR_DEG
records = [
    # --- dry pre-storm period ---
    (7, 1,  700,  80, 45, 0.00, 6, 270),
    (7, 1,  800,  85, 40, 0.00, 6, 270),
    (7, 1,  900,  88, 38, 0.00, 6, 265),
    # --- precipitation onset: 0.15 in/hr = 3.81 mm/hr ---
    (7, 1, 1000,  75, 70, 0.15, 4, 240),
    (7, 1, 1100,  72, 80, 0.15, 4, 235),
    (7, 1, 1200,  70, 85, 0.15, 3, 230),
    (7, 1, 1300,  70, 88, 0.15, 3, 230),
    # --- rain ends, start drying ---
    (7, 1, 1400,  78, 60, 0.00, 7, 265),
]

lines = ["# Synthetic FARSITE .wtr file – wtr_rain_wetting regtest\n",
         "# MONTH DAY HOUR TEMP_F RH_PCT PRECIP_IN WSPD_MPH WDIR_DEG\n"]
for m, d, h, T, rh, p, ws, wd in records:
    lines.append(f"{m:2d} {d:2d} {h:4d}  {T:3d}  {rh:3d}  {p:.2f}  {ws:2d}  {wd:3d}\n")

out = Path("fire.wtr")
out.write_text("".join(lines))
print(f"Created {out}  ({len(records)} hourly records, "
      f"{sum(1 for *_, p, _, _ in records if p > 0)} rainy hours)")
