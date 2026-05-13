# WTR Rain-Wetting Regtest
# Tests: FARSITE .wtr weather file with non-zero precipitation records driving
#        the rain-wetting model (precipitation_moisture.H + fuel_moisture_scheduler.H).
#
# Setup: run create_wtr.py first (CMake fixture) to generate fire.wtr
#
# Checks exercised:
#   * load_wtr_weather() parses non-zero PRECIP_IN records correctly.
#   * get_precip_at_time() returns positive rates during hours 10:00–13:00.
#   * precipitation_moisture wetting activates when rain rate > threshold:
#     dead-fuel moisture climbs toward M_sat = 1.20 during the rain window.
#   * Moisture dries back after rain ends (hour 14:00+) under the diurnal
#     T/RH schedule.
#   * Solver completes without aborting when wtr_file, diurnal_moisture, and
#     precipitation wetting parameters are all active simultaneously.

# Grid & domain (500 m x 500 m, UTM Zone 11N – Southern California reference)
n_cell_x = 32
n_cell_y = 32
max_grid_size = 16
prob_lo_x = 330000.0
prob_lo_y = 3775000.0
prob_hi_x = 330500.0
prob_hi_y = 3775500.0

# Simulate 7 hours to span the dry → rainy → clearing transition
final_time = 25200.0
cfl = 0.5
plot_int = 20
reinit_int = -1

# Light initial wind; WTR schedule will override this each step
u_x = 2.7
u_y = 0.0

# Ignition at domain centre
source_type = sphere
center_x = 330250.0
center_y = 3775250.0
sphere_radius = 30.0

# Fuel: FM4 chaparral with dry initial moisture (pre-storm)
rothermel.fuel_model = FM4
rothermel.M_f = 0.08

# WTR weather file provides time-varying wind + T/RH + precipitation
wtr_file = fire.wtr
# Sim start aligns with the first WTR record (07:00 on 1 July)
wtr_start_year  = 2000
wtr_start_month = 7
wtr_start_day   = 1
wtr_start_hour  = 700

# Diurnal moisture model – driven by T/RH from the WTR file
diurnal_moisture.enable = 1
diurnal_moisture.T_min  = 21.1
diurnal_moisture.T_max  = 31.1
diurnal_moisture.RH_min = 38.0
diurnal_moisture.RH_max = 88.0

# Precipitation wetting: wets dead fuels when rain rate > threshold
# The WTR file supplies 0.15 in/hr = 3.81 mm/hr during hours 10:00–13:00,
# well above the default 0.25 mm/hr threshold.
diurnal_moisture.precip_rain_rate_mm_hr  = 3.81
diurnal_moisture.precip_threshold_mm_hr  = 0.25
diurnal_moisture.M_sat = 1.20

# Level-set propagation
propagation_method = levelset
