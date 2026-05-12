# WTR Diurnal Weather Regtest
# Tests: FARSITE .wtr weather file parser with time-varying wind and diurnal
#        temperature / RH cycle driving Nelson (2000) dead-fuel moisture updates.
#
# Setup: run create_wtr.py first (CMake fixture) to generate fire.wtr
#
# Checks exercised:
#   * load_wtr_weather() parses all 13 hourly records without aborting.
#   * WindDirSchedule is populated and wind speed/direction changes over time.
#   * Diurnal moisture conditioning: dead-fuel moisture rises during the cooler
#     morning records and falls back during the hot afternoon records.
#   * No precipitation records – moisture changes are purely thermally driven.

# Grid & domain (500 m x 500 m, UTM Zone 11N – Southern California reference)
n_cell_x = 32
n_cell_y = 32
max_grid_size = 16
prob_lo_x = 330000.0
prob_lo_y = 3775000.0
prob_hi_x = 330500.0
prob_hi_y = 3775500.0

# Simulate 6 hours to span several WTR records and observe moisture variation
final_time = 21600.0
cfl = 0.5
plot_int = 20
reinit_int = -1

# Ignition at domain centre
source_type = sphere
center_x = 330250.0
center_y = 3775250.0
sphere_radius = 30.0

# Fuel: FM4 chaparral with initial moderate moisture
rothermel.fuel_model = FM4
rothermel.M_f = 0.10

# WTR weather file provides time-varying wind + diurnal T/RH for moisture
wtr_file = fire.wtr
# Sim start aligns with the first WTR record (07:00 on 1 July)
wtr_start_year  = 2000
wtr_start_month = 7
wtr_start_day   = 1
wtr_start_hour  = 700

# Enable diurnal moisture model so that the WTR T/RH records drive the
# Nelson (2000) fuel moisture update each time step
diurnal_moisture.enable = 1
diurnal_moisture.T_min  = 22.2
diurnal_moisture.T_max  = 37.8
diurnal_moisture.RH_min = 15.0
diurnal_moisture.RH_max = 60.0

# Level-set propagation
propagation_method = levelset
