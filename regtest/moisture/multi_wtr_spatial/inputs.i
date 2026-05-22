# Multi-Station Weather Spatial Interpolation Regtest
# Tests: Multiple weather stations with IDW interpolation of T/RH creating
#        spatially-varying fuel moisture content via apply_multi_wtr_TRH_to_spatial()
#
# Setup: run create_stations.py first (CMake fixture) to generate:
#        - stations.csv (3 stations with different T/RH)
#        - station1.wtr, station2.wtr, station3.wtr
#
# Checks exercised:
#   * load_multi_wtr_stations() parses station list and loads all .wtr files
#   * apply_multi_wtr_TRH_to_spatial() interpolates T and RH using IDW
#   * Spatial T/RH gradient creates varying EMC across the domain
#   * Fire spreads faster in hot/dry regions (SW) vs cool/wet regions (N)
#
# Expected behavior:
#   * SW corner (near station 1): T~40°C, RH~15%, low moisture, fast spread
#   * SE corner (near station 2): T~30°C, RH~30%, moderate moisture
#   * North (near station 3):     T~20°C, RH~60%, high moisture, slow spread

# Grid & domain (1000 m x 1000 m, UTM Zone 11N)
n_cell_x = 64
n_cell_y = 64
max_grid_size = 16
prob_lo_x = 330000.0
prob_lo_y = 3775000.0
prob_hi_x = 331000.0
prob_hi_y = 3776000.0

# Simulate 2 hours to observe spatial spread variation
final_time = 7200.0
cfl = 0.5
plot_int = 20
reinit_int = -1

# Ignition at domain centre (halfway between all stations)
source_type = sphere
center_x = 330500.0
center_y = 3775500.0
sphere_radius = 50.0

# Fuel: FM4 chaparral
rothermel.fuel_model = FM4
# Initial moisture will be overridden by spatially-varying EMC from multi_wtr
rothermel.M_f = 0.10

# Multi-station weather file provides spatially-varying wind, T, and RH
multi_wtr_file = stations.csv
multi_wtr_idw_power = 2.0

# Sim start aligns with the first WTR record (07:00 on 1 July)
wtr_start_year  = 2000
wtr_start_month = 7
wtr_start_day   = 1
wtr_start_hour  = 700

# Enable diurnal moisture model (required for spatial T/RH to affect EMC)
diurnal_moisture.enable = 1

# Optional: Enable solar radiation to see shade-adjusted EMC with spatial T/RH
# solar_radiation.enable = 1
# solar_radiation.latitude  = 34.0
# solar_radiation.longitude = -118.0
# solar_radiation.year  = 2000
# solar_radiation.month = 7
# solar_radiation.day   = 1
# solar_radiation.sim_start_hour = 7.0
# solar_radiation.timezone_offset = -8.0
# solar_radiation.solar_heating_C = 17.0

# Level-set propagation
propagation_method = levelset

# Output fields to verify spatial T/RH interpolation
plot_vars = phi arrival_time fireline_intensity spatial_moisture_d1 spatial_moisture_d10 spatial_moisture_d100
