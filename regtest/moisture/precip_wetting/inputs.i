# Precipitation Wetting Test
# Tests: rain-driven wetting of dead fuels on top of diurnal moisture cycle
# Dead fuel moisture rises toward saturation (M_sat) while rain rate > threshold

# Grid & domain (500 m x 500 m, UTM Zone 11N)
n_cell_x = 64
n_cell_y = 64
n_cell_z = 1
max_grid_size = 32
prob_lo_x = 330000.0
prob_lo_y = 3775000.0
prob_lo_z = 0.0
prob_hi_x = 330500.0
prob_hi_y = 3775500.0
prob_hi_z = 1.0

# Simulate 3 hours; moisture should rise toward saturation
final_time = 10800.0
cfl = 0.5
plot_int = 20
reinit_int = -1

# Light wind
u_x = 2.0
u_y = 0.0
u_z = 0.0

# Spherical ignition at domain centre
source_type = sphere
center_x = 330250.0
center_y = 3775250.0
center_z = 0.5
sphere_radius = 30.0

# Fuel: FM4 chaparral with moderate initial moisture
rothermel.fuel_model = FM4
rothermel.M_f = 0.12

# Diurnal moisture cycle (required by precipitation wetting)
diurnal_moisture.enable = 1
diurnal_moisture.T_min  = 12.0
diurnal_moisture.T_max  = 22.0
diurnal_moisture.RH_min = 50.0
diurnal_moisture.RH_max = 90.0

# Precipitation wetting: 2 mm/hr steady rain wets fuels toward saturation
diurnal_moisture.precip_rain_rate_mm_hr = 2.0
diurnal_moisture.precip_threshold_mm_hr = 0.25
diurnal_moisture.M_sat = 1.20

# Level-set propagation
propagation_method = levelset
