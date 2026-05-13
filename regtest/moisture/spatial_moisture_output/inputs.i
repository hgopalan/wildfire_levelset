# Spatial Fuel Moisture Output Regtest
# Tests: feature #10 – spatial_moisture_mf (moisture_d1 through moisture_lw)
#        is written to every AMReX plotfile.
# Also tests: diurnal moisture + FMD moisture schedule interaction with plotfile output.

# Grid & domain (500 m x 500 m, UTM Zone 11N)
n_cell_x = 32
n_cell_y = 32
max_grid_size = 16
prob_lo_x = 330000.0
prob_lo_y = 3775000.0
prob_hi_x = 330500.0
prob_hi_y = 3775500.0

# Simulate 2 hours
final_time = 7200.0
cfl = 0.5
plot_int = 10
reinit_int = -1

# Moderate wind
u_x = 3.0
u_y = 0.0
u_z = 0.0

# Spherical ignition at domain centre
source_type = sphere
center_x = 330250.0
center_y = 3775250.0
center_z = 0.5
sphere_radius = 30.0

# Fuel: FM4 chaparral
rothermel.fuel_model = FM4
rothermel.M_f = 0.10

# Diurnal moisture (drives spatial_moisture_mf)
diurnal_moisture.enable = 1
diurnal_moisture.T_min  = 15.0
diurnal_moisture.T_max  = 30.0
diurnal_moisture.RH_min = 25.0
diurnal_moisture.RH_max = 80.0

# Level-set propagation
propagation_method = levelset
