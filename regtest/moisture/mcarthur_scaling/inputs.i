# McArthur Moisture Scaling Regression Test
# Tests: McArthur-style temperature/RH-dependent moisture response time scaling
# Reference: McArthur (1967) Australian fire danger rating system

# Grid & domain (1000 m x 1000 m, UTM Zone 11N)
n_cell_x = 64
n_cell_y = 64
max_grid_size = 32
prob_lo_x = 330000.0
prob_lo_y = 3775000.0
prob_hi_x = 331000.0
prob_hi_y = 3776000.0

# Time & output (simulate 4 hours to observe moisture response)
nsteps = 200
cfl = 0.5
plot_int = 20
reinit_int = -1

# Wind: moderate eastward
u_x = 5.0
u_y = 0.0

# Spherical ignition at western portion
source_type = sphere
center_x = 330200.0
center_y = 3775500.0
sphere_radius = 30.0

# Fuel: FM4 chaparral
rothermel.fuel_model = FM4
rothermel.M_f = 0.08

# Enable precipitation moisture model with McArthur scaling
# McArthur scaling adjusts moisture response time based on temperature and RH
precipitation_moisture.enable = 1
precipitation_moisture.use_mcarthur_scaling = 1
precipitation_moisture.initial_rain_mm = 10.0    # recent 10mm rainfall
precipitation_moisture.tau_1hr = 3600.0          # base 1-hr response time [s]
precipitation_moisture.tau_10hr = 14400.0        # base 10-hr response time [s]
precipitation_moisture.tau_100hr = 86400.0       # base 100-hr response time [s]

# Set temperature and RH for McArthur scaling test
# McArthur formula adjusts drying rate based on T and RH
precipitation_moisture.T_celsius = 25.0
precipitation_moisture.RH_percent = 30.0

# Level-set propagation
propagation_method = levelset
