# Time-Dependent Wind Test (2-D build only)
# Tests: temporally-varying wind field loaded from numbered CSV snapshots
# At t=0: loads time_wind.csv; at t=wind_time_spacing: loads time_wind_1.csv
# Requires: cmake -DLEVELSET_DIM_2D=ON
# Data files: time_wind.csv, time_wind_1.csv (11x11 grid, 0–1000 m)

# Grid & domain (1000 m x 1000 m)
n_cell_x = 64
n_cell_y = 64
n_cell_z = 1
max_grid_size = 32
prob_lo_x = 0.0
prob_lo_y = 0.0
prob_lo_z = 0.0
prob_hi_x = 1000.0
prob_hi_y = 1000.0
prob_hi_z = 1.0

# Simulate two wind-time intervals
final_time = 120.0
cfl = 0.5
plot_int = 10
reinit_int = -1

# Time-dependent wind: wind field switches at t=60 s
velocity_file = time_wind.csv
use_time_dependent_wind = 1
wind_time_spacing = 60.0

# Spherical ignition at domain centre
source_type = sphere
center_x = 200.0
center_y = 500.0
center_z = 0.5
sphere_radius = 40.0

# Fuel: FM4 chaparral
rothermel.fuel_model = FM4
rothermel.M_f = 0.08

# Level-set propagation
propagation_method = levelset
