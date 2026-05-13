# Compact Wind Direction Schedule Test
# Tests: time-varying wind direction and speed from a compact CSV schedule
# Data file: wind_schedule.csv (time_s, speed_ms, dir_deg)

# Grid & domain (500 m x 500 m, UTM Zone 11N)
n_cell_x = 64
n_cell_y = 64
max_grid_size = 32
prob_lo_x = 330000.0
prob_lo_y = 3775000.0
prob_hi_x = 330500.0
prob_hi_y = 3775500.0

# Simulate the full schedule (schedule covers ~2 h)
final_time = 7200.0
cfl = 0.5
plot_int = 20
reinit_int = -1

# Initial wind (overridden each step by schedule)
u_x = 3.0
u_y = 0.0

# Spherical ignition at domain centre
source_type = sphere
center_x = 330250.0
center_y = 3775250.0
sphere_radius = 30.0

# Fuel: FM4 chaparral
rothermel.fuel_model = FM4
rothermel.M_f = 0.08

# Compact wind direction schedule (overrides u_x/u_y each step)
wind_dir_schedule_file = wind_schedule.csv

# Level-set propagation
propagation_method = levelset
