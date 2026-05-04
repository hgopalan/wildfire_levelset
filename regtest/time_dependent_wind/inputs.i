# Time-Dependent Wind Field Test
# Tests: Level-set fire spread with temporally-varying wind field
# Wind rotates from east (0 s) to northeast (3600 s) via linear interpolation

# Grid & domain (1 km x 1 km, flat terrain)
n_cell_x = 50
n_cell_y = 50
prob_lo_x = 0.0
prob_lo_y = 0.0
prob_hi_x = 1000.0
prob_hi_y = 1000.0

# Time & output
final_time = 750.0
cfl = 0.5
plot_int = 10

# Time-dependent wind field
# Snapshots: time_wind.csv (t=0) and time_wind_1.csv (t=3600 s)
velocity_file = time_wind.csv
use_time_dependent_wind = 1
wind_time_spacing = 3600.0

# Initial source - line fire at western edge
source_type = box
box_xmin = 100.0
box_xmax = 130.0
box_ymin = 400.0
box_ymax = 600.0

# Reinitialization
reinit_int = 10

# Rothermel model with FM4 chaparral fuel
rothermel.fuel_model = FM4
rothermel.M_f = 0.08

# FARSITE ellipse parameters (used when propagation_method = farsite)
farsite.use_anderson_LW = 1
farsite.phi_threshold = 0.1

# Level set control
propagation_method = levelset
