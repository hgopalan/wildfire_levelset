# Anderson L/W Ratio Test
# Tests: Dynamic ellipse elongation based on wind speed (Anderson 1983)

# Domain: UTM Zone 11N, Santa Monica Mountains, Southern California (100m x 100m)
# Grid & domain
n_cell_x = 64    # split from n_cell=64 to allow non-square UTM domain in future
n_cell_y = 64
max_grid_size = 32
prob_lo_x = 330000.0
prob_lo_y = 3775000.0
prob_lo_z = 0.0
prob_hi_x = 330100.0
prob_hi_y = 3775100.0
prob_hi_z = 1.0

# Time & output
final_time = 120.0
cfl = 0.5
plot_int = 10

# Velocity (moderate wind - ~10 mph)
u_x = 0.4
u_y = 0.0
u_z = 0.0

# Initial source - box
source_type = box
box_xmin = 330045.0
box_xmax = 330055.0
box_ymin = 3775045.0
box_ymax = 3775055.0
box_zmin = 0.45
box_zmax = 0.55

# Reinitialization
reinit_int = -1

# FARSITE with Anderson L/W
farsite.enable = 1
farsite.use_anderson_LW = 1
farsite.phi_threshold = 0.1

# Rothermel model
rothermel.fuel_model = FM4
rothermel.M_f = 0.08

# Level set control
skip_levelset = 0
