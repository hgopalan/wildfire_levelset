# External Terrain and Wind Test
# Tests: Gaussian hill terrain with spatially-varying wind field

# Grid & domain (1 km x 1 km)
n_cell_x = 100
n_cell_y = 100
prob_lo_x = 0.0
prob_lo_y = 0.0
prob_hi_x = 1000.0
prob_hi_y = 1000.0

# Time & output
nsteps = 200
cfl = 0.5
plot_int = 20

# Velocity from file
velocity_file = gaussian_hill_wind.csv

# Initial source - line fire at base of hill
source_type = box
box_xmin = 100.0
box_xmax = 120.0
box_ymin = 200.0
box_ymax = 800.0

# Reinitialization
reinit_int = 20

# Rothermel model with terrain from file
rothermel.fuel_model = FM4
rothermel.M_f = 0.08
rothermel.terrain_file = gaussian_hill_terrain.csv

# FARSITE model
farsite.enable = 1
farsite.use_anderson_LW = 1
farsite.phi_threshold = 0.1

# Level set control
skip_levelset = 0
