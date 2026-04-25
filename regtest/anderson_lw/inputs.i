# Anderson L/W Ratio Test
# Tests: Dynamic ellipse elongation based on wind speed (Anderson 1983)

# Grid & domain
n_cell = 64
max_grid_size = 32
prob_lo_x = 0.0
prob_lo_y = 0.0
prob_lo_z = 0.0
prob_hi_x = 1.0
prob_hi_y = 1.0
prob_hi_z = 1.0

# Time & output
nsteps = 100
cfl = 0.5
plot_int = 10

# Velocity (moderate wind - ~10 mph)
u_x = 0.4
u_y = 0.0
u_z = 0.0

# Initial source - box
source_type = box
box_xmin = 0.45
box_xmax = 0.55
box_ymin = 0.45
box_ymax = 0.55
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
