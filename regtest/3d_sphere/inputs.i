# 3D Fire Spread Test
# Tests: Full 3D fire spread with FARSITE and Rothermel

# Grid & domain
n_cell_x = 48
n_cell_y = 48
n_cell_z = 48
max_grid_size = 24
prob_lo_x = 0.0
prob_lo_y = 0.0
prob_lo_z = 0.0
prob_hi_x = 1.0
prob_hi_y = 1.0
prob_hi_z = 1.0

# Time & output
nsteps = 100
cfl = 0.5
plot_int = 20

# Velocity (3D wind)
u_x = 0.25
u_y = 0.15
u_z = 0.05

# Initial source - sphere
source_type = sphere
sphere_center_x = 0.3
sphere_center_y = 0.3
sphere_center_z = 0.3
sphere_radius = 0.12

# Reinitialization
reinit_int = 20
reinit_iters = 20
reinit_dtau = 0.5

# FARSITE model
farsite.enable = 1
farsite.use_anderson_LW = 1
farsite.phi_threshold = 0.1

# Rothermel model
rothermel.fuel_model = FM3
rothermel.M_f = 0.07
rothermel.slope_x = 0.1
rothermel.slope_y = 0.0

# Level set control
skip_levelset = 0
