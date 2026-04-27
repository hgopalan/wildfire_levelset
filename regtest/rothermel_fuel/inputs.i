# Rothermel Fuel Model Test
# Tests: Different fuel types from NFFL database

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

# Velocity (moderate wind)
u_x = 0.2
u_y = 0.1
u_z = 0.0

# Initial source - sphere
source_type = sphere
sphere_center_x = 0.5
sphere_center_y = 0.5
sphere_center_z = 0.5
sphere_radius = 0.1

# Reinitialization
reinit_int = 20

# Rothermel model - Short Grass (FM1)
rothermel.fuel_model = FM1
rothermel.M_f = 0.06

# FARSITE disabled for this test
farsite.enable = 0

# Level set control
skip_levelset = 0
