# Basic Level-Set Advection Test
# Tests: Simple sphere advection with constant velocity

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
plot_int = 20

# Velocity (constant)
u_x = 0.25
u_y = 0.1
u_z = 0.0

# Initial source - sphere
source_type = sphere
sphere_center_x = 0.3
sphere_center_y = 0.5
sphere_center_z = 0.5
sphere_radius = 0.15

# Reinitialization
reinit_int = 20
reinit_iters = 20
reinit_dtau = 0.5

# Level set control
skip_levelset = 0
farsite.enable = 0
