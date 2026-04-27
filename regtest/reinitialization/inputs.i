# Level-Set Reinitialization Test
# Tests: Periodic reinitialization to maintain signed distance function

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
nsteps = 200
cfl = 0.5
plot_int = 10

# Velocity (diagonal)
u_x = 0.3
u_y = 0.2
u_z = 0.1

# Initial source - sphere
source_type = sphere
sphere_center_x = 0.25
sphere_center_y = 0.25
sphere_center_z = 0.25
sphere_radius = 0.15

# Reinitialization (aggressive)
reinit_int = 5
reinit_iters = 30
reinit_dtau = 0.3

# FARSITE disabled
farsite.enable = 0

# Level set control
skip_levelset = 0
