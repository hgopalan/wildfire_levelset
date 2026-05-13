# Basic Level-Set Advection Test
# Tests: fundamental level-set advection with constant velocity

# Grid & domain (unit cube)
n_cell = 64
max_grid_size = 32
prob_lo_x = 0.0
prob_lo_y = 0.0
prob_hi_x = 1.0
prob_hi_y = 1.0

# Time & output
nsteps = 100
cfl = 0.5
plot_int = 10
reinit_int = 10

# Constant wind field
u_x = 0.25
u_y = 0.0

# Spherical ignition at domain centre
source_type = sphere
center_x = 0.5
center_y = 0.5
sphere_radius = 0.1

# Level-set propagation (no fire spread model — pure advection)
propagation_method = levelset
