# Elliptical SDF Initial Condition Test
# Tests: ellipsoidal signed-distance-function ignition shape

# Grid & domain (unit cube)
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
reinit_int = 20

# Constant wind field
u_x = 0.25
u_y = 0.0
u_z = 0.0

# Elliptical initial condition (off-centre to show advection)
source_type = ellipse
ellipse_center_x = 0.4
ellipse_center_y = 0.5
ellipse_center_z = 0.5
ellipse_radius_x = 0.25
ellipse_radius_y = 0.15
ellipse_radius_z = 0.10

# Level-set propagation
propagation_method = levelset
