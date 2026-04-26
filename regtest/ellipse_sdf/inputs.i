# Elliptical SDF Test
# Tests: Elliptical initial condition with signed distance function

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
u_x = 0.2
u_y = 0.15
u_z = 0.0

# Initial source - ellipse
source_type = ellipse
ellipse_center_x = 0.4
ellipse_center_y = 0.5
ellipse_center_z = 0.5
ellipse_radius_x = 0.25
ellipse_radius_y = 0.15
ellipse_radius_z = 0.10

# Reinitialization
reinit_int = 20
reinit_iters = 20
reinit_dtau = 0.5

# Level set control
skip_levelset = 0
farsite.enable = 0
