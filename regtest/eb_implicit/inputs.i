# EB Implicit Function Test
# Tests: Embedded boundary implicit function for initial condition

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

# Initial source - EB implicit function (ellipsoid)
source_type = eb
eb_type = ellipsoid
eb_param1 = 0.3    # center_x
eb_param2 = 0.5    # center_y
eb_param3 = 0.5    # center_z
eb_param4 = 0.20   # radius_x
eb_param5 = 0.15   # radius_y
eb_param6 = 0.10   # radius_z

# Reinitialization
reinit_int = 20
reinit_iters = 20
reinit_dtau = 0.5

# Level set control
skip_levelset = 0
farsite.enable = 0
