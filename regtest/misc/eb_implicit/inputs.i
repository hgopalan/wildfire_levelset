# EB Implicit Function Test
# Tests: embedded-boundary implicit function as initial fire perimeter
# Supported EB types: plane, cylinder, sphere, ellipsoid
# This test uses an ellipsoid centred off-origin

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
reinit_int = 20

# Constant wind field
u_x = 0.25
u_y = 0.1

# EB ellipsoid initial condition
# phi = 0 on the ellipsoid surface; phi < 0 inside (burning region)
source_type = eb
eb_type = ellipsoid
eb_param1 = 0.3     # center_x
eb_param2 = 0.5     # center_y
eb_param3 = 0.5     # center_z
eb_param4 = 0.20    # radius_x
eb_param5 = 0.15    # radius_y
eb_param6 = 0.10    # radius_z

# Fuel defaults
rothermel.fuel_model = FM4
rothermel.M_f = 0.08

# Level-set propagation
propagation_method = levelset
