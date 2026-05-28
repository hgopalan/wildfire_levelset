# Regression test for Python fire solver API - basic fire spread
# Tests initialization, time-stepping, and state extraction

# Grid setup
n_cell_x = 32
n_cell_y = 32
max_grid = 16

# Domain (400m x 400m)
plo_x = 329800.0
plo_y = 3774800.0
phi_x = 330200.0
phi_y = 3775200.0

# Constant wind
ux = 4.0
uy = 0.5

# Fuel model (Anderson fuel model 1 - short grass)
fuel_model = 1

# Simple circular ignition at center
ignition.type = circle
ignition.x0 = 330000.0
ignition.y0 = 3775000.0
ignition.radius = 25.0

# Time integration
cfl = 0.5
nsteps = 20
dt_plot = 50.0

# Propagation method
propagation_method = levelset

# Output control
plot_file = plt
verbose = 1
