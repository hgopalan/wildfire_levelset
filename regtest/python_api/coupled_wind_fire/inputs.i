# Regression test for coupled Python wind-fire simulation
# Tests 3D wind field updates from an external wind solver

# Grid setup
n_cell_x = 24
n_cell_y = 24
max_grid = 16

# Domain (300m x 300m)
plo_x = 329850.0
plo_y = 3774850.0
phi_x = 330150.0
phi_y = 3775150.0

# Initial wind (will be overridden by Python)
ux = 3.0
uy = 0.0

# Fuel model (Anderson fuel model 2 - timber/grass)
fuel_model = 2

# Circular ignition
ignition.type = circle
ignition.x0 = 330000.0
ignition.y0 = 3775000.0
ignition.radius = 20.0

# Time integration
cfl = 0.4
nsteps = 15
dt_plot = 30.0

# Propagation method
propagation_method = levelset

# Output control
plot_file = plt
verbose = 1
