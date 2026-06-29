# Fire solver inputs for two-way coupling test

# Grid configuration
n_cell_x = 32
n_cell_y = 32

# Domain bounds (must match wind solver)
prob_lo_x = 0.0
prob_lo_y = 0.0
prob_hi_x = 1000.0
prob_hi_y = 1000.0

# AMR configuration
max_grid = 32

# Ignition setup
source_type = sphere
center_x = 250.0
center_y = 250.0
sphere_radius = 50.0

# Time control
cfl = 0.5
nsteps = 20
max_time = 600.0

# Propagation
propagation_method = farsite

# Fuel model
rothermel.model_number = 1

# Output
plot_interval = 1000000
write_plotfile = 0
