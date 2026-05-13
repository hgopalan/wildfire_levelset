# Wilson (1988) Single Ellipse (Rear-Focus) Shape Test
# Tests: ellipse with fire origin at rear focus

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
reinit_int = -1

# Constant wind field
u_x = 0.25
u_y = 0.0

# Box (line fire) ignition
source_type = box
box_xmin = 0.05
box_xmax = 0.12
box_ymin = 0.3
box_ymax = 0.7
box_zmin = 0.0
box_zmax = 1.0

# Wilson (1988) single ellipse with origin at rear focus
farsite.use_anderson_LW = 1
farsite.phi_threshold = 0.1
farsite.fire_shape_model = wilson

rothermel.fuel_model = FM4
rothermel.M_f = 0.08

propagation_method = farsite
