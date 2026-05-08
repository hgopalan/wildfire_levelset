# LANDFIRE / FARSITE Landscape Test
# Tests: FARSITE fire spread with a binary LCP landscape file
# Setup: run create_landscape.py first (CMake fixture)
#   Downloads LANDFIRE data or falls back to a synthetic 11x11 chaparral landscape
#   Output: landscape.lcp (300 m x 300 m, 30 m grid, UTM Zone 11N)

# Grid & domain (300 m x 300 m, UTM Zone 11N, matches synthetic landscape)
n_cell_x = 64
n_cell_y = 64
n_cell_z = 1
max_grid_size = 32
prob_lo_x = 330000.0
prob_lo_y = 3775000.0
prob_lo_z = 0.0
prob_hi_x = 330300.0
prob_hi_y = 3775300.0
prob_hi_z = 1.0

# Time & output
nsteps = 80
cfl = 0.5
plot_int = 10
reinit_int = -1

# Wind: 5 m/s eastward (typical Santa Ana conditions)
u_x = 5.0
u_y = 0.0
u_z = 0.0

# Box ignition at western edge
source_type = box
box_xmin = 330010.0
box_xmax = 330040.0
box_ymin = 3775100.0
box_ymax = 3775200.0
box_zmin = 0.0
box_zmax = 1.0

# Landscape file provides per-cell elevation, slope, aspect, and fuel model
rothermel.landscape_file = landscape.lcp
rothermel.landscape_fuel_type = 13   # Anderson FBFM13 system
rothermel.M_f = 0.08

# FARSITE ellipse with Anderson L/W (spatially-varying L/W from landscape)
farsite.use_anderson_LW = 1
farsite.phi_threshold = 0.1
farsite.fire_shape_model = richards

# FARSITE propagation
propagation_method = farsite
