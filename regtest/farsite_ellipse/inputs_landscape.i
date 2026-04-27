# FARSITE with Landscape File Test
# Tests: Richards' (1990) elliptical fire expansion with landscape file
# This test uses a FARSITE landscape file (.lcp) with terrain data
# The landscape file contains elevation, slope, aspect, and fuel model data

# Grid & domain (100m x 100m to match landscape file)
n_cell_x = 50
n_cell_y = 50
max_grid_size = 32
prob_lo_x = 0.0
prob_lo_y = 0.0
prob_lo_z = 0.0
prob_hi_x = 100.0
prob_hi_y = 100.0
prob_hi_z = 1.0

# Time & output
nsteps = 50
cfl = 0.5
plot_int = 10

# Velocity (constant wind from west, typical Santa Ana conditions)
u_x = 5.0
u_y = 0.0
u_z = 0.0

# Initial source - line fire at western edge
source_type = box
box_xmin = 10.0
box_xmax = 15.0
box_ymin = 30.0
box_ymax = 70.0
box_zmin = 0.0
box_zmax = 1.0

# Reinitialization
reinit_int = -1

# FARSITE ellipse model (Richards 1990)
farsite.enable = 1
farsite.use_anderson_LW = 1
farsite.phi_threshold = 0.1

# Rothermel model - use landscape file for terrain
# When landscape file is specified, terrain_file slopes/elevation are ignored
rothermel.fuel_model = FM4
rothermel.M_f = 0.08
rothermel.landscape_file = socal_chaparral_landscape.lcp

# Level set control
skip_levelset = 1
