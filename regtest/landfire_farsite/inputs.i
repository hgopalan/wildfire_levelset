# FARSITE with Auto-Downloaded LANDFIRE Landscape
# Tests: FARSITE elliptical fire expansion with landscape from LANDFIRE data
# (or synthetic fallback).  The landscape file (landscape.lcp) is generated
# by create_landscape.py before this test runs.
#
# Domain: 300m x 300m to cover the synthetic 11x11 grid at 30m spacing

# Grid & domain
n_cell_x = 30
n_cell_y = 30
max_grid_size = 32
prob_lo_x = 0.0
prob_lo_y = 0.0
prob_lo_z = 0.0
prob_hi_x = 300.0
prob_hi_y = 300.0
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
box_xmin = 20.0
box_xmax = 30.0
box_ymin = 80.0
box_ymax = 220.0
box_zmin = 0.0
box_zmax = 1.0

# Reinitialization
reinit_int = -1

# FARSITE ellipse model (Richards 1990) with Anderson L/W ratio
farsite.enable = 1
farsite.use_anderson_LW = 1
farsite.phi_threshold = 0.1

# Rothermel model - landscape file provides terrain and fuel data
rothermel.fuel_model = FM4
rothermel.M_f = 0.08
rothermel.landscape_file = landscape.lcp

# Level set control
skip_levelset = 0
