# FARSITE Ellipse Model Test
# Tests: Richards' (1990) elliptical fire expansion with fixed coefficients

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
plot_int = 10

# Velocity (constant wind)
u_x = 0.3
u_y = 0.0
u_z = 0.0

# Initial source - box (line fire)
source_type = box
box_xmin = 0.4
box_xmax = 0.6
box_ymin = 0.4
box_ymax = 0.6
box_zmin = 0.4
box_zmax = 0.6

# Reinitialization
reinit_int = -1

# FARSITE ellipse model (Richards 1990)
farsite.enable = 1
farsite.use_anderson_LW = 0
farsite.length_to_width_ratio = 3.0
farsite.phi_threshold = 0.1
farsite.coeff_a = 1.0
farsite.coeff_b = 0.4
farsite.coeff_c = 0.2

# Rothermel model (use default chaparral)
rothermel.fuel_model = FM4

# Level set control
skip_levelset = 0
