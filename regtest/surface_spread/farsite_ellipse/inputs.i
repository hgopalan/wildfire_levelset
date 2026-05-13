# FARSITE Elliptical Fire Spread Test
# Tests: Richards (1990) ellipse model with fixed L/W and explicit coefficients

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

# Box (line fire) ignition at western edge
source_type = box
box_xmin = 0.05
box_xmax = 0.12
box_ymin = 0.3
box_ymax = 0.7
box_zmin = 0.0
box_zmax = 1.0

# FARSITE ellipse model (Richards 1990)
farsite.use_anderson_LW = 0
farsite.length_to_width_ratio = 3.0
farsite.phi_threshold = 0.1
farsite.coeff_a = 1.0
farsite.coeff_b = 0.4
farsite.coeff_c = 0.2
farsite.fire_shape_model = richards

# Rothermel fuel (FM4 chaparral defaults)
rothermel.fuel_model = FM4
rothermel.M_f = 0.08

# FARSITE propagation
propagation_method = farsite
