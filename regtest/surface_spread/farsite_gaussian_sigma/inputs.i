# FARSITE Gaussian Sigma Spread-Point Smoothing Test
# Tests: farsite.gaussian_sigma > 0 replaces single-cell phi stamping with
#        an SDF disk of radius sigma, giving a smooth fire-front boundary
#        analogous to fire_points_file / fire_gaussian_sigma at initialisation.

# Grid & domain (unit square, 20 m cells)
n_cell_x = 50
n_cell_y = 50
max_grid_size = 50
prob_lo_x = 0.0
prob_lo_y = 0.0
prob_hi_x = 1000.0
prob_hi_y = 1000.0

# Time & output
nsteps = 20
cfl = 0.5
plot_int = 20
reinit_int = -1

# Constant wind field (5 m/s in x-direction)
u_x = 5.0
u_y = 0.0

# Box (line fire) ignition at western edge
source_type = box
box_xmin = 50.0
box_xmax = 100.0
box_ymin = 300.0
box_ymax = 700.0
box_zmin = 0.0
box_zmax = 1000.0

# FARSITE ellipse model (Richards 1990) with Gaussian smoothing enabled.
# farsite.gaussian_sigma > 0: propagated fire-front points are stamped as
# SDF disks of this radius [m] rather than single cells.
farsite.use_anderson_LW = 0
farsite.length_to_width_ratio = 3.0
farsite.phi_threshold = 0.1
farsite.coeff_a = 1.0
farsite.coeff_b = 0.4
farsite.coeff_c = 0.2
farsite.fire_shape_model = richards
farsite.gaussian_sigma = 40.0

# Rothermel fuel (FM4 chaparral defaults)
rothermel.fuel_model = FM4
rothermel.M_f = 0.08

# FARSITE propagation
propagation_method = farsite
