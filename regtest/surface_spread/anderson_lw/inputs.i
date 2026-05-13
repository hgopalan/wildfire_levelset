# Anderson (1983) Dynamic L/W Ratio Test
# Tests: wind-speed-dependent ellipse elongation
# Formula: L/W = 0.936*exp(0.2566*U) + 0.461*exp(-0.1548*U) - 0.397

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

# Wind ~10 mph = 4.47 m/s → expected L/W ≈ 2.5
u_x = 4.47
u_y = 0.0
u_z = 0.0

# Spherical ignition
source_type = sphere
center_x = 0.5
center_y = 0.5
center_z = 0.5
sphere_radius = 0.1

# FARSITE with Anderson L/W
farsite.use_anderson_LW = 1
farsite.phi_threshold = 0.1
farsite.fire_shape_model = richards

# Rothermel fuel defaults
rothermel.fuel_model = FM4
rothermel.M_f = 0.08

# FARSITE propagation
propagation_method = farsite
