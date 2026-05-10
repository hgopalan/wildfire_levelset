# 2D Fire Spread Test
# Tests: 2D level-set advection, 2D FARSITE ellipse, terrain slope, Rothermel

# Grid & domain (unit square)
n_cell = 64
max_grid_size = 32
prob_lo_x = 0.0
prob_lo_y = 0.0
prob_hi_x = 1.0
prob_hi_y = 1.0

# Time & output
nsteps = 100
cfl = 0.5
plot_int = 20
reinit_int = 20

# 2-D wind vector
u_x = 0.25
u_y = 0.15

# Spherical ignition (off-centre to show asymmetric spread)
source_type = sphere
center_x = 0.3
center_y = 0.3
sphere_radius = 0.12

# Fuel: FM3 Tall Grass with upslope terrain
rothermel.fuel_model = FM3
rothermel.M_f = 0.08
rothermel.slope_x = 0.1    # 10% upslope in x-direction (tan of slope angle)

# FARSITE ellipse with Anderson L/W
farsite.use_anderson_LW = 1
farsite.phi_threshold = 0.1
farsite.fire_shape_model = richards

# FARSITE propagation
propagation_method = farsite
