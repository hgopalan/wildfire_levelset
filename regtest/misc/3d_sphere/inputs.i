# 3-D Fire Spread Test (3-D build only)
# Tests: 3-D level-set advection, 3-D FARSITE ellipse, terrain slope, Rothermel
# Requires: default 3-D cmake build (cmake -S . -B build)

# Grid & domain (unit cube, reduced for 3-D efficiency)
n_cell = 48
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
plot_int = 20
reinit_int = 20

# 3-D wind vector
u_x = 0.25
u_y = 0.15
u_z = 0.05

# Spherical ignition (off-centre to show 3-D spread asymmetry)
source_type = sphere
center_x = 0.3
center_y = 0.3
center_z = 0.3
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
