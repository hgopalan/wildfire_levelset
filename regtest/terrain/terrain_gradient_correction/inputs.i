# Terrain Gradient Correction Test (2-D build only)
# Tests: terrain-corrected level set gradient (arc-length spacings)
# Steep Gaussian hill forces a large slope factor so the corrected
# gradient magnitude visibly differs from the flat-Cartesian result.
# Requires: cmake -DLEVELSET_DIM_2D=ON
# Data file: steep_gaussian_terrain.csv

# Grid & domain (1000 m x 1000 m, UTM Zone 11N)
n_cell_x = 64
n_cell_y = 64
max_grid_size = 32
prob_lo_x = 330000.0
prob_lo_y = 3775000.0
prob_hi_x = 331000.0
prob_hi_y = 3776000.0

# Time & output
nsteps = 100
cfl    = 0.5
plot_int = 25
reinit_int = 20

# Constant westerly wind (5 m/s)
wind_speed = 5.0
wind_angle = 270.0

# Point ignition at western edge of the domain
source_type = sphere
sphere_x    = 330100.0
sphere_y    = 3775500.0
sphere_z    = 0.0
sphere_r    = 50.0

# Steep terrain (Gaussian hill H=200m, σ=100m)
# At r=σ the slope magnitude is ~1.21 (≈50°), giving a 57% larger
# effective grid spacing so the surface gradient differs significantly
# from the flat-Cartesian value.
rothermel.terrain_file = steep_gaussian_terrain.csv
rothermel.fuel_model   = FM4
rothermel.M_f          = 0.08

# Level-set propagation (exercises godunov_norm_grad_phi with terrain slopes)
propagation_method = levelset
