# Terrain Wind Test (2-D build only)
# Tests: Gaussian-hill terrain + spatially-varying wind, terrain-slope ROS correction
# Requires: cmake -DLEVELSET_DIM_2D=ON
# Data files: gaussian_hill_terrain.csv, gaussian_hill_wind.csv

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
cfl = 0.5
plot_int = 10
reinit_int = -1

# Spatially-varying wind from file (Gaussian hill speedup)
velocity_file = gaussian_hill_wind.csv

# Box ignition at western edge
source_type = box
box_xmin = 330010.0
box_xmax = 330060.0
box_ymin = 3775350.0
box_ymax = 3775650.0
box_zmin = 0.0
box_zmax = 1.0

# Terrain from file (Gaussian hill, 100 m peak)
rothermel.terrain_file = gaussian_hill_terrain.csv
rothermel.fuel_model = FM4
rothermel.M_f = 0.08

# FARSITE ellipse with Anderson L/W (terrain effects on shape)
farsite.use_anderson_LW = 1
farsite.phi_threshold = 0.1
farsite.fire_shape_model = richards

# FARSITE propagation
propagation_method = farsite
