# FARSITE Wind Stream Simulation Test
# Tests: FARSITE-style empirical wind deflection using terrain curvature
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
nsteps = 80
cfl = 0.5
plot_int = 10
reinit_int = -1

# Spatially-varying wind from file
velocity_file = ../windninja_ridge_canyon/gaussian_hill_wind.csv

# Spherical ignition at western side
source_type = sphere
center_x = 330250.0
center_y = 3775500.0
sphere_radius = 40.0

# Terrain (Gaussian hill - reuse from windninja test)
rothermel.terrain_file = ../windninja_ridge_canyon/gaussian_hill_terrain.csv
rothermel.fuel_model = FM4
rothermel.M_f = 0.08

# FARSITE wind stream simulation (Option 8)
# Computes terrain curvature to identify ridges/valleys
# Modifies both wind speed and direction
wind_terrain.model = farsite_wind
wind_terrain.k_ridge_farsite = 1.5    # ridge speed-up coefficient
wind_terrain.k_shelter = 0.6           # lee-side sheltering coefficient
wind_terrain.k_valley = 0.8            # valley channeling coefficient
wind_terrain.k_deflection = 0.3        # wind direction deflection coefficient
wind_terrain.min_curvature = 0.0001    # minimum curvature threshold

# Level-set propagation
propagation_method = levelset
