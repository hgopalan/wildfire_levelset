# WindNinja Ridge/Canyon Wind-Terrain Test
# Tests: WindNinja-style empirical ridge speed-up and canyon channeling
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
velocity_file = gaussian_hill_wind.csv

# Spherical ignition at western side
source_type = sphere
center_x = 330250.0
center_y = 3775500.0
center_z = 0.5
sphere_radius = 40.0

# Terrain (Gaussian hill)
rothermel.terrain_file = gaussian_hill_terrain.csv
rothermel.fuel_model = FM4
rothermel.M_f = 0.08

# WindNinja ridge/canyon wind-terrain model
# Ridge  (wind upslope):   U_eff = U * (1 + k_ridge * tan_phi * alignment)
# Canyon (wind downslope): U_eff = U * (1 + k_canyon_wn * tan_phi * |alignment|)
wind_terrain.model = windninja_ridge_canyon
wind_terrain.k_ridge = 1.0
wind_terrain.k_canyon_wn = 0.5

# Level-set propagation
propagation_method = levelset
