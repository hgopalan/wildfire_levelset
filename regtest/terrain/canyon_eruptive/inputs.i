# Canyon Eruptive Wind-Terrain Test
# Tests: Viegas' (2004) canyon eruptive acceleration model by dynamically coupling front orientation to terrain-induced draft
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
nsteps = 10
cfl = 0.5
plot_int = 10
reinit_int = -1

# Spatially-varying wind from file
velocity_file = gaussian_hill_wind.csv

# Spherical ignition at western side
source_type = sphere
center_x = 330250.0
center_y = 3775500.0
sphere_radius = 40.0

# Terrain (Gaussian hill)
rothermel.terrain_file = gaussian_hill_terrain.csv
rothermel.fuel_model = FM4
rothermel.M_f = 0.08

# Viegas eruptive model parameters
viegas.enable = 1
viegas.T_a = 300.0
viegas.T_f = 1000.0
viegas.tan_phi_c = 0.4

# Canyon Eruptive model
wind_terrain.model = canyon_eruptive

# Level-set propagation
propagation_method = levelset
