# Balbi + Viegas + Heat Flux Test
# Tests: Balbi (2009) physical spread model, Viegas (2004) eruptive diagnostics,
#        and heat-flux-induced wind correction on Gaussian hill terrain
# Data files: gaussian_hill_terrain.csv, gaussian_hill_wind.csv

# Grid & domain (1000 m x 1000 m, UTM Zone 11N)
n_cell_x = 64
n_cell_y = 64
n_cell_z = 1
max_grid_size = 32
prob_lo_x = 330000.0
prob_lo_y = 3775000.0
prob_lo_z = 0.0
prob_hi_x = 331000.0
prob_hi_y = 3776000.0
prob_hi_z = 1.0

# Time & output
nsteps = 80
cfl = 0.5
plot_int = 10
reinit_int = -1

# Spatially-varying wind from file
velocity_file = gaussian_hill_wind.csv

# Spherical ignition near western slope
source_type = sphere
center_x = 330250.0
center_y = 3775500.0
center_z = 0.5
sphere_radius = 40.0

# Terrain (Gaussian hill)
rothermel.terrain_file = gaussian_hill_terrain.csv
rothermel.fuel_model = FM4
rothermel.M_f = 0.08

# Balbi (2009) physical fire spread model
fire_spread_model = balbi
balbi.T_a = 300.0
balbi.T_f = 1000.0
balbi.T_i = 600.0

# Viegas (2004) eruptive fire diagnostics
viegas.enable = 1
viegas.a_V = 1.83
viegas.tan_phi_c = 0.4

# Heat flux wind correction (plume buoyancy)
heat_flux.value = 1000.0
heat_flux.enable_upward = 1
heat_flux.k_upward = 1.0
heat_flux.ref_height = 10.0

# Level-set propagation
propagation_method = levelset
