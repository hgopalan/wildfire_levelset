# Cell Size Effects Correction Test - With Correction
# FARSITE propagation with cell size correction enabled
# Cell size: 20m, Reference size: 30m
# Correction exponent: 0.1
# Expected correction factor: (30/20)^0.1 ≈ 1.0414 (+4.14%)

# Grid & domain (400 m × 400 m, 20 m cells)
n_cell_x = 20
n_cell_y = 20
max_grid_size = 20
prob_lo_x = 0.0
prob_lo_y = 0.0
prob_hi_x = 400.0
prob_hi_y = 400.0

# Time & output
nsteps = 20
cfl = 0.5
plot_int = 5
reinit_int = -1

# Wind
u_x = 5.0
u_y = 1.0

# Spherical ignition at domain centre
source_type = sphere
center_x = 200.0
center_y = 200.0
sphere_radius = 20.0

# Fuel: FM2 Timber (grass understory)
rothermel.fuel_model = FM2
rothermel.M_f = 0.08

# FARSITE propagation
propagation_method = farsite

# Cell size correction ENABLED
cellsize.enable = 1
cellsize.dx_ref = 30.0
cellsize.correction_exponent = 0.1

# Fire statistics CSV
fire_stats_file = fire_corrected.csv
