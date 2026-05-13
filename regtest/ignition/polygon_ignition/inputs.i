# Polygon Ignition Test
# Tests: rasterisation of a closed CSV polygon as initial fire perimeter
# Data file: ignition_polygon.csv (~20 m x 15 m polygon, UTM Zone 11N)

# Grid & domain (100 m x 100 m, UTM Zone 11N)
n_cell_x = 64
n_cell_y = 64
max_grid_size = 32
prob_lo_x = 330000.0
prob_lo_y = 3775000.0
prob_hi_x = 330100.0
prob_hi_y = 3775100.0

# Time & output
nsteps = 80
cfl = 0.5
plot_int = 10
reinit_int = 20

# Wind: 3 m/s eastward
u_x = 3.0
u_y = 0.0
u_z = 0.0

# Polygon ignition from CSV file
source_type = polygon
fire_polygon_file = ignition_polygon.csv

# Fuel: FM4 chaparral
rothermel.fuel_model = FM4
rothermel.M_f = 0.08

# FARSITE ellipse propagation
farsite.use_anderson_LW = 1
farsite.phi_threshold = 0.1

propagation_method = farsite
