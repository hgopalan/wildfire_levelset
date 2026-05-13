# Polyline Ignition Test
# Tests: rasterisation of a polyline (line fire) from a CSV vertex file
# Data file: ignition_line.csv (east-west line, UTM Zone 11N)

# Grid & domain (100 m x 50 m, UTM Zone 11N)
n_cell_x = 64
n_cell_y = 32
max_grid_size = 32
prob_lo_x = 330000.0
prob_lo_y = 3775000.0
prob_hi_x = 330100.0
prob_hi_y = 3775050.0

# Time & output
nsteps = 80
cfl = 0.5
plot_int = 10
reinit_int = 20

# Wind: 4 m/s northward
u_x = 0.0
u_y = 4.0
u_z = 0.0

# Polyline ignition from CSV file
source_type = polyline
fire_polygon_file = ignition_line.csv
polyline_width = 4.0    # half-width [m] burned on each side of line

# Fuel: FM4 chaparral
rothermel.fuel_model = FM4
rothermel.M_f = 0.08

# FARSITE ellipse propagation
farsite.use_anderson_LW = 1
farsite.phi_threshold = 0.1

propagation_method = farsite
