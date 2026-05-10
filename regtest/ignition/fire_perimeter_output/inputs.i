# Fire Perimeter Output Test
# Tests: perimeter_NNNN.csv, perimeter_NNNN.geojson, and fire_stats.csv output
# A separate validate.py script checks the written files after the run.

# Grid & domain (500 m x 500 m, UTM Zone 11N)
n_cell_x = 64
n_cell_y = 64
max_grid_size = 32
prob_lo_x = 330000.0
prob_lo_y = 3775000.0
prob_hi_x = 330500.0
prob_hi_y = 3775500.0

# Time & output – write several plotfiles so validate.py finds perimeter files
nsteps = 50
cfl = 0.5
plot_int = 10
reinit_int = -1

# Moderate wind
u_x = 4.0
u_y = 0.0

# Spherical ignition at domain centre
source_type = sphere
center_x = 330250.0
center_y = 3775250.0
sphere_radius = 30.0

# Fuel: FM4 chaparral
rothermel.fuel_model = FM4
rothermel.M_f = 0.08

# Enable all perimeter and stats output (validate.py checks these)
write_perimeter_csv     = 1
write_perimeter_geojson = 1
fire_stats_file         = fire_stats.csv

# Level-set propagation
propagation_method = levelset
