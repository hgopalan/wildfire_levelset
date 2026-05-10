# MTT (Minimum Travel Time) Propagation Test
# Tests: Dijkstra fast-marching pre-computes arrival times; phi = arrival - t
# Reference: Finney (2002) USDA Forest Service Research Paper RMRS-RP-XX

# Grid & domain (500 m x 500 m, UTM Zone 11N)
n_cell_x = 64
n_cell_y = 64
max_grid_size = 32
prob_lo_x = 330000.0
prob_lo_y = 3775000.0
prob_hi_x = 330500.0
prob_hi_y = 3775500.0

# Simulate 600 s
final_time = 600.0
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

# MTT propagation: pre-computes arrival_time field at startup
# phi(i,j,k) = arrival_time(i,j,k) - current_time at each step
propagation_method = mtt
