# Rothermel Fuel Model Test
# Tests: Rothermel (1972) fire spread with FM1 (Short Grass) fuel

# Grid & domain (unit cube)
n_cell = 64
max_grid_size = 32
prob_lo_x = 0.0
prob_lo_y = 0.0
prob_hi_x = 1.0
prob_hi_y = 1.0

# Time & output
nsteps = 100
cfl = 0.5
plot_int = 10
reinit_int = 20

# Moderate wind
u_x = 0.25
u_y = 0.1
u_z = 0.0

# Spherical ignition at domain centre
source_type = sphere
center_x = 0.5
center_y = 0.5
center_z = 0.5
sphere_radius = 0.1

# Fuel: FM1 Short Grass with low moisture (fire-prone conditions)
rothermel.fuel_model = FM1
rothermel.M_f = 0.06

# Level-set propagation
propagation_method = levelset
