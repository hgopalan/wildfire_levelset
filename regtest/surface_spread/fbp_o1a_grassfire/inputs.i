# Canadian FBP O1a Grassfire Spread Test
# Tests: FBP O1a (matted grass) rate of spread model for open grassland fire

# Grid & domain (600 m x 600 m, UTM Zone 11N Southern California reference)
n_cell = 64
max_grid_size = 32
prob_lo_x = 330000.0
prob_lo_y = 3775000.0
prob_hi_x = 330600.0
prob_hi_y = 3775600.0

# Simulation: 20 min grassfire spread
final_time = 1200.0
cfl = 0.5
plot_int = 20
reinit_int = -1

# Wind: 5 m/s eastward (18 km/h), typical afternoon grassland wind
u_x = 5.0
u_y = 0.0
u_z = 0.0

# Point ignition near domain centre
source_type = sphere
center_x = 330300.0
center_y = 3775300.0
center_z = 0.5
sphere_radius = 15.0

# Canadian FBP O1a (matted grass) spread model
fire_spread_model = fbp_o1a
fbp.fuel_type = o1a
fbp.curing = 80.0     # percent curing [%] – late summer, heavily cured
fbp.moisture = 8.0    # dead fine fuel moisture [%]

# Level-set propagation
propagation_method = levelset
