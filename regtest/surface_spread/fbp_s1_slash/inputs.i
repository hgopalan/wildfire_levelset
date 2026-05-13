# Canadian FBP S1 Slash Spread Test
# Tests: FBP S1 (Jack or Lodgepole Pine Slash) rate of spread model

# Grid & domain (600 m x 600 m, UTM Zone 11N Southern California reference)
n_cell = 64
max_grid_size = 32
prob_lo_x = 330000.0
prob_lo_y = 3775000.0
prob_hi_x = 330600.0
prob_hi_y = 3775600.0

# Simulation: 30 min slash fire spread
final_time = 1800.0
cfl = 0.5
plot_int = 20
reinit_int = -1

# Wind: 4 m/s eastward (14.4 km/h)
u_x = 4.0
u_y = 0.0

# Point ignition near domain centre
source_type = sphere
center_x = 330300.0
center_y = 3775300.0
sphere_radius = 15.0

# Canadian FBP S1 (Jack or Lodgepole Pine Slash) spread model
fire_spread_model = fbp_s1
fbp.fuel_type = s1
fbp.moisture = 10.0   # dead fine fuel moisture [%]

# Level-set propagation
propagation_method = levelset
