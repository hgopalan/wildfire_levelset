# Cheney & Gould (1995 / 1998) Grassland Fire Spread Test
# Tests: Empirical grassland ROS model with moisture and curing corrections

# Domain: 1 km x 1 km flat grassland patch (UTM-like coordinates in metres)
n_cell = 64
max_grid_size = 32
prob_lo_x = 0.0
prob_lo_y = 0.0
prob_lo_z = 0.0
prob_hi_x = 1000.0
prob_hi_y = 1000.0
prob_hi_z = 1.0

# Time & output
final_time = 600.0
cfl = 0.5
plot_int = 20

# Wind: 10 m/s from west (= 36 km/h → well into the high-wind regime of the model)
u_x = 10.0
u_y = 0.0
u_z = 0.0

# Initial fire source: sphere at domain centre, radius 20 m
source_type = sphere
sphere_center_x = 500.0
sphere_center_y = 500.0
sphere_center_z = 0.5
sphere_radius = 20.0

# Reinitialization
reinit_int = 20

# Cheney & Gould grassland spread model
fire_spread_model = cheney_gould
cheney_gould.moisture = 8.0   # 8% dead fine fuel moisture
cheney_gould.curing   = 1.0   # fully cured grassland

# Level-set propagation
propagation_method = levelset
