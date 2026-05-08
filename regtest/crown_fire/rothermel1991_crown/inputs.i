# Rothermel (1991) Active Crown Fire ROS Test
# Tests: R_crown = 3.34 × R_surface crown fire spread override

# Grid & domain (600 m x 600 m, UTM Zone 11N Southern California reference)
n_cell = 64
max_grid_size = 32
prob_lo_x = 330000.0
prob_lo_y = 3775000.0
prob_lo_z = 0.0
prob_hi_x = 330600.0
prob_hi_y = 3775600.0
prob_hi_z = 1.0

# Simulation: 30 min crown fire spread
final_time = 1800.0
cfl = 0.5
plot_int = 20
reinit_int = -1

# Wind: 7 m/s eastward (25.2 km/h), strong fire-weather conditions
u_x = 7.0
u_y = 0.0
u_z = 0.0

# Point ignition near domain centre
source_type = sphere
center_x = 330300.0
center_y = 3775300.0
center_z = 0.5
sphere_radius = 20.0

# Rothermel surface spread model (base for crown multiplier)
fire_spread_model = rothermel
rothermel.fuel_model = FM8   # compact timber litter – prone to crown fire
rothermel.M_f = 0.07

# Crown fire model: Rothermel (1991) multiplier
crown.enable = 1
crown.CBH = 3.0      # canopy base height [m]
crown.CBD = 0.20     # canopy bulk density [kg/m3]
crown.FMC = 90.0     # foliar moisture content [%]
crown.use_rothermel1991_crown = 1   # R_crown = 3.34 × R_surface

# Level-set propagation
propagation_method = levelset
