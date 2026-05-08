# Cruz, Alexander & Wakimoto (2005) Crown Fire Spread Test
# Tests: algebraic crown fire ROS for Sierra Nevada mixed-conifer stand
# R = 11.02 * U10^0.90 * CBD^0.19 * exp(-0.17 * MC10)

# Grid & domain (2000 m x 2000 m, UTM Zone 11N)
n_cell_x = 64
n_cell_y = 64
n_cell_z = 1
max_grid_size = 32
prob_lo_x = 330000.0
prob_lo_y = 3775000.0
prob_lo_z = 0.0
prob_hi_x = 332000.0
prob_hi_y = 3777000.0
prob_hi_z = 1.0

# Simulation: 30 min crown fire spread
final_time = 1800.0
cfl = 0.5
plot_int = 20
reinit_int = -1

# Wind: 6 m/s eastward (21.6 km/h), typical Diablo/Santa Ana conditions
u_x = 6.0
u_y = 0.0
u_z = 0.0

# Point ignition at domain centre
source_type = sphere
center_x = 331000.0
center_y = 3776000.0
center_z = 0.5
sphere_radius = 30.0

# Cruz et al. (2005) crown fire spread model
fire_spread_model = cruz_crown
cruz_crown.CBD  = 0.15    # canopy bulk density [kg/m3] – dense Sierra Nevada mixed-conifer
cruz_crown.MC10 = 8.0     # 10-h fuel moisture [%] – dry summer conditions

# Level-set propagation
propagation_method = levelset
