# Cheney & Gould (1995/1998) Grassland Fire Spread Test
# Tests: piecewise-linear wind formula for fully cured grassland
# Expected head-fire ROS at 10 m/s: ~ 2.69 m/s (wind > 5 km/h branch)

# Grid & domain (1000 m x 1000 m, UTM Zone 11N)
n_cell_x = 64
n_cell_y = 64
n_cell_z = 1
max_grid_size = 32
prob_lo_x = 330000.0
prob_lo_y = 3775000.0
prob_lo_z = 0.0
prob_hi_x = 331000.0
prob_hi_y = 3776000.0
prob_hi_z = 1.0

# Time & output
nsteps = 100
cfl = 0.5
plot_int = 20
reinit_int = 20

# Wind: 10 m/s eastward (36 km/h 10-m open wind)
u_x = 10.0
u_y = 0.0
u_z = 0.0

# Spherical ignition at domain centre
source_type = sphere
center_x = 330500.0
center_y = 3775500.0
center_z = 0.5
sphere_radius = 40.0

# Cheney & Gould grassland fire spread model
fire_spread_model = cheney_gould
cheney_gould.moisture = 8.0   # dead fine fuel moisture [%]
cheney_gould.curing   = 1.0   # fully cured grassland

# Level-set propagation
propagation_method = levelset
