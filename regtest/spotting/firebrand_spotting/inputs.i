# Firebrand Spotting Test
# Tests: stochastic firebrand spotting model (probability-based)

# Grid & domain (500 m x 500 m, UTM Zone 11N)
n_cell_x = 64
n_cell_y = 64
max_grid_size = 32
prob_lo_x = 330000.0
prob_lo_y = 3775000.0
prob_hi_x = 330500.0
prob_hi_y = 3775500.0

# Time & output
nsteps = 100
cfl = 0.5
plot_int = 10
reinit_int = -1

# Wind: 5 m/s eastward (drives spotting downwind)
u_x = 5.0
u_y = 0.0

# Spherical ignition at western portion of domain
source_type = sphere
center_x = 330100.0
center_y = 3775250.0
sphere_radius = 25.0

# Fuel: FM4 chaparral
rothermel.fuel_model = FM4
rothermel.M_f = 0.07

# Stochastic firebrand spotting model
spotting.enable = 1
spotting.P_base = 0.02
spotting.k_wind = 0.3
spotting.I_critical = 800.0
spotting.d_mean = 30.0
spotting.d_sigma = 0.5
spotting.distance_model = lognormal
spotting.lateral_spread_angle = 15.0
spotting.spot_radius = 10.0
spotting.check_interval = 5
spotting.random_seed = 42

# Level-set propagation
propagation_method = levelset
