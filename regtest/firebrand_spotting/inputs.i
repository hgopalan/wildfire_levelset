# Firebrand Spotting Model Test
# Tests: Probability-based firebrand spotting with FARSITE ellipse model

# Domain: UTM Zone 11N, Southern California (100m x 100m)
# Grid & domain
n_cell = 64
max_grid_size = 32
prob_lo_x = 330000.0
prob_lo_y = 3775000.0
prob_lo_z = 0.0
prob_hi_x = 330100.0
prob_hi_y = 3775100.0
prob_hi_z = 1.0

# Time & output
nsteps = 50
cfl = 0.5
plot_int = 5

# Velocity (moderate constant wind in x-direction)
u_x = 0.4
u_y = 0.0
u_z = 0.0

# Initial source - small sphere (point ignition)
source_type = sphere
center_x = 330030.0
center_y = 3775050.0
center_z = 0.5
sphere_radius = 10.0

# Reinitialization
reinit_int = -1

# FARSITE ellipse model (Richards 1990)
# Must be enabled for spotting to work
farsite.enable = 1
farsite.use_anderson_LW = 0
farsite.length_to_width_ratio = 3.0
farsite.phi_threshold = 0.1
farsite.coeff_a = 1.0
farsite.coeff_b = 0.5
farsite.coeff_c = 0.2

# Firebrand spotting model
spotting.enable = 1
spotting.P_base = 0.05
spotting.k_wind = 0.3
spotting.I_critical = 800.0
spotting.d_mean = 15.0
spotting.d_sigma = 40.0
spotting.d_lambda = 0.08    # decay rate [1/m]: was 8.0 in 0-1 domain → 8.0/100 = 0.08 /m (inverse scaling)
spotting.distance_model = lognormal
spotting.lateral_spread_angle = 20.0
spotting.spot_radius = 3.0
spotting.random_seed = 12345
spotting.check_interval = 3

# Rothermel model (use chaparral fuel type)
rothermel.fuel_model = FM4

# Level set control - skip to use FARSITE only
skip_levelset = 1
