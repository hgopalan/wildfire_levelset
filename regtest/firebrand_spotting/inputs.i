# Firebrand Spotting Model Test
# Tests: Probability-based firebrand spotting with FARSITE ellipse model

# Grid & domain
n_cell = 64
max_grid_size = 32
prob_lo_x = 0.0
prob_lo_y = 0.0
prob_lo_z = 0.0
prob_hi_x = 1.0
prob_hi_y = 1.0
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
center_x = 0.3
center_y = 0.5
center_z = 0.5
sphere_radius = 0.1

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
spotting.d_mean = 0.15
spotting.d_sigma = 0.4
spotting.d_lambda = 8.0
spotting.distance_model = lognormal
spotting.lateral_spread_angle = 20.0
spotting.spot_radius = 0.03
spotting.random_seed = 12345
spotting.check_interval = 3

# Rothermel model (use chaparral fuel type)
rothermel.fuel_model = FM4

# Level set control - skip to use FARSITE only
skip_levelset = 1
