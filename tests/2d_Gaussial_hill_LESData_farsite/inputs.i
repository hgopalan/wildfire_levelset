# Grid & domain
n_cell_x = 250
n_cell_y = 250
prob_lo_x = 0
prob_hi_x = 5000.0
prob_lo_y = 0.0  
prob_hi_y = 5000.0

# Time & output
nsteps = 100
plot_int = 10

# Velocity
velocity_file = "turbulent_field_2d.csv"

# Fire 
source_type=box
box_xmin = 2500
box_xmax = 2550 
box_ymin = 1000
box_ymax = 4000

# FARSITE — required for spotting to work
farsite.enable = 1
farsite.use_anderson_LW = 1
farsite.phi_threshold = 0.1
farsite.coeff_a = 1.0
farsite.coeff_b = 0.5
farsite.coeff_c = 0.2
farsite.use_bulk_fuel_consumption = 1

# Firebrand spotting model
spotting.enable = 0
spotting.P_base = 0.50
spotting.k_wind = 1.0
spotting.I_critical = 100.0
spotting.d_mean = 500.0
spotting.d_sigma = 0.6
spotting.d_lambda = 0.005
spotting.distance_model = lognormal
spotting.lateral_spread_angle = 45.0
spotting.spot_radius = 50.0
spotting.random_seed = 12345
spotting.check_interval = 1

# Rothermel model (use chaparral fuel type)
rothermel.fuel_model = FM4
rothermel.terrain_file = gaussian_hill_topography.csv

# ---------------- Level Set Reinitialization ----------------
reinit_int   = -1
skip_levelset = 1
