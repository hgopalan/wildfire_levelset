# Grid & domain
n_cell_x = 250
n_cell_y = 250
prob_lo_x = 0
prob_hi_x = 5000.0
prob_lo_y = 0.0  
prob_hi_y = 5000.0

# Time & output
nsteps = 10
plot_int = 100

# Velocity
velocity_file = "turbulent_field_2d.csv"

# Fire 
source_type=box
box_xmin = 2500
box_xmax = 2550 
box_ymin = 1000
box_ymax = 4000

# Firebrand spotting model
spotting.enable = 1
spotting.P_base = 0.15          # increased from 0.05
spotting.k_wind = 0.5           # increased from 0.3
spotting.I_critical = 400.0     # decreased from 800.0
spotting.d_mean = 0.15
spotting.d_sigma = 0.4
spotting.d_lambda = 8.0
spotting.distance_model = lognormal
spotting.lateral_spread_angle = 30.0  # increased from 20.0
spotting.spot_radius = 25
spotting.random_seed = 12345
spotting.check_interval = 1     # decreased from 3

# Rothermel model (use chaparral fuel type)
rothermel.fuel_model = FM4

# ---------------- Level Set Reinitialization ----------------
reinit_int   = -1
farsite.phi_threshold = 0.0
farsite.use_anderson_LW = 1
rothermel.terrain_file = gaussian_hill_topography.csv
skip_levelset = 1
