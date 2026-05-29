# Python API Integration Demo: 2026 Enhancement Features
# Demonstrates comprehensive use of the Python wildfire solver API
# Tests: Multi-feature simulation with 2026 enhancements via Python bindings

# Grid & domain (2 km x 2 km domain for demonstration)
n_cell_x = 64
n_cell_y = 64
max_grid_size = 32
prob_lo_x = 330000.0
prob_lo_y = 3775000.0
prob_hi_x = 332000.0
prob_hi_y = 3777000.0

# Moderate simulation time for API demo
nsteps = 50
cfl = 0.5
plot_int = 10
reinit_int = 10

# Wind (moderate westerly)
u_x = 5.0
u_y = 1.0

# Spherical ignition
source_type = sphere
center_x = 330500.0
center_y = 3776000.0
sphere_radius = 50.0

# Terrain (simple slope for demo)
rothermel.terrain_file = simple_slope.csv

# Fuel model
rothermel.fuel_model = FM4
rothermel.M_f = 0.08

# FARSITE propagation
propagation_method = levelset
farsite.enable = 1
farsite.use_anderson_lw = 1

# McArthur moisture scaling
precipitation_moisture.enable = 1
precipitation_moisture.use_mcarthur_scaling = 1
precipitation_moisture.T_celsius = 28.0
precipitation_moisture.RH_percent = 25.0

# FMC phenology (sinusoidal)
fmc_phenology.model = sinusoidal
fmc_phenology.FMC_mean = 100.0
fmc_phenology.FMC_amplitude = 30.0
fmc_phenology.peak_doy = 150
fmc_phenology.start_doy = 200

# Ember accumulation
ember_accumulation.enable = 1
ember_accumulation.k_decay = 0.00167
ember_accumulation.rho_threshold = 10.0

# Albini spotting
albini_spotting.enable = 1
albini_spotting.P_base = 0.015
albini_spotting.I_B_min = 100.0

# Radiation preheating with slope tilt
radiation_preheating.enable = 1
radiation_preheating.use_slope_tilt = 1

# Periodic wind gusts
turb_wind.enable = 1
turb_wind.gust_amplitude = 0.30
turb_wind.gust_period = 900.0

# Fire stats
fire_stats_file = "fire_stats.csv"
