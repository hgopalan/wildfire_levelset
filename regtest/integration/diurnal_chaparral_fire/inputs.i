# Integration Scenario A: Complete Diurnal Wildfire
# Purpose: Demonstrate realistic operational forecast scenario with multiple features
# Location: Southern California chaparral with terrain
# Duration: 24-hour simulation

# Grid & domain (2 km x 2 km, UTM Zone 11N - Southern California)
n_cell_x = 128
n_cell_y = 128
max_grid_size = 64
prob_lo_x = 330000.0
prob_lo_y = 3775000.0
prob_hi_x = 332000.0
prob_hi_y = 3777000.0

# Time & output (24-hour simulation with 30-min intervals)
nsteps = 400
cfl = 0.5
plot_int = 20
reinit_int = 10

# Base wind: moderate westerly (will be modulated by gusts)
u_x = 4.0
u_y = 0.0

# Spherical ignition at western edge (ignition at 10:00 local time)
source_type = sphere
center_x = 330400.0
center_y = 3776000.0
sphere_radius = 40.0

# Terrain: Gaussian hill for slope effects
rothermel.terrain_file = terrain_hill.csv

# Fuel: FM4 Southern California chaparral
rothermel.fuel_model = FM4
rothermel.M_f = 0.07   # initial moisture

# FARSITE elliptical expansion with Anderson dynamic L/W
farsite.enable = 1
farsite.use_anderson_lw = 1
farsite.use_bulk_fuel_consumption = 1

# FMC sinusoidal phenology (summer peak conditions)
fmc_phenology.model = sinusoidal
fmc_phenology.FMC_mean = 100.0
fmc_phenology.FMC_amplitude = 35.0
fmc_phenology.peak_doy = 150
fmc_phenology.start_doy = 200    # mid-July

# McArthur moisture scaling with diurnal temperature variation
precipitation_moisture.enable = 1
precipitation_moisture.use_mcarthur_scaling = 1
precipitation_moisture.initial_rain_mm = 0.0
precipitation_moisture.tau_1hr = 3600.0
precipitation_moisture.T_celsius = 28.0      # afternoon temperature
precipitation_moisture.RH_percent = 25.0     # low afternoon humidity

# Periodic wind gusts (afternoon thermal turbulence)
# Peak gusts at ~14:00-16:00 (4-6 hours into simulation)
turb_wind.enable = 1
turb_wind.gust_amplitude = 0.35              # 35% variation
turb_wind.gust_period = 900.0                # 15-min gust cycle
turb_wind.gust_phase = 0.0

# Albini spotting with ember accumulation
albini_spotting.enable = 1
albini_spotting.terminal_velocity = 1.0
albini_spotting.P_base = 0.015
albini_spotting.I_B_min = 100.0
albini_spotting.spot_radius = 12.0
albini_spotting.check_interval = 5
albini_spotting.n_traj_steps = 100
albini_spotting.random_seed = 456

ember_accumulation.enable = 1
ember_accumulation.k_decay = 0.00167         # 10-min burnout
ember_accumulation.rho_threshold = 10.0
ember_accumulation.k_ignition = 1.0
ember_accumulation.dt_check = 60.0
ember_accumulation.use_moisture_damping = 1

# Radiation preheating with slope-dependent flame tilt
radiation_preheating.enable = 1
radiation_preheating.use_slope_tilt = 1
radiation_preheating.I_min = 150.0
radiation_preheating.view_factor_model = simple

flame_tilt.T_a = 301.0    # 28°C ambient
flame_tilt.T_f = 1000.0

# Burn period gate (10:00–20:00 local time)
# Fire can only spread during daylight hours
burn_period.enable = 1
burn_period.start_hour = 10.0
burn_period.end_hour = 20.0
burn_period.sim_start_hour = 10.0

# Simulation date/time (July 19, 2024)
sim_datetime.year = 2024
sim_datetime.month = 7
sim_datetime.day = 19

# Perimeter output
write_perimeter_csv = 1
write_perimeter_geojson = 1
fire_stats_file = fire_stats.csv

# FARSITE propagation for elliptical spread
propagation_method = farsite
