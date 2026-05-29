# Integration Scenario C: Spotting Cascade with Terrain
# Purpose: Demonstrate multi-generation spotting cascade across complex terrain
# Location: Mountainous terrain with ridge-valley system
# Duration: 4-hour extreme fire weather scenario

# Grid & domain (4 km x 4 km, UTM Zone 11N - mountain terrain)
n_cell_x = 128
n_cell_y = 128
max_grid_size = 64
prob_lo_x = 340000.0
prob_lo_y = 4200000.0
prob_hi_x = 344000.0
prob_hi_y = 4204000.0

# Time & output (4-hour extreme fire behavior simulation)
final_time = 14400.0
cfl = 0.5
plot_int = 25
reinit_int = 10

# Wind: strong mountain wind with terrain channeling
u_x = 12.0
u_y = 3.0

# Spherical ignition at southwestern corner (windward ridge)
source_type = sphere
center_x = 340600.0
center_y = 4200800.0
sphere_radius = 60.0

# Terrain: Complex ridge-valley system
rothermel.terrain_file = ridge_valley_terrain.csv

# Fuel: FM4 chaparral (dense shrubland)
rothermel.fuel_model = FM4
rothermel.M_f = 0.05   # extremely dry conditions

# FARSITE elliptical expansion
farsite.enable = 1
farsite.use_anderson_lw = 1

# Albini spotting with ember accumulation
# High spotting activity due to extreme conditions
albini_spotting.enable = 1
albini_spotting.terminal_velocity = 1.2
albini_spotting.P_base = 0.025              # increased probability
albini_spotting.I_B_min = 80.0              # lower threshold (more spots)
albini_spotting.spot_radius = 15.0
albini_spotting.check_interval = 3          # frequent checks
albini_spotting.n_traj_steps = 150          # longer trajectories
albini_spotting.random_seed = 789
albini_spotting.max_loft_height = 1200.0    # tall plume
albini_spotting.wind_reduction = 0.5        # moderate wind reduction

# Ember accumulation with multi-generation ignition
ember_accumulation.enable = 1
ember_accumulation.k_decay = 0.001          # 16-min burnout (longer-lived embers)
ember_accumulation.rho_threshold = 8.0      # lower threshold (easier ignition)
ember_accumulation.k_ignition = 1.5         # faster ignition
ember_accumulation.dt_check = 45.0          # frequent ignition checks
ember_accumulation.use_moisture_damping = 1
ember_accumulation.allow_cascade = 1        # enable multi-generation spotting

# Spotting cascade tracking
# Records generations of spot fires
spotting_diagnostics.enable = 1
spotting_diagnostics.track_generations = 1
spotting_diagnostics.max_generations = 5

# Slope-dependent flame tilt for enhanced upslope spread
radiation_preheating.enable = 1
radiation_preheating.use_slope_tilt = 1
radiation_preheating.I_min = 120.0
radiation_preheating.view_factor_model = simple

flame_tilt.T_a = 308.0    # 35°C extreme heat
flame_tilt.T_f = 1100.0

# Periodic wind gusts (extreme fire weather)
turb_wind.enable = 1
turb_wind.gust_amplitude = 0.45              # 45% variation (extreme gusts)
turb_wind.gust_period = 600.0                # 10-min gust cycle
turb_wind.gust_phase = 0.0

# McArthur moisture scaling
precipitation_moisture.enable = 1
precipitation_moisture.use_mcarthur_scaling = 1
precipitation_moisture.initial_rain_mm = 0.0
precipitation_moisture.tau_1hr = 2400.0
precipitation_moisture.T_celsius = 35.0      # extreme temperature
precipitation_moisture.RH_percent = 15.0     # very low humidity

# Bulk fuel consumption
farsite.use_bulk_fuel_consumption = 1
fuel_consumption.enable = 1
fuel_consumption.tau_consumption = 900.0   # 15-min consumption

# Simulation date/time (September 15, 2024 - extreme fire weather)
sim_datetime.year = 2024
sim_datetime.month = 9
sim_datetime.day = 15
sim_datetime.hour = 12.0

# Fire behavior diagnostics
fire_stats_file = "fire_stats.csv"
write_perimeter_csv = 1
write_perimeter_geojson = 1
spotting_log_file = "spotting_events.csv"

# Propagation method
propagation_method = levelset
