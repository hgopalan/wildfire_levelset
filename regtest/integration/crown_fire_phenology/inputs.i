# Integration Scenario B: Crown Fire with Phenology
# Purpose: Demonstrate active crown fire behavior modulated by seasonal FMC phenology
# Location: Sierra Nevada mixed-conifer forest with seasonal moisture patterns
# Duration: Spring greenup period (May-June transition)

# Grid & domain (3 km x 3 km, UTM Zone 11N - Sierra Nevada)
n_cell_x = 96
n_cell_y = 96
max_grid_size = 48
prob_lo_x = 325000.0
prob_lo_y = 4180000.0
prob_hi_x = 328000.0
prob_hi_y = 4183000.0

# Time & output (2-hour active crown fire simulation)
final_time = 7200.0
cfl = 0.5
plot_int = 20
reinit_int = 10

# Wind: moderate afternoon Santa Ana-like wind
u_x = 8.0
u_y = 2.0

# Spherical ignition at western edge
source_type = sphere
center_x = 325500.0
center_y = 4181500.0
sphere_radius = 50.0

# Terrain: moderate slope (10% grade eastward)
rothermel.terrain_file = slope_terrain.csv

# Surface fuel: FM10 Timber litter and understory
rothermel.fuel_model = FM10
rothermel.M_f = 0.06   # dry understory conditions

# Cruz crown fire model
fire_spread_model = cruz_crown
cruz_crown.CBD  = 0.18    # high canopy bulk density (dense mixed-conifer)
cruz_crown.MC10 = 6.0     # dry 10-h fuel moisture

# GDD-based FMC phenology (spring greenup)
# Simulates transition from dormant to active growth
fmc_phenology.model = gdd
fmc_phenology.FMC_dormant = 85.0
fmc_phenology.FMC_peak = 130.0
fmc_phenology.GDD_threshold = 100.0
fmc_phenology.GDD_current = 250.0     # mid-greenup (late May)
fmc_phenology.GDD_peak = 800.0
fmc_phenology.T_base = 5.0

# Van Wagner crown fire initiation
crown_fire.enable = 1
crown_fire.CBH = 3.0        # 3 m canopy base height
crown_fire.FMC = 95.0       # foliar moisture from phenology
crown_fire.use_phenology_fmc = 1

# Rothermel1991 crown fire spread
rothermel1991_crown.enable = 1
rothermel1991_crown.use_active_crown = 1

# McArthur moisture scaling
precipitation_moisture.enable = 1
precipitation_moisture.use_mcarthur_scaling = 1
precipitation_moisture.initial_rain_mm = 0.0
precipitation_moisture.tau_1hr = 4200.0
precipitation_moisture.T_celsius = 22.0      # spring temperature
precipitation_moisture.RH_percent = 30.0     # moderate spring humidity

# Radiation preheating
radiation_preheating.enable = 1
radiation_preheating.I_min = 200.0
radiation_preheating.view_factor_model = simple

flame_tilt.T_a = 295.0    # 22°C ambient
flame_tilt.T_f = 1200.0   # hotter crown fire flames

# Bulk fuel consumption
farsite.use_bulk_fuel_consumption = 1
fuel_consumption.enable = 1
fuel_consumption.tau_consumption = 600.0   # 10-min consumption time

# Simulation date/time (May 25, 2025 - spring greenup period)
sim_datetime.year = 2025
sim_datetime.month = 5
sim_datetime.day = 25
sim_datetime.hour = 14.0

# Fire behavior diagnostics
fire_stats_file = "fire_stats.csv"
write_perimeter_csv = 1
write_perimeter_geojson = 0

# Propagation method
propagation_method = levelset
