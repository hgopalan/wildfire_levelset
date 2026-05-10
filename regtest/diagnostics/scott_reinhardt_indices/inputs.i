# Scott & Reinhardt (2001) Torching/Crowning Index Diagnostics Test
# Tests: fire ecology diagnostics including TI and CI proxy ratios,
#        scorch height, probability of ignition, and tree mortality
# Reference: Scott & Reinhardt (2001) USDA Forest Service RMRS-RP-29

# Grid & domain (500 m x 500 m, UTM Zone 11N)
n_cell_x = 64
n_cell_y = 64
max_grid_size = 32
prob_lo_x = 330000.0
prob_lo_y = 3775000.0
prob_hi_x = 330500.0
prob_hi_y = 3775500.0

# Time & output
nsteps = 60
cfl = 0.5
plot_int = 10
reinit_int = -1

# Strong wind to generate high fireline intensity and trigger crown activity
u_x = 8.0
u_y = 0.0

# Spherical ignition at domain centre
source_type = sphere
center_x = 330250.0
center_y = 3775250.0
sphere_radius = 30.0

# Fuel: FM10 Timber Litter (higher intensity than grass fuels)
rothermel.fuel_model = FM10
rothermel.M_f = 0.07

# Crown parameters for ecology diagnostics (Van Wagner 1977)
# I_o = 0.010 * CBH * (460 + 25.9 * FMC) = 0.010 * 2.5 * (460 + 25.9*100) = ~76 kW/m
crown.CBH = 2.5       # canopy base height [m]
crown.CBD = 0.15      # canopy bulk density [kg/m3]
crown.FMC = 100.0     # foliar moisture content [%]

# Fire ecology diagnostics parameters
fire_ecology.T_a_C = 35.0           # hot summer day [°C]
fire_ecology.solar_heating_F = 30.0 # solar fuel-temperature increment [°F]
fire_ecology.tree_height = 12.0     # stand height [m] for mortality calc

# Level-set propagation
propagation_method = levelset
