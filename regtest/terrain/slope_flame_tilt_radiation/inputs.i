# Slope-Dependent Flame Tilt for Radiation Test
# Tests: Radiation preheating with slope-dependent flame tilt angle
# Feature: Enhanced upslope fire spread via flame tilt coupling

# Grid & domain (1000 m x 1000 m, UTM Zone 11N)
n_cell_x = 64
n_cell_y = 64
max_grid_size = 32
prob_lo_x = 330000.0
prob_lo_y = 3775000.0
prob_hi_x = 331000.0
prob_hi_y = 3776000.0

# Time & output
nsteps = 120
cfl = 0.5
plot_int = 12
reinit_int = -1

# Wind: moderate eastward
u_x = 4.0
u_y = 0.0

# Spherical ignition at base of slope
source_type = sphere
center_x = 330200.0
center_y = 3775500.0
sphere_radius = 25.0

# Fuel: FM4 chaparral (sensitive to radiation preheating)
rothermel.fuel_model = FM4
rothermel.M_f = 0.08

# Terrain: create simple slope file
# 30° upslope in +x direction
rothermel.terrain_file = slope_terrain.csv

# Radiation preheating with slope-dependent flame tilt
radiation_preheating.enable = 1
radiation_preheating.use_slope_tilt = 1     # key feature: slope enhances flame tilt
radiation_preheating.I_min = 100.0          # min intensity for preheating [kW/m]
radiation_preheating.view_factor_model = simple

# Flame tilt parameters
flame_tilt.T_a = 300.0     # ambient temperature [K]
flame_tilt.T_f = 1000.0    # flame temperature [K]

# Level-set propagation
propagation_method = levelset
