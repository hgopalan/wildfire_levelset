# Crown Fire Initiation Test (Van Wagner 1977)
# Tests: Van Wagner's threshold model for crown fire initiation

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
plot_int = 10

# Velocity (moderate wind to drive fire)
u_x = 0.3
u_y = 0.1
u_z = 0.0

# Initial source - sphere
source_type = sphere
sphere_center_x = 0.5
sphere_center_y = 0.5
sphere_center_z = 0.5
sphere_radius = 0.1

# Reinitialization
reinit_int = 20

# Rothermel model - Timber with understory (FM10) - high intensity fuel
rothermel.fuel_model = FM10
rothermel.M_f = 0.08

# FARSITE configuration
farsite.enable = 1

# Level set control (use FARSITE mode)
skip_levelset = 1

# Van Wagner Crown Fire Initiation Model
crown.enable = 1
crown.CBH = 4.0              # Canopy base height [m] - typical conifer forest
crown.CBD = 0.15             # Canopy bulk density [kg/m³] - moderate density
crown.FMC = 100.0            # Foliar moisture content [%] - typical live fuel
crown.crown_fraction_weight = 1.0
crown.use_metric_units = 1   # Use metric units (m, kW/m)
