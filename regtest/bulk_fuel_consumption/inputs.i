# Bulk Fuel Consumption Fraction Model Test
# Tests: Fuel consumption fraction model based on fire intensity

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
nsteps = 80
cfl = 0.5
plot_int = 10

# Velocity (moderate wind)
u_x = 0.25
u_y = 0.1
u_z = 0.0

# Initial source - sphere
source_type = sphere
sphere_center_x = 0.5
sphere_center_y = 0.5
sphere_center_z = 0.5
sphere_radius = 0.12

# Reinitialization
reinit_int = 20

# Rothermel model - Use chaparral (FM4) for moderate to high intensity
rothermel.fuel_model = FM4
rothermel.M_f = 0.07

# FARSITE configuration
farsite.enable = 1
farsite.use_anderson_LW = 1

# Enable bulk fuel consumption model
farsite.use_bulk_fuel_consumption = 1
farsite.tau_residence = 60.0        # Residence time [seconds]
farsite.f_consumed_min = 0.5        # Minimum consumption fraction (fast fires)
farsite.f_consumed_max = 0.9        # Maximum consumption fraction (intense fires)

# Level set control (use FARSITE mode)
skip_levelset = 1
