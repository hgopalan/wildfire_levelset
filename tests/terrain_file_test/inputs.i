# Grid & domain
n_cell_x = 64
n_cell_y = 64
prob_lo_x = 0
prob_hi_x = 1000.0
prob_lo_y = 0.0  
prob_hi_y = 1000.0

# Time & output
nsteps = 100
cfl    = 0.25
plot_int = 20

# Velocity
u_x = 2.0
u_y = 0.0
u_z = 0.0

# Sphere SDF
source_type=sphere 
center_x = 500.0
center_y = 500.0
sphere_radius = 50.0

# Level Set Reinitialization
reinit_int = 20

# Rothermel fire spread model with terrain file
rothermel.fuel_model = FM4
rothermel.M_f = 0.08
rothermel.terrain_file = simple_slope_terrain.txt

# FARSITE ellipse model
farsite.enable = 1
farsite.use_anderson_LW = 1
farsite.phi_threshold = 5.0
