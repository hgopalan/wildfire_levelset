# Grid & domain
n_cell_x = 500
n_cell_y = 500
prob_lo_x = 0
prob_hi_x = 5000.0
prob_lo_y = 0.0  
prob_hi_y = 5000.0


# Time & output
nsteps = 5000
cfl    = 0.25
plot_int = 5000

# Velocity
u_x = 1.0
u_y = 0.00
u_z = 0.00
velocity_file = "turbulent_field_2d.csv"

# Sphere SDF
source_type=sphere 
center_x = 2500.0
center_y = 2500.0
sphere_radius   = 50.0


# ---------------- Level Set Reinitialization ----------------
reinit_int   = -1
farsite.phi_threshold = 5.0
farsite.use_anderson_LW = 1
rothermel.terrain_file = gaussian_hill_topography.csv
