# Grid & domain
n_cell_x = 64
n_cell_y = 32
prob_lo_x = 0.0  
prob_hi_x = 40.0
prob_lo_y = 0.0  
prob_hi_y = 20.0

# Time & output
nsteps = 10
cfl = 0.25
plot_int = 5

# Velocity from file
velocity_file = velocity_data.txt

# Sphere SDF
source_type = sphere 
center_x = 10.0
center_y = 10.0
sphere_radius = 5.0

# Level Set Reinitialization
reinit_int = 5
