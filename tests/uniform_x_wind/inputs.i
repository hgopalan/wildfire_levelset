# Grid & domain
n_cell_x = 128
n_cell_y = 64
n_cell_z = 64
prob_lo_x = -20.0  
prob_hi_x = 40.0
prob_lo_y = -20.0  
prob_hi_y = 20.0
prob_lo_z = -20.0  
prob_hi_z = 20.0
amr_enable_negative_phi_refine = 0
amr_regrid_int = 25
amr_refine_ratio = 2
amr_max_refinements = 1

# Time & output
nsteps = 200
cfl    = 0.25
plot_int = 5

# Velocity
u_x = 1.0
u_y = 0.00
u_z = 0.00

# Sphere SDF
source_type=sphere 
center_x = 0.0
center_y = 0.0
center_z = 0.0
sphere_radius   = 5.0


# ---------------- Level Set Reinitialization ----------------
reinit_int   = 5
