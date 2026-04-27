# Grid & domain
n_cell_x = 250
n_cell_y = 250
prob_lo_x = 0
prob_hi_x = 5000.0
prob_lo_y = 0.0  
prob_hi_y = 5000.0

# Time & output
nsteps = 1000
cfl = 0.25
plot_int = 100

# Fire 
source_type=box
box_xmin = 100
box_xmax = 120 
box_ymin = 1000
box_ymax = 4000


# ---------------- Level Set Reinitialization ----------------
reinit_int   = -1
rothermel.terrain_file = gaussian_hill_topography.csv
skip_levelset = 0
farsite.enable = 1

