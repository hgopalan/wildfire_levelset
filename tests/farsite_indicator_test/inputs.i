# Test case for FARSITE indicator initialization
# When skip_levelset=1 and farsite.enable=1, phi is initialized as:
#   phi = 1 inside the fire region
#   phi = 0 outside the fire region
# This is different from normal mode where phi is a signed distance function

# Grid & domain
n_cell_x = 64
n_cell_y = 64
n_cell_z = 64
prob_lo_x = 0.0
prob_hi_x = 1000.0
prob_lo_y = 0.0
prob_hi_y = 1000.0
prob_lo_z = 0.0
prob_hi_z = 1000.0

# Time & output
nsteps = 100
cfl = 0.5
plot_int = 10

# Velocity field
u_x = 2.0
u_y = 0.0
u_z = 0.0

# Initial fire geometry (sphere)
source_type = sphere
center_x = 500.0
center_y = 500.0
center_z = 500.0
sphere_radius = 100.0

# FARSITE mode: skip level set evolution and use indicator initialization
skip_levelset = 1
farsite.enable = 1
farsite.length_to_width_ratio = 3.0
farsite.phi_threshold = 0.5
farsite.coeff_a = 1.0
farsite.coeff_b = 0.5
farsite.coeff_c = 0.2

# Rothermel parameters (using default chaparral fuel)
rothermel.M_f = 0.08

# Disable reinitialization (not needed when level set is skipped)
reinit_int = -1
