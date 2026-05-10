# Level-Set Reinitialization Test
# Tests: aggressive reinitialization preserves |grad phi| = 1

# Grid & domain (unit cube)
n_cell = 64
max_grid_size = 32
prob_lo_x = 0.0
prob_lo_y = 0.0
prob_hi_x = 1.0
prob_hi_y = 1.0

# Time & output
nsteps = 100
cfl = 0.5
plot_int = 10

# Aggressive reinitialization every 5 steps
reinit_int = 5

# Constant wind field
u_x = 0.25
u_y = 0.0

# Spherical ignition
source_type = sphere
center_x = 0.5
center_y = 0.5
sphere_radius = 0.15

# Rothermel fuel defaults
rothermel.fuel_model = FM4
rothermel.M_f = 0.08

# Level-set propagation
propagation_method = levelset
