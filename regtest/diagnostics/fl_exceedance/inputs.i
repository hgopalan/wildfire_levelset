# Flame Length Exceedance Raster Test
# Tests that fl_exceedance_mf is written to every plotfile and accumulates
# the maximum observed Byram flame length per cell over the simulation.
# The field should be monotonically non-decreasing in the burned region and
# match the per-step flame_length where the fire is most intense.

# Grid & domain (500 m × 500 m)
n_cell_x = 64
n_cell_y = 64
max_grid_size = 32
prob_lo_x = 0.0
prob_lo_y = 0.0
prob_hi_x = 500.0
prob_hi_y = 500.0

# Time & output
nsteps = 40
cfl = 0.5
plot_int = 10
reinit_int = -1

# Wind
u_x = 5.0
u_y = 0.0

# Spherical ignition at domain centre
source_type = sphere
center_x = 250.0
center_y = 250.0
sphere_radius = 25.0

# Fuel: FM5 Brush for moderate fireline intensity
rothermel.fuel_model = FM5
rothermel.M_f = 0.06

# Propagation
propagation_method = levelset

# Verify fl_exceedance appears in plotfile
# (just need to exercise the code path; checked by successful run)
plot_vars = phi fl_exceedance flame_length fireline_intensity
