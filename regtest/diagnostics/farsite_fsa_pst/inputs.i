# FARSITE FSA + PST Output Test
# Tests that the Fire Spread Atlas (.fsa) and Post-processing Statistics (.pst)
# files are written correctly during a FARSITE-mode simulation.
# Verifies:
#   - fire.fsa is created and non-empty after the run
#   - fire.pst is created, has a header row, and numeric data rows

# Grid & domain (400 m × 400 m)
n_cell_x = 64
n_cell_y = 64
max_grid_size = 32
prob_lo_x = 0.0
prob_lo_y = 0.0
prob_hi_x = 400.0
prob_hi_y = 400.0

# Time & output
nsteps = 30
cfl = 0.5
plot_int = 10
reinit_int = -1

# Wind
u_x = 4.0
u_y = 1.0

# Spherical ignition at domain centre
source_type = sphere
center_x = 200.0
center_y = 200.0
sphere_radius = 20.0

# Fuel: FM2 Timber (grass understory)
rothermel.fuel_model = FM2
rothermel.M_f = 0.08

# FARSITE propagation
propagation_method = farsite

# Fire statistics CSV
fire_stats_file = fire_stats.csv

# FARSITE Fire Spread Atlas output
fsa_file = fire.fsa

# FARSITE Post-processing Statistics output
pst_file = fire.pst
