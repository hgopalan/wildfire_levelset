# Smoke Plume-Rise Model Regression Test (Briggs 1965)
# Tests that the smoke_plume.enable = 1 option produces non-zero plume rise
# heights in fire-front cells and that the model integrates correctly with
# the plotfile output pipeline (plume_rise_m field is written).

# Grid & domain (1 km x 1 km)
n_cell_x = 64
n_cell_y = 64
max_grid_size = 32
prob_lo_x = 0.0
prob_lo_y = 0.0
prob_hi_x = 1000.0
prob_hi_y = 1000.0

# Time: short run, 20 steps
nsteps = 20
cfl = 0.5
plot_int = 10
reinit_int = -1

# Moderate wind (5 m/s from west)
u_x = 5.0
u_y = 0.0

# Circular ignition at domain centre
source_type = sphere
center_x = 500.0
center_y = 500.0
sphere_radius = 50.0

# Fuel: FM4 chaparral (moderate intensity, generates measurable plume rise)
rothermel.fuel_model = FM4
rothermel.M_f = 0.08

# ---- Smoke plume-rise model (Briggs 1965) ----
# Using default atmospheric conditions (T_a=303.15 K, rho_a=1.20 kg/m3)
smoke_plume.enable = 1
smoke_plume.T_a    = 303.15
smoke_plume.rho_a  = 1.20
smoke_plume.Cp_a   = 1005.0

# Level-set propagation
propagation_method = levelset
