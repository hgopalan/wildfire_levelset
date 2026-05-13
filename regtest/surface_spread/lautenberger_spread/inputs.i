# Lautenberger (2013) Physics-Based Spread Test
# Tests: Lautenberger physics-based fire spread model

# Grid & domain (600 m x 600 m, UTM Zone 11N Southern California reference)
n_cell = 64
max_grid_size = 32
prob_lo_x = 330000.0
prob_lo_y = 3775000.0
prob_hi_x = 330600.0
prob_hi_y = 3775600.0

# Simulation: 30 min chaparral spread
final_time = 1800.0
cfl = 0.5
plot_int = 20
reinit_int = -1

# Wind: 6 m/s eastward (21.6 km/h), typical Santa Ana-like conditions
u_x = 6.0
u_y = 0.0
u_z = 0.0

# Point ignition near domain centre
source_type = sphere
center_x = 330300.0
center_y = 3775300.0
center_z = 0.5
sphere_radius = 15.0

# Lautenberger (2013) physics-based fire spread model
fire_spread_model = lautenberger

# Rothermel fuel parameters (used by Lautenberger for fuel moisture)
rothermel.fuel_model = FM4
rothermel.M_f = 0.08

# Lautenberger model coefficients
lautenberger.A_L = 1.05e-5   # pre-factor [m^2/s]
lautenberger.B_L = 2.5       # moisture sensitivity [-]
lautenberger.C_L = 0.45      # wind correction [(m/s)^-1]
lautenberger.D_L = 0.50      # slope correction [-]

# Level-set propagation
propagation_method = levelset
