# Scott & Reinhardt (2001) Full TI/CI Diagnostic Test
# Tests: full bisection-based Torching Index / Crowning Index computation
# and output of torching_index_kmh and crowning_index_kmh plotfile fields

# Grid & domain (600 m x 600 m, UTM Zone 11N Southern California reference)
n_cell = 64
max_grid_size = 32
prob_lo_x = 330000.0
prob_lo_y = 3775000.0
prob_hi_x = 330600.0
prob_hi_y = 3775600.0

# Short run: we just need the initial diagnostics
final_time = 300.0
cfl = 0.5
plot_int = 10
reinit_int = -1

# Moderate wind
u_x = 3.0
u_y = 0.0

# Small ignition at domain centre
source_type = sphere
center_x = 330300.0
center_y = 3775300.0
sphere_radius = 10.0

# Surface spread model
fire_spread_model = rothermel
rothermel.fuel_model = FM8
rothermel.M_f = 0.08

# Crown canopy parameters (needed for TI threshold I_o and active crown R'_SA)
crown.enable = 1
crown.CBH = 4.0      # canopy base height [m]
crown.CBD = 0.15     # canopy bulk density [kg/m3]
crown.FMC = 100.0    # foliar moisture content [%]

# Enable full Scott & Reinhardt bisection-based TI/CI
scott_reinhardt_full.enable = 1

# Level-set propagation
propagation_method = levelset
