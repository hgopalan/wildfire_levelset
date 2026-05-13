# Post-fire Re-entry Spotting Fuel Adjustment Regression Test
# Tests that fuel_depletion.adjust_spotting_reentry = 1 causes firebrand spots
# landing in previously-burned fuel-depleted areas to be suppressed.
#
# This test runs a simulation with:
#   - Fuel depletion enabled (fuel_depletion.enable = 1)
#   - Re-entry spotting adjustment enabled
#   - Albini spotting enabled
# The test just verifies the simulation runs without error (smoke test).
# The validate.py script checks the log output for the expected messages.

# Grid & domain (500 m x 500 m)
n_cell_x = 64
n_cell_y = 64
max_grid_size = 32
prob_lo_x = 0.0
prob_lo_y = 0.0
prob_hi_x = 500.0
prob_hi_y = 500.0

# Time
nsteps = 30
cfl = 0.5
plot_int = 10
reinit_int = -1

# Wind: strong enough to drive spotting
u_x = 8.0
u_y = 0.0

# Ignition
source_type = sphere
center_x = 100.0
center_y = 250.0
sphere_radius = 25.0

# Fuel: FM4 chaparral (high-intensity, drives spotting)
rothermel.fuel_model = FM4
rothermel.M_f = 0.06

# ---- Fuel depletion: track per-cell fuel load ----
fuel_depletion.enable              = 1
fuel_depletion.tau_burnout         = 600.0
fuel_depletion.couple_to_ros       = 0
fuel_depletion.adjust_spotting_reentry   = 1
fuel_depletion.spotting_fuel_threshold   = 0.05

# ---- Albini spotting ----
albini_spotting.enable          = 1
albini_spotting.P_base          = 0.5
albini_spotting.P_catch         = 1.0
albini_spotting.I_B_min         = 10.0
albini_spotting.check_interval  = 5
albini_spotting.random_seed     = 42
albini_spotting.spot_radius     = 15.0

# Level-set propagation
propagation_method = levelset
