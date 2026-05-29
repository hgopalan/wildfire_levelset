# ============================================================================
# Regression Test: Duff Moisture and Smoldering Combustion
# ============================================================================
# Tests duff moisture tracking and smoldering combustion model
# ============================================================================

# Domain
geometry.coord_sys   = 0
geometry.prob_lo     = 330000.0 3775000.0 0.0
geometry.prob_hi     = 330500.0 3775500.0 0.0
amr.n_cell           = 32 32 1

# Time (longer for smoldering)
max_step             = 100
stop_time            = 1200.0
levelset.cfl         = 0.3

# Fire spread model
fire_spread_model    = rothermel

# Simple ignition
ignition.enable      = 1
ignition.type        = point
ignition.x           = 330250.0
ignition.y           = 3775250.0
ignition.time        = 0.0

# Fuel model: timber with duff
rothermel.use_anderson13 = 1
rothermel.anderson13_code = 8

# Light wind
wind.enable          = 1
wind.u_mag           = 2.0
wind.u_dir_deg       = 0.0

# Fuel moisture
fuel_moisture.d1     = 0.08
fuel_moisture.d10    = 0.10
fuel_moisture.d100   = 0.12

# Output
output.plot_int      = 20
output.plot_file     = plt
