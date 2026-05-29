# ============================================================================
# Regression Test: Canadian FWI System
# ============================================================================
# Tests the Canadian Forest Fire Weather Index System implementation
# ============================================================================

# Domain
geometry.coord_sys   = 0
geometry.prob_lo     = 330000.0 3775000.0 0.0
geometry.prob_hi     = 330500.0 3775500.0 0.0
amr.n_cell           = 32 32 1

# Time
max_step             = 20
stop_time            = 180.0
levelset.cfl         = 0.3

# Fire spread model
fire_spread_model    = rothermel

# Simple ignition
ignition.enable      = 1
ignition.type        = point
ignition.x           = 330250.0
ignition.y           = 3775250.0
ignition.time        = 0.0

# Fuel model
rothermel.use_anderson13 = 1
rothermel.anderson13_code = 4

# Wind
wind.enable          = 1
wind.u_mag           = 3.0
wind.u_dir_deg       = 45.0

# Fuel moisture
fuel_moisture.d1     = 0.06
fuel_moisture.d10    = 0.08
fuel_moisture.d100   = 0.10

# Output
output.plot_int      = 20
output.plot_file     = plt
