# ============================================================================
# Regression Test: Two-Fuel Model Blending
# ============================================================================
# Tests two-fuel model blending for heterogeneous landscapes
# ============================================================================

# Domain
geometry.coord_sys   = 0
geometry.prob_lo     = 330000.0 3775000.0 0.0
geometry.prob_hi     = 330500.0 3775500.0 0.0
amr.n_cell           = 32 32 1

# Time
max_step             = 30
stop_time            = 240.0
levelset.cfl         = 0.3

# Fire spread model
fire_spread_model    = rothermel

# Simple ignition
ignition.enable      = 1
ignition.type        = point
ignition.x           = 330250.0
ignition.y           = 3775250.0
ignition.time        = 0.0

# Fuel model: grass
rothermel.use_anderson13 = 1
rothermel.anderson13_code = 2

# Wind
wind.enable          = 1
wind.u_mag           = 4.0
wind.u_dir_deg       = 0.0

# Fuel moisture
fuel_moisture.d1     = 0.05
fuel_moisture.d10    = 0.07
fuel_moisture.d100   = 0.09

# Output
output.plot_int      = 30
output.plot_file     = plt
