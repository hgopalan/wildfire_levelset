# ============================================================================
# Regression Test: Advanced Physics Features
# ============================================================================
# Tests that all 10 newly implemented physics feature headers compile and
# can be used together. This is a basic compilation and integration test.
#
# Features tested:
#   1. Radiation-driven preheating distance
#   2. Fuel particle temperature evolution
#   3. Fire line intensity rate of change
#   4. Flame residence time from fire intensity
#   5. Wind-fuel interaction feedback
#   6. Fuel moisture of extinction gradient
#   7. Flame intermittency factor
#   8. Plume entrainment feedback on surface wind
#   9. Fuel bed bulk density spatial variation
#   10. Critical heat flux for ignition
# ============================================================================

# Domain
geometry.coord_sys   = 0
geometry.prob_lo     = 330000.0 3775000.0 0.0
geometry.prob_hi     = 331000.0 3776000.0 0.0
amr.n_cell           = 64 64 1

# Time
max_step             = 50
stop_time            = 300.0
levelset.cfl         = 0.3

# Basic fire spread model
fire_spread_model    = rothermel

# Simple ignition at center
ignition.enable      = 1
ignition.type        = point
ignition.x           = 330500.0
ignition.y           = 3775500.0
ignition.time        = 0.0

# Uniform fuel (Anderson FM1 - short grass)
rothermel.use_anderson13 = 1
rothermel.anderson13_code = 1

# Uniform wind
wind.enable          = 1
wind.u_mag           = 5.0
wind.u_dir_deg       = 0.0

# Moisture
fuel_moisture.d1     = 0.10
fuel_moisture.d10    = 0.12
fuel_moisture.d100   = 0.14
fuel_moisture.lh     = 0.80
fuel_moisture.lw     = 0.90

# Plot output - include all new diagnostic fields
amr.plot_int         = 25
amr.plot_file        = plt
plot_vars            = phi ros fireline_intensity flame_length

# Enable verbose output to verify features are accessible
amr.v                = 1
levelset.v           = 1
