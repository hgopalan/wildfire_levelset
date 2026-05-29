# ============================================================================
# Regression Test: Fire Behavior Diagnostics
# ============================================================================
# Tests the five new fire behavior diagnostic fields:
#   1. Byram convective number (Nc)
#   2. Flame tilt angle
#   3. Packing ratio (β and β/β_opt)
#   4. Flame front depth
#   5. McArthur Forest Fire Danger Index (FFDI)
#
# This test verifies that all new diagnostic fields are computed and output
# correctly to plotfiles without affecting fire spread behavior.
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

# Fuel model: Anderson FM4 (chaparral)
rothermel.use_anderson13 = 1
rothermel.anderson13_code = 4

# Uniform moderate wind (5 m/s from west)
wind.enable          = 1
wind.u_mag           = 5.0
wind.u_dir_deg       = 0.0

# Fuel moisture
fuel_moisture.d1     = 0.08
fuel_moisture.d10    = 0.10
fuel_moisture.d100   = 0.12
fuel_moisture.lh     = 0.90
fuel_moisture.lw     = 1.20

# Plot output - include all new diagnostic fields
amr.plot_int         = 25
amr.plot_file        = plt
plot_vars            = phi R fireline_intensity flame_length \
                       convective_number flame_tilt packing_ratio \
                       relative_packing flame_depth mcarthur_ffdi

# Enable verbose output to verify diagnostics are computed
amr.v                = 1
levelset.v           = 1
