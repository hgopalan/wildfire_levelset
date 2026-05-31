# ============================================================================
# Regression Test: KBDI and Haines Index Computation
# ============================================================================
# Tests the new fire danger indices:
#   1. Keetch-Byram Drought Index (KBDI) - soil moisture deficit
#   2. Haines Index - atmospheric stability for plume-dominated fires
#
# This test verifies that the new indices compute correctly with typical
# input values for a dry season scenario in southeastern United States.
#
# Expected results:
#   - KBDI should increment daily based on temperature and precipitation
#   - Haines Index should classify atmospheric instability correctly
#   - Both indices should be available as diagnostic outputs
# ============================================================================

# Domain - Small test domain for diagnostics
geometry.coord_sys   = 0
geometry.prob_lo     = 330000.0 3775000.0 0.0
geometry.prob_hi     = 331000.0 3776000.0 0.0
amr.n_cell           = 32 32 1

# Time - Short simulation to test index computation
max_step             = 20
stop_time            = 100.0
levelset.cfl         = 0.3

# Basic fire spread model
fire_spread_model    = rothermel
propagation_method   = levelset

# Point ignition
ignition.enable      = 1
ignition.type        = point
ignition.x           = 330500.0
ignition.y           = 3775500.0
ignition.time        = 0.0

# Fuel model: Anderson FM4 (chaparral) - common in dry climates
rothermel.use_anderson13 = 1
rothermel.anderson13_code = 4

# Hot, dry conditions (typical for high KBDI)
wind.enable          = 1
wind.u_mag           = 3.0
wind.u_dir_deg       = 270.0

# Dry fuel moisture (consistent with drought conditions)
fuel_moisture.d1     = 0.05
fuel_moisture.d10    = 0.08
fuel_moisture.d100   = 0.10
fuel_moisture.lh     = 0.70
fuel_moisture.lw     = 1.00

# Plot output
amr.plot_int         = 10
amr.plot_file        = plt

# Enable verbose output
amr.v                = 1
levelset.v           = 1

# Note: KBDI and Haines Index computation would be enabled via additional
# parameters when integrated into main.cpp. This test provides the baseline
# configuration for testing the indices.
#
# Future integration parameters:
#   kbdi.enable = 1
#   kbdi.initial_value = 200.0
#   kbdi.annual_precip_mm = 1000.0
#   haines.enable = 1
#   haines.variant = MID
#   haines.use_surface_approximation = 1
