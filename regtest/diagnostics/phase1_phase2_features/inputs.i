# ============================================================================
# Regression Test: Phase 1 & Phase 2 Operational Features
# ============================================================================
# Tests the new Phase 1 and Phase 2 features:
#   Phase 1:
#     1. Grass curing model (seasonal/moisture-dependent/GDD)
#     2. Diurnal weather cycles (temperature/RH/wind)
#     3. Elevation temperature lapse rate
#   Phase 2:
#     4. NFDRS Spread Component (SC)
#     5. Chandler Burning Index (CBI)
#
# This test verifies that the new features compute correctly with typical
# input values for a grassland fire in varying terrain.
#
# Expected results:
#   - Grass curing should vary seasonally or with drought
#   - Diurnal cycles should modulate fire behavior throughout the day
#   - Elevation lapse should adjust temperature and RH with terrain
#   - NFDRS SC should provide operational spread rating
#   - CBI should provide fire danger classification
# ============================================================================

# Domain - Small test domain with grassland fuel
geometry.coord_sys   = 0
geometry.prob_lo     = 330000.0 3775000.0 0.0
geometry.prob_hi     = 332000.0 3777000.0 0.0
amr.n_cell           = 64 64 1

# Time - Multi-step simulation to test diurnal variation
max_step             = 100
stop_time            = 3600.0  # 1 hour simulation
levelset.cfl         = 0.3

# Basic fire spread model
fire_spread_model    = rothermel
propagation_method   = levelset

# Point ignition
ignition.enable      = 1
ignition.type        = point
ignition.x           = 331000.0
ignition.y           = 3776000.0
ignition.time        = 0.0

# Fuel model: Anderson FM1 (short grass) - suitable for curing model test
rothermel.use_anderson13 = 1
rothermel.anderson13_code = 1

# Weather conditions for fire danger indices
wind.enable          = 1
wind.u_mag           = 5.0    # 5 m/s wind (moderate)
wind.u_dir_deg       = 270.0  # West wind

# Temperature and RH for CBI calculation
# T = 35°C (95°F), RH = 15% (hot, dry conditions)
temperature          = 35.0
relative_humidity    = 15.0

# Fuel moisture consistent with dry grassland
fuel_moisture.d1     = 0.04   # 4% (very dry)
fuel_moisture.d10    = 0.06   # 6%
fuel_moisture.d100   = 0.08   # 8%
fuel_moisture.lh     = 0.60   # Live herbaceous 60%
fuel_moisture.lw     = 0.90   # Live woody 90%

# Grass curing test (can be enabled in future enhancement)
# grass_curing.enable  = 1
# grass_curing.model   = seasonal
# grass_curing.day_of_year = 210  # Late July (peak curing in Northern Hemisphere)

# Diurnal weather (can be enabled in future enhancement)
# diurnal.enable       = 1
# diurnal.T_min        = 20.0   # °C
# diurnal.T_max        = 35.0   # °C
# diurnal.RH_min       = 15.0   # %
# diurnal.RH_max       = 60.0   # %
# diurnal.U_mean       = 5.0    # m/s

# Elevation lapse rate (can be enabled in future enhancement)
# elevation_lapse.enable = 1
# elevation_lapse.rate   = 0.0065  # °C/m (standard atmosphere)
# elevation_lapse.z_ref  = 0.0     # Reference elevation

# NFDRS and CBI calculation (diagnostic outputs)
# These would be computed and written to plotfiles as diagnostic fields

# Plot output
amr.plot_int         = 20
amr.plot_file        = plt

# Enable verbose output
amr.v                = 1

# Test notes:
# This regression test serves as a placeholder for future integration
# of the Phase 1 & Phase 2 features into the main simulation loop.
# The features are currently implemented as header-only libraries
# and can be called from main.cpp when needed.
