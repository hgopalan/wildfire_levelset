# Example: Complete Feature Integration Test
# This inputs file demonstrates all 10 newly integrated wildfire features
# Run with: ./levelset inputs_features_demo.i

# ============================================================================
# Basic Setup
# ============================================================================
n_cell_x = 64
n_cell_y = 64
max_grid_size = 32

prob_lo_x = 0.0
prob_lo_y = 0.0
prob_hi_x = 1000.0
prob_hi_y = 1000.0

reinit_int = 20
nsteps = 100
cfl = 0.5
plot_int = 20

# ============================================================================
# Velocity (constant wind in x-direction)
# ============================================================================
u_x = 5.0    # 5 m/s wind
u_y = 0.0
u_z = 0.0

# ============================================================================
# Initial Ignition (sphere)
# ============================================================================
source_type = sphere
center_x = 0.1
center_y = 0.5
sphere_radius = 0.05

# ============================================================================
# Rothermel Fire Spread Model
# ============================================================================
rothermel.fuel_model = "Model 4"   # Southern California chaparral

# Feature 1: Fuel Continuity Factor
# Default 1.0 (continuous fuels); reduce for patchy/discontinuous fuels
rothermel.fuel_continuity = 0.9    # 90% fuel coverage

# ============================================================================
# FARSITE Level-Set Propagation
# ============================================================================
propagation_method = levelset
farsite.length_to_width_ratio = 3.0

# ============================================================================
# Feature 3: Crown Fraction Burned Diagnostic
# ============================================================================
crown_fraction.enable = 1          # Compute CFB diagnostic

# ============================================================================
# Feature 4: Effective Wind Speed (Combined Wind + Slope)
# ============================================================================
effective_wind.enable = 1          # Requires terrain file for slopes
# Note: Only active if landscape_file is configured

# ============================================================================
# Feature 5: Thomas Flame Length Model (Alternative to Byram)
# ============================================================================
# Choose between "byram" (default) and "thomas"
flame_length_model.model = byram    # Change to "thomas" to use Thomas 1963

# ============================================================================
# Feature 6: Fuel Boundary Smoothing
# ============================================================================
# Eliminates unrealistic ROS discontinuities at fuel model boundaries
fuel_boundary.enable = 0           # Requires landscape file with fuel data
fuel_boundary.transition_cells = 2.5

# ============================================================================
# Feature 7: CSIRO Grassfire Acceleration (Non-Equilibrium Growth)
# ============================================================================
# Models acceleration of fire spread during initial phase
grassfire_accel.enable = 0         # Set to 1 for grassland scenarios
grassfire_accel.t_accel = 600.0    # Time constant [seconds] = 10 minutes

# ============================================================================
# Feature 8: Burnout Time Separation (Flaming vs Smoldering)
# ============================================================================
# Splits residence time into flaming and smoldering phases by fuel type
burnout_separation.enable = 1
burnout_separation.flaming_fraction_fine = 0.70      # Fine fuels: 70% flaming
burnout_separation.flaming_fraction_medium = 0.40    # Medium: 40% flaming
burnout_separation.flaming_fraction_heavy = 0.20     # Heavy: 20% flaming
burnout_separation.flaming_fraction_duff = 0.10      # Duff: 10% flaming

# ============================================================================
# Feature 9: Simard Moisture Model (Exponential Time-Lag)
# ============================================================================
# Alternative moisture update using exponential approach to equilibrium
simard_moisture.enable = 0         # Not yet integrated into main loop
simard_moisture.tau_1hr = 1.0      # 1-hour fuel lag [hours]
simard_moisture.tau_10hr = 10.0    # 10-hour fuel lag [hours]
simard_moisture.tau_100hr = 100.0  # 100-hour fuel lag [hours]

# ============================================================================
# Feature 10: Post-Frontal Smoldering (Residual Heat Release)
# ============================================================================
# Tracks residual combustion after flame front passage
post_frontal.enable = 1            # Enable post-frontal tracking
post_frontal.tau_fine = 1800.0     # Fine fuels decay [seconds] = 30 min
post_frontal.tau_medium = 3600.0   # Medium fuels decay = 1 hour
post_frontal.tau_heavy = 7200.0    # Heavy fuels decay = 2 hours
post_frontal.tau_duff = 21600.0    # Duff decay = 6 hours

# ============================================================================
# Output Configuration
# ============================================================================
plot_int = 20
plot_dir = plt

# Optional: Select specific variables to output (reduces file size)
# plot_vars = phi R arrival_time flame_length crown_fraction_burned
#             effective_wind_speed burnout_flaming_time residual_heat_release

# ============================================================================
# Fire Report Output
# ============================================================================
fire_stats_file = fire_stats.csv
write_perimeter_geojson = 0

# ============================================================================
# Crown Fire Model (Optional)
# ============================================================================
crown.enable = 1
crown.CBH = 3.0          # Crown base height [m]
crown.CBD = 0.15         # Canopy bulk density [kg/m³]
crown.FMC = 100.0        # Foliar moisture content [%]

# ============================================================================
# Viegas Eruptive Fire Diagnostics (Optional)
# ============================================================================
viegas.enable = 0
