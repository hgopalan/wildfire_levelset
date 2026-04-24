# ============================================================================
# Gaussian Hill Test Case with Anderson L/W Ratio
# ============================================================================
# This test case simulates fire spread on a Gaussian-shaped hill using the
# Anderson (1983) L/W ratio formulation for elliptical fire spread.
#
# Terrain: Simplified representation of a Gaussian hill using constant slope
#          In a full implementation, slope would vary spatially
# Wind: Uniform wind field to demonstrate Anderson L/W ratio response
# Fire Model: FARSITE with Anderson (1983) dynamic L/W ratio based on wind speed
# ============================================================================

# -------------------- Grid & Domain --------------------
n_cell_x = 200
n_cell_y = 200
n_cell_z = 64

# Domain size: 1000m x 1000m horizontal, 200m vertical
prob_lo_x = 0.0
prob_hi_x = 1000.0
prob_lo_y = 0.0
prob_hi_y = 1000.0
prob_lo_z = 0.0
prob_hi_z = 200.0

max_grid_size = 32

# AMR disabled for this test
amr_enable_negative_phi_refine = 0

# -------------------- Time & Output --------------------
nsteps = 500
cfl    = 0.4
plot_int = 25

# -------------------- Velocity Field --------------------
# Moderate wind: 5 m/s in x-direction (~11 mph)
# This will produce Anderson L/W ratio of approximately 2.3
u_x = 5.0
u_y = 0.0
u_z = 0.0

# -------------------- Initial Fire Source --------------------
source_type = sphere
# Fire starts at the base of the hill (center of domain)
center_x = 500.0
center_y = 500.0
center_z = 10.0
sphere_radius = 25.0

# -------------------- Rothermel Fire Model --------------------
# Using Anderson fuel model 4 (Chaparral)
# This is a common fuel type for hillside fires
rothermel.fuel_model = FM4

# Gaussian hill terrain parameters
# Simplified constant slope approximation
# For a Gaussian hill with height H and characteristic width σ:
# Maximum slope ≈ H/(σ√e) ≈ 0.6 * H/σ
# Using H=100m, σ=300m gives max slope ≈ 0.2 (11.3 degrees)
rothermel.slope_x = 0.15    # tan(8.5°) - moderate uphill slope in x-direction
rothermel.slope_y = 0.0     # flat in y-direction

# Note: In a full Gaussian hill implementation, slope would be computed as:
# z(x,y) = H * exp(-((x-x0)²+(y-y0)²)/(2σ²))
# slope_x = ∂z/∂x = -H(x-x0)/σ² * exp(-r²/(2σ²))
# slope_y = ∂z/∂y = -H(y-y0)/σ² * exp(-r²/(2σ²))
# where r² = (x-x0)² + (y-y0)²

# -------------------- FARSITE Ellipse Model --------------------
farsite.enable = 1

# Enable Anderson (1983) L/W ratio based on wind speed
farsite.use_anderson_LW = 1

# Fire front detection threshold
farsite.phi_threshold = 0.5

# Note: When use_anderson_LW = 1, the coefficients a, b, c are computed
# dynamically from wind speed using Anderson (1983) formula.
# The following fixed coefficients are only used when use_anderson_LW = 0:
farsite.coeff_a = 1.0
farsite.coeff_b = 0.5
farsite.coeff_c = 0.2

# -------------------- Level Set Reinitialization --------------------
reinit_int = 10
