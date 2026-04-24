# ============================================================================
# Gaussian Hill Test Case with Fixed Richards' Coefficients
# ============================================================================
# This test case is identical to gaussian_hill_anderson except it uses
# fixed Richards' coefficients instead of Anderson L/W ratio.
# Use this for comparison to see the difference between dynamic and fixed L/W.
# ============================================================================

# -------------------- Grid & Domain --------------------
n_cell_x = 200
n_cell_y = 200
n_cell_z = 64

prob_lo_x = 0.0
prob_hi_x = 1000.0
prob_lo_y = 0.0
prob_hi_y = 1000.0
prob_lo_z = 0.0
prob_hi_z = 200.0

max_grid_size = 32
amr_enable_negative_phi_refine = 0

# -------------------- Time & Output --------------------
nsteps = 500
cfl    = 0.4
plot_int = 25

# -------------------- Velocity Field --------------------
# Same wind as Anderson test: 5 m/s
u_x = 5.0
u_y = 0.0
u_z = 0.0

# -------------------- Initial Fire Source --------------------
source_type = sphere
center_x = 500.0
center_y = 500.0
center_z = 10.0
sphere_radius = 25.0

# -------------------- Rothermel Fire Model --------------------
rothermel.fuel_model = FM4

# Same terrain slope as Anderson test
rothermel.slope_x = 0.15
rothermel.slope_y = 0.0

# -------------------- FARSITE Ellipse Model --------------------
farsite.enable = 1

# Use FIXED Richards' coefficients (not Anderson L/W)
farsite.use_anderson_LW = 0

farsite.phi_threshold = 0.5

# Fixed Richards' coefficients
# These correspond to L/W ratio of approximately 1.2
# (much less elongated than Anderson's ~2.3 for this wind speed)
farsite.coeff_a = 1.0
farsite.coeff_b = 0.5
farsite.coeff_c = 0.2

# -------------------- Level Set Reinitialization --------------------
reinit_int = 10
