# Wind Adjustment Factor – BehavePlus linear formula
# Tests: rothermel.use_waf = 1 with waf_formula = "behaviorplus"
#        and the exponential canopy attenuation parameter waf_canopy_alpha.
#
# BehavePlus linear WAF for open/shrub fuels:
#   WAF = 0.36 + 0.004 × h_in   [h_in = fuel depth in inches]
# For FM4 (h = 6 ft = 72 in): WAF = 0.36 + 0.288 = 0.648.
# For FM1 (h = 1 ft = 12 in): WAF = 0.36 + 0.048 = 0.408.
#
# The waf_canopy_alpha parameter controls the exponential sheltering decay under
# closed canopies (only active when a landscape file provides canopy cover and
# height data).  This test exercises the formula on open fuel only.
#
# Domain: 500 m × 500 m (UTM Zone 11N, Southern California reference)
# Requires: 2-D build  (cmake -DLEVELSET_DIM_2D=ON)

# Grid & domain
n_cell_x = 64
n_cell_y = 64
max_grid_size = 32
prob_lo_x = 330000.0
prob_lo_y = 3775000.0
prob_hi_x = 330500.0
prob_hi_y = 3775500.0

# Time & output
nsteps = 60
cfl = 0.5
plot_int = 20
reinit_int = -1

# Wind: 5 m/s eastward (20-ft reference; WAF will reduce to midflame)
u_x = 5.0
u_y = 0.0

# Spherical ignition at domain centre
source_type = sphere
center_x = 330250.0
center_y = 3775250.0
sphere_radius = 30.0

# Fuel: FM4 chaparral (h = 6 ft = 72 in) – moderate moisture
rothermel.fuel_model = FM4
rothermel.M_f = 0.08

# BehavePlus WAF
rothermel.use_waf          = 1
rothermel.waf_formula      = behaviorplus  # linear: WAF = 0.36 + 0.004*h_in
rothermel.waf_canopy_alpha = 1.5          # canopy attenuation α_c (default)
rothermel.use_wind_limit   = 1            # cap wind speed at MEWS

# Level-set propagation
propagation_method = levelset
