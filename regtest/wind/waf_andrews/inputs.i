# Wind Adjustment Factor – Andrews (2018) logarithmic formula
# Tests: rothermel.use_waf = 1 with waf_formula = "andrews" (default)
#        and the Maximum Effective Wind Speed (MEWS) limit.
#
# The Andrews logarithmic WAF converts 20-ft open wind to midflame height:
#   WAF = 1.83 / ln((20 + 0.36*h) / (0.13*h))   [h in ft]
# For FM4 (h = 6 ft): WAF ≈ 0.36.  Wind input 5 m/s → midflame ≈ 1.8 m/s.
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

# Fuel: FM4 chaparral (h = 6 ft) – moderate moisture
rothermel.fuel_model = FM4
rothermel.M_f = 0.08

# Andrews (2018) wind adjustments
rothermel.use_waf        = 1
rothermel.waf_formula    = andrews   # logarithmic Albini & Baughman (1979) WAF
rothermel.use_wind_limit = 1         # cap wind speed at MEWS (Rothermel 1972)

# Level-set propagation
propagation_method = levelset
