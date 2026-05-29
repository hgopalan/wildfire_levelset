# Scott & Reinhardt (2001) Crown Fire Surface Area (CFSA) Test
# Tests: Crown fire surface area calculation based on canopy structure
# Reference: Scott & Reinhardt (2001) USDA Forest Service RMRS-RP-29

# Grid & domain (400 m x 400 m, UTM Zone 11N)
n_cell_x = 64
n_cell_y = 64
max_grid_size = 32
prob_lo_x = 330000.0
prob_lo_y = 3775000.0
prob_hi_x = 330400.0
prob_hi_y = 3775400.0

# Time & output
nsteps = 50
cfl = 0.5
plot_int = 10
reinit_int = -1

# Moderate wind to generate crown fire
u_x = 6.0
u_y = 0.0

# Spherical ignition at domain centre
source_type = sphere
center_x = 330200.0
center_y = 3775200.0
sphere_radius = 25.0

# Fuel: FM10 Timber Litter with crown parameters
rothermel.fuel_model = FM10
rothermel.M_f = 0.08

# Crown parameters for active crown fire
# These values should produce significant CFSA
crown.enable = 1
crown.CBH = 3.0        # canopy base height [m]
crown.CBD = 0.20       # canopy bulk density [kg/m3] (moderate-high)
crown.FMC = 100.0      # foliar moisture content [%]
crown.CH = 15.0        # canopy height [m]

# CFSA-specific parameters
# k_sa controls the surface area coefficient [m²/kg]
# Higher CBD and canopy depth → higher CFSA
# Expected CFSA (uncapped) = CBD × (CH - CBH) × k_sa
#                          = 0.20 × (15.0 - 3.0) × 6.0 = 14.4
# Expected CFSA (capped)   = min(14.4, max_cfsa) = min(14.4, 3.0) = 3.0
crown.cfsa_k_sa = 6.0
crown.cfsa_max = 3.0

# Level-set propagation with crown fire
propagation_method = levelset
crown.method = rothermel1991  # Active crown fire: R_crown = 3.34 × R_surface
