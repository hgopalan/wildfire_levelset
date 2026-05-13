# Crown Fire Initiation Test
# Tests: Van Wagner (1977) crown fire initiation with Rothermel surface fire

# Grid & domain (unit cube)
n_cell = 64
max_grid_size = 32
prob_lo_x = 0.0
prob_lo_y = 0.0
prob_hi_x = 1.0
prob_hi_y = 1.0

# Time & output
nsteps = 50
cfl = 0.5
plot_int = 10
reinit_int = 20

# Moderate wind
u_x = 0.3
u_y = 0.1

# Spherical ignition at domain centre
source_type = sphere
center_x = 0.5
center_y = 0.5
sphere_radius = 0.1

# Fuel: FM10 Timber Litter and Understory; moderate moisture
rothermel.fuel_model = FM10
rothermel.M_f = 0.08

# Van Wagner (1977) crown fire initiation
crown.enable = 1
crown.CBH = 4.0       # canopy base height [m]
crown.CBD = 0.15      # canopy bulk density [kg/m3]
crown.FMC = 100.0     # foliar moisture content [%]
crown.crown_fraction_weight = 1.0
crown.use_metric_units = 1

# Level-set propagation
propagation_method = levelset
