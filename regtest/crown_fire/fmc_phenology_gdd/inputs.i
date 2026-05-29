# FMC Growing Degree Day (GDD) Phenology Test
# Tests: GDD-based foliar moisture content phenology model
# Feature: Temperature-driven greenup affects crown fire behavior

# Grid & domain (1000 m x 1000 m, UTM Zone 11N)
n_cell_x = 64
n_cell_y = 64
max_grid_size = 32
prob_lo_x = 330000.0
prob_lo_y = 3775000.0
prob_hi_x = 331000.0
prob_hi_y = 3776000.0

# Time & output (simulate 1 hour)
nsteps = 100
cfl = 0.5
plot_int = 10
reinit_int = -1

# Wind: moderate eastward
u_x = 6.0
u_y = 0.0

# Spherical ignition
source_type = sphere
center_x = 330300.0
center_y = 3775500.0
sphere_radius = 30.0

# Fuel: FM4 chaparral with crown layer
rothermel.fuel_model = FM4
rothermel.M_f = 0.06

# Crown fire parameters
crown_fire.enable = 1
crown_fire.CBD = 0.15
crown_fire.CBH = 2.0
crown_fire.h_c = 8.0
crown_fire.M_lw = 1.0

# Growing Degree Day (GDD) phenology model
# FMC increases as GDD accumulates above base temperature
fmc_phenology.model = gdd
fmc_phenology.FMC_dormant = 80.0        # dormant/winter FMC [%]
fmc_phenology.FMC_peak = 140.0          # peak greenup FMC [%]
fmc_phenology.GDD_start = 50.0          # GDD threshold for greenup start
fmc_phenology.GDD_peak = 400.0          # GDD for full greenup
fmc_phenology.T_base = 5.0              # base temperature [°C]
fmc_phenology.start_doy = 100           # simulation starts at day 100 (early April)
fmc_phenology.mean_T = 15.0             # mean daily temperature [°C]

# At DOY 100 with mean_T = 15°C:
# Daily GDD = max(0, 15 - 5) = 10
# Total GDD = 100 days × 10 = 1000 (exceeds GDD_peak = 400)
# Expected FMC = FMC_peak = 140% (full greenup)

# Level-set propagation
propagation_method = levelset
