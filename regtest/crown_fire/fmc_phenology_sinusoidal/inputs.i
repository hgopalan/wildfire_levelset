# FMC Sinusoidal Phenology Test
# Tests: Sinusoidal foliar moisture content phenology model
# Feature: Seasonal FMC variation affects crown fire initiation

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
rothermel.M_f = 0.06   # dry surface fuel for high intensity

# Crown fire parameters
crown_fire.enable = 1
crown_fire.CBD = 0.15      # canopy bulk density [kg/m³]
crown_fire.CBH = 2.0       # canopy base height [m]
crown_fire.h_c = 8.0       # canopy height [m]
crown_fire.M_lw = 1.0      # live woody moisture (will be overridden by phenology)

# Sinusoidal FMC phenology model
# FMC varies seasonally: FMC(t) = FMC_mean + FMC_amplitude × sin(2π(t - t_peak)/T_year)
fmc_phenology.model = sinusoidal
fmc_phenology.FMC_mean = 100.0          # mean FMC [%]
fmc_phenology.FMC_amplitude = 40.0      # seasonal variation ±40%
fmc_phenology.peak_doy = 150            # day of year for peak FMC (late May)
fmc_phenology.start_doy = 200           # simulation starts at day 200 (mid-July)

# Expected FMC ≈ 70% at mid-summer (DOY 200)
# Lower FMC → lower crown fire initiation threshold

# Level-set propagation
propagation_method = levelset
