# Multi-Layer Canopy Wind Profile Test
# Tests: Vertical wind profile through multi-layer canopy
# References: Massman (1997), Inoue (1963)

# Grid & domain (300 m x 300 m, UTM Zone 11N)
n_cell_x = 48
n_cell_y = 48
max_grid_size = 24
prob_lo_x = 330000.0
prob_lo_y = 3775000.0
prob_hi_x = 330300.0
prob_hi_y = 3775300.0

# Time & output
nsteps = 40
cfl = 0.5
plot_int = 10
reinit_int = -1

# Reference wind speed at 10m height
# Wind will be attenuated in canopy layers
u_x = 8.0
u_y = 0.0

# Spherical ignition at domain centre
source_type = sphere
center_x = 330150.0
center_y = 3775150.0
sphere_radius = 20.0

# Fuel: FM10 Timber Litter
rothermel.fuel_model = FM10
rothermel.M_f = 0.10

# Dense canopy for significant wind attenuation
crown.enable = 1
crown.CBH = 4.0        # canopy base height [m]
crown.CBD = 0.25       # canopy bulk density [kg/m3] (dense)
crown.FMC = 110.0      # foliar moisture content [%]
crown.CH = 20.0        # canopy height [m]

# Multi-layer canopy wind parameters
canopy_wind.enable = 1
canopy_wind.n_layers = 5             # 5 vertical layers
canopy_wind.alpha = 2.5              # attenuation coefficient (moderate density)
canopy_wind.LAI = 4.0                # leaf area index [m²/m²]
canopy_wind.d_ratio = 0.7            # displacement height ratio
canopy_wind.z0_ratio = 0.1           # roughness length ratio

# Reference height for input wind speed
canopy_wind.z_ref = 10.0             # [m] - wind specified at 10m

# Level-set propagation
propagation_method = levelset
