# Fine Fuel Moisture Time-Lag Differential Equations Test
# Tests: Physically-based fuel moisture evolution with time-lag response
# References: Nelson (2000), Viney (1991)

# Grid & domain (250 m x 250 m, UTM Zone 11N)
n_cell_x = 40
n_cell_y = 40
max_grid_size = 20
prob_lo_x = 330000.0
prob_lo_y = 3775000.0
prob_hi_x = 330250.0
prob_hi_y = 3775250.0

# Time & output (longer simulation to observe moisture dynamics)
nsteps = 100
cfl = 0.5
plot_int = 20
reinit_int = -1

# Moderate wind
u_x = 5.0
u_y = 0.0

# Spherical ignition at domain centre
source_type = sphere
center_x = 330125.0
center_y = 3775125.0
sphere_radius = 15.0

# Fuel: FM1 Short Grass (responds quickly to moisture changes)
rothermel.fuel_model = FM1
rothermel.M_f = 0.12  # Initial moisture: 12%

# Fuel moisture time-lag differential equation parameters
fuel_moisture_de.enable = 1
fuel_moisture_de.method = timelag_de  # Use differential equation solver

# Weather conditions (changing during simulation)
# Temperature affects drying rate
# RH affects equilibrium moisture content
fuel_moisture_de.T_ambient = 25.0    # Temperature [°C]
fuel_moisture_de.RH = 40.0           # Relative humidity [%]

# Time-lag constants for different size classes [hours]
fuel_moisture_de.tau_1hr = 1.0       # 1-hour fuels
fuel_moisture_de.tau_10hr = 10.0     # 10-hour fuels
fuel_moisture_de.tau_100hr = 100.0   # 100-hour fuels

# Initial moisture contents [fraction]
fuel_moisture_de.M_1hr_init = 0.12   # 12%
fuel_moisture_de.M_10hr_init = 0.15  # 15%
fuel_moisture_de.M_100hr_init = 0.18 # 18%

# Precipitation wetting (can add rainfall event during simulation)
fuel_moisture_de.rain_rate = 0.0     # [mm/h] - no rain initially

# Temperature correction for time-lag
fuel_moisture_de.temp_correction = 1  # Enable temperature-dependent drying

# Hysteresis (adsorption vs desorption curves)
fuel_moisture_de.use_hysteresis = 1   # Enable hysteresis

# Level-set propagation
propagation_method = levelset
