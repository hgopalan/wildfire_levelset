# Periodic Wind Gust Factor Test
# Tests: Sinusoidal periodic wind modulation
# Feature: Periodic gust factor with adjustable amplitude and period

# Grid & domain (1000 m x 1000 m, UTM Zone 11N)
n_cell_x = 64
n_cell_y = 64
max_grid_size = 32
prob_lo_x = 330000.0
prob_lo_y = 3775000.0
prob_hi_x = 331000.0
prob_hi_y = 3776000.0

# Time & output (simulate 2 hours to see multiple gust cycles)
nsteps = 200
cfl = 0.5
plot_int = 20
reinit_int = -1

# Base wind: moderate eastward (will be modulated by gusts)
u_x = 5.0
u_y = 0.0

# Spherical ignition at center
source_type = sphere
center_x = 330500.0
center_y = 3775500.0
sphere_radius = 30.0

# Fuel: FM4 chaparral
rothermel.fuel_model = FM4
rothermel.M_f = 0.08

# Periodic gust factor
# Wind is modulated as: V(t) = V_base × (1 + A × sin(2π × t / T))
# where A = amplitude (0.3 = 30% variation), T = period
turb_wind.enable = 1
turb_wind.gust_amplitude = 0.4      # 40% amplitude: wind varies from 0.6× to 1.4× base
turb_wind.gust_period = 600.0       # 10-min period (600 s)
turb_wind.gust_phase = 0.0          # initial phase [radians]

# Level-set propagation
propagation_method = levelset
