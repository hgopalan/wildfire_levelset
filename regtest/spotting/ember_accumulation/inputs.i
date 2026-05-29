# Ember Accumulation and Probabilistic Ignition Test
# Tests: Ember density accumulation with decay and probabilistic ignition
# Features: Albini spotting + ember accumulation tracking

# Grid & domain (1000 m x 1000 m, UTM Zone 11N)
n_cell_x = 64
n_cell_y = 64
max_grid_size = 32
prob_lo_x = 330000.0
prob_lo_y = 3775000.0
prob_hi_x = 331000.0
prob_hi_y = 3776000.0

# Time & output (2 hour simulation)
nsteps = 150
cfl = 0.5
plot_int = 15
reinit_int = -1

# Wind: strong eastward for long-distance spotting
u_x = 8.0
u_y = 0.0

# Spherical ignition at western portion
source_type = sphere
center_x = 330150.0
center_y = 3775500.0
sphere_radius = 35.0

# Fuel: FM4 chaparral (high intensity for strong lofting)
rothermel.fuel_model = FM4
rothermel.M_f = 0.06   # dry fuel for high fire intensity

# Albini spotting to generate landing flux
albini_spotting.enable = 1
albini_spotting.terminal_velocity = 1.2     # firebrand fall speed [m/s]
albini_spotting.P_base = 0.02               # higher launch probability
albini_spotting.I_B_min = 50.0              # minimum Byram intensity [kW/m]
albini_spotting.spot_radius = 10.0          # ignition zone radius [m]
albini_spotting.check_interval = 3
albini_spotting.n_traj_steps = 100
albini_spotting.random_seed = 123

# Ember accumulation tracking with decay and probabilistic ignition
ember_accumulation.enable = 1
ember_accumulation.k_decay = 0.00167        # decay rate [1/s] (10-min burnout: 1/600 s)
ember_accumulation.rho_threshold = 8.0      # min density for ignition [embers/m²]
ember_accumulation.k_ignition = 1.2         # ignition rate coefficient
ember_accumulation.dt_check = 30.0          # check interval [s]
ember_accumulation.use_moisture_damping = 1 # enable fuel moisture effect on ignition

# Level-set propagation
propagation_method = levelset
