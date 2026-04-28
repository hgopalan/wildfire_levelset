# Albini (1983) Spotting Model with 2-D Trajectory Integration Test
#
# A 1 km x 1 km domain with a moderate 2 m/s wind blowing in the x-direction.
# A box ignition at the left edge starts the fire; FARSITE ellipse spread is
# used to evolve the fire front.  The Albini spotting model is enabled and
# will loft firebrands using Byram-derived fire line intensity, then integrate
# their 2-D horizontal trajectory through the wind field.
#
# Expected behaviour:
#   - Firebrand spot fires appear ahead of the main fire front in the
#     downwind direction (~100-200 m ahead).
#   - Plot files contain albini_Hz, albini_count, albini_dist, albini_active
#     diagnostic fields.
#
# Physics check (FM4 defaults, u=2 m/s):
#   R ≈ 88 ft/min, I_B ≈ 9355 kW/m, H_z ≈ 257 m
#   Flight time = 257 m / 5 m/s ≈ 51 s  →  landing ≈ 2 m/s * 51 s ≈ 103 m

# Grid & domain (1 km x 1 km, dx = 10 m)
n_cell_x = 100
n_cell_y = 100
max_grid_size = 50
prob_lo_x = 0.0
prob_lo_y = 0.0
prob_hi_x = 1000.0
prob_hi_y = 1000.0

# Time & output
nsteps = 100
cfl = 0.5
plot_int = 25
reinit_int = -1

# Constant wind: 2 m/s in x, 0.5 m/s in y (slight cross-wind)
u_x = 2.0
u_y = 0.5
u_z = 0.0

# Initial source: box ignition on the left edge
source_type = box
box_xmin = 50.0
box_xmax = 100.0
box_ymin = 300.0
box_ymax = 700.0

# Rothermel fuel model (NFFL Model 4 - chaparral, standard moisture)
rothermel.fuel_model = FM4
rothermel.M_f = 0.08

# FARSITE ellipse model (Richards 1990) with Anderson L/W
farsite.enable = 1
farsite.use_anderson_LW = 1
farsite.phi_threshold = 0.1

# Skip level set; use FARSITE ellipse only
skip_levelset = 1

# Albini (1983) spotting model with 2-D trajectory
albini_spotting.enable = 1
albini_spotting.terminal_velocity = 5.0   # firebrand descent speed [m/s]
albini_spotting.P_base = 0.05             # max launch probability per fire-front cell
albini_spotting.I_B_min = 500.0           # min fire line intensity to generate firebrands [kW/m]
albini_spotting.spot_radius = 15.0        # radius of new spot-fire ignition zone [m]
albini_spotting.random_seed = 42          # reproducible RNG
albini_spotting.check_interval = 5        # run spotting every 5 steps
albini_spotting.n_traj_steps = 200        # forward-Euler trajectory sub-steps
