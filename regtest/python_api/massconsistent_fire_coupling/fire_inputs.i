# Fire solver inputs for Rothermel model coupling
# This file demonstrates Rothermel spread model for use with massconsistent_amr

# Problem setup
max_level = 0

# Domain
geometry.prob_lo = 0.0 0.0
geometry.prob_hi = 1000.0 1000.0
amr.n_cell = 32 32

# Grid refinement (disable for 2D)
amr.max_level = 0

# Spread model: Use Rothermel for two-way wind coupling
fire.spread_model = rothermel
fire.fuel_model = 8  # Timber litter with some underbrush (Anderson 13)

# Ignition configuration
fire.ignition.type = circle
fire.ignition.x_center = 500.0
fire.ignition.y_center = 500.0
fire.ignition.radius = 25.0
fire.ignition.time = 0.0

# Initial conditions
fire.fuel_moisture = 10.0 12.0 15.0 80.0  # 1hr, 10hr, 100hr, live herb (%)

# Wind configuration
# Wind is provided via update_wind() or update_wind_3d() calls from Python
# This is a baseline/default wind for if not coupled
fire.wind.speed = 5.0  # m/s
fire.wind.direction = 270.0  # from west (degrees)

# Terrain (flat for this example)
fire.terrain.type = flat
fire.terrain.zmin = 0.0
fire.terrain.zmax = 0.0

# Output
amr.plot_int = 10
amr.checkpoint_int = -1

# Time stepping
fire.cfl = 0.8
fire.initial_dt = 0.1
fire.dt_multiplier = 1.1
fire.dt_max = 1.0

# Simulation control
max_time = 3600.0  # 1 hour
nsteps = 36000  # Safety limit

# Diagnostic/debug output
amr.verbose = 1
