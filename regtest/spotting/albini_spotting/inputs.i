# Albini (1983) Firebrand Spotting Test
# Tests: physics-based Albini spotting with 2-D trajectory integration
# Reference: Albini (1983) USDA Forest Service Research Paper INT-309

# Grid & domain (500 m x 500 m, UTM Zone 11N)
n_cell_x = 64
n_cell_y = 64
max_grid_size = 32
prob_lo_x = 330000.0
prob_lo_y = 3775000.0
prob_hi_x = 330500.0
prob_hi_y = 3775500.0

# Time & output
nsteps = 100
cfl = 0.5
plot_int = 10
reinit_int = -1

# Wind: 5 m/s eastward
u_x = 5.0
u_y = 0.0
u_z = 0.0

# Spherical ignition at western portion
source_type = sphere
center_x = 330100.0
center_y = 3775250.0
center_z = 0.5
sphere_radius = 25.0

# Fuel: FM4 chaparral (high heat content drives strong lofting)
rothermel.fuel_model = FM4
rothermel.M_f = 0.07

# Albini (1983) physics-based firebrand spotting
# Lofting height: H_z = 12.2 * I_B^(1/3) [m]
# Trajectory: 2-D forward-Euler with 100 sub-steps
albini_spotting.enable = 1
albini_spotting.terminal_velocity = 1.0     # firebrand fall speed [m/s]
albini_spotting.P_base = 0.01               # launch probability per front cell
albini_spotting.I_B_min = 20.0             # minimum Byram intensity [kW/m]
albini_spotting.spot_radius = 8.0          # ignition zone radius [m]
albini_spotting.check_interval = 5
albini_spotting.n_traj_steps = 100
albini_spotting.random_seed = 42

# Level-set propagation
propagation_method = levelset
