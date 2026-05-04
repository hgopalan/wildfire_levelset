# Basic Level-Set Advection Test
# Tests: Simple sphere advection with constant velocity

# Grid & domain
n_cell = 64
max_grid_size = 32
prob_lo_x = 0.0
prob_lo_y = 0.0
prob_lo_z = 0.0
prob_hi_x = 1.0
prob_hi_y = 1.0
prob_hi_z = 1.0

# Time & output
final_time = 2.0
cfl = 0.5
plot_int = 20

# Velocity (constant)
u_x = 0.25
u_y = 0.1
u_z = 0.0

# Initial source - sphere
source_type = sphere
sphere_center_x = 0.3
sphere_center_y = 0.5
sphere_center_z = 0.5
sphere_radius = 0.15

# Reinitialization
reinit_int = 20
reinit_iters = 20
reinit_dtau = 0.5

# Level set control
# Level set control
propagation_method = levelset


# ---- Dynamic fire points (optional) ----
# When this file is present on disk at the start of any time step, the ignition
# points listed in it (CSV format: X Y [Z]) are merged into phi and the file is
# renamed to <name>.applied so it is not re-applied the next step.  Drop a new
# file under the same name to trigger additional ignitions at any later time.
#
# dynamic_fire_points_file = new_ignitions.csv
# fire_gaussian_sigma = -1.0   # <= 0 means auto (3 cells)

# ---- Checkpoint / restart (optional) ----
# Write a checkpoint directory every chk_int steps.
# chk_int = 50
#
# Restart from a previously written checkpoint directory.
# restart_chkfile = chk0050
