# Simple test: verify level 0 calculations produce finite phi values.
#
# Run with: mpirun -n 4 ./levelset tests/regrid_finite_check/inputs.i
# (or any number of MPI ranks >= 1)
#
# Verification: the simulation completes without crashing and the final
# phi_min reported in stdout is a finite number (not -inf or nan).

# Small grid to keep the test fast
n_cell_x = 32
n_cell_y = 32
n_cell_z = 32
prob_lo_x = -10.0
prob_hi_x =  10.0
prob_lo_y = -10.0
prob_hi_y =  10.0
prob_lo_z = -10.0
prob_hi_z =  10.0

# Short run: enough steps to verify basic functionality
nsteps   = 20
cfl      = 0.1
plot_int = -1

# Constant x-wind
u_x = 0.5
u_y = 0.0
u_z = 0.0

# Sphere SDF centred at origin with radius 3 (negative phi region exists)
source_type   = sphere
center_x      = 0.0
center_y      = 0.0
center_z      = 0.0
sphere_radius = 3.0

# Reinitialization disabled to keep test simple
reinit_int = -1
