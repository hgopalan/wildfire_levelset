# Regression test: verify regrid_negative_phi produces finite phi after regridding.
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

# Enable AMR regridding — this exercises the fail-safe paths
amr_enable_negative_phi_refine = 1
amr_regrid_int      = 5
amr_refine_ratio    = 2
amr_max_refinements = 1

# Short run: enough steps to trigger at least one regrid
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

# Reinitialization disabled so we isolate the regrid path
reinit_int = -1
