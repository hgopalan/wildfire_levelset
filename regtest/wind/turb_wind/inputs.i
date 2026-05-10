# Turbulent Wind Test (2-D build only)
# Tests: spectral-noise turbulent wind perturbation model
# The spectral_noise model uses Random Fourier Features with OU temporal evolution.
# Requires: cmake -DLEVELSET_DIM_2D=ON

# Grid & domain (500 m x 500 m)
n_cell_x = 64
n_cell_y = 64
max_grid_size = 32
prob_lo_x = 0.0
prob_lo_y = 0.0
prob_hi_x = 500.0
prob_hi_y = 500.0

# Time & output
nsteps = 80
cfl = 0.5
plot_int = 10
reinit_int = -1

# Mean wind: 4 m/s eastward
u_x = 4.0
u_y = 0.0

# Spherical ignition at domain centre
source_type = sphere
center_x = 250.0
center_y = 250.0
sphere_radius = 30.0

# Fuel: FM4 chaparral
rothermel.fuel_model = FM4
rothermel.M_f = 0.08

# Spectral-noise turbulent wind perturbation
# OU temporal correlation with spatially-correlated Fourier modes
turb_wind.model = spectral_noise
turb_wind.theta = 0.1       # OU reversion rate [1/s]
turb_wind.sigma = 0.5       # perturbation std dev [m/s]
turb_wind.L_c   = 50.0      # spatial correlation length [m] (required for spectral_noise)
turb_wind.N_modes = 32
turb_wind.random_seed = 7

# Level-set propagation
propagation_method = levelset
