# Bulk Fuel Consumption Fraction Test
# Tests: sigmoid fuel consumption model (f_c varies from f_min to f_max with intensity)

# Grid & domain (unit cube)
n_cell = 64
max_grid_size = 32
prob_lo_x = 0.0
prob_lo_y = 0.0
prob_hi_x = 1.0
prob_hi_y = 1.0

# Time & output
nsteps = 80
cfl = 0.5
plot_int = 10
reinit_int = -1

# Moderate wind
u_x = 0.25
u_y = 0.1

# Spherical ignition
source_type = sphere
center_x = 0.5
center_y = 0.5
sphere_radius = 0.12

# Fuel: FM4 chaparral with dry conditions for higher intensity
rothermel.fuel_model = FM4
rothermel.M_f = 0.07

# Bulk fuel consumption fraction model
# f_c = f_min + (f_max - f_min) * sigmoid(I_norm - 1)
farsite.use_bulk_fuel_consumption = 1
farsite.tau_residence = 60.0        # post-frontal combustion residence time [s]
farsite.f_consumed_max = 0.90       # max fraction (slow, high-intensity fires)
farsite.f_consumed_min = 0.50       # min fraction (fast, low-intensity fires)

# FARSITE propagation with Anderson L/W
farsite.use_anderson_LW = 1
farsite.phi_threshold = 0.1
propagation_method = farsite
