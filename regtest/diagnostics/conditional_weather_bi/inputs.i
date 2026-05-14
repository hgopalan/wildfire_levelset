# Conditional Weather BI / SC Trigger Test
# Tests that conditional_weather_trigger = "bi" (NFDRS Burning Index) selects
# the correct weather scenario from the conditional weather table.
# Verifies:
#   - Simulation completes without error when trigger = "bi"
#   - Dead fuel moistures are updated from the conditional weather table each step
#   - Wind components are updated from the table

# Grid & domain (500 m × 500 m)
n_cell_x = 64
n_cell_y = 64
max_grid_size = 32
prob_lo_x = 0.0
prob_lo_y = 0.0
prob_hi_x = 500.0
prob_hi_y = 500.0

# Time & output
nsteps = 40
cfl = 0.5
plot_int = 20
reinit_int = -1

# Wind (initial; may be overridden by conditional weather table)
u_x = 3.0
u_y = 0.0

# Spherical ignition at domain centre
source_type = sphere
center_x = 250.0
center_y = 250.0
sphere_radius = 25.0

# Fuel: FM4 chaparral – produces measurable BI and ERC
rothermel.fuel_model = FM4
rothermel.M_f = 0.08

# Fire ecology required to compute BI (ecology component 8)
# (enabled by default; fire_ecology parameters are at defaults)

# Conditional weather: trigger on NFDRS Burning Index
conditional_weather_file    = cond_weather.csv
conditional_weather_trigger = bi

# Fire statistics
fire_stats_file = fire_stats.csv

# Level-set propagation
propagation_method = levelset
