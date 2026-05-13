# FMC Seasonal Schedule Test
# Tests: foliar moisture content updated via built-in parametric phenological curve
# Van Wagner (1977) crown fire initiation uses time-varying FMC

# Grid & domain (500 m x 500 m, UTM Zone 11N)
n_cell_x = 64
n_cell_y = 64
max_grid_size = 32
prob_lo_x = 330000.0
prob_lo_y = 3775000.0
prob_hi_x = 330500.0
prob_hi_y = 3775500.0

# Simulation: 600 s (~10 min) to observe FMC variation
final_time = 600.0
cfl = 0.5
plot_int = 10
reinit_int = 20

# Moderate wind
u_x = 4.0
u_y = 0.0

# Spherical ignition at domain centre
source_type = sphere
center_x = 330250.0
center_y = 3775250.0
sphere_radius = 30.0

# Fuel: chaparral
rothermel.fuel_model = FM4
rothermel.M_f = 0.08

# Crown fire (Van Wagner 1977) – uses FMC from the schedule
crown.enable = 1
crown.CBH = 3.5
crown.CBD = 0.12
crown.FMC = 100.0       # initial FMC [%]; overridden each step by fmc_schedule
crown.use_metric_units = 1

# FMC seasonal schedule: FARSITE phenological curve starting on day-of-year 200
# (mid-July: near peak green-up in chaparral)
fmc_schedule.enable = 1
fmc_schedule.use_farsite_curve = 1
fmc_schedule.start_doy = 200
fmc_schedule.spring_start = 90
fmc_schedule.summer_peak = 150
fmc_schedule.fall_start  = 240
fmc_schedule.fall_end    = 300
fmc_schedule.fmc_min = 85.0
fmc_schedule.fmc_max = 140.0

# Level-set propagation
propagation_method = levelset
