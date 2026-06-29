# Wind solver inputs for one-way coupling test
# Simple flat terrain for fire coupling tests

dx = 31.25
dy = 31.25
dz = 30.0
domain_height = 300.0

# Reference wind
U_ref = 10.0
V_ref = 0.0
z_ref = 10.0

# Roughness
z0 = 0.1

# Solver settings
alpha_h = 1.0
alpha_v = 1.0
tol_rel = 1.e-8
max_iter = 200
