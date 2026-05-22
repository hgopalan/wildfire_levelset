# Fire Acceleration Regression Test

## Purpose

Validates the FARSITE temporal acceleration model (McAlpine & Wakimoto 1991) implementation including:

1. **Size-based model** (Catchpole et al. 1992): Fire acceleration based on current fire size
2. **FARSITE temporal model** (McAlpine & Wakimoto 1991): Per-cell temporal acceleration with VanWagner's equation
3. **Wind-onset time-lag**: Fire ROS response to sudden wind changes

## Test Scenarios

### Test 1: Size-Based Model (`inputs.size_based`)
- Single point ignition
- Constant wind (5 m/s)
- Acceleration enabled with L_acc = 50 m
- Expected: Fire ROS increases as fire grows from ~0 to full ROS

### Test 2: FARSITE Temporal - Point Ignition (`inputs.temporal_point`)
- Single point ignition
- Constant wind (5 m/s)
- Temporal acceleration with A_point = 0.115 1/min
- Expected: Exponential ramp-up: R(t) = R_E × (1 - exp(-0.115 × t/60))

### Test 3: Wind-Onset Time-Lag (`inputs.wind_lag`)
- Single point ignition
- Wind schedule: sudden increase at t=300s
- Wind-lag enabled with tau_wind = 60.0 s
- Expected: ROS ramps up gradually after wind change, not instantly

## References

- McAlpine, R.S. & Wakimoto, R.H. (1991). The acceleration of fire from point source to equilibrium spread. Forest Science, 37(5), 1314–1337.
- Alexander, M.E., Stocks, B.J. & Lawson, B.D. (1992). Fire behavior in Black Spruce-lichen woodland. Info. Rep. NOR-X-310. USDA/CFS.
- Catchpole, E.A., de Mestre, N.J. & Gill, A.M. (1992). Intensity of fire at its perimeter. Australian Journal of Ecology, 17(1), 1–4.
- Finney, M.A. (1998/2004). FARSITE: Fire Area Simulator. USDA Forest Service RMRS-RP-4.
