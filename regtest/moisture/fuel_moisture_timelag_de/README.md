# Fine Fuel Moisture Time-Lag Differential Equations Test

## Purpose

Tests the physically-based fuel moisture dynamics using time-lag differential equations (Nelson 2000, Viney 1991).

## Physics

Fuel moisture responds to equilibrium moisture content (EMC) with a time lag that depends on fuel size class. The governing equation is:

```
dM/dt = (M_e - M) / τ + P(t)
```

where:
- M = fuel moisture content [fraction]
- M_e = equilibrium moisture content [fraction]
- τ = time-lag constant [hours]
- P(t) = precipitation wetting term [1/hour]

### Equilibrium Moisture Content

EMC depends on relative humidity with separate curves for adsorption (wetting) and desorption (drying):

**Adsorption (Nelson 2000):**
```
M_e = 0.03229 + 0.2810×H + 0.4093×H² - 1.3560×H³ + 1.6596×H⁴
```

**Desorption (Nelson 2000):**
```
M_e = 0.05800 + 0.1985×H + 0.6250×H² - 1.1830×H³ + 1.0570×H⁴
```

where H = RH/100

### Test Conditions

- Initial moisture: 12% (dry)
- RH = 40% → EMC ≈ 8-10% (desorption)
- Temperature = 25°C (warm → faster drying)
- Time-lag constants: 1hr, 10hr, 100hr

### Expected Results

Over time, moisture should:
1. 1-hr fuels: equilibrate quickly (within ~3-5 hours)
2. 10-hr fuels: equilibrate moderately (within ~30-50 hours)  
3. 100-hr fuels: equilibrate slowly (within ~300-500 hours)

Final equilibrium moisture ≈ 8-10% for RH = 40%.

## References

- Nelson, R.M. (2000). "Prediction of diurnal change in 10-h fuel stick moisture content." Canadian Journal of Forest Research, 30, 1071-1087.
- Viney, N.R. (1991). "A review of fine fuel moisture modelling." International Journal of Wildland Fire, 1(4), 215-234.
- Van Wagner, C.E. (1972). "Equilibrium moisture contents of some fine forest fuels in eastern Canada." Canadian Forestry Service Report 7.

## Validation

- Verify exponential approach to equilibrium moisture
- Check 1-hr fuels respond much faster than 100-hr fuels
- Confirm hysteresis (different adsorption/desorption curves)
- Validate temperature correction accelerates drying
- Test precipitation wetting when rain_rate > 0
