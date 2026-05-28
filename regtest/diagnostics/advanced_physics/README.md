Advanced Physics Features Regression Test
==========================================

This test verifies that all 10 newly implemented physics feature headers
can be included and compiled successfully. It runs a simple fire simulation
with the basic Rothermel model.

Features Tested
---------------

1. Radiation-driven preheating distance (radiation_preheating.H)
2. Fuel particle temperature evolution (fuel_temperature.H)
3. Fire line intensity rate of change (intensity_rate_of_change.H)
4. Flame residence time from fire intensity (intensity_residence_time.H)
5. Wind-fuel interaction feedback (wind_fuel_interaction.H)
6. Fuel moisture of extinction gradient (moisture_of_extinction.H)
7. Flame intermittency factor (flame_intermittency.H)
8. Plume entrainment feedback (plume_momentum_feedback.H)
9. Fuel bed bulk density spatial variation (fuel_loading_variation.H)
10. Critical heat flux for ignition (critical_heat_flux.H)

Test Configuration
------------------

- Simple 1 km x 1 km domain with 64x64 cells
- Point ignition at center
- Uniform short grass fuel (Anderson FM1)
- Constant 5 m/s easterly wind
- 300 second simulation (50 time steps)

Expected Behavior
-----------------

The test should:
- Compile successfully with all new headers included
- Run without errors
- Produce basic fire spread output
- Generate plotfiles at t=0s, 150s, 300s

Success Criteria
----------------

- Exit code 0
- No compilation warnings from new headers
- At least 2 plotfiles generated
- Fire spreads eastward (downwind direction)

Notes
-----

This is a basic compilation and integration test. Individual feature
functionality tests would require more specific input configurations
to exercise each physics model separately.
