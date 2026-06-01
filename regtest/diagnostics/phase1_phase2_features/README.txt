Regression Test: Phase 1 & Phase 2 Operational Features
========================================================

This test validates the new operational fire behavior and fire danger features
added in the Phase 1 and Phase 2 implementation:

Phase 1 Features:
-----------------
1. Grass Curing Model (grass_curing_model.H)
   - Seasonal sinusoidal curing cycle
   - Moisture-dependent curing (KBDI-driven)
   - Growing degree day (GDD) phenology model

2. Diurnal Weather Cycles (diurnal_weather.H)
   - Daily temperature variation (sinusoidal)
   - Daily RH variation (anti-phase with temperature)
   - Daily wind speed variation

3. Elevation Lapse Rate (elevation_lapse_rate.H)
   - Temperature decrease with elevation (~6.5°C/km)
   - RH increase with elevation
   - Barometric pressure calculation

Phase 2 Features:
-----------------
4. NFDRS Spread Component (nfdrs_spread_component.H)
   - U.S. National Fire Danger Rating System spread rating
   - Wind, slope, and fuel moisture effects
   - Burning Index calculation (SC × ERC / 10)

5. Chandler Burning Index (chandler_burning_index.H)
   - NWS fire danger index
   - Red Flag Warning criteria
   - Fosberg Fire Weather Index variant

Test Scenario:
--------------
- Fuel type: Anderson FM1 (short grass)
- Weather: Hot, dry, windy (35°C, 15% RH, 5 m/s wind)
- Terrain: Flat domain (2 km × 2 km)
- Duration: 1 hour simulation
- Expected: High CBI (~90+), high NFDRS SC (~40+)

Integration Status:
-------------------
The header files are complete and GPU-compatible. Full integration into
the main simulation loop (main.cpp) requires:
1. Adding input parameters to parse_inputs.H
2. Creating MultiFabs for spatial fields
3. Calling update functions in time-stepping loop
4. Writing diagnostic outputs to plotfiles

This regression test serves as a template for future full integration.

References:
-----------
See docs/new_features.rst for complete mathematical formulations and usage.
See docs/references.rst for literature citations.
