Fire Behavior Diagnostics Regression Test
===========================================

This test verifies that the five new fire behavior diagnostic fields are computed
correctly:

1. **Byram convective number (Nc)**: Classifies fire as wind-driven vs. plume-dominated
2. **Flame tilt angle**: Flame angle from vertical due to wind/buoyancy
3. **Packing ratio**: Fuel bed compaction metrics (β and β/β_opt)
4. **Flame front depth**: Spatial thickness of burning zone (R × τ_res)
5. **McArthur FFDI**: Australian fire danger index

Test Configuration
------------------
- Domain: 1 km × 1 km (64×64 grid)
- Fuel: Anderson FM4 (chaparral)
- Wind: 5 m/s from west
- Ignition: Point source at domain center
- Duration: 300 seconds (50 steps with CFL=0.3)

Expected Behavior
-----------------
All five new diagnostic fields should be computed and written to plotfiles without
errors. The values should be physically reasonable:

- Convective number: Typically 0.1–10 (wind-driven to plume-dominated)
- Flame tilt: 0–90° (0° = vertical, higher with stronger wind)
- Packing ratio: ~0.001–0.1 for natural fuels
- Relative packing: ~0.5–2.0 (optimal near 1.0)
- Flame depth: ~1–50 m depending on ROS and residence time
- McArthur FFDI: ~5–30 for moderate conditions

Pass Criteria
-------------
1. Simulation completes without errors
2. All diagnostic fields are present in plotfiles
3. Values are non-negative and finite
4. Diagnostic outputs do not affect fire spread (phi, R fields unchanged)
