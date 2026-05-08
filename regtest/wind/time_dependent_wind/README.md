# Time-Dependent Wind Field Test

## Purpose
Tests fire spread with a wind field that changes direction over time. The simulation uses linear temporal interpolation between two wind snapshots to update the velocity field every time step, exercising the `update_time_dependent_velocity` code path.

## Configuration
- **Domain**: 1000 m × 1000 m, flat terrain
- **Grid**: 50 × 50 cells (20 m resolution)
- **Wind**:
  - Snapshot 0 (`time_wind.csv`, t = 0 s): uniform 5 m/s eastward (u = 5.0, v = 0.0)
  - Snapshot 1 (`time_wind_1.csv`, t = 3600 s): uniform 5 m/s northeast (u ≈ 3.54, v ≈ 3.54)
  - `wind_time_spacing = 3600.0` s – both snapshots cover the full 50-step run
- **Initial condition**: Line fire box at x = 100–130 m, y = 400–600 m (western edge)
- **Time steps**: 50
- **Output interval**: Every 10 steps
- **Fuel**: FM4 Chaparral (NFFL Model 4)
- **FARSITE**: Enabled with Anderson (1983) dynamic L/W ratio

## Wind Model
Two uniform wind snapshots are linearly interpolated in time:
```
u(t) = u_0 + (t / 3600) * (u_1 - u_0)
v(t) = v_0 + (t / 3600) * (v_1 - v_0)
```
where `(u_0, v_0) = (5.0, 0.0)` m/s and `(u_1, v_1) ≈ (3.54, 3.54)` m/s.

## Expected Behaviour
- Fire starts as a line at the western edge of the domain.
- Wind initially drives the fire eastward.
- Over the simulation the wind gradually rotates toward the northeast,
  causing the fire front to shift its direction of maximum spread.
- The elliptical FARSITE spread pattern reflects the changing wind direction.

## Files
- `inputs.i`: Input parameters
- `time_wind.csv`: Initial wind snapshot (t = 0 s), 11 × 11 uniform grid
- `time_wind_1.csv`: Second wind snapshot (t = 3600 s), 11 × 11 uniform grid

## Run Command
```bash
cd regtest/time_dependent_wind
../../build/levelset inputs.i
```

## Notes
- This test requires a 2D build (`-DLEVELSET_DIM_2D=ON`).
- The `update_time_dependent_velocity` function loads wind snapshots on demand
  and interpolates spatially via inverse-distance weighting and temporally via
  linear interpolation.
- The wind field files contain 121 uniformly-distributed sample points that are
  sufficient for IDW interpolation over the 1 km × 1 km domain.
