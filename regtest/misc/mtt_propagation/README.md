# MTT (Minimum Travel Time) Propagation Regression Test

## Purpose

Tests the **Minimum Travel Time (MTT)** propagation method
(`propagation_method = mtt`).

MTT pre-computes the fire arrival time for every grid cell once at startup
using a Dijkstra fast-marching sweep over the selected ROS field.  The
level-set field phi is then updated at every time step as:

    phi(i,j,k) = arrival_time(i,j,k) - current_time

giving phi < 0 (burned), phi = 0 (fire front), phi > 0 (unburned).

## Configuration

| Parameter | Value |
|-----------|-------|
| Domain | 500 m × 500 m × 500 m |
| Grid | 64³ cells (≈ 7.8 m resolution) |
| Wind | 3 m/s eastward |
| Fuel | FM4 (Southern California chaparral) |
| Propagation | `mtt` |
| Simulation time | 600 s |

## Expected Behaviour

- Fire grows roughly elliptically downwind (eastward).
- Unlike the level-set method, the fire front does not need reinitialisation.
- Arrival-time field (`arrival_time`) increases monotonically from the
  ignition point.
- The phi field is a smooth signed-distance-like field derived analytically
  from the arrival times.

## Run Command

```bash
./build/levelset regtest/mtt_propagation/inputs.i
```

## Reference

Finney, M.A. (2002). "Fire growth using minimum travel time methods."
*Canadian Journal of Forest Research*, 32(8), 1420–1424.
https://doi.org/10.1139/x02-068
