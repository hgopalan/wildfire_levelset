#!/usr/bin/env python3
"""
create_terrain.py – Generate a synthetic ridge-canyon-ridge terrain for the
solar_horizon_shading regression test.

Writes ``canyon_terrain.csv`` (X Y Z format) to the current working directory.

Terrain geometry
----------------
The terrain models a N-S running canyon flanked by E and W ridges:

  W ridge (x ≈ 330 000 m UTM):  elevation ~200 m
  Canyon floor (x ≈ 330 500 m): elevation  ~50 m
  E ridge (x ≈ 331 000 m UTM):  elevation ~200 m

A gentle N-S gradient (20 m across the domain) ensures that slope and aspect
fields are non-trivial everywhere.

Topographic horizon from the canyon floor
-----------------------------------------
Horizontal distance to each ridge crest = 500 m.
Elevation difference = 200 − 50 = 150 m.
Nominal E/W horizon angle ≈ atan2(150, 500) ≈ 16.7°.

At 06:30 AM PDT on July 1, 2024 (lat 34.1° N) the solar elevation is roughly
10°, which is below the ridge horizon.  The FARSITE topographic horizon scan
therefore shadows the canyon floor, whereas local surface-normal shading
(cos_i > 0 on the nearly flat valley floor) would leave it unshaded.
"""

import math
import os
import sys

# ---------------------------------------------------------------------------
# Domain parameters – must match inputs.i
# ---------------------------------------------------------------------------
PROB_LO_X = 330_000.0
PROB_HI_X = 331_000.0
PROB_LO_Y = 3_775_000.0
PROB_HI_Y = 3_776_000.0

NX = 21   # terrain points in x
NY = 21   # terrain points in y

# Canyon geometry
RIDGE_ELEV_M = 200.0   # ridge-crest elevation [m]
FLOOR_ELEV_M =  50.0   # canyon-floor elevation [m]
NS_GRADIENT_M =  20.0  # elevation rise from S to N edge [m]


def _canyon_elevation(x: float, y: float) -> float:
    """Return terrain elevation [m] at UTM coordinates (x, y).

    Uses a Gaussian cross-valley profile that is high at the E/W edges and
    low at the domain x-centre, plus a linear N-S gradient.
    """
    cx      = 0.5 * (PROB_LO_X + PROB_HI_X)   # domain x-centre
    half_w  = 0.5 * (PROB_HI_X - PROB_LO_X)   # half-width

    # Normalised distance from x-centre in [0, 1]
    norm_x = abs(x - cx) / half_w

    # Gaussian ridge profile: maximum at edges (norm_x = 1), minimum at centre
    # Sigma chosen so the ridge is sharp near the walls.
    sigma = 0.35
    ridge_frac = math.exp(-((1.0 - norm_x) ** 2) / (2.0 * sigma ** 2))
    elev = FLOOR_ELEV_M + (RIDGE_ELEV_M - FLOOR_ELEV_M) * ridge_frac

    # Linear N-S gradient
    norm_y = (y - PROB_LO_Y) / (PROB_HI_Y - PROB_LO_Y)
    elev += NS_GRADIENT_M * norm_y

    return elev


def main(argv=None):
    out_path = os.path.join(os.getcwd(), "canyon_terrain.csv")

    dx = (PROB_HI_X - PROB_LO_X) / (NX - 1)
    dy = (PROB_HI_Y - PROB_LO_Y) / (NY - 1)

    rows = []
    for j in range(NY):
        for i in range(NX):
            x = PROB_LO_X + i * dx
            y = PROB_LO_Y + j * dy
            z = _canyon_elevation(x, y)
            rows.append((x, y, z))

    hz_angle = math.degrees(
        math.atan2(RIDGE_ELEV_M - FLOOR_ELEV_M,
                   0.5 * (PROB_HI_X - PROB_LO_X))
    )

    with open(out_path, "w") as fh:
        fh.write("# Synthetic canyon terrain for solar_horizon_shading regtest\n")
        fh.write("# Format: X[m]  Y[m]  Z[m]  (UTM Zone 11N)\n")
        fh.write(f"# Grid: {NX}×{NY} = {NX * NY} points, "
                 f"dx={dx:.1f} m, dy={dy:.1f} m\n")
        fh.write(f"# W/E ridges at ~{RIDGE_ELEV_M:.0f} m, "
                 f"canyon floor at ~{FLOOR_ELEV_M:.0f} m\n")
        fh.write(f"# Nominal E/W horizon angle from floor: {hz_angle:.1f} deg\n")
        for (x, y, z) in rows:
            fh.write(f"{x:.2f}  {y:.2f}  {z:.2f}\n")

    print(f"Wrote {len(rows)} terrain points → '{out_path}'")
    print(f"Nominal E/W ridge horizon angle from canyon floor: {hz_angle:.1f}°")
    print("At 06:30 AM PDT on Jul 1, 2024 (lat 34.1°N) solar elevation ≈ 10°,")
    print(f"which is below the {hz_angle:.1f}° horizon → canyon floor shaded by horizon scan.")


if __name__ == "__main__":
    main(sys.argv[1:])
