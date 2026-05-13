#!/usr/bin/env python3
"""
create_detections.py
====================
Generates fire_detections.csv for the satellite_assimilation regression test.

The CSV contains three synthetic active-fire detections (in lon/lat/confidence
format) that map into the simulation domain [0, 20000] × [0, 20000] m when the
UTM projection parameters match those in inputs.i:

    utm_zone             = 10  (central meridian 123°W)
    utm_northern         = 1
    prob_lo_easting_m    = 570000.0 m
    prob_lo_northing_m   = 4100000.0 m

The three detections are placed at:
    Detection 1: sim (3000,  5000) m  → (easting 573000, northing 4105000)
    Detection 2: sim (10000, 10000) m → (easting 580000, northing 4110000)
    Detection 3: sim (16000, 14000) m → (easting 586000, northing 4114000)

We compute the approximate lon/lat from UTM using the inverse projection so
the values round-trip correctly through the forward UTM projection in
satellite_assimilation.H.

Usage
-----
    python3 create_detections.py

Output
------
    fire_detections.csv   (lon_deg, lat_deg, confidence_pct)
"""

import math

# UTM inverse projection (approximate, accurate to < 1 m within a zone)
def utm_to_lonlat(easting, northing, zone, northern=True):
    """Return (lon_deg, lat_deg) from UTM easting/northing."""
    # WGS-84 constants
    a  = 6378137.0
    f  = 1.0 / 298.257223563
    k0 = 0.9996
    E0 = 500000.0
    N0 = 0.0 if northern else 10000000.0

    e2   = 2*f - f**2
    ep2  = e2 / (1 - e2)
    n    = f / (2 - f)
    A    = a / (1 + n) * (1 + n**2/4 + n**4/64)

    xi   = (northing - N0) / (k0 * A)
    eta  = (easting  - E0) / (k0 * A)

    # Series coefficients (Karney 2011, 3rd-order)
    b1 =  n/2   - 2*n**2/3 + 37*n**3/96
    b2 =          n**2/48  + n**3/15
    b3 =                     17*n**3/480

    xi_p  = xi  - b1*math.sin(2*xi)*math.cosh(2*eta) \
                - b2*math.sin(4*xi)*math.cosh(4*eta) \
                - b3*math.sin(6*xi)*math.cosh(6*eta)
    eta_p = eta - b1*math.cos(2*xi)*math.sinh(2*eta) \
                - b2*math.cos(4*xi)*math.sinh(4*eta) \
                - b3*math.cos(6*xi)*math.sinh(6*eta)

    chi  = math.asin(math.sin(xi_p) / math.cosh(eta_p))
    lon0 = math.radians((zone - 1)*6 - 180 + 3)

    dlon = math.atan2(math.sinh(eta_p), math.cos(xi_p))

    # Iterate conformal lat → geodetic lat
    e = math.sqrt(e2)
    lat = chi
    for _ in range(3):
        sin_lat = math.sin(lat)
        lat = 2*math.atan(
            math.tan(math.pi/4 + chi/2) *
            ((1 + e*sin_lat) / (1 - e*sin_lat))**(e/2)
        ) - math.pi/2

    lon_deg = math.degrees(lon0 + dlon)
    lat_deg = math.degrees(lat)
    return lon_deg, lat_deg


# Parameters from inputs.i
zone        = 10
northern    = True
elo         = 570000.0   # prob_lo_easting_m
nlo         = 4100000.0  # prob_lo_northing_m

# Simulation-domain coordinates of the three detections [m]
detections_sim = [
    (3000,  5000,  85),   # sim_x, sim_y, confidence%
    (10000, 10000, 75),
    (16000, 14000, 90),
]

rows = []
for (sx, sy, conf) in detections_sim:
    east   = elo + sx
    north  = nlo + sy
    lon, lat = utm_to_lonlat(east, north, zone, northern)
    rows.append((lon, lat, conf))

output_file = "fire_detections.csv"
with open(output_file, "w") as f:
    f.write("# lon_deg,lat_deg,confidence_pct\n")
    f.write("# Synthetic VIIRS-style fire detections for satellite_assimilation regtest\n")
    for (lon, lat, conf) in rows:
        f.write(f"{lon:.6f},{lat:.6f},{conf}\n")

print(f"Written {len(rows)} detections to {output_file}")
for (lon, lat, conf) in rows:
    print(f"  lon={lon:.6f}  lat={lat:.6f}  conf={conf}%")
