#!/usr/bin/env python3
"""
behavior_matrix.py – BehavePlus-style fuel condition fire behavior matrix.

Computes Rothermel (1972) fire behavior metrics across a grid of wind speeds
and dead fuel moisture values for a chosen fuel model (Anderson 13 or
Scott & Burgan 40).

Outputs a CSV matrix suitable for calibration, sensitivity analysis, and
field briefings — analogous to the tabular output produced by BehavePlus.

Metrics computed per (wind, moisture) combination
--------------------------------------------------
  R_ros   – Rothermel rate of spread          [m/min]
  I_R     – Reaction intensity                [kW/m²]
  I_B     – Byram fireline intensity          [kW/m]
  L_f     – Byram flame length               [m]
  phi_w   – Wind factor                      [-]
  phi_s   – Slope factor (at user slope)     [-]

Requirements
------------
  None (pure Python).  Optional ``numpy`` speeds up batch computation.
  Optional ``matplotlib`` for heatmap plots.

Usage
-----
  # Matrix for Anderson FM 4 (chaparral), winds 0–10 m/s, moisture 4–20%
  python3 tools/behavior_matrix.py \\
      --fuel-model 4 --fuel-system 13 \\
      --wind-min 0 --wind-max 10 --wind-steps 11 \\
      --moisture-min 0.04 --moisture-max 0.20 --moisture-steps 9 \\
      --slope 0.3 \\
      --out fm4_matrix.csv

  # Heatmap plot as well
  python3 tools/behavior_matrix.py --fuel-model 4 --plot

  # Scott & Burgan fuel model 145 (tall grass)
  python3 tools/behavior_matrix.py --fuel-model 145 --fuel-system 40 \\
      --out fm145_matrix.csv

References
----------
  Rothermel, R.C. (1972). A mathematical model for predicting fire spread
    in wildland fuels. USDA Forest Service Research Paper INT-115.
  Andrews, P.L. (2018). The Rothermel surface fire spread model and
    associated developments. USDA Forest Service General Technical Report
    RMRS-GTR-371.
  BehavePlus: https://www.firelab.org/project/behaveplusfiremodeling
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from typing import Dict, List, Optional, Tuple

# ---------------------------------------------------------------------------
# Fuel model database (Anderson 13 + Scott & Burgan 40, SI-converted)
# All fuel loads in lb/ft², SAV in ft⁻¹, depth in ft.
# ---------------------------------------------------------------------------

# Anderson (1982) FBFM13 fuel models
# Format: (code, name, w0[lb/ft²], sigma[ft⁻¹], delta[ft], M_x[-])
_FBFM13: Dict[int, Tuple] = {
    1:  (1,  "Short grass",           0.034, 3500, 1.0, 0.12),
    2:  (2,  "Timber grass/shrub",    0.092, 2784, 1.0, 0.15),
    3:  (3,  "Tall grass",            0.138, 1500, 2.5, 0.25),
    4:  (4,  "Chaparral",             0.230, 1739, 6.0, 0.20),
    5:  (5,  "Brush",                 0.046, 1683, 2.0, 0.20),
    6:  (6,  "Dormant brush",         0.069, 1564, 2.5, 0.25),
    7:  (7,  "Southern rough",        0.052, 1552, 2.5, 0.40),
    8:  (8,  "Compact timber litter", 0.069, 1889, 0.2, 0.30),
    9:  (9,  "Hardwood litter",       0.134, 2484, 0.2, 0.25),
    10: (10, "Timber (understory)",   0.138, 1764, 1.0, 0.25),
    11: (11, "Light slash",           0.069, 1182, 1.0, 0.15),
    12: (12, "Medium slash",          0.184, 1145, 2.3, 0.20),
    13: (13, "Heavy slash",           0.322, 1159, 3.0, 0.25),
}

# Scott & Burgan (2005) FBFM40 – representative subset
# (w0 = total oven-dry surface fuel load lb/ft²; sigma = char. SAV ft⁻¹)
_FBFM40: Dict[int, Tuple] = {
    # NB (non-burnable)
    91: (91,  "Urban/Developed",          0.000, 100,  0.0, 0.15),
    92: (92,  "Snow/Ice",                 0.000, 100,  0.0, 0.15),
    93: (93,  "Agriculture",              0.000, 100,  0.0, 0.15),
    98: (98,  "Open Water",               0.000, 100,  0.0, 0.15),
    99: (99,  "Bare Ground",              0.000, 100,  0.0, 0.15),
    # GR (grass)
    101:(101, "GR1 Short sparse dry grass",    0.011, 2200, 0.4, 0.15),
    102:(102, "GR2 Low grow dry grass",        0.046, 2000, 1.0, 0.15),
    103:(103, "GR3 Low grass",                 0.057, 1500, 2.0, 0.30),
    104:(104, "GR4 Mod. load dry grass",       0.087, 2000, 2.0, 0.15),
    105:(105, "GR5 Low load moist grass",      0.092, 1800, 1.5, 0.40),
    106:(106, "GR6 Mod. load humid grass",     0.115, 2200, 1.5, 0.40),
    107:(107, "GR7 High load dry grass",       0.230, 2000, 3.0, 0.15),
    108:(108, "GR8 High load moist grass",     0.299, 1500, 4.0, 0.30),
    109:(109, "GR9 Very high load dry grass",  0.414, 1800, 5.0, 0.40),
    # GS (grass-shrub)
    121:(121, "GS1 Low load dry grass-shrub",  0.057, 2000, 0.9, 0.15),
    122:(122, "GS2 Mod. load dry grass-shrub", 0.092, 1800, 1.5, 0.15),
    123:(123, "GS3 Mod. load humid grass-shrub",0.138,1800, 1.8, 0.40),
    124:(124, "GS4 High load humid grass-shrub",0.207,1800, 2.1, 0.40),
    # SH (shrub)
    141:(141, "SH1 Low load dry shrub",        0.046, 1600, 1.0, 0.15),
    142:(142, "SH2 Mod. load dry shrub",       0.161, 1600, 1.0, 0.15),
    143:(143, "SH3 Mod. load humid shrub",     0.069, 1600, 2.4, 0.40),
    144:(144, "SH4 Low load humid shrub",      0.092, 2000, 3.0, 0.30),
    145:(145, "SH5 High load dry shrub",       0.299, 750,  6.0, 0.15),
    146:(146, "SH6 Low load humid shrub",      0.115, 1600, 2.0, 0.30),
    147:(147, "SH7 Very high load dry shrub",  0.460, 750,  6.0, 0.15),
    148:(148, "SH8 High load humid shrub",     0.184, 750,  3.0, 0.40),
    149:(149, "SH9 Very high load humid shrub",0.345, 750,  4.4, 0.40),
    # TU (timber-understory)
    161:(161, "TU1 Low load dry grass-shrub-tim",0.092,2000,0.6,0.20),
    162:(162, "TU2 Mod. load humid grass-shrub", 0.057,2000,1.0,0.30),
    163:(163, "TU3 Mod. load humid grass-shrub", 0.184,1500,1.3,0.30),
    164:(164, "TU4 Dwarf conifer shrub",         0.069,2300,0.5,0.12),
    165:(165, "TU5 High load conifer litter",    0.322,1500,1.0,0.25),
    # TL (timber litter)
    181:(181, "TL1 Low load compact conifer litter",  0.069,2000,0.2,0.30),
    182:(182, "TL2 Low load broadleaf litter",         0.069,2000,0.2,0.25),
    183:(183, "TL3 Mod. load conifer litter",          0.115,2000,0.3,0.20),
    184:(184, "TL4 Small downed logs",                 0.115,2000,0.4,0.25),
    185:(185, "TL5 High load conifer litter",          0.161,1500,0.6,0.25),
    186:(186, "TL6 High load broadleaf litter",        0.184,2000,0.3,0.25),
    187:(187, "TL7 Large downed logs",                 0.230,2000,0.4,0.25),
    188:(188, "TL8 Long-needle litter",                0.299,1800,0.3,0.35),
    189:(189, "TL9 Very high load broadleaf litter",   0.414,1600,0.6,0.35),
    # SB (slash-blowdown)
    201:(201, "SB1 Low load activity fuel",            0.069,2000,1.0,0.25),
    202:(202, "SB2 Mod. load activity/wind",           0.207,2000,1.0,0.25),
    203:(203, "SB3 High load activity/wind/blowdown",  0.414,2000,1.2,0.25),
    204:(204, "SB4 High load blowdown",                0.460,2000,2.7,0.25),
}


def _get_fuel(code: int, system: str) -> Optional[Tuple]:
    db = _FBFM13 if system == "13" else _FBFM40
    return db.get(code)


# ---------------------------------------------------------------------------
# Rothermel (1972) single-class model – pure Python
# ---------------------------------------------------------------------------

def _rothermel(
    w0: float,       # oven-dry fuel load [lb/ft²]
    sigma: float,    # surface-area-to-volume ratio [ft⁻¹]
    delta: float,    # fuel bed depth [ft]
    M_f: float,      # fuel moisture [fraction]
    M_x: float,      # moisture of extinction [fraction]
    h_heat: float,   # heat content [BTU/lb]  default 8000
    S_T: float,      # total mineral content  default 0.0555
    S_e: float,      # effective mineral content default 0.010
    rho_p: float,    # particle density [lb/ft³]  default 32.0
    U_ftmin: float,  # mid-flame wind speed [ft/min]
    slope_tan: float,  # tan(slope angle)
) -> Dict[str, float]:
    """Rothermel (1972) single-class surface fire spread model.

    Returns a dict with keys: R0, R_ros, I_R, I_B, L_f, phi_w, phi_s.
    All outputs are in the Rothermel SI-mixed system:
      R [ft/min], I_R [BTU/ft²/min], I_B [BTU/ft/s], L_f [ft].
    """
    # Avoid division by zero
    if w0 <= 0.0 or sigma <= 0.0:
        return dict(R0=0.0, R_ros=0.0, I_R=0.0, I_B=0.0, L_f=0.0,
                    phi_w=0.0, phi_s=0.0)

    # Bulk density
    rho_b = w0 / delta
    beta = rho_b / rho_p
    beta_op = 3.348 * sigma**(-0.8189)

    # Optimum reaction velocity [1/min]
    Gamma_max = (sigma**1.5) / (495.0 + 0.0594 * sigma**1.5)
    A = 133.0 / sigma**0.7913
    Gamma_prime = Gamma_max * (beta / beta_op)**A * math.exp(A * (1.0 - beta / beta_op))

    # Moisture damping
    if M_x <= 0.0:
        eta_M = 0.0
    else:
        r_M = min(M_f / M_x, 1.0)
        eta_M = 1.0 - 2.59 * r_M + 5.11 * r_M**2 - 3.52 * r_M**3

    # Mineral damping
    eta_s = 0.174 * S_e**(-0.19)

    # Net fuel load
    w_n = w0 * (1.0 - S_T)

    # Reaction intensity [BTU/ft²/min]
    I_R = Gamma_prime * w_n * h_heat * eta_M * eta_s
    I_R = max(I_R, 0.0)

    # Propagating flux ratio
    xi = math.exp((0.792 + 0.681 * sigma**0.5) * (beta + 0.1)) / (192.0 + 0.2595 * sigma)

    # Heat of preignition [BTU/lb]
    eps_h = math.exp(-138.0 / sigma)
    Q_ig = 250.0 + 1116.0 * M_f

    # No-wind, no-slope ROS [ft/min]
    denom = rho_b * eps_h * Q_ig
    if denom <= 0.0:
        return dict(R0=0.0, R_ros=0.0, I_R=I_R, I_B=0.0, L_f=0.0,
                    phi_w=0.0, phi_s=0.0)
    R0 = I_R * xi / denom

    # Wind factor phi_w
    C = 7.47 * math.exp(-0.133 * sigma**0.55)
    B = 0.02526 * sigma**0.54
    E = 0.715 * math.exp(-3.59e-4 * sigma)
    beta_ratio = beta / beta_op
    phi_w = C * (U_ftmin**B) * beta_ratio**(-E) if U_ftmin > 0.0 else 0.0

    # Slope factor phi_s
    phi_s = 5.275 * beta**(-0.3) * slope_tan**2 if slope_tan > 0.0 else 0.0

    # Final ROS [ft/min]
    R_ros = R0 * (1.0 + phi_w + phi_s)

    # Byram fireline intensity [BTU/ft/s]
    # I_B = H * w_a * R  (Byram 1959; H in BTU/lb, w_a in lb/ft², R in ft/s)
    I_B = h_heat * w_n * (R_ros / 60.0)

    # Flame length [ft]  (Byram 1959)
    L_f_ft = 0.45 * I_B**0.46 if I_B > 0.0 else 0.0

    return dict(
        R0=R0,
        R_ros=R_ros,
        I_R=I_R,
        I_B=I_B,
        L_f=L_f_ft,
        phi_w=phi_w,
        phi_s=phi_s,
    )


# ---------------------------------------------------------------------------
# Unit conversions
# ---------------------------------------------------------------------------

_FT_MIN_TO_M_MIN = 0.3048
_BTU_FT2_MIN_TO_KW_M2 = 0.18941        # 1 BTU/ft²/min = 0.18941 kW/m²
_BTU_FT_S_TO_KW_M = 0.34879            # 1 BTU/ft/s = 0.34879 kW/m
_FT_TO_M = 0.3048
_M_S_TO_FT_MIN = 196.85


# ---------------------------------------------------------------------------
# Matrix computation
# ---------------------------------------------------------------------------

def compute_matrix(
    fuel_code: int,
    fuel_system: str,
    wind_speeds_ms: List[float],   # [m/s]
    moistures: List[float],         # [fraction]
    slope_tan: float = 0.0,
    h_heat: float = 8000.0,
    S_T: float = 0.0555,
    S_e: float = 0.010,
    rho_p: float = 32.0,
) -> List[Dict]:
    """Compute Rothermel fire behavior matrix.

    Returns a list of dicts, one per (wind, moisture) combination.
    """
    fuel = _get_fuel(fuel_code, fuel_system)
    if fuel is None:
        raise ValueError(
            f"Fuel model {fuel_code} not found in FBFM{fuel_system} database."
        )
    _code, name, w0, sigma, delta, M_x = fuel

    rows = []
    for U_ms in wind_speeds_ms:
        for M_f in moistures:
            U_ftmin = U_ms * _M_S_TO_FT_MIN
            res = _rothermel(w0, sigma, delta, M_f, M_x, h_heat, S_T, S_e, rho_p,
                             U_ftmin, slope_tan)
            rows.append({
                "fuel_code":     fuel_code,
                "fuel_name":     name,
                "wind_m_s":      round(U_ms, 3),
                "moisture_pct":  round(M_f * 100.0, 1),
                "slope_pct":     round(math.degrees(math.atan(slope_tan)), 1),
                "R_ros_m_min":   round(res["R_ros"] * _FT_MIN_TO_M_MIN, 4),
                "I_R_kW_m2":     round(res["I_R"] * _BTU_FT2_MIN_TO_KW_M2, 3),
                "I_B_kW_m":      round(res["I_B"] * _BTU_FT_S_TO_KW_M, 3),
                "L_f_m":         round(res["L_f"] * _FT_TO_M, 3),
                "phi_w":         round(res["phi_w"], 4),
                "phi_s":         round(res["phi_s"], 4),
            })
    return rows


# ---------------------------------------------------------------------------
# CSV writer
# ---------------------------------------------------------------------------

def write_matrix_csv(rows: List[Dict], out_path: str) -> None:
    if not rows:
        print("(no rows to write)")
        return
    fieldnames = list(rows[0].keys())
    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} rows → {out_path}")


# ---------------------------------------------------------------------------
# ASCII table
# ---------------------------------------------------------------------------

def print_matrix_table(rows: List[Dict], metric: str = "R_ros_m_min") -> None:
    """Print a pivot table: rows = wind speeds, columns = moisture values."""
    if not rows:
        return

    winds = sorted({r["wind_m_s"] for r in rows})
    moist = sorted({r["moisture_pct"] for r in rows})

    label_map = {
        "R_ros_m_min": "ROS [m/min]",
        "I_B_kW_m":    "IB [kW/m]",
        "L_f_m":       "L_f [m]",
        "I_R_kW_m2":   "IR [kW/m²]",
    }
    title = label_map.get(metric, metric)

    # Build lookup
    lookup = {(r["wind_m_s"], r["moisture_pct"]): r[metric] for r in rows}

    col_w = 9
    header = f"{'Wind [m/s]':>12s}" + "".join(f"  MC={m:4.1f}%" for m in moist)
    print(f"\n{title}")
    print("-" * len(header))
    print(header)
    print("-" * len(header))
    for U in winds:
        row_str = f"{U:12.2f}"
        for m in moist:
            val = lookup.get((U, m), float("nan"))
            row_str += f"  {val:>9.3f}"
        print(row_str)
    print()


# ---------------------------------------------------------------------------
# Optional heatmap
# ---------------------------------------------------------------------------

def make_heatmaps(rows: List[Dict], out_path: str) -> None:
    """Save heatmaps of ROS, IB, L_f as a PNG figure."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        print("WARNING: matplotlib/numpy not available – skipping heatmap.", file=sys.stderr)
        return

    winds = sorted({r["wind_m_s"] for r in rows})
    moist = sorted({r["moisture_pct"] for r in rows})

    metrics = [
        ("R_ros_m_min", "ROS [m/min]",  "Reds"),
        ("I_B_kW_m",    "IB [kW/m]",   "YlOrRd"),
        ("L_f_m",       "L_f [m]",     "OrRd"),
    ]

    fig, axes = plt.subplots(1, len(metrics), figsize=(5 * len(metrics), 4))
    if len(metrics) == 1:
        axes = [axes]

    for ax, (metric, label, cmap) in zip(axes, metrics):
        data = np.array([[
            next((r[metric] for r in rows if r["wind_m_s"]==U and r["moisture_pct"]==m), 0)
            for m in moist] for U in winds])
        im = ax.imshow(data, aspect="auto", origin="lower", cmap=cmap,
                       extent=[min(moist)-0.5, max(moist)+0.5,
                               min(winds)-0.25, max(winds)+0.25])
        ax.set_xlabel("Dead fuel moisture [%]")
        ax.set_ylabel("Wind speed [m/s]")
        ax.set_title(label)
        fig.colorbar(im, ax=ax, fraction=0.04)

    fuel_name = rows[0].get("fuel_name", "") if rows else ""
    fig.suptitle(f"Fire behavior matrix – FM {rows[0]['fuel_code']} {fuel_name}" if rows else "")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"Saved heatmap → {out_path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv=None):
    parser = argparse.ArgumentParser(
        description="BehavePlus-style Rothermel fire behavior matrix.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--fuel-model",  type=int, default=4,
                        help="Fuel model code (default: 4)")
    parser.add_argument("--fuel-system", default="13", choices=["13", "40"],
                        help="Fuel model system: 13 (FBFM13) or 40 (FBFM40) (default: 13)")
    parser.add_argument("--wind-min",    type=float, default=0.0,
                        help="Minimum wind speed [m/s] (default: 0.0)")
    parser.add_argument("--wind-max",    type=float, default=10.0,
                        help="Maximum wind speed [m/s] (default: 10.0)")
    parser.add_argument("--wind-steps",  type=int,   default=11,
                        help="Number of wind speed steps (default: 11)")
    parser.add_argument("--moisture-min", type=float, default=0.04,
                        help="Min dead fuel moisture fraction (default: 0.04)")
    parser.add_argument("--moisture-max", type=float, default=0.30,
                        help="Max dead fuel moisture fraction (default: 0.30)")
    parser.add_argument("--moisture-steps", type=int, default=9,
                        help="Number of moisture steps (default: 9)")
    parser.add_argument("--slope",       type=float, default=0.0,
                        help="Slope [tan(angle)] for phi_s (default: 0.0)")
    parser.add_argument("--out",         default="behavior_matrix.csv",
                        help="Output CSV path (default: behavior_matrix.csv)")
    parser.add_argument("--plot",        action="store_true",
                        help="Save heatmap PNG (requires matplotlib + numpy)")
    parser.add_argument("--plot-out",    default="behavior_matrix_heatmap.png",
                        help="Heatmap output path (default: behavior_matrix_heatmap.png)")
    parser.add_argument("--list-fuels",  action="store_true",
                        help="Print available fuel model codes and exit.")

    args = parser.parse_args(argv)

    if args.list_fuels:
        db = _FBFM13 if args.fuel_system == "13" else _FBFM40
        print(f"\nFBFM{args.fuel_system} fuel models:")
        print(f"  {'Code':>6}  Name")
        print(f"  {'----':>6}  ----")
        for code, entry in sorted(db.items()):
            print(f"  {code:>6}  {entry[1]}")
        return

    # Build parameter grids
    def linspace(lo, hi, n):
        if n == 1:
            return [lo]
        return [lo + (hi - lo) * i / (n - 1) for i in range(n)]

    winds   = linspace(args.wind_min, args.wind_max, args.wind_steps)
    moist   = linspace(args.moisture_min, args.moisture_max, args.moisture_steps)

    fuel = _get_fuel(args.fuel_model, args.fuel_system)
    if fuel is None:
        print(f"ERROR: fuel model {args.fuel_model} not found in FBFM{args.fuel_system}.",
              file=sys.stderr)
        sys.exit(1)
    print(f"Fuel model: FM{args.fuel_model} – {fuel[1]}")
    print(f"Wind speeds: {len(winds)} steps  [{args.wind_min:.2f} – {args.wind_max:.2f} m/s]")
    print(f"Moisture:    {len(moist)} steps  [{args.moisture_min*100:.1f}% – {args.moisture_max*100:.1f}%]")
    if args.slope > 0:
        print(f"Slope:       {math.degrees(math.atan(args.slope)):.1f}°  (tan = {args.slope:.3f})")

    rows = compute_matrix(args.fuel_model, args.fuel_system, winds, moist, args.slope)

    print_matrix_table(rows, "R_ros_m_min")
    print_matrix_table(rows, "I_B_kW_m")

    write_matrix_csv(rows, args.out)

    if args.plot:
        make_heatmaps(rows, args.plot_out)


if __name__ == "__main__":
    main()
