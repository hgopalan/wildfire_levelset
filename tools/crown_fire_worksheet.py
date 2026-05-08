#!/usr/bin/env python3
"""
crown_fire_worksheet.py – Van Wagner (1977) crown fire assessment worksheet.

Chains the Rothermel (1972) surface fire model with Van Wagner's (1977) crown
fire initiation and spread models to produce a BehavePlus-style crown fire
worksheet.  Given stand and weather inputs it computes:

  1. Surface ROS (Rothermel 1972)            [m/min]
  2. Byram fireline intensity (Byram 1959)   [kW/m]
  3. Van Wagner critical initiation intensity [kW/m]
  4. Crown fire status (surface / passive / active)
  5. Active crown fire ROS (Van Wagner 1977) [m/min]
  6. Crown fraction of total spread          [-]
  7. Total (surface + crown) ROS             [m/min]
  8. Crowning index (critical wind for active crown fire)  [m/s]

The worksheet can be swept over a range of wind speeds and/or moisture values
to produce a tabular output identical in structure to the BehavePlus *Crown
Fire* module output.

Background
----------
Van Wagner's (1977) crown fire initiation condition:

    I_B ≥ I_o    where  I_o [kW/m] = 0.010 × CBH × (460 + 25.9 × FMC)

    CBH – canopy base height [m]
    FMC – foliar moisture content [%]

Critical active crown fire ROS (Van Wagner 1977):

    R'_SA [m/min] = 3.0 / CBD    (CBD in kg/m³)

Crowning index (CI) – the mid-flame wind speed at which the surface fire
is intense enough to initiate crown fire (Andrews 2018 approximation):

    CI is computed by iterating / inverting the Rothermel ROS to find U
    such that I_B(U, MC) = I_o.

Requirements
------------
  None (pure Python).  Optional ``matplotlib`` for heatmap plots.

Usage
-----
  # Default conditions: FM4 chaparral, 4 m canopy
  python3 tools/crown_fire_worksheet.py

  # Custom stand and weather parameters
  python3 tools/crown_fire_worksheet.py \\
      --fuel-model 10 --fuel-system 13 \\
      --cbh 5.0 --cbd 0.15 --fmc 100 \\
      --wind-min 0 --wind-max 15 --wind-steps 16 \\
      --moisture 0.08 \\
      --out crown_worksheet.csv

  # Sweep moisture as well
  python3 tools/crown_fire_worksheet.py \\
      --fuel-model 4 \\
      --moisture-min 0.04 --moisture-max 0.20 --moisture-steps 5 \\
      --wind-min 0 --wind-max 12 --wind-steps 7

  # Save heatmap (crown status)
  python3 tools/crown_fire_worksheet.py --plot

References
----------
  Van Wagner, C.E. (1977). Conditions for the start and spread of crown fire.
    Canadian Journal of Forest Research, 7(1), 23–34.
  Byram, G.M. (1959). Combustion of forest fuels.  In: Davis, K.P. (ed.)
    Forest Fire: Control and Use.  McGraw-Hill, New York, pp. 61–89.
  Rothermel, R.C. (1972). A mathematical model for predicting fire spread in
    wildland fuels.  USDA Forest Service Research Paper INT-115.
  Scott, J.H. & Reinhardt, E.D. (2001). Assessing crown fire potential by
    linking models of surface and crown fire behavior.  USDA Forest Service
    Research Paper RMRS-RP-29.
  Andrews, P.L. (2018). The Rothermel surface fire spread model and
    associated developments.  USDA Forest Service GTR RMRS-GTR-371.
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from typing import Dict, List, Optional, Tuple


# ---------------------------------------------------------------------------
# Rothermel (1972) single-class surface fire model
# (Same implementation as in behavior_matrix.py – kept self-contained.)
# ---------------------------------------------------------------------------

def _rothermel(
    w0: float, sigma: float, delta: float,
    M_f: float, M_x: float,
    h_heat: float = 8000.0, S_T: float = 0.0555,
    S_e: float = 0.010, rho_p: float = 32.0,
    U_ftmin: float = 0.0, slope_tan: float = 0.0,
) -> Dict[str, float]:
    """Rothermel (1972) single-class surface fire spread model.

    All inputs are in the Rothermel BG (British-Gravitational) system:
      w0 [lb/ft²], sigma [ft⁻¹], delta [ft], M_f/M_x [-], h_heat [BTU/lb],
      U_ftmin [ft/min], slope_tan [-].

    Returns a dict with R_ros [ft/min], I_R [BTU/ft²/min],
    I_B [BTU/ft/s], L_f [ft], phi_w [-], phi_s [-].
    """
    if w0 <= 0.0 or sigma <= 0.0:
        return dict(R0=0.0, R_ros=0.0, I_R=0.0, I_B=0.0, L_f=0.0,
                    phi_w=0.0, phi_s=0.0)

    rho_b = w0 / delta
    beta  = rho_b / rho_p
    beta_op = 3.348 * sigma ** (-0.8189)

    Gamma_max  = (sigma ** 1.5) / (495.0 + 0.0594 * sigma ** 1.5)
    A          = 133.0 / sigma ** 0.7913
    Gamma_p    = Gamma_max * (beta / beta_op) ** A * math.exp(A * (1.0 - beta / beta_op))

    r_M    = min(M_f / M_x, 1.0) if M_x > 0.0 else 0.0
    eta_M  = max(0.0, 1.0 - 2.59 * r_M + 5.11 * r_M ** 2 - 3.52 * r_M ** 3)
    eta_s  = 0.174 * S_e ** (-0.19)
    w_n    = w0 * (1.0 - S_T)
    I_R    = max(0.0, Gamma_p * w_n * h_heat * eta_M * eta_s)

    xi    = math.exp((0.792 + 0.681 * sigma ** 0.5) * (beta + 0.1)) / (192.0 + 0.2595 * sigma)
    eps_h = math.exp(-138.0 / sigma)
    Q_ig  = 250.0 + 1116.0 * M_f

    denom = rho_b * eps_h * Q_ig
    if denom <= 0.0:
        return dict(R0=0.0, R_ros=0.0, I_R=I_R, I_B=0.0, L_f=0.0,
                    phi_w=0.0, phi_s=0.0)
    R0 = I_R * xi / denom

    C  = 7.47 * math.exp(-0.133 * sigma ** 0.55)
    B  = 0.02526 * sigma ** 0.54
    E  = 0.715 * math.exp(-3.59e-4 * sigma)
    phi_w = C * (U_ftmin ** B) * (beta / beta_op) ** (-E) if U_ftmin > 0.0 else 0.0
    phi_s = 5.275 * beta ** (-0.3) * slope_tan ** 2 if slope_tan > 0.0 else 0.0
    R_ros = R0 * (1.0 + phi_w + phi_s)

    I_B   = h_heat * w_n * (R_ros / 60.0)
    L_f   = 0.45 * I_B ** 0.46 if I_B > 0.0 else 0.0

    return dict(R0=R0, R_ros=R_ros, I_R=I_R, I_B=I_B, L_f=L_f,
                phi_w=phi_w, phi_s=phi_s)


# ---------------------------------------------------------------------------
# Fuel model database (abbreviated – same subset as behavior_matrix.py)
# ---------------------------------------------------------------------------

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

_FBFM40: Dict[int, Tuple] = {
    101:(101, "GR1 Short sparse dry grass",    0.011, 2200, 0.4, 0.15),
    102:(102, "GR2 Low grow dry grass",        0.046, 2000, 1.0, 0.15),
    103:(103, "GR3 Low grass",                 0.057, 1500, 2.0, 0.30),
    104:(104, "GR4 Mod. load dry grass",       0.087, 2000, 2.0, 0.15),
    105:(105, "GR5 Low load moist grass",      0.092, 1800, 1.5, 0.40),
    106:(106, "GR6 Mod. load humid grass",     0.115, 2200, 1.5, 0.40),
    107:(107, "GR7 High load dry grass",       0.230, 2000, 3.0, 0.15),
    108:(108, "GR8 High load moist grass",     0.299, 1500, 4.0, 0.30),
    109:(109, "GR9 Very high load dry grass",  0.414, 1800, 5.0, 0.40),
    141:(141, "SH1 Low load dry shrub",        0.046, 1600, 1.0, 0.15),
    142:(142, "SH2 Mod. load dry shrub",       0.161, 1600, 1.0, 0.15),
    145:(145, "SH5 High load dry shrub",       0.299,  750, 6.0, 0.15),
    161:(161, "TU1 Low load dry grass-shrub",  0.092, 2000, 0.6, 0.20),
    165:(165, "TU5 High load conifer litter",  0.322, 1500, 1.0, 0.25),
    181:(181, "TL1 Low load conifer litter",   0.069, 2000, 0.2, 0.30),
    185:(185, "TL5 High load conifer litter",  0.161, 1500, 0.6, 0.25),
}

_M_S_TO_FT_MIN    = 196.85
_BTU_FT_S_TO_KW_M = 0.34879
_FT_MIN_TO_M_MIN  = 0.3048


# ---------------------------------------------------------------------------
# Van Wagner (1977) crown fire empirical constants
# ---------------------------------------------------------------------------
# Van Wagner, C.E. (1977). Conditions for the start and spread of crown fire.
#   Canadian Journal of Forest Research, 7(1), 23–34.
#
# Critical initiation intensity (Van Wagner 1977, Eq. 1):
#   I_o [kW/m] = _VW_I0_COEFF × CBH × (_VW_I0_CONST + _VW_I0_FMC × FMC%)
_VW_I0_COEFF  = 0.010   # Van Wagner (1977) Eq. 1 leading coefficient
_VW_I0_CONST  = 460.0   # Van Wagner (1977) Eq. 1 intercept [°C·s/m]
_VW_I0_FMC    = 25.9    # Van Wagner (1977) Eq. 1 FMC coefficient

# Critical active crown ROS (Van Wagner 1977, Eq. 2):
#   R'_SA [m/min] = _VW_CROWN_ROS_NUM / CBD [kg/m³]
_VW_CROWN_ROS_NUM  = 3.0   # Van Wagner (1977) Eq. 2 numerator [m²/min·kg]

# Moisture reduction factor bounds for active crown ROS
# (Scott & Reinhardt 2001, following Van Wagner 1977 discussion):
#   m_factor = 1 - (FMC - _VW_MF_FMC_REF) / _VW_MF_FMC_RANGE
#   clamped to [_VW_MF_MIN, 1.0]
_VW_MF_FMC_REF   = 100.0   # reference foliar moisture [%] (Scott & Reinhardt 2001)
_VW_MF_FMC_RANGE = 200.0   # FMC denominator [%]
_VW_MF_MIN       = 0.3     # minimum moisture reduction factor (Scott & Reinhardt 2001)


def _get_fuel(code: int, system: str) -> Optional[Tuple]:
    db = _FBFM13 if system == "13" else _FBFM40
    return db.get(code)


# ---------------------------------------------------------------------------
# Van Wagner (1977) crown fire functions
# ---------------------------------------------------------------------------

def critical_crown_intensity(CBH: float, FMC: float) -> float:
    """Compute Van Wagner (1977) critical surface fire intensity for crown initiation.

    Uses Eq. 1 from Van Wagner (1977):

        I_o [kW/m] = _VW_I0_COEFF × CBH × (_VW_I0_CONST + _VW_I0_FMC × FMC%)

    Parameters
    ----------
    CBH : float
        Canopy base height [m].
    FMC : float
        Foliar moisture content [%].

    Returns
    -------
    float
        Critical fireline intensity I_o [kW/m].
    """
    CBH = max(CBH, 0.1)
    FMC = max(FMC, 50.0)
    return _VW_I0_COEFF * CBH * (_VW_I0_CONST + _VW_I0_FMC * FMC)


def active_crown_ros(CBD: float, FMC: float = 100.0) -> float:
    """Compute Van Wagner (1977) active crown fire rate of spread.

    Uses Eq. 2 from Van Wagner (1977):

        R'_SA [m/min] = _VW_CROWN_ROS_NUM / CBD

    A moisture reduction factor (Scott & Reinhardt 2001) is applied:

        m_factor = 1 - (FMC - _VW_MF_FMC_REF) / _VW_MF_FMC_RANGE
        m_factor = clamp(m_factor, _VW_MF_MIN, 1.0)

    Parameters
    ----------
    CBD : float
        Canopy bulk density [kg/m³].
    FMC : float
        Foliar moisture content [%] (used for moisture reduction factor).

    Returns
    -------
    float
        Active crown fire ROS R'_SA [m/min].
    """
    CBD = max(CBD, 0.01)
    R_sa = _VW_CROWN_ROS_NUM / CBD
    m_factor = 1.0 - (FMC - _VW_MF_FMC_REF) / _VW_MF_FMC_RANGE
    m_factor = max(_VW_MF_MIN, min(1.0, m_factor))
    return R_sa * m_factor


def crown_status_label(I_B_kW: float, R_ros_m_min: float,
                       I_o: float, R_sa: float) -> str:
    """Return crown fire activity label.

    Parameters
    ----------
    I_B_kW : float
        Byram fireline intensity [kW/m].
    R_ros_m_min : float
        Surface fire ROS [m/min].
    I_o : float
        Critical crown initiation intensity [kW/m].
    R_sa : float
        Critical active crown ROS [m/min].

    Returns
    -------
    str
        One of ``"surface"``, ``"passive"``, or ``"active"``.
    """
    if I_B_kW < I_o:
        return "surface"
    return "active" if R_ros_m_min >= R_sa else "passive"


def crowning_index(
    w0: float, sigma: float, delta: float,
    M_f: float, M_x: float,
    I_o: float,
    slope_tan: float = 0.0,
    U_max_ms: float = 50.0,
    n_iter: int = 60,
) -> float:
    """Estimate the crowning-index wind speed [m/s].

    The crowning index (CI) is the mid-flame wind speed at which the
    surface fire just reaches the critical crown initiation intensity
    I_o (Van Wagner 1977).  It is found by binary search.

    Parameters
    ----------
    w0, sigma, delta, M_f, M_x : float
        Rothermel fuel parameters (BG units).
    I_o : float
        Critical crown initiation intensity [kW/m].
    slope_tan : float
        Terrain slope as tan(angle).
    U_max_ms : float
        Upper bound for the binary search [m/s] (default: 50 m/s).
    n_iter : int
        Number of binary-search iterations (default: 60).

    Returns
    -------
    float
        Crowning-index wind speed [m/s], or ``U_max_ms`` if not reachable.
    """
    # Convert I_o to BTU/ft/s for comparison inside Rothermel
    I_o_btu = I_o / _BTU_FT_S_TO_KW_M

    lo, hi = 0.0, U_max_ms
    # Check if CI is achievable at all
    res_hi = _rothermel(w0, sigma, delta, M_f, M_x,
                        U_ftmin=hi * _M_S_TO_FT_MIN, slope_tan=slope_tan)
    if res_hi["I_B"] < I_o_btu:
        return U_max_ms  # CI not reachable within search range

    for _ in range(n_iter):
        mid = (lo + hi) / 2.0
        res = _rothermel(w0, sigma, delta, M_f, M_x,
                         U_ftmin=mid * _M_S_TO_FT_MIN, slope_tan=slope_tan)
        if res["I_B"] >= I_o_btu:
            hi = mid
        else:
            lo = mid
    return round((lo + hi) / 2.0, 3)


# ---------------------------------------------------------------------------
# Worksheet computation
# ---------------------------------------------------------------------------

def compute_crown_worksheet(
    fuel_code: int,
    fuel_system: str,
    wind_speeds_ms: List[float],
    moistures: List[float],
    CBH: float = 4.0,
    CBD: float = 0.15,
    FMC: float = 100.0,
    slope_tan: float = 0.0,
) -> List[Dict]:
    """Compute the crown fire assessment worksheet.

    For each (wind, moisture) combination, returns a dict containing:
      * Surface ROS, fireline intensity, flame length
      * Crown initiation threshold and status
      * Active crown ROS and crowning index

    Parameters
    ----------
    fuel_code : int
        Fuel model code.
    fuel_system : str
        ``"13"`` (FBFM13) or ``"40"`` (FBFM40).
    wind_speeds_ms : list[float]
        Wind speed values [m/s].
    moistures : list[float]
        Dead fuel moisture values [fraction, e.g. 0.08 = 8%].
    CBH : float
        Canopy base height [m] (default: 4.0).
    CBD : float
        Canopy bulk density [kg/m³] (default: 0.15).
    FMC : float
        Foliar moisture content [%] (default: 100).
    slope_tan : float
        Terrain slope as tan(angle) (default: 0.0).

    Returns
    -------
    list[dict]
        One dict per (wind, moisture) pair.
    """
    fuel = _get_fuel(fuel_code, fuel_system)
    if fuel is None:
        raise ValueError(
            f"Fuel model {fuel_code} not found in FBFM{fuel_system} database."
        )
    _code, name, w0, sigma, delta, M_x = fuel

    I_o   = critical_crown_intensity(CBH, FMC)
    R_sa  = active_crown_ros(CBD, FMC)        # m/min

    rows = []
    for U_ms in wind_speeds_ms:
        # Crowning index is independent of moisture, but depends on wind range,
        # so compute it once per wind-moisture pair with the given moisture.
        for M_f in moistures:
            U_ftmin = U_ms * _M_S_TO_FT_MIN
            res  = _rothermel(w0, sigma, delta, M_f, M_x, U_ftmin=U_ftmin,
                               slope_tan=slope_tan)
            R_surf    = res["R_ros"] * _FT_MIN_TO_M_MIN          # m/min
            I_B_kW    = res["I_B"]   * _BTU_FT_S_TO_KW_M        # kW/m
            L_f_m     = res["L_f"]   * 0.3048                    # m
            status    = crown_status_label(I_B_kW, R_surf, I_o, R_sa)
            R_crown   = active_crown_ros(CBD, FMC)               # m/min
            R_total   = max(R_surf, R_crown) if status == "active" else R_surf
            cf        = (R_crown / (R_crown + R_surf)
                         if status == "active" and (R_crown + R_surf) > 0 else 0.0)
            CI        = crowning_index(w0, sigma, delta, M_f, M_x, I_o, slope_tan)

            rows.append({
                "fuel_code":      fuel_code,
                "fuel_name":      name,
                "wind_m_s":       round(U_ms, 3),
                "moisture_pct":   round(M_f * 100.0, 1),
                "slope_deg":      round(math.degrees(math.atan(slope_tan)), 1),
                "CBH_m":          round(CBH, 2),
                "CBD_kgm3":       round(CBD, 4),
                "FMC_pct":        round(FMC, 1),
                "I_o_kW_m":       round(I_o, 2),
                "R_sa_m_min":     round(R_sa, 3),
                "R_surf_m_min":   round(R_surf, 4),
                "I_B_kW_m":       round(I_B_kW, 3),
                "L_f_m":          round(L_f_m, 3),
                "crown_status":   status,
                "R_crown_m_min":  round(R_crown, 3),
                "crown_fraction": round(cf, 4),
                "R_total_m_min":  round(R_total, 4),
                "CI_m_s":         round(CI, 3),
            })
    return rows


# ---------------------------------------------------------------------------
# ASCII worksheet printer
# ---------------------------------------------------------------------------

def print_crown_worksheet(rows: List[Dict]) -> None:
    """Print the crown fire assessment worksheet as an ASCII table.

    Displays the crown fire status for each (wind, moisture) pair in a
    condensed, BehavePlus-style layout.

    Parameters
    ----------
    rows : list[dict]
        Output of :func:`compute_crown_worksheet`.
    """
    if not rows:
        return

    r0 = rows[0]
    print(f"\nCrown fire assessment worksheet")
    print(f"  Fuel: FM{r0['fuel_code']} {r0['fuel_name']}")
    print(f"  CBH={r0['CBH_m']} m   CBD={r0['CBD_kgm3']} kg/m³   FMC={r0['FMC_pct']}%")
    print(f"  I_o (critical initiation intensity) = {r0['I_o_kW_m']:.1f} kW/m")
    print(f"  R'_SA (critical active crown ROS)   = {r0['R_sa_m_min']:.2f} m/min")
    print()

    winds = sorted({r["wind_m_s"] for r in rows})
    moist = sorted({r["moisture_pct"] for r in rows})

    # Header
    col_w = 10
    hdr = f"{'Wind [m/s]':>12s}" + "".join(f"  MC={m:4.1f}%" for m in moist)

    for label, key in [
        ("ROS surf [m/min]",   "R_surf_m_min"),
        ("I_B [kW/m]",         "I_B_kW_m"),
        ("Crown status",        "crown_status"),
        ("ROS total [m/min]",  "R_total_m_min"),
        ("CI [m/s]",           "CI_m_s"),
    ]:
        lookup = {(r["wind_m_s"], r["moisture_pct"]): r[key] for r in rows}
        print(f"{label}")
        print("-" * len(hdr))
        print(hdr)
        print("-" * len(hdr))
        for U in winds:
            row_str = f"{U:12.2f}"
            for m in moist:
                val = lookup.get((U, m), "")
                if isinstance(val, float):
                    row_str += f"  {val:>9.3f}"
                else:
                    row_str += f"  {str(val):>9s}"
            print(row_str)
        print()


# ---------------------------------------------------------------------------
# CSV output
# ---------------------------------------------------------------------------

def write_crown_csv(rows: List[Dict], out_path: str) -> None:
    """Write the crown fire worksheet rows to a CSV file.

    Parameters
    ----------
    rows : list[dict]
        Output of :func:`compute_crown_worksheet`.
    out_path : str
        Destination file path.
    """
    if not rows:
        return
    fieldnames = list(rows[0].keys())
    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} rows → {out_path}")


# ---------------------------------------------------------------------------
# Optional heatmap
# ---------------------------------------------------------------------------

def make_crown_heatmap(rows: List[Dict], out_path: str) -> None:
    """Save a heatmap of crown fire status and total ROS as a PNG figure.

    Requires ``matplotlib`` and ``numpy``.

    Parameters
    ----------
    rows : list[dict]
        Output of :func:`compute_crown_worksheet`.
    out_path : str
        Destination PNG file path.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        print("WARNING: matplotlib/numpy not available – skipping heatmap.",
              file=sys.stderr)
        return

    winds = sorted({r["wind_m_s"] for r in rows})
    moist = sorted({r["moisture_pct"] for r in rows})
    status_map = {"surface": 0, "passive": 1, "active": 2}

    def grid(key, transform=None):
        lu = {(r["wind_m_s"], r["moisture_pct"]): r[key] for r in rows}
        arr = np.array([[lu.get((U, m), 0) for m in moist] for U in winds],
                       dtype=float)
        return arr if transform is None else transform(arr)

    status_arr = grid("crown_status",
                      lambda _: np.array(
                          [[status_map.get(
                              next((r["crown_status"] for r in rows
                                    if r["wind_m_s"]==U and r["moisture_pct"]==m), "surface"),
                              0)
                            for m in moist] for U in winds], dtype=float))
    ros_arr = grid("R_total_m_min")

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    im0 = axes[0].imshow(status_arr, aspect="auto", origin="lower",
                          cmap="RdYlGn_r", vmin=0, vmax=2,
                          extent=[min(moist)-0.5, max(moist)+0.5,
                                  min(winds)-0.25, max(winds)+0.25])
    axes[0].set_xlabel("Dead Fuel Moisture [%]")
    axes[0].set_ylabel("Wind Speed [m/s]")
    axes[0].set_title("Crown Fire Status  (0=surface, 1=passive, 2=active)")
    cb0 = fig.colorbar(im0, ax=axes[0], ticks=[0, 1, 2])
    cb0.ax.set_yticklabels(["surface", "passive", "active"])

    im1 = axes[1].imshow(ros_arr, aspect="auto", origin="lower",
                          cmap="Reds",
                          extent=[min(moist)-0.5, max(moist)+0.5,
                                  min(winds)-0.25, max(winds)+0.25])
    axes[1].set_xlabel("Dead Fuel Moisture [%]")
    axes[1].set_ylabel("Wind Speed [m/s]")
    axes[1].set_title("Total ROS [m/min]")
    fig.colorbar(im1, ax=axes[1], label="m/min")

    r0 = rows[0]
    fig.suptitle(
        f"Crown fire worksheet – FM{r0['fuel_code']} {r0['fuel_name']}  "
        f"CBH={r0['CBH_m']} m  CBD={r0['CBD_kgm3']} kg/m³  FMC={r0['FMC_pct']}%"
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"Saved heatmap → {out_path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv=None):
    parser = argparse.ArgumentParser(
        description=(
            "Van Wagner (1977) crown fire assessment worksheet.\n\n"
            "Chains Rothermel (1972) surface ROS with Van Wagner's crown fire\n"
            "initiation and spread thresholds to classify each (wind, moisture)\n"
            "combination as surface, passive, or active crown fire and compute\n"
            "the crowning index (critical wind speed)."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    # --- fuel ---
    parser.add_argument("--fuel-model",   type=int,   default=4,
                        help="Fuel model code (default: 4)")
    parser.add_argument("--fuel-system",  default="13", choices=["13", "40"],
                        help="Fuel model system: 13 or 40 (default: 13)")
    # --- canopy stand ---
    parser.add_argument("--cbh",          type=float, default=4.0,
                        help="Canopy base height [m] (default: 4.0)")
    parser.add_argument("--cbd",          type=float, default=0.15,
                        help="Canopy bulk density [kg/m³] (default: 0.15)")
    parser.add_argument("--fmc",          type=float, default=100.0,
                        help="Foliar moisture content [%%] (default: 100)")
    # --- weather ---
    parser.add_argument("--wind-min",     type=float, default=0.0,
                        help="Minimum wind speed [m/s] (default: 0.0)")
    parser.add_argument("--wind-max",     type=float, default=12.0,
                        help="Maximum wind speed [m/s] (default: 12.0)")
    parser.add_argument("--wind-steps",   type=int,   default=13,
                        help="Number of wind speed steps (default: 13)")
    parser.add_argument("--moisture",     type=float, default=None,
                        help="Single dead fuel moisture fraction (default: 0.08)."
                             " Overrides --moisture-min/max/steps.")
    parser.add_argument("--moisture-min", type=float, default=0.04,
                        help="Min dead fuel moisture fraction (default: 0.04)")
    parser.add_argument("--moisture-max", type=float, default=0.20,
                        help="Max dead fuel moisture fraction (default: 0.20)")
    parser.add_argument("--moisture-steps", type=int, default=5,
                        help="Number of moisture steps (default: 5)")
    parser.add_argument("--slope",        type=float, default=0.0,
                        help="Terrain slope [tan(angle)] (default: 0.0)")
    # --- output ---
    parser.add_argument("--out",          default="crown_fire_worksheet.csv",
                        help="Output CSV path (default: crown_fire_worksheet.csv)")
    parser.add_argument("--plot",         action="store_true",
                        help="Save heatmap PNG (requires matplotlib + numpy)")
    parser.add_argument("--plot-out",     default="crown_fire_worksheet_heatmap.png",
                        help="Heatmap PNG path (default: crown_fire_worksheet_heatmap.png)")

    args = parser.parse_args(argv)

    def linspace(lo, hi, n):
        if n == 1:
            return [lo]
        return [lo + (hi - lo) * i / (n - 1) for i in range(n)]

    winds = linspace(args.wind_min, args.wind_max, args.wind_steps)
    if args.moisture is not None:
        moist = [args.moisture]
    else:
        moist = linspace(args.moisture_min, args.moisture_max, args.moisture_steps)

    fuel = _get_fuel(args.fuel_model, args.fuel_system)
    if fuel is None:
        print(f"ERROR: fuel model {args.fuel_model} not found in FBFM{args.fuel_system}.",
              file=sys.stderr)
        sys.exit(1)

    I_o = critical_crown_intensity(args.cbh, args.fmc)
    R_sa = active_crown_ros(args.cbd, args.fmc)

    print(f"Fuel model:  FM{args.fuel_model} – {fuel[1]}  (FBFM{args.fuel_system})")
    print(f"Canopy:      CBH={args.cbh} m   CBD={args.cbd} kg/m³   FMC={args.fmc}%")
    print(f"Critical I_o = {I_o:.1f} kW/m")
    print(f"Critical R'_SA = {R_sa:.2f} m/min")
    print(f"Wind speeds: {len(winds)} steps  [{args.wind_min:.1f} – {args.wind_max:.1f} m/s]")
    print(f"Moisture:    {len(moist)} steps  [{min(moist)*100:.1f}% – {max(moist)*100:.1f}%]")

    rows = compute_crown_worksheet(
        args.fuel_model, args.fuel_system, winds, moist,
        CBH=args.cbh, CBD=args.cbd, FMC=args.fmc, slope_tan=args.slope,
    )

    print_crown_worksheet(rows)
    write_crown_csv(rows, args.out)

    if args.plot:
        make_crown_heatmap(rows, args.plot_out)


if __name__ == "__main__":
    main()
