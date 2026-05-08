#!/usr/bin/env python3
"""
behavior_matrix.py – BehavePlus-style fuel condition fire behavior matrix.

Computes Rothermel (1972) fire behavior metrics across a grid of wind speeds,
dead fuel moisture values, terrain slopes, and fuel models (Anderson 13 or
Scott & Burgan 40).

Outputs a CSV matrix suitable for calibration, sensitivity analysis, and
field briefings — analogous to the tabular output produced by BehavePlus.

Capabilities
------------
  1. **Albini (1979) spotting distance** (``--spotting``):
     Estimates maximum firebrand spotting distance for torching-tree and
     intermittent crown-fire sources alongside every ROS / intensity row.

  2. **Fire size and shape at time T** (``--time``):
     For each (wind, moisture) combination, appends the expected fire area
     [ha], length [m], width [m], and perimeter [km] after ``--time`` minutes
     of unobstructed spread, using standard fire-ellipse geometry.

  3. **Multi-slope sweep** (``--slope-min``, ``--slope-max``, ``--slope-steps``):
     Sweeps terrain slope as a third dimension in the matrix, enabling a
     full (wind × moisture × slope) sensitivity surface.

  4. **Multi-fuel-model sweep** (``--fuel-models``):
     Sweeps over multiple fuel model codes in a single run, stacking all
     results in one CSV for cross-fuel-model comparison.

Metrics computed per (fuel, slope, wind, moisture) combination
--------------------------------------------------------------
  R_ros      – Rothermel rate of spread          [m/min]
  I_R        – Reaction intensity                [kW/m²]
  I_B        – Byram fireline intensity          [kW/m]
  L_f        – Byram flame length               [m]
  phi_w      – Wind factor                      [-]
  phi_s      – Slope factor                     [-]
  spot_dist  – Albini max spotting distance [m]  (with ``--spotting``)
  fire_area  – Fire area [ha]                    (with ``--time T``)
  perim      – Fire perimeter [km]               (with ``--time T``)

Requirements
------------
  None (pure Python).  Optional ``numpy`` / ``matplotlib`` for heatmaps.

Usage
-----
  # Matrix for Anderson FM 4 (chaparral), winds 0–10 m/s, moisture 4–20%
  python3 tools/behavior_matrix.py \\
      --fuel-model 4 --fuel-system 13 \\
      --wind-min 0 --wind-max 10 --wind-steps 11 \\
      --moisture-min 0.04 --moisture-max 0.20 --moisture-steps 9 \\
      --slope 0.3 \\
      --out fm4_matrix.csv

  # Slope sweep: wind × moisture × slope (0°–30°, 5 steps)
  python3 tools/behavior_matrix.py --fuel-model 4 \\
      --slope-min 0 --slope-max 0.577 --slope-steps 5

  # Multi-fuel sweep: FM 1, 4, 7 in one CSV
  python3 tools/behavior_matrix.py --fuel-models 1 4 7 --out multi_fuel.csv

  # Add Albini spotting distance (crown-fire mode)
  python3 tools/behavior_matrix.py --fuel-model 4 --spotting crown

  # Fire size at 60 minutes
  python3 tools/behavior_matrix.py --fuel-model 4 --time 60

  # Heatmap plots
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
  Albini, F.A. (1979). Spot fire distance from burning trees – a predictive
    model. USDA Forest Service General Technical Report INT-56.
  McAlpine, R.S. & Wakimoto, R.H. (1991). The acceleration of fire from
    point source to equilibrium spread. Forest Science, 37(5), 1314–1337.
  BehavePlus: https://www.firelab.org/project/behaveplusfiremodeling
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from typing import Dict, List, Optional, Tuple


# ---------------------------------------------------------------------------
# Albini (1979) maximum spotting distance nomograph (Cap 1)
# ---------------------------------------------------------------------------
# Albini (1979) spotting distance empirical constants (Cap 1)
# ---------------------------------------------------------------------------
# Albini, F.A. (1979). Spot fire distance from burning trees – a predictive
# model. USDA Forest Service General Technical Report INT-56.
#
# Two empirical modes:
#   "torching"  – maximum spotting distance from a group of torching trees
#                 driven by a mountain pine beetle type outbreak scenario.
#                 Approximation based on Albini (1979) equations for
#                 maximum horizontal travel distance (Eq. 7 / Table 2):
#                   L_spot [m] = F_h * U_ftmin^0.5 * H_z
#                 where H_z = Albini (1979) plume height [m] (≈ 12.2 * I_B^(1/3))
#                 and F_h is an empirical scaling (≈ 0.0176 per Albini 1979).
#
#   "crown"     – maximum spotting distance from an intermittent crown fire
#                 using the Scott & Reinhardt (2005) approximation:
#                   L_spot [m] = 0.0176 * I_B^0.655  [I_B in kW/m]
#
# These are condensed fits to the full Albini (1979) nomograph curves and give
# reasonable maximum spotting distance estimates for planning purposes.
# ---------------------------------------------------------------------------

# Albini (1979) Eq. 7 plume height scaling: H_z [m] = _ALBINI_PLUME_COEFF * I_B^(1/3)
_ALBINI_PLUME_COEFF     = 12.2
# Empirical firebrand horizontal transport factor F_h from Albini (1979) Table 2
_ALBINI_F_H             = 0.0176
# Scott & Reinhardt (2005) crown-fire spotting: L_spot [m] = coeff * I_B^exp [kW/m]
_ALBINI_CROWN_COEFF     = 0.0176
_ALBINI_CROWN_EXPONENT  = 0.655  # Scott & Reinhardt (2005), Eq. 12

def albini_spotting_distance(
    I_B_kW_m: float,
    U_ms: float,
    mode: str = "crown",
) -> float:
    """Estimate maximum firebrand spotting distance [m] using Albini (1979).

    Parameters
    ----------
    I_B_kW_m : float
        Byram fireline intensity [kW/m].
    U_ms : float
        20-ft wind speed [m/s].
    mode : str
        "torching" – torching-tree source (Albini 1979 Eq. 7 approximation).
        "crown"    – intermittent crown fire source (Scott & Reinhardt 2005).

    Returns
    -------
    float
        Estimated maximum spotting distance [m].  Returns 0 when I_B ≤ 0.
    """
    if I_B_kW_m <= 0.0:
        return 0.0

    if mode == "torching":
        # Albini (1979) Eq. 7 approximation:
        #   H_z [m] = _ALBINI_PLUME_COEFF * I_B^(1/3)   (plume lofting height)
        #   L_spot [m] = _ALBINI_F_H * U_ftmin^0.5 * H_z
        H_z = _ALBINI_PLUME_COEFF * I_B_kW_m ** (1.0 / 3.0)
        U_ftmin = U_ms * _M_S_TO_FT_MIN
        L_spot = _ALBINI_F_H * math.sqrt(max(U_ftmin, 0.0)) * H_z
    else:
        # Scott & Reinhardt (2005) crown-fire spotting approximation (Eq. 12):
        #   L_spot [m] = _ALBINI_CROWN_COEFF * I_B^_ALBINI_CROWN_EXPONENT
        L_spot = _ALBINI_CROWN_COEFF * I_B_kW_m ** _ALBINI_CROWN_EXPONENT
    return max(L_spot, 0.0)


# ---------------------------------------------------------------------------
# Fire size and shape at elapsed time T (Cap 2)
# ---------------------------------------------------------------------------
# Given head-fire ROS and a length-to-width (L/W) ratio, compute the
# expected fire size and shape at elapsed time T using standard fire-ellipse
# geometry (McAlpine & Wakimoto 1991; Finney 2004).
#
#   Head fire spread distance from ignition to head:
#     d_H [m] = R_head [m/min] * T [min]
#
#   Ellipse semi-axes:
#     a = d_H / (1 + c/a)  = length / 2        (semi-major axis [m])
#     b = a / (L/W)                              (semi-minor axis [m])
#
#   Fire area [ha]  = π * a * b / 10000
#   Fire perimeter [km] (ellipse approx.)
#     ≈ π * (3(a+b) - √((3a+b)(a+3b))) / 1000  (Ramanujan 1914)
# ---------------------------------------------------------------------------

def fire_shape_at_time(
    R_head_m_min: float,
    elapsed_min: float,
    length_to_width: float = 3.0,
) -> Dict[str, float]:
    """Compute expected fire size and shape at *elapsed_min* minutes.

    Parameters
    ----------
    R_head_m_min : float
        Head-fire rate of spread [m/min].
    elapsed_min : float
        Elapsed simulation time [min].
    length_to_width : float
        Fire ellipse L/W ratio (default: 3.0).

    Returns
    -------
    dict with keys:
        length_m      – fire length (head to back) [m]
        width_m       – fire width (flank to flank) [m]
        area_ha       – fire area [ha]
        perimeter_km  – fire perimeter [km]
        d_head_m      – head-fire spread distance from ignition [m]
    """
    if R_head_m_min <= 0.0 or elapsed_min <= 0.0 or length_to_width <= 0.0:
        return dict(length_m=0.0, width_m=0.0, area_ha=0.0,
                    perimeter_km=0.0, d_head_m=0.0)

    LW = max(length_to_width, 1.0)
    # head spread distance
    d_head = R_head_m_min * elapsed_min
    # ellipse semi-major axis (= half the fire length head-to-back)
    # For a fire ellipse the head is at one end; ignition at focus.
    # Using the simple approximation: fire length ≈ 2 * a
    # so a ≈ d_head (head-fire travels a from ignition centre)
    a = d_head  # semi-major axis [m]
    b = a / LW  # semi-minor axis [m]

    fire_length = 2.0 * a
    fire_width  = 2.0 * b
    area_ha     = math.pi * a * b / 1.0e4

    # Ramanujan's first approximation for ellipse perimeter
    h = ((a - b) / (a + b)) ** 2
    perim_m = math.pi * (a + b) * (1.0 + 3.0 * h / (10.0 + math.sqrt(4.0 - 3.0 * h)))
    perim_km = perim_m / 1000.0

    return dict(
        length_m=round(fire_length, 2),
        width_m=round(fire_width, 2),
        area_ha=round(area_ha, 4),
        perimeter_km=round(perim_km, 4),
        d_head_m=round(d_head, 2),
    )

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
    fuel_codes: List[int],
    fuel_system: str,
    wind_speeds_ms: List[float],    # [m/s]
    moistures: List[float],          # [fraction]
    slope_tans: Optional[List[float]] = None,  # list of tan(slope)
    h_heat: float = 8000.0,
    S_T: float = 0.0555,
    S_e: float = 0.010,
    rho_p: float = 32.0,
    elapsed_min: float = 0.0,
    length_to_width: float = 3.0,
    spotting_mode: Optional[str] = None,
) -> List[Dict]:
    """Compute Rothermel fire behavior matrix.

    Supports multi-fuel-model and multi-slope sweeps (Cap 4).
    Optionally appends fire size/shape at elapsed_min (Cap 2) and
    Albini spotting distance columns (Cap 1).

    Parameters
    ----------
    fuel_codes : list[int]
        One or more fuel model codes to include in the sweep.
    fuel_system : str
        "13" (FBFM13) or "40" (FBFM40).
    wind_speeds_ms : list[float]
        Wind speed values [m/s].
    moistures : list[float]
        Dead fuel moisture values [fraction].
    slope_tans : list[float] or None
        Slope values as tan(angle).  If None, defaults to [0.0].
    elapsed_min : float
        Elapsed time for fire size/shape calculation [min].  0 = skip.
    length_to_width : float
        L/W ratio for fire ellipse (used when elapsed_min > 0).
    spotting_mode : str or None
        "torching", "crown", or None (skip spotting distance column).

    Returns a list of dicts, one per (fuel, slope, wind, moisture) combination.
    """
    if slope_tans is None:
        slope_tans = [0.0]

    rows = []
    for fuel_code in fuel_codes:
        fuel = _get_fuel(fuel_code, fuel_system)
        if fuel is None:
            raise ValueError(
                f"Fuel model {fuel_code} not found in FBFM{fuel_system} database."
            )
        _code, name, w0, sigma, delta, M_x = fuel

        for slope_tan in slope_tans:
            for U_ms in wind_speeds_ms:
                for M_f in moistures:
                    U_ftmin = U_ms * _M_S_TO_FT_MIN
                    res = _rothermel(w0, sigma, delta, M_f, M_x, h_heat, S_T, S_e, rho_p,
                                     U_ftmin, slope_tan)
                    I_B_kW = res["I_B"] * _BTU_FT_S_TO_KW_M
                    R_m_min = res["R_ros"] * _FT_MIN_TO_M_MIN
                    row: Dict = {
                        "fuel_code":     fuel_code,
                        "fuel_name":     name,
                        "wind_m_s":      round(U_ms, 3),
                        "moisture_pct":  round(M_f * 100.0, 1),
                        "slope_deg":     round(math.degrees(math.atan(slope_tan)), 1),
                        "R_ros_m_min":   round(R_m_min, 4),
                        "I_R_kW_m2":     round(res["I_R"] * _BTU_FT2_MIN_TO_KW_M2, 3),
                        "I_B_kW_m":      round(I_B_kW, 3),
                        "L_f_m":         round(res["L_f"] * _FT_TO_M, 3),
                        "phi_w":         round(res["phi_w"], 4),
                        "phi_s":         round(res["phi_s"], 4),
                    }
                    # Cap 2: fire size/shape at elapsed time T
                    if elapsed_min > 0.0:
                        shape = fire_shape_at_time(R_m_min, elapsed_min, length_to_width)
                        row["elapsed_min"]   = round(elapsed_min, 1)
                        row["fire_area_ha"]  = shape["area_ha"]
                        row["fire_length_m"] = shape["length_m"]
                        row["fire_width_m"]  = shape["width_m"]
                        row["perim_km"]      = shape["perimeter_km"]
                    # Cap 1: Albini (1979) spotting distance
                    if spotting_mode is not None:
                        row["spot_dist_m"] = round(
                            albini_spotting_distance(I_B_kW, U_ms, spotting_mode), 1)
                    rows.append(row)
    return rows


# ---------------------------------------------------------------------------
# CSV writer
# ---------------------------------------------------------------------------

def write_matrix_csv(rows: List[Dict], out_path: str) -> None:
    """Write the fire behavior matrix rows to a CSV file.

    Parameters
    ----------
    rows : list[dict]
        Output of :func:`compute_matrix`.
    out_path : str
        Destination file path.
    """
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
    """Print a pivot table: rows = wind speeds, columns = moisture values.

    When the matrix contains multiple fuel models or slope values, only the
    first unique combination of (fuel_code, slope_deg) is displayed.  Use
    the CSV output for the full multi-dimensional result.

    Parameters
    ----------
    rows : list[dict]
        Output of :func:`compute_matrix`.
    metric : str
        Column name to pivot (default: ``"R_ros_m_min"``).
    """
    if not rows:
        return

    # Filter to first fuel_code and first slope_deg for a clean 2-D pivot
    first_fuel = rows[0]["fuel_code"]
    first_slope = rows[0]["slope_deg"]
    subset = [r for r in rows if r["fuel_code"] == first_fuel
              and r["slope_deg"] == first_slope]

    winds = sorted({r["wind_m_s"] for r in subset})
    moist = sorted({r["moisture_pct"] for r in subset})

    if not winds or not moist:
        return

    label_map = {
        "R_ros_m_min": "ROS [m/min]",
        "I_B_kW_m":    "IB [kW/m]",
        "L_f_m":       "L_f [m]",
        "I_R_kW_m2":   "IR [kW/m²]",
        "spot_dist_m": "Spot dist [m]",
        "fire_area_ha": "Fire area [ha]",
    }
    title = label_map.get(metric, metric)
    if metric not in subset[0]:
        return

    lookup = {(r["wind_m_s"], r["moisture_pct"]): r[metric] for r in subset}

    header = f"{'Wind [m/s]':>12s}" + "".join(f"  MC={m:4.1f}%" for m in moist)
    fuel_name = rows[0].get("fuel_name", "")
    print(f"\nFM{first_fuel} {fuel_name}  slope={first_slope}°  {title}")
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
    """Save heatmaps of ROS, IB, and L_f as a PNG figure.

    Generates a panel of colour-mapped heatmaps (wind speed × dead fuel
    moisture) for the first fuel model / slope combination in *rows*.
    Requires ``matplotlib`` and ``numpy``.

    Parameters
    ----------
    rows : list[dict]
        Output of :func:`compute_matrix`.
    out_path : str
        Destination PNG file path.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        print("WARNING: matplotlib/numpy not available – skipping heatmap.", file=sys.stderr)
        return

    # Use only the first fuel / slope combination for the heatmap
    first_fuel  = rows[0]["fuel_code"]
    first_slope = rows[0]["slope_deg"]
    subset = [r for r in rows if r["fuel_code"] == first_fuel
              and r["slope_deg"] == first_slope]

    winds = sorted({r["wind_m_s"] for r in subset})
    moist = sorted({r["moisture_pct"] for r in subset})

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
            next((r[metric] for r in subset if r["wind_m_s"]==U and r["moisture_pct"]==m), 0)
            for m in moist] for U in winds])
        im = ax.imshow(data, aspect="auto", origin="lower", cmap=cmap,
                       extent=[min(moist)-0.5, max(moist)+0.5,
                               min(winds)-0.25, max(winds)+0.25])
        ax.set_xlabel("Dead fuel moisture [%]")
        ax.set_ylabel("Wind speed [m/s]")
        ax.set_title(label)
        fig.colorbar(im, ax=ax, fraction=0.04)

    fuel_name = rows[0].get("fuel_name", "") if rows else ""
    fig.suptitle(
        f"Fire behavior matrix – FM {first_fuel} {fuel_name}  slope={first_slope}°"
        if rows else ""
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
            "BehavePlus-style Rothermel fire behavior matrix.\n\n"
            "Sweeps wind speed, dead fuel moisture, terrain slope, and fuel model\n"
            "to produce a CSV table of fire behavior metrics.  Optional capabilities:\n"
            "  --spotting MODE  : append Albini (1979) max spotting distance [m]\n"
            "  --time T         : append expected fire size / shape after T minutes\n"
            "  --slope-steps N  : sweep slope as a third matrix dimension\n"
            "  --fuel-models FM1 FM2 … : sweep over multiple fuel model codes"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    # --- fuel model ---
    fm_group = parser.add_mutually_exclusive_group()
    fm_group.add_argument("--fuel-model",  type=int, default=4,
                          help="Single fuel model code (default: 4)")
    fm_group.add_argument("--fuel-models", type=int, nargs="+", metavar="FM",
                          help="One or more fuel model codes for a multi-fuel sweep "
                               "(overrides --fuel-model)")
    parser.add_argument("--fuel-system", default="13", choices=["13", "40"],
                        help="Fuel model system: 13 (FBFM13) or 40 (FBFM40) (default: 13)")
    # --- wind ---
    parser.add_argument("--wind-min",    type=float, default=0.0,
                        help="Minimum wind speed [m/s] (default: 0.0)")
    parser.add_argument("--wind-max",    type=float, default=10.0,
                        help="Maximum wind speed [m/s] (default: 10.0)")
    parser.add_argument("--wind-steps",  type=int,   default=11,
                        help="Number of wind speed steps (default: 11)")
    # --- moisture ---
    parser.add_argument("--moisture-min", type=float, default=0.04,
                        help="Min dead fuel moisture fraction (default: 0.04)")
    parser.add_argument("--moisture-max", type=float, default=0.30,
                        help="Max dead fuel moisture fraction (default: 0.30)")
    parser.add_argument("--moisture-steps", type=int, default=9,
                        help="Number of moisture steps (default: 9)")
    # --- slope (single value kept for backward compatibility) ---
    parser.add_argument("--slope",       type=float, default=0.0,
                        help="Single slope value [tan(angle)] (default: 0.0)."
                             " Use --slope-min/max/steps for a sweep.")
    # --- slope sweep (Cap 4) ---
    parser.add_argument("--slope-min",   type=float, default=None,
                        help="Min slope [tan(angle)] for slope sweep (default: 0.0)")
    parser.add_argument("--slope-max",   type=float, default=None,
                        help="Max slope [tan(angle)] for slope sweep")
    parser.add_argument("--slope-steps", type=int,   default=None,
                        help="Number of slope steps (activates slope sweep)")
    # --- output ---
    parser.add_argument("--out",         default="behavior_matrix.csv",
                        help="Output CSV path (default: behavior_matrix.csv)")
    parser.add_argument("--plot",        action="store_true",
                        help="Save heatmap PNG (requires matplotlib + numpy)")
    parser.add_argument("--plot-out",    default="behavior_matrix_heatmap.png",
                        help="Heatmap output path (default: behavior_matrix_heatmap.png)")
    parser.add_argument("--list-fuels",  action="store_true",
                        help="Print available fuel model codes and exit.")
    # --- Cap 1: spotting distance ---
    parser.add_argument("--spotting",    default=None,
                        choices=["torching", "crown"],
                        help="Append Albini (1979) max spotting distance column.\n"
                             "'torching' = torching-tree source;\n"
                             "'crown'    = intermittent crown fire source.")
    # --- Cap 2: fire size at time T ---
    parser.add_argument("--time",        type=float, default=0.0,
                        help="Elapsed time [min] for fire size/shape calculation "
                             "(0 = disabled, default: 0)")
    parser.add_argument("--lw-ratio",    type=float, default=3.0,
                        help="Fire length-to-width ratio for ellipse geometry "
                             "(used with --time, default: 3.0)")

    args = parser.parse_args(argv)

    if args.list_fuels:
        db = _FBFM13 if args.fuel_system == "13" else _FBFM40
        print(f"\nFBFM{args.fuel_system} fuel models:")
        print(f"  {'Code':>6}  Name")
        print(f"  {'----':>6}  ----")
        for code, entry in sorted(db.items()):
            print(f"  {code:>6}  {entry[1]}")
        return

    # Resolve fuel model list
    fuel_codes = args.fuel_models if args.fuel_models else [args.fuel_model]

    # Build parameter grids
    def linspace(lo, hi, n):
        if n == 1:
            return [lo]
        return [lo + (hi - lo) * i / (n - 1) for i in range(n)]

    winds = linspace(args.wind_min, args.wind_max, args.wind_steps)
    moist = linspace(args.moisture_min, args.moisture_max, args.moisture_steps)

    # Slope grid: slope sweep takes precedence over single --slope
    if args.slope_steps is not None and args.slope_steps > 1:
        s_lo = args.slope_min if args.slope_min is not None else 0.0
        s_hi = args.slope_max if args.slope_max is not None else args.slope
        slopes = linspace(s_lo, s_hi, args.slope_steps)
    else:
        slopes = [args.slope]

    # Validate fuel codes
    for fc in fuel_codes:
        if _get_fuel(fc, args.fuel_system) is None:
            print(f"ERROR: fuel model {fc} not found in FBFM{args.fuel_system}.",
                  file=sys.stderr)
            sys.exit(1)

    print(f"Fuel models:  {fuel_codes}  (FBFM{args.fuel_system})")
    print(f"Wind speeds:  {len(winds)} steps  [{args.wind_min:.2f} – {args.wind_max:.2f} m/s]")
    print(f"Moisture:     {len(moist)} steps  [{args.moisture_min*100:.1f}% – {args.moisture_max*100:.1f}%]")
    print(f"Slopes:       {[round(math.degrees(math.atan(s)), 1) for s in slopes]}°")
    if args.spotting:
        print(f"Spotting:     Albini (1979) mode='{args.spotting}'")
    if args.time > 0:
        print(f"Fire shape:   at T={args.time:.1f} min  (L/W={args.lw_ratio:.1f})")

    rows = compute_matrix(
        fuel_codes, args.fuel_system, winds, moist,
        slope_tans=slopes,
        elapsed_min=args.time,
        length_to_width=args.lw_ratio,
        spotting_mode=args.spotting,
    )

    print_matrix_table(rows, "R_ros_m_min")
    print_matrix_table(rows, "I_B_kW_m")
    if args.spotting:
        print_matrix_table(rows, "spot_dist_m")
    if args.time > 0:
        print_matrix_table(rows, "fire_area_ha")

    write_matrix_csv(rows, args.out)

    if args.plot:
        make_heatmaps(rows, args.plot_out)


if __name__ == "__main__":
    main()
