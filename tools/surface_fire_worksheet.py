#!/usr/bin/env python3
"""
surface_fire_worksheet.py – BehavePlus-style comprehensive surface fire
behavior worksheet.

Computes the full Rothermel (1972) surface fire behavior chain for a single
point (or for a batch of rows from a CSV) and reports every intermediate and
final quantity in a human-readable worksheet analogous to the BehavePlus
*Surface Fire* module:

  • Fuel-bed physical properties (loading, SAV, bulk density, packing ratio)
  • Reaction intensity  I_R   [kW/m²]
  • Propagating flux ratio  ξ  [-]
  • Wind factor  φ_w           [-]
  • Slope factor  φ_s          [-]
  • No-wind / no-slope ROS  R₀ [m/min]
  • Head-fire ROS  R           [m/min]
  • Byram fireline intensity  I_B  [kW/m]
  • Byram flame length  L_f       [m]
  • Fire area, length, width, perimeter at elapsed time T  (optional)
  • Albini (1979) max spotting distance (optional)
  • Anderson (1970) probability of ignition  P_i  (optional)
  • Crowning potential via Van Wagner (1977) criterion (optional)

Output is printed as a formatted ASCII worksheet and optionally saved to CSV.

Requirements
------------
  None (pure Python).  Optional ``numpy`` / ``matplotlib`` for heatmap plots.

Usage
-----
  # Default FM4 fuel, 8% moisture, 10 m/s wind, flat terrain
  python3 tools/surface_fire_worksheet.py

  # Custom fuel and conditions
  python3 tools/surface_fire_worksheet.py \\
      --fuel-model 10 --fuel-system 13 \\
      --moisture 0.10 \\
      --wind 8.0 \\
      --slope 0.20 \\
      --out fm10_worksheet.csv

  # Scott & Burgan FM 145, crown fire assessment
  python3 tools/surface_fire_worksheet.py \\
      --fuel-model 145 --fuel-system 40 \\
      --moisture 0.08 \\
      --wind 12.0 \\
      --cbh 4.0 --cbd 0.12 --fmc 100 \\
      --crown

  # Spotting distance estimate (crown-fire firebrand source)
  python3 tools/surface_fire_worksheet.py \\
      --fuel-model 4 --wind 10.0 --moisture 0.06 --spotting crown

  # Fire size after 60 min
  python3 tools/surface_fire_worksheet.py \\
      --fuel-model 4 --wind 10.0 --moisture 0.06 --time 60

  # Batch mode: one row per CSV record
  python3 tools/surface_fire_worksheet.py \\
      --batch conditions.csv --out results.csv

Batch CSV format (header required)::

    fuel_model,fuel_system,moisture,wind_ms,slope_tan,
    4,13,0.08,10.0,0.0
    4,13,0.10,8.0,0.1

References
----------
  Rothermel, R.C. (1972). A mathematical model for predicting fire spread in
    wildland fuels. USDA Forest Service Research Paper INT-115.
  Andrews, P.L. (2018). The Rothermel surface fire spread model and associated
    developments. USDA Forest Service GTR RMRS-GTR-371.
  Byram, G.M. (1959). Combustion of forest fuels. In Davis (ed.) Forest Fire.
  Anderson, H.E. (1970). Forest fuel ignitibility. Fire Technology 6(4):312-319.
  Albini, F.A. (1979). Spot fire distance from burning trees. USDA FS INT-56.
  Van Wagner, C.E. (1977). Conditions for the start and spread of crown fire.
    Canadian Journal of Forest Research 7(1):23-34.
  BehavePlus: https://www.firelab.org/project/behaveplusfiremodeling
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from typing import Dict, List, Optional, Tuple


# ---------------------------------------------------------------------------
# Unit-conversion constants
# ---------------------------------------------------------------------------
_M_S_TO_FT_MIN    = 196.85      # m/s → ft/min
_BTU_FT_S_TO_KW_M = 0.34879     # BTU/(ft·s) → kW/m
_BTU_FT2_MIN_TO_KW_M2 = 0.18942 # BTU/(ft²·min) → kW/m²
_FT_MIN_TO_M_MIN  = 0.30480     # ft/min → m/min
_FT2_TO_M2        = 0.09290     # ft² → m²
_LB_FT2_TO_KG_M2  = 4.8824      # lb/ft² → kg/m²


# ---------------------------------------------------------------------------
# Fuel model database – Anderson 13 (FBFM13) and Scott & Burgan 40 (FBFM40)
# ---------------------------------------------------------------------------
# Tuple: (code, name, w0 [lb/ft²], sigma [ft⁻¹], delta [ft], M_x [fraction])
# This is the same condensed subset used in behavior_matrix.py / crown_fire_worksheet.py;
# only a representative sample of FBFM40 is listed for brevity.

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


def _get_fuel(code: int, system: str) -> Optional[Tuple]:
    """Return fuel model tuple (code, name, w0, sigma, delta, M_x) or None."""
    db = _FBFM40 if system == "40" else _FBFM13
    return db.get(code)


# ---------------------------------------------------------------------------
# Core Rothermel (1972) model — full-detail single-class implementation
# ---------------------------------------------------------------------------

def rothermel_full(
    w0:        float,  # oven-dry fuel load [lb/ft²]
    sigma:     float,  # surface-area-to-volume [ft⁻¹]
    delta:     float,  # fuel bed depth [ft]
    M_f:       float,  # fuel moisture [fraction]
    M_x:       float,  # moisture of extinction [fraction]
    h_heat:    float = 8000.0,  # heat content [BTU/lb]
    S_T:       float = 0.0555,  # total mineral content [fraction]
    S_e:       float = 0.010,   # effective mineral content [fraction]
    rho_p:     float = 32.0,    # particle density [lb/ft³]
    U_ms:      float = 0.0,     # midflame wind speed [m/s]
    slope_tan: float = 0.0,     # terrain slope as tan(angle)
) -> Dict[str, float]:
    """Full Rothermel (1972) single-class surface fire calculation.

    All inputs in British-Gravitational units (Rothermel 1972 convention);
    outputs converted to SI where labelled.

    Returns a dict with every intermediate and final quantity.
    """
    result: Dict[str, float] = {}

    if w0 <= 0.0 or sigma <= 0.0 or delta <= 0.0:
        # Non-burnable — return zeroed result
        for k in ("rho_b", "beta", "beta_op", "Gamma_max", "Gamma_p", "A",
                  "r_M", "eta_M", "eta_s", "w_n", "I_R", "xi", "eps_h",
                  "Q_ig", "R0", "phi_w", "phi_s", "R_ros",
                  "I_B", "L_f",
                  "rho_b_SI", "I_R_SI", "R0_SI", "R_ros_SI", "I_B_SI"):
            result[k] = 0.0
        return result

    # --- Packing ratio β ---
    rho_b    = w0 / delta                      # bulk density [lb/ft³]
    beta     = rho_b / rho_p                   # packing ratio [-]
    beta_op  = 3.348 * sigma ** (-0.8189)      # optimum packing ratio
    result["rho_b"]   = rho_b
    result["beta"]    = beta
    result["beta_op"] = beta_op

    # --- Optimum reaction velocity Γ' ---
    Gamma_max = (sigma ** 1.5) / (495.0 + 0.0594 * sigma ** 1.5)
    A         = 133.0 / sigma ** 0.7913
    Gamma_p   = Gamma_max * (beta / beta_op) ** A * math.exp(A * (1.0 - beta / beta_op))
    result["Gamma_max"] = Gamma_max
    result["Gamma_p"]   = Gamma_p
    result["A"]         = A

    # --- Moisture damping η_M ---
    r_M   = min(M_f / M_x, 1.0) if M_x > 0.0 else 0.0
    eta_M = max(0.0, 1.0 - 2.59*r_M + 5.11*r_M**2 - 3.52*r_M**3)
    result["r_M"]   = r_M
    result["eta_M"] = eta_M

    # --- Mineral damping η_s ---
    eta_s = 0.174 * S_e ** (-0.19)
    result["eta_s"] = eta_s

    # --- Net fuel loading ---
    w_n = w0 * (1.0 - S_T)
    result["w_n"] = w_n

    # --- Reaction intensity I_R [BTU/ft²/min] ---
    I_R = max(0.0, Gamma_p * w_n * h_heat * eta_M * eta_s)
    result["I_R"] = I_R

    # --- Propagating flux ratio ξ ---
    xi = math.exp((0.792 + 0.681 * sigma**0.5) * (beta + 0.1)) / (192.0 + 0.2595 * sigma)
    result["xi"] = xi

    # --- Heat of pre-ignition Q_ig [BTU/lb] ---
    eps_h = math.exp(-138.0 / sigma)
    Q_ig  = 250.0 + 1116.0 * M_f
    result["eps_h"] = eps_h
    result["Q_ig"]  = Q_ig

    # --- No-wind / no-slope ROS R₀ [ft/min] ---
    denom = rho_b * eps_h * Q_ig
    R0 = (I_R * xi / denom) if denom > 0.0 else 0.0
    result["R0"] = R0

    # --- Wind factor φ_w ---
    U_ftmin = U_ms * _M_S_TO_FT_MIN
    C = 7.47 * math.exp(-0.133 * sigma**0.55)
    B = 0.02526 * sigma**0.54
    E = 0.715  * math.exp(-3.59e-4 * sigma)
    phi_w = C * (U_ftmin ** B) * (beta / beta_op) ** (-E) if U_ftmin > 0.0 else 0.0
    result["phi_w"] = phi_w
    result["U_ftmin"] = U_ftmin
    result["C"] = C
    result["B"] = B
    result["E"] = E

    # --- Slope factor φ_s ---
    phi_s = 5.275 * beta**(-0.3) * slope_tan**2 if slope_tan > 0.0 else 0.0
    result["phi_s"] = phi_s

    # --- Head-fire ROS R [ft/min] ---
    R_ros = R0 * (1.0 + phi_w + phi_s)
    result["R_ros"] = R_ros

    # --- Byram fireline intensity I_B [BTU/(ft·s)] ---
    I_B  = h_heat * w_n * (R_ros / 60.0)
    result["I_B"] = I_B

    # --- Byram flame length L_f [ft] ---
    L_f = 0.45 * I_B**0.46 if I_B > 0.0 else 0.0
    result["L_f"] = L_f

    # --- SI conversions ---
    result["rho_b_SI"]  = rho_b * 16.0185                 # lb/ft³ → kg/m³
    result["I_R_SI"]    = I_R  * _BTU_FT2_MIN_TO_KW_M2   # BTU/(ft²·min) → kW/m²
    result["R0_SI"]     = R0   * _FT_MIN_TO_M_MIN         # ft/min → m/min
    result["R_ros_SI"]  = R_ros * _FT_MIN_TO_M_MIN        # ft/min → m/min
    result["I_B_SI"]    = I_B  * _BTU_FT_S_TO_KW_M       # BTU/(ft·s) → kW/m
    result["L_f_SI"]    = L_f  * 0.3048                   # ft → m

    return result


# ---------------------------------------------------------------------------
# Fire size at elapsed time T (standard fire-ellipse geometry)
# ---------------------------------------------------------------------------

def fire_size_at_time(R_head_m_min: float, elapsed_min: float,
                      LW: float = 3.0) -> Dict[str, float]:
    """Compute expected fire size and shape at *elapsed_min* minutes.

    Uses the McAlpine & Wakimoto (1991) / Finney (2004) fire-ellipse model:
      head spread distance d_H = R_head * T
      semi-major axis a = d_H / (1 + c/a)  ≈ d_H / (1 + 1/LW)  [simplified]
      semi-minor axis b = a / LW
    """
    if R_head_m_min <= 0.0 or elapsed_min <= 0.0 or LW <= 0.0:
        return dict(length_m=0.0, width_m=0.0, area_ha=0.0,
                    perimeter_km=0.0, d_head_m=0.0)
    d_H = R_head_m_min * elapsed_min
    a   = d_H / (1.0 + 1.0 / max(LW, 1.0))   # semi-major (head + back)
    b   = a / max(LW, 1.0)                    # semi-minor (flank)
    length_m = 2.0 * a
    width_m  = 2.0 * b
    area_ha  = math.pi * a * b / 1.0e4
    # Ramanujan perimeter approximation
    h     = ((a - b) / (a + b)) ** 2
    perim = math.pi * (a + b) * (1.0 + 3.0*h / (10.0 + math.sqrt(4.0 - 3.0*h)))
    return dict(length_m=length_m, width_m=width_m, area_ha=area_ha,
                perimeter_km=perim / 1000.0, d_head_m=d_H)


# ---------------------------------------------------------------------------
# Albini (1979) spotting distance estimate
# ---------------------------------------------------------------------------

_ALBINI_PLUME_COEFF = 12.2
_ALBINI_F_H         = 0.0176
_ALBINI_CROWN_COEFF = 0.0176
_ALBINI_CROWN_EXP   = 0.655   # Scott & Reinhardt (2005) Eq. 12

def albini_spotting(I_B_kW_m: float, U_ms: float, mode: str = "crown") -> float:
    """Albini (1979) max spotting distance [m].

    mode: "crown"    – Scott & Reinhardt (2005) crown-fire source
          "torching" – Albini (1979) plume-lofting approximation
    """
    if I_B_kW_m <= 0.0:
        return 0.0
    if mode == "torching":
        H_z    = _ALBINI_PLUME_COEFF * I_B_kW_m ** (1.0 / 3.0)
        U_ftmin = U_ms * _M_S_TO_FT_MIN
        return max(0.0, _ALBINI_F_H * math.sqrt(max(U_ftmin, 0.0)) * H_z)
    return max(0.0, _ALBINI_CROWN_COEFF * I_B_kW_m ** _ALBINI_CROWN_EXP)


# ---------------------------------------------------------------------------
# Anderson (1970) probability of ignition
# ---------------------------------------------------------------------------

def prob_ignition(T_ambient_C: float, MC_pct: float,
                  solar_heating_F: float = 25.0) -> float:
    """Anderson (1970) probability of ignition P_i.

    T_ambient_C     – ambient air temperature [°C]
    MC_pct          – 1-hr dead fuel moisture [%]
    solar_heating_F – solar temperature increment [°F] (default: 25 °F)
    """
    T_a_F    = T_ambient_C * 9.0 / 5.0 + 32.0
    T_fuel_F = T_a_F + solar_heating_F
    if T_fuel_F <= 0.0:
        return 0.0
    raw = 0.000048 * (T_fuel_F ** 1.4) * math.exp(-0.07 * MC_pct)
    return min(1.0, max(0.0, raw))


# ---------------------------------------------------------------------------
# Van Wagner (1977) crown fire initiation check
# ---------------------------------------------------------------------------

def van_wagner_crown(I_B_kW_m: float, CBH: float, FMC: float,
                     CBD: float) -> Dict[str, float]:
    """Van Wagner (1977) crown fire initiation and active-crown check.

    Parameters
    ----------
    I_B_kW_m : Byram fireline intensity [kW/m]
    CBH      : canopy base height [m]
    FMC      : foliar moisture content [%]
    CBD      : canopy bulk density [kg/m³]

    Returns a dict with:
      I_o        – critical initiation intensity [kW/m]
      initiates  – True when I_B >= I_o
      R_prime_SA – critical active-crown ROS [m/s]
    """
    I_o       = 0.010 * max(CBH, 0.01) * (460.0 + 25.9 * FMC)
    R_prime   = 3.0 / max(CBD, 0.01) / 60.0   # [m/s]
    initiates = I_B_kW_m >= I_o
    return dict(I_o=I_o, initiates=initiates, R_prime_SA=R_prime)


# ---------------------------------------------------------------------------
# Worksheet formatter
# ---------------------------------------------------------------------------

def print_worksheet(fuel: Tuple, r: Dict[str, float],
                    U_ms: float, slope_tan: float, M_f: float,
                    fire_size: Optional[Dict] = None,
                    spot_dist: Optional[float] = None,
                    p_ignition: Optional[float] = None,
                    crown_result: Optional[Dict] = None) -> None:
    """Print a BehavePlus-style ASCII worksheet."""
    code, name = fuel[0], fuel[1]
    w0, sigma, delta, M_x = fuel[2], fuel[3], fuel[4], fuel[5]

    hdr = "=" * 68
    print(hdr)
    print(f"  SURFACE FIRE BEHAVIOR WORKSHEET")
    print(f"  Rothermel (1972)  |  Andrews (2018 update)")
    print(hdr)
    print(f"\n  Fuel model  : FM{code} – {name}")
    print(f"  Fuel system : {'Scott & Burgan 40' if code >= 100 else 'Anderson 13'}")
    print(f"  Moisture    : {M_f*100:.1f} %  (M_x = {M_x*100:.0f} %)")
    print(f"  Wind speed  : {U_ms:.2f} m/s  ({U_ms * _M_S_TO_FT_MIN:.0f} ft/min)")
    print(f"  Slope       : {math.degrees(math.atan(slope_tan)):.1f}°  "
          f"(tan = {slope_tan:.3f})")

    print(f"\n  ── Fuel-Bed Properties ──────────────────────────────────────")
    print(f"    w0     (oven-dry load)      = {w0:.4f} lb/ft²"
          f"  = {w0*_LB_FT2_TO_KG_M2:.3f} kg/m²")
    print(f"    σ      (SAV ratio)          = {sigma:.0f} ft⁻¹")
    print(f"    δ      (fuel bed depth)     = {delta:.2f} ft"
          f"  = {delta*0.3048:.3f} m")
    print(f"    ρ_b    (bulk density)       = {r['rho_b']:.4f} lb/ft³"
          f"  = {r['rho_b_SI']:.3f} kg/m³")
    print(f"    β      (packing ratio)      = {r['beta']:.5f}")
    print(f"    β_op   (optimum packing)    = {r['beta_op']:.5f}")

    print(f"\n  ── Reaction Intensity ───────────────────────────────────────")
    print(f"    Γ'     (opt. rxn velocity)  = {r['Gamma_p']:.4f} lb/(ft²·min)")
    print(f"    η_M    (moisture damping)   = {r['eta_M']:.4f}")
    print(f"    η_s    (mineral damping)    = {r['eta_s']:.4f}")
    print(f"    I_R    (reaction intensity) = {r['I_R']:.2f} BTU/(ft²·min)"
          f"  = {r['I_R_SI']:.3f} kW/m²")

    print(f"\n  ── Wind & Slope Factors ─────────────────────────────────────")
    print(f"    ξ      (prop. flux ratio)   = {r['xi']:.5f}")
    print(f"    Q_ig   (heat of preignition)= {r['Q_ig']:.1f} BTU/lb")
    print(f"    R₀     (no-wind/slope ROS)  = {r['R0']:.4f} ft/min"
          f"  = {r['R0_SI']:.4f} m/min")
    print(f"    φ_w    (wind factor)        = {r['phi_w']:.4f}")
    print(f"    φ_s    (slope factor)       = {r['phi_s']:.4f}")

    print(f"\n  ── Fire Behavior ────────────────────────────────────────────")
    print(f"    R      (head-fire ROS)      = {r['R_ros']:.4f} ft/min"
          f"  = {r['R_ros_SI']:.4f} m/min"
          f"  = {r['R_ros_SI']*100.0/60.0:.3f} cm/s")
    print(f"    I_B    (fireline intensity) = {r['I_B']:.2f} BTU/(ft·s)"
          f"  = {r['I_B_SI']:.2f} kW/m")
    print(f"    L_f    (flame length)       = {r['L_f']:.2f} ft"
          f"  = {r['L_f_SI']:.2f} m")

    if fire_size is not None:
        print(f"\n  ── Fire Size at Time T ──────────────────────────────────────")
        print(f"    Head spread distance        = {fire_size['d_head_m']:.1f} m")
        print(f"    Fire length (head–back)     = {fire_size['length_m']:.1f} m")
        print(f"    Fire width  (flank–flank)   = {fire_size['width_m']:.1f} m")
        print(f"    Fire area                   = {fire_size['area_ha']:.4f} ha")
        print(f"    Fire perimeter              = {fire_size['perimeter_km']:.4f} km")

    if spot_dist is not None:
        print(f"\n  ── Spotting Distance ────────────────────────────────────────")
        print(f"    Albini (1979) max spotting  = {spot_dist:.0f} m"
              f"  = {spot_dist/1000.0:.3f} km")

    if p_ignition is not None:
        print(f"\n  ── Ignition Probability ─────────────────────────────────────")
        print(f"    P_i  (Anderson 1970)        = {p_ignition:.3f}  "
              f"({p_ignition*100.0:.1f} %)")

    if crown_result is not None:
        print(f"\n  ── Crown Fire Assessment (Van Wagner 1977) ──────────────────")
        print(f"    I_o  (initiation threshold) = {crown_result['I_o']:.1f} kW/m")
        print(f"    I_B  (computed)             = {r['I_B_SI']:.1f} kW/m")
        status = ("INITIATION EXPECTED" if crown_result["initiates"]
                  else "surface fire only")
        print(f"    Crown fire status           = {status}")
        print(f"    R'_SA (active-crown ROS)    = {crown_result['R_prime_SA']*60.0:.3f} m/min")

    print(f"\n{hdr}\n")


# ---------------------------------------------------------------------------
# Single-condition computation
# ---------------------------------------------------------------------------

def run_one(
    fuel_code:      int,
    fuel_system:    str,
    M_f:            float,
    U_ms:           float,
    slope_tan:      float,
    elapsed_min:    float = 0.0,
    LW:             float = 3.0,
    spotting_mode:  Optional[str] = None,
    T_ambient_C:    float = 25.0,
    solar_heating_F: float = 25.0,
    crown:          bool = False,
    CBH:            float = 4.0,
    CBD:            float = 0.15,
    FMC:            float = 100.0,
    h_heat:         float = 8000.0,
    S_T:            float = 0.0555,
    S_e:            float = 0.010,
    rho_p:          float = 32.0,
) -> Dict[str, float]:
    """Run the worksheet for one set of conditions; return result dict."""
    fuel = _get_fuel(fuel_code, fuel_system)
    if fuel is None:
        raise ValueError(f"Fuel model {fuel_code} not found in "
                         f"{'FBFM40' if fuel_system == '40' else 'FBFM13'} database")
    w0, sigma, delta, M_x = fuel[2], fuel[3], fuel[4], fuel[5]

    r = rothermel_full(w0, sigma, delta, M_f, M_x, h_heat, S_T, S_e, rho_p,
                       U_ms, slope_tan)

    fire_size   = fire_size_at_time(r["R_ros_SI"], elapsed_min, LW) if elapsed_min > 0 else None
    spot_dist   = albini_spotting(r["I_B_SI"], U_ms, spotting_mode) if spotting_mode else None
    p_ign       = prob_ignition(T_ambient_C, M_f * 100.0, solar_heating_F)
    crown_res   = van_wagner_crown(r["I_B_SI"], CBH, FMC, CBD) if crown else None

    return dict(
        fuel_code=fuel_code, fuel_name=fuel[1], fuel_system=fuel_system,
        M_f_pct=M_f * 100.0, M_x_pct=M_x * 100.0,
        U_ms=U_ms, slope_tan=slope_tan,
        rho_b_kgm3=r["rho_b_SI"], beta=r["beta"], beta_op=r["beta_op"],
        I_R_kW_m2=r["I_R_SI"], xi=r["xi"], eta_M=r["eta_M"], eta_s=r["eta_s"],
        phi_w=r["phi_w"], phi_s=r["phi_s"],
        R0_m_min=r["R0_SI"], R_ros_m_min=r["R_ros_SI"],
        I_B_kW_m=r["I_B_SI"], L_f_m=r["L_f_SI"],
        P_ignition=p_ign,
        fire_area_ha=fire_size["area_ha"] if fire_size else float("nan"),
        fire_perimeter_km=fire_size["perimeter_km"] if fire_size else float("nan"),
        spot_dist_m=spot_dist if spot_dist is not None else float("nan"),
        crown_initiates=int(crown_res["initiates"]) if crown_res else -1,
        I_o_kW_m=crown_res["I_o"] if crown_res else float("nan"),
        R_prime_SA_m_min=crown_res["R_prime_SA"] * 60.0 if crown_res else float("nan"),
        _r=r, _fuel=fuel, _fire_size=fire_size,
        _spot=spot_dist, _p_ign=p_ign, _crown=crown_res,
    )


# ---------------------------------------------------------------------------
# Batch mode
# ---------------------------------------------------------------------------

def run_batch(csv_path: str, args) -> List[Dict[str, float]]:
    """Run worksheet for every row in a CSV file; return list of result dicts."""
    results = []
    with open(csv_path, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            code    = int(row.get("fuel_model",  args.fuel_model))
            system  = row.get("fuel_system",     args.fuel_system)
            M_f     = float(row.get("moisture",  args.moisture))
            U_ms    = float(row.get("wind_ms",   args.wind))
            slope_t = float(row.get("slope_tan", args.slope))
            try:
                res = run_one(code, system, M_f, U_ms, slope_t,
                              elapsed_min=args.time, LW=args.lw,
                              spotting_mode=args.spotting,
                              T_ambient_C=args.temperature,
                              solar_heating_F=args.solar_heating,
                              crown=args.crown,
                              CBH=args.cbh, CBD=args.cbd, FMC=args.fmc)
            except ValueError as e:
                print(f"WARNING: {e} – skipping row", file=sys.stderr)
                continue
            results.append(res)
    return results


# ---------------------------------------------------------------------------
# CSV output
# ---------------------------------------------------------------------------

_CSV_FIELDS = [
    "fuel_code", "fuel_name", "fuel_system", "M_f_pct", "M_x_pct",
    "U_ms", "slope_tan", "rho_b_kgm3", "beta", "beta_op",
    "I_R_kW_m2", "xi", "eta_M", "eta_s", "phi_w", "phi_s",
    "R0_m_min", "R_ros_m_min", "I_B_kW_m", "L_f_m", "P_ignition",
    "fire_area_ha", "fire_perimeter_km", "spot_dist_m",
    "crown_initiates", "I_o_kW_m", "R_prime_SA_m_min",
]

def write_csv(results: List[Dict], out_path: str) -> None:
    """Write result dicts to CSV."""
    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=_CSV_FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(results)
    print(f"Wrote {len(results)} row(s) → {out_path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="BehavePlus-style surface fire behavior worksheet (Rothermel 1972).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("--fuel-model",  type=int,   default=4,       metavar="N",
                   help="Fuel model code (default: 4)")
    p.add_argument("--fuel-system", choices=["13", "40"], default="13",
                   help="Fuel model system: '13' (Anderson) or '40' (Scott&Burgan) (default: 13)")
    p.add_argument("--moisture",    type=float, default=0.08,    metavar="FRAC",
                   help="Dead fuel moisture fraction (default: 0.08 = 8%)")
    p.add_argument("--wind",        type=float, default=5.0,     metavar="M/S",
                   help="Midflame wind speed [m/s] (default: 5.0)")
    p.add_argument("--slope",       type=float, default=0.0,     metavar="TAN",
                   help="Terrain slope as tan(angle) (default: 0.0 = flat)")
    p.add_argument("--time",        type=float, default=0.0,     metavar="MIN",
                   help="Elapsed time [min] for fire-size computation (0 = skip)")
    p.add_argument("--lw",          type=float, default=3.0,
                   help="Fire length-to-width ratio for fire-size computation (default: 3.0)")
    p.add_argument("--spotting",    choices=["crown", "torching"], default=None,
                   help="Compute Albini (1979) max spotting distance (disabled by default)")
    p.add_argument("--temperature", type=float, default=25.0,   metavar="C",
                   help="Ambient temperature [°C] for P_i calculation (default: 25.0)")
    p.add_argument("--solar-heating", type=float, default=25.0, metavar="F",
                   help="Solar temperature increment [°F] for P_i (default: 25.0)")
    p.add_argument("--crown",       action="store_true",
                   help="Include Van Wagner (1977) crown fire assessment")
    p.add_argument("--cbh",         type=float, default=4.0,
                   help="Canopy base height [m] for crown assessment (default: 4.0)")
    p.add_argument("--cbd",         type=float, default=0.15,
                   help="Canopy bulk density [kg/m³] for crown assessment (default: 0.15)")
    p.add_argument("--fmc",         type=float, default=100.0,
                   help="Foliar moisture content [%] for crown assessment (default: 100)")
    p.add_argument("--batch",       default=None, metavar="CSV",
                   help="Batch mode: CSV file of conditions (header required)")
    p.add_argument("--out",         default=None, metavar="FILE",
                   help="Output CSV file (optional)")
    p.add_argument("--quiet",       action="store_true",
                   help="Suppress ASCII worksheet (useful in batch mode)")
    return p


def main(argv=None) -> None:
    args = _build_parser().parse_args(argv)

    if args.batch:
        results = run_batch(args.batch, args)
        if not args.quiet:
            for res in results:
                print_worksheet(res["_fuel"], res["_r"],
                                res["U_ms"], res["slope_tan"], res["M_f_pct"] / 100.0,
                                fire_size=res["_fire_size"],
                                spot_dist=res["_spot"],
                                p_ignition=res["_p_ign"],
                                crown_result=res["_crown"])
        if args.out:
            write_csv(results, args.out)
        return

    # Single-condition mode
    res = run_one(
        args.fuel_model, args.fuel_system, args.moisture, args.wind, args.slope,
        elapsed_min=args.time, LW=args.lw, spotting_mode=args.spotting,
        T_ambient_C=args.temperature, solar_heating_F=args.solar_heating,
        crown=args.crown, CBH=args.cbh, CBD=args.cbd, FMC=args.fmc,
    )

    if not args.quiet:
        print_worksheet(res["_fuel"], res["_r"],
                        args.wind, args.slope, args.moisture,
                        fire_size=res["_fire_size"],
                        spot_dist=res["_spot"],
                        p_ignition=res["_p_ign"],
                        crown_result=res["_crown"])

    if args.out:
        write_csv([res], args.out)


if __name__ == "__main__":
    main()
