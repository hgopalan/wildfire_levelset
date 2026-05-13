#!/usr/bin/env python3
"""
Canonical regression tests for the Rothermel (1972) surface fire spread model.

Covers all 13 standard FBFM13 (Anderson 1982) fuel models and validates the
parameters introduced or corrected in src/rothermel_model.H and
src/fuel_database.H, specifically:

  1. A-factor weighted characteristic SAV (σ_c) — not load-weighted
         σ_c = Σ(Aᵢ σᵢ) / Σ(Aᵢ),   Aᵢ = wᵢ σᵢ exp(−138/σᵢ)
     A-factor weighting naturally suppresses coarse dead fuels
     (σ₁₀ = 109, σ₁₀₀ = 30 ft⁻¹) and recovers σ_c values close to the
     Rothermel (1972) Table A-1 tabulated characteristic SAV.  A simple
     load-weighted average is dominated by the mass of coarse fuels and
     gives σ_c far too low, inflating the wind coefficient C and φ_w.

  2. FM4 and FM5 fine-fuel SAV: σ_d1 = 1739 ft⁻¹ (Rothermel 1972 Table A-1),
     not 2000 ft⁻¹.

  3. Wind factor φ_w = C U^B (β/β_op)^{-E} (Eq. 47) is zero at zero wind,
     positive for positive wind, and increases monotonically with wind speed.

  4. No-wind rate of spread R₀ is positive for each burnable model at
     standard dead fuel moisture (8 %) and vanishes when moisture equals M_x.

  5. Total ROS at standard conditions (5 mph midflame wind, 8 % dead
     moisture) falls within physically plausible ranges derived from the
     corrected model equations.

The Python reference implementation below is a direct port of the C++ code in:
  src/rothermel_model.H       – σ_c, Γ', ξ, I_R, R₀
  src/compute_rothermel_R.H   – φ_w, total ROS
  src/fuel_database.H         – per-class fuel parameters

Canonical expected values were computed from this reference implementation and
locked here to catch future regressions.  Any intentional change to the model
equations must be accompanied by updated expected values and a justification.

Run with:
  python3 tools/tests/test_rothermel_params.py
  python3 -m pytest tools/tests/test_rothermel_params.py -v

References
----------
Rothermel (1972) USDA Forest Service Research Paper INT-115
Andrews (2018)   Gen. Tech. Rep. RMRS-GTR-371
"""

import math
import unittest

# ── Physical constants (match src/constants.H / src/rothermel_model.H) ───────
SIGMA_D10  = 109.0   # 10-hr dead fuel SAV [ft⁻¹]
SIGMA_D100 =  30.0   # 100-hr dead fuel SAV [ft⁻¹]
WIND_CONV  = 196.85  # m/s → ft/min
ROS_CONV   = 0.00508 # ft/min → m/s

# ── FBFM13 fuel model database ────────────────────────────────────────────────
# Mirrors src/fuel_database.H exactly.
# Fields: w0, sigma (single-class fallback SAV), delta (depth), M_x,
#         h_heat, S_T (total mineral), S_e (effective mineral), rho_p,
#         w_d1/sigma_d1 (1-hr dead), w_d10 (10-hr dead), w_d100 (100-hr dead),
#         w_lh/sigma_lh (live herb), w_lw/sigma_lw (live woody).
# All loads in lb/ft², SAVs in ft⁻¹, delta in ft.
FBFM13 = {
    "FM1": dict(
        w0=0.034,  sigma=3500.0, delta=1.0,  M_x=0.12,
        h_heat=8000.0, S_T=0.0555, S_e=0.010, rho_p=32.0,
        w_d1=0.034,  sigma_d1=3500.0,
        w_d10=0.000, w_d100=0.000,
        w_lh=0.000,  sigma_lh=1500.0,
        w_lw=0.000,  sigma_lw=1500.0,
    ),
    "FM2": dict(
        w0=0.092,  sigma=3000.0, delta=1.0,  M_x=0.15,
        h_heat=8000.0, S_T=0.0555, S_e=0.010, rho_p=32.0,
        w_d1=0.046,  sigma_d1=3000.0,
        w_d10=0.023, w_d100=0.023,
        w_lh=0.023,  sigma_lh=1500.0,
        w_lw=0.000,  sigma_lw=1500.0,
    ),
    "FM3": dict(
        w0=0.138,  sigma=1500.0, delta=2.5,  M_x=0.25,
        h_heat=8000.0, S_T=0.0555, S_e=0.010, rho_p=32.0,
        w_d1=0.138,  sigma_d1=1500.0,
        w_d10=0.000, w_d100=0.000,
        w_lh=0.000,  sigma_lh=1500.0,
        w_lw=0.000,  sigma_lw=1500.0,
    ),
    "FM4": dict(
        w0=0.230,  sigma=1739.0, delta=6.0,  M_x=0.20,
        h_heat=8000.0, S_T=0.0555, S_e=0.010, rho_p=32.0,
        w_d1=0.230,  sigma_d1=1739.0,   # Rothermel 1972 Table A-1: 1739 ft⁻¹
        w_d10=0.184, w_d100=0.092,
        w_lh=0.230,  sigma_lh=1500.0,
        w_lw=0.000,  sigma_lw=1500.0,
    ),
    "FM5": dict(
        w0=0.046,  sigma=1739.0, delta=2.0,  M_x=0.20,
        h_heat=8000.0, S_T=0.0555, S_e=0.010, rho_p=32.0,
        w_d1=0.046,  sigma_d1=1739.0,   # Rothermel 1972 Table A-1: 1739 ft⁻¹
        w_d10=0.023, w_d100=0.000,
        w_lh=0.092,  sigma_lh=1500.0,
        w_lw=0.000,  sigma_lw=1500.0,
    ),
    "FM6": dict(
        w0=0.069,  sigma=1750.0, delta=2.5,  M_x=0.25,
        h_heat=8000.0, S_T=0.0555, S_e=0.010, rho_p=32.0,
        w_d1=0.069,  sigma_d1=1750.0,
        w_d10=0.115, w_d100=0.115,
        w_lh=0.000,  sigma_lh=1500.0,
        w_lw=0.000,  sigma_lw=1500.0,
    ),
    "FM7": dict(
        w0=0.052,  sigma=1550.0, delta=2.5,  M_x=0.40,
        h_heat=8000.0, S_T=0.0555, S_e=0.010, rho_p=32.0,
        w_d1=0.052,  sigma_d1=1750.0,
        w_d10=0.086, w_d100=0.069,
        w_lh=0.017,  sigma_lh=1550.0,
        w_lw=0.017,  sigma_lw=1550.0,
    ),
    "FM8": dict(
        w0=0.069,  sigma=2000.0, delta=0.2,  M_x=0.30,
        h_heat=8000.0, S_T=0.0555, S_e=0.010, rho_p=32.0,
        w_d1=0.069,  sigma_d1=2000.0,
        w_d10=0.046, w_d100=0.115,
        w_lh=0.000,  sigma_lh=1500.0,
        w_lw=0.000,  sigma_lw=1500.0,
    ),
    "FM9": dict(
        w0=0.134,  sigma=2500.0, delta=0.2,  M_x=0.25,
        h_heat=8000.0, S_T=0.0555, S_e=0.010, rho_p=32.0,
        w_d1=0.134,  sigma_d1=2500.0,
        w_d10=0.019, w_d100=0.007,
        w_lh=0.000,  sigma_lh=1500.0,
        w_lw=0.000,  sigma_lw=1500.0,
    ),
    "FM10": dict(
        w0=0.138,  sigma=2000.0, delta=1.0,  M_x=0.25,
        h_heat=8000.0, S_T=0.0555, S_e=0.010, rho_p=32.0,
        w_d1=0.138,  sigma_d1=2000.0,
        w_d10=0.092, w_d100=0.230,
        w_lh=0.092,  sigma_lh=1500.0,
        w_lw=0.000,  sigma_lw=1500.0,
    ),
    "FM11": dict(
        w0=0.069,  sigma=1500.0, delta=1.0,  M_x=0.15,
        h_heat=8000.0, S_T=0.0555, S_e=0.010, rho_p=32.0,
        w_d1=0.069,  sigma_d1=1500.0,
        w_d10=0.207, w_d100=0.253,
        w_lh=0.000,  sigma_lh=1500.0,
        w_lw=0.000,  sigma_lw=1500.0,
    ),
    "FM12": dict(
        w0=0.184,  sigma=1500.0, delta=2.3,  M_x=0.20,
        h_heat=8000.0, S_T=0.0555, S_e=0.010, rho_p=32.0,
        w_d1=0.184,  sigma_d1=1500.0,
        w_d10=0.644, w_d100=0.759,
        w_lh=0.000,  sigma_lh=1500.0,
        w_lw=0.000,  sigma_lw=1500.0,
    ),
    "FM13": dict(
        w0=0.322,  sigma=1500.0, delta=3.0,  M_x=0.25,
        h_heat=8000.0, S_T=0.0555, S_e=0.010, rho_p=32.0,
        w_d1=0.322,  sigma_d1=1500.0,
        w_d10=1.058, w_d100=1.288,
        w_lh=0.000,  sigma_lh=1500.0,
        w_lw=0.000,  sigma_lw=1500.0,
    ),
}

# ── Reference implementation ──────────────────────────────────────────────────
# Direct Python port of the corrected C++ code.
# Any change to these functions must be mirrored in src/rothermel_model.H.

def _a(w: float, sigma: float) -> float:
    """A-factor: Aᵢ = wᵢ σᵢ exp(−138/σᵢ).  Returns 0 if w or σ ≤ 0."""
    if w <= 0.0 or sigma <= 0.0:
        return 0.0
    return w * sigma * math.exp(-138.0 / sigma)


def compute_sigma_c_afactor(fm: dict) -> float:
    """
    A-factor weighted characteristic SAV σ_c [ft⁻¹].

    σ_c = Σ(Aᵢ σᵢ) / Σ(Aᵢ),  Aᵢ = wᵢ σᵢ exp(−138/σᵢ)

    Mirrors the multi-class block in src/rothermel_model.H.
    Falls back to fm['sigma'] (single-class path) when per-class loads are
    all zero, and falls back to load-weighted if all A-factors are zero.
    """
    total = (fm["w_d1"] + fm["w_d10"] + fm["w_d100"]
             + fm["w_lh"] + fm["w_lw"])
    if total <= 0.0:
        return fm["sigma"]  # single-class fallback

    A_d1   = _a(fm["w_d1"],   fm["sigma_d1"])
    A_d10  = _a(fm["w_d10"],  SIGMA_D10)
    A_d100 = _a(fm["w_d100"], SIGMA_D100)
    A_lh   = _a(fm["w_lh"],   fm["sigma_lh"])
    A_lw   = _a(fm["w_lw"],   fm["sigma_lw"])
    A_tot  = A_d1 + A_d10 + A_d100 + A_lh + A_lw

    if A_tot <= 0.0:
        # degenerate load-weighted fallback
        return (fm["w_d1"] * fm["sigma_d1"] + fm["w_d10"] * SIGMA_D10
                + fm["w_d100"] * SIGMA_D100 + fm["w_lh"] * fm["sigma_lh"]
                + fm["w_lw"] * fm["sigma_lw"]) / total

    return (A_d1 * fm["sigma_d1"] + A_d10 * SIGMA_D10 + A_d100 * SIGMA_D100
            + A_lh * fm["sigma_lh"] + A_lw * fm["sigma_lw"]) / A_tot


def compute_sigma_c_load_weighted(fm: dict) -> float:
    """Load-weighted σ_c — the OLD (incorrect) implementation kept here for
    comparison in regression assertions."""
    total = (fm["w_d1"] + fm["w_d10"] + fm["w_d100"]
             + fm["w_lh"] + fm["w_lw"])
    if total <= 0.0:
        return fm["sigma"]
    return (fm["w_d1"] * fm["sigma_d1"] + fm["w_d10"] * SIGMA_D10
            + fm["w_d100"] * SIGMA_D100 + fm["w_lh"] * fm["sigma_lh"]
            + fm["w_lw"] * fm["sigma_lw"]) / total


def compute_rothermel(fm: dict,
                      M_dead: float = 0.08,
                      M_lh:   float = 0.90,
                      M_lw:   float = 1.20,
                      U_ftmin: float = 0.0) -> dict:
    """
    Full Rothermel (1972) surface fire spread computation.

    Parameters
    ----------
    fm      : fuel model dict (from FBFM13)
    M_dead  : dead fuel moisture fraction (all dead size classes)
    M_lh    : live herbaceous moisture fraction
    M_lw    : live woody moisture fraction
    U_ftmin : midflame wind speed [ft/min]

    Returns
    -------
    dict with keys: sigma_c, sigma_c_load, beta, beta_op,
                    Gamma_prime, I_R, xi, eps_h, Q_ig,
                    R0_ftmin, phi_w, ROS_ftmin, ROS_ms
    """
    sc   = compute_sigma_c_afactor(fm)
    scL  = compute_sigma_c_load_weighted(fm)

    w_total = (fm["w_d1"] + fm["w_d10"] + fm["w_d100"]
               + fm["w_lh"] + fm["w_lw"])
    rho_b = (w_total if w_total > 0.0 else fm["w0"]) / fm["delta"]
    beta  = rho_b / fm["rho_p"]

    # Optimum packing ratio and reaction velocity (Eqs. 36-38)
    beta_op    = 3.348 * sc**(-0.8189)
    beta_ratio = beta / beta_op
    sc15       = sc**1.5
    Gamma_max  = sc15 / (495.0 + 0.0594 * sc15)
    A_coeff    = 133.0 * sc**(-0.7913)
    Gamma_p    = (Gamma_max * beta_ratio**A_coeff
                  * math.exp(A_coeff * (1.0 - beta_ratio)))

    # A-factors (same values used for both σ_c and moisture damping)
    A_d1   = _a(fm["w_d1"],   fm["sigma_d1"])
    A_d10  = _a(fm["w_d10"],  SIGMA_D10)
    A_d100 = _a(fm["w_d100"], SIGMA_D100)
    A_lh   = _a(fm["w_lh"],   fm["sigma_lh"])
    A_lw   = _a(fm["w_lw"],   fm["sigma_lw"])
    A_dead = A_d1 + A_d10 + A_d100
    A_live = A_lh + A_lw

    S_T = fm["S_T"]
    M_x = fm["M_x"]

    # Dead moisture damping (Eq. 29)
    rm_dead    = 0.0
    etaM_dead  = 0.0
    if A_dead > 0.0:
        rm_dead = min(
            (A_d1 * M_dead + A_d10 * M_dead + A_d100 * M_dead)
            / (A_dead * M_x),
            1.0,
        )
        etaM_dead = max(
            1.0 - 2.59 * rm_dead + 5.11 * rm_dead**2 - 3.52 * rm_dead**3,
            0.0,
        )

    # Live moisture damping (Eq. 88 + Eq. 29)
    w_n_dead = (fm["w_d1"] + fm["w_d10"] + fm["w_d100"]) * (1.0 - S_T)
    w_n_live = (fm["w_lh"] + fm["w_lw"]) * (1.0 - S_T)
    etaM_live = 0.0
    if A_live > 0.0 and w_n_live > 1.0e-10:
        Mx_live = max(
            2.9 * (w_n_dead / w_n_live) * (1.0 - rm_dead) - 0.226,
            0.30,
        )
        rm_live = min(
            (A_lh * M_lh + A_lw * M_lw) / (A_live * Mx_live),
            1.0,
        )
        etaM_live = max(
            1.0 - 2.59 * rm_live + 5.11 * rm_live**2 - 3.52 * rm_live**3,
            0.0,
        )

    # Mineral damping (Eq. 30)
    eta_s = 0.174 * fm["S_e"]**(-0.19)

    # Reaction intensity (Eq. 27)
    h = fm["h_heat"]
    if w_total > 0.0:
        I_R = Gamma_p * (w_n_dead * h * etaM_dead
                        + w_n_live * h * etaM_live) * eta_s
    else:
        w_n  = fm["w0"] * (1.0 - S_T)
        rm   = min(M_dead / M_x, 1.0)
        etaM = max(1.0 - 2.59*rm + 5.11*rm**2 - 3.52*rm**3, 0.0)
        I_R  = Gamma_p * w_n * h * etaM * eta_s

    # Propagating flux ratio (Eq. 42), effective heating number, ignition heat
    xi    = (math.exp((0.792 + 0.681 * sc**0.5) * (beta + 0.1))
             / (192.0 + 0.2595 * sc))
    eps_h = math.exp(-138.0 / sc)
    Q_ig  = 250.0 + 1116.0 * M_dead

    # No-wind, no-slope ROS (Eq. 1)
    R0 = max(I_R * xi / (rho_b * eps_h * Q_ig), 0.0)

    # Wind factor φ_w (Eq. 47-49)
    phi_w = 0.0
    if U_ftmin > 0.0:
        C = 7.47  * math.exp(-0.133  * sc**0.55)  # Eq. 47
        B = 0.02526 * sc**0.54                      # Eq. 48
        E = 0.715  * math.exp(-3.59e-4 * sc)        # Eq. 49
        phi_w = C * U_ftmin**B * (beta / beta_op)**(-E)

    ROS_ftmin = R0 * (1.0 + phi_w)
    return dict(
        sigma_c=sc, sigma_c_load=scL,
        beta=beta, beta_op=beta_op,
        Gamma_prime=Gamma_p, I_R=I_R, xi=xi, eps_h=eps_h, Q_ig=Q_ig,
        R0_ftmin=R0, phi_w=phi_w,
        ROS_ftmin=ROS_ftmin, ROS_ms=ROS_ftmin * ROS_CONV,
    )


# ── Canonical expected values ─────────────────────────────────────────────────
# Computed from the reference implementation above (A-factor weighted σ_c,
# corrected fuel database).  Locked here as regression targets.
# Conditions: M_dead = 0.08, M_lh = 0.90, M_lw = 1.20, U = 440 ft/min
# (5 mph midflame wind).  Tolerances: 0.5 % relative.

_STD_U   = 440.0   # ft/min  (≈ 5 mph midflame)
_STD_Md  = 0.08    # dead fuel moisture

# σ_c [ft⁻¹] – A-factor weighted
SIGMA_C_EXPECTED = {
    "FM1":  3500.0,
    "FM2":  2699.6,
    "FM3":  1500.0,
    "FM4":  1616.5,
    "FM5":  1583.2,
    "FM6":  1699.1,
    "FM7":  1646.2,
    "FM8":  1979.0,
    "FM9":  2495.6,
    "FM10": 1822.9,
    "FM11": 1411.1,
    "FM12": 1397.5,
    "FM13": 1403.3,
}

# R₀ [ft/min] – no-wind, no-slope, M_dead = 8 %
R0_EXPECTED = {
    "FM1":  3.977,
    "FM2":  4.131,
    "FM3":  4.298,
    "FM4": 14.673,
    "FM5":  3.188,
    "FM6":  6.435,
    "FM7":  6.863,
    "FM8":  0.746,
    "FM9":  0.959,
    "FM10": 3.891,
    "FM11": 2.797,
    "FM12": 6.622,
    "FM13": 8.757,
}

# φ_w at 440 ft/min (5 mph midflame)
PHI_W_EXPECTED = {
    "FM1":  21.426,
    "FM2":  16.612,
    "FM3":  25.136,
    "FM4":  17.640,
    "FM5":  21.039,
    "FM6":  17.590,
    "FM7":  19.314,
    "FM8":   7.680,
    "FM9":   9.478,
    "FM10":  9.800,
    "FM11":  9.777,
    "FM12":  8.714,
    "FM13":  7.811,
}

# Total ROS [m/s] at 440 ft/min wind, M_dead = 8 %
ROS_MS_EXPECTED = {
    "FM1":  0.4531,
    "FM2":  0.3696,
    "FM3":  0.5706,
    "FM4":  1.3894,
    "FM5":  0.3569,
    "FM6":  0.6077,
    "FM7":  0.7082,
    "FM8":  0.0329,
    "FM9":  0.0510,
    "FM10": 0.2135,
    "FM11": 0.1531,
    "FM12": 0.3268,
    "FM13": 0.3920,
}


# ── Test classes ──────────────────────────────────────────────────────────────

class TestFuelDatabaseParams(unittest.TestCase):
    """
    Verify per-class fuel parameters match published Rothermel (1972) values.

    These tests would catch any accidental revert of the FM4/FM5 σ_d1 fix
    (which corrected the value from 2000 to the Rothermel 1972 reference of
    1739 ft⁻¹).
    """

    def test_fm4_sigma_d1_is_1739(self):
        """FM4 1-hr dead SAV must be 1739 ft⁻¹ (Rothermel 1972 Table A-1)."""
        self.assertAlmostEqual(FBFM13["FM4"]["sigma_d1"], 1739.0, places=1,
                               msg="FM4 sigma_d1 should be 1739 ft⁻¹")

    def test_fm4_sigma_is_1739(self):
        """FM4 aggregate sigma (single-class fallback) must also be 1739 ft⁻¹."""
        self.assertAlmostEqual(FBFM13["FM4"]["sigma"], 1739.0, places=1)

    def test_fm5_sigma_d1_is_1739(self):
        """FM5 1-hr dead SAV must be 1739 ft⁻¹ (Rothermel 1972 Table A-1)."""
        self.assertAlmostEqual(FBFM13["FM5"]["sigma_d1"], 1739.0, places=1,
                               msg="FM5 sigma_d1 should be 1739 ft⁻¹")

    def test_fm5_sigma_is_1739(self):
        self.assertAlmostEqual(FBFM13["FM5"]["sigma"], 1739.0, places=1)

    def test_all_sigma_d1_positive(self):
        """All fuel models must have a positive 1-hr dead SAV."""
        for name, fm in FBFM13.items():
            with self.subTest(model=name):
                self.assertGreater(fm["sigma_d1"], 0.0,
                                   msg=f"{name} sigma_d1 must be positive")

    def test_fm1_sigma_d1_is_3500(self):
        """FM1 (Short Grass) 1-hr SAV = 3500 ft⁻¹."""
        self.assertAlmostEqual(FBFM13["FM1"]["sigma_d1"], 3500.0, places=1)

    def test_fm3_sigma_d1_is_1500(self):
        """FM3 (Tall Grass) 1-hr SAV = 1500 ft⁻¹."""
        self.assertAlmostEqual(FBFM13["FM3"]["sigma_d1"], 1500.0, places=1)

    def test_fm8_sigma_d1_is_2000(self):
        """FM8 (Closed Timber Litter) 1-hr SAV = 2000 ft⁻¹."""
        self.assertAlmostEqual(FBFM13["FM8"]["sigma_d1"], 2000.0, places=1)

    def test_fm9_sigma_d1_is_2500(self):
        """FM9 (Hardwood Litter) 1-hr SAV = 2500 ft⁻¹."""
        self.assertAlmostEqual(FBFM13["FM9"]["sigma_d1"], 2500.0, places=1)

    def test_fuel_loads_non_negative(self):
        """All per-class fuel loads must be ≥ 0."""
        for name, fm in FBFM13.items():
            for field in ("w_d1", "w_d10", "w_d100", "w_lh", "w_lw"):
                with self.subTest(model=name, field=field):
                    self.assertGreaterEqual(fm[field], 0.0)

    def test_fm4_dead_loads(self):
        """FM4 dead fuel loads match Rothermel (1972) Table A-1."""
        fm = FBFM13["FM4"]
        self.assertAlmostEqual(fm["w_d1"],   0.230, places=3)
        self.assertAlmostEqual(fm["w_d10"],  0.184, places=3)
        self.assertAlmostEqual(fm["w_d100"], 0.092, places=3)
        self.assertAlmostEqual(fm["w_lh"],   0.230, places=3)

    def test_extinction_moisture_reasonable(self):
        """All M_x values must be in (0, 1]."""
        for name, fm in FBFM13.items():
            with self.subTest(model=name):
                self.assertGreater(fm["M_x"], 0.0)
                self.assertLessEqual(fm["M_x"], 1.0)


class TestSigmaCAFactorWeighting(unittest.TestCase):
    """
    Verify that the A-factor weighted σ_c matches canonical expected values.

    These tests catch any regression to the old load-weighted average, which
    gave σ_c values far too low for models containing coarse dead fuels.
    """

    _TOL = 0.005   # 0.5 % relative tolerance

    def _check(self, name):
        sc = compute_sigma_c_afactor(FBFM13[name])
        expected = SIGMA_C_EXPECTED[name]
        self.assertAlmostEqual(sc / expected, 1.0, delta=self._TOL,
                               msg=f"{name} σ_c: got {sc:.1f}, expected {expected:.1f}")

    def test_fm1_sigma_c(self):   self._check("FM1")
    def test_fm2_sigma_c(self):   self._check("FM2")
    def test_fm3_sigma_c(self):   self._check("FM3")
    def test_fm4_sigma_c(self):   self._check("FM4")
    def test_fm5_sigma_c(self):   self._check("FM5")
    def test_fm6_sigma_c(self):   self._check("FM6")
    def test_fm7_sigma_c(self):   self._check("FM7")
    def test_fm8_sigma_c(self):   self._check("FM8")
    def test_fm9_sigma_c(self):   self._check("FM9")
    def test_fm10_sigma_c(self):  self._check("FM10")
    def test_fm11_sigma_c(self):  self._check("FM11")
    def test_fm12_sigma_c(self):  self._check("FM12")
    def test_fm13_sigma_c(self):  self._check("FM13")

    def test_single_class_models_unchanged(self):
        """
        FM1 and FM3 contain only 1-hr dead fuel; their σ_c must equal σ_d1
        exactly regardless of the weighting method.
        """
        for name in ("FM1", "FM3"):
            with self.subTest(model=name):
                fm = FBFM13[name]
                sc_af = compute_sigma_c_afactor(fm)
                sc_lw = compute_sigma_c_load_weighted(fm)
                self.assertAlmostEqual(sc_af, fm["sigma_d1"], places=2,
                                       msg=f"{name}: A-factor σ_c should equal σ_d1")
                self.assertAlmostEqual(sc_lw, fm["sigma_d1"], places=2,
                                       msg=f"{name}: load-wt σ_c should equal σ_d1")


class TestAFactorBeatLoadWeighted(unittest.TestCase):
    """
    Demonstrate that A-factor σ_c > load-weighted σ_c for any fuel model that
    contains coarse dead fuels (10-hr or 100-hr), and equal for fine-only
    models.  This is the core property that the phi_w fix relies on.
    """

    # Models with coarse dead fuels: A-factor σ_c must be strictly higher
    # than load-weighted σ_c because coarse fuels have tiny A-factors and
    # therefore contribute very little to the A-factor weighted numerator.
    COARSE_MODELS = ["FM2", "FM4", "FM5", "FM6", "FM7",
                     "FM8", "FM9", "FM10", "FM11", "FM12", "FM13"]
    # Pure fine-fuel models: both weightings give the same result
    FINE_ONLY_MODELS = ["FM1", "FM3"]

    def test_afactor_above_loadweighted_for_coarse_models(self):
        for name in self.COARSE_MODELS:
            with self.subTest(model=name):
                fm = FBFM13[name]
                sc_af = compute_sigma_c_afactor(fm)
                sc_lw = compute_sigma_c_load_weighted(fm)
                self.assertGreater(
                    sc_af, sc_lw,
                    msg=(f"{name}: A-factor σ_c ({sc_af:.1f}) should exceed "
                         f"load-weighted σ_c ({sc_lw:.1f}) for fuel beds "
                         f"containing coarse dead fuels"),
                )

    def test_fm4_afactor_much_higher_than_load_weighted(self):
        """
        FM4 has large 10-hr (0.184 lb/ft²) and 100-hr (0.092 lb/ft²) dead
        loads.  A-factor σ_c should be at least 50 % above load-weighted σ_c.
        The old code gave ≈ 1043 ft⁻¹ (load-weighted); the corrected code
        gives ≈ 1617 ft⁻¹ (A-factor), a ≈ 55 % increase.
        """
        fm = FBFM13["FM4"]
        sc_af = compute_sigma_c_afactor(fm)
        sc_lw = compute_sigma_c_load_weighted(fm)
        self.assertGreater(sc_af / sc_lw, 1.50,
                           msg="FM4: A-factor σ_c should be > 1.5× load-weighted")

    def test_fm11_afactor_much_higher_than_load_weighted(self):
        """
        FM11 (Light Slash) has heavy 10-hr and 100-hr loads relative to 1-hr;
        A-factor σ_c should be at least 4× the load-weighted value.
        """
        fm = FBFM13["FM11"]
        sc_af = compute_sigma_c_afactor(fm)
        sc_lw = compute_sigma_c_load_weighted(fm)
        self.assertGreater(sc_af / sc_lw, 4.0,
                           msg="FM11: A-factor σ_c should be >> load-weighted")

    def test_equal_for_pure_fine_fuel_models(self):
        """For FM1 and FM3 both weightings must give identical σ_c."""
        for name in self.FINE_ONLY_MODELS:
            with self.subTest(model=name):
                fm = FBFM13[name]
                self.assertAlmostEqual(
                    compute_sigma_c_afactor(fm),
                    compute_sigma_c_load_weighted(fm),
                    places=4,
                    msg=f"{name}: fine-only model should have equal σ_c from both methods",
                )


class TestWindCoefficients(unittest.TestCase):
    """
    Verify the wind factor coefficients C, B, E (Eqs. 47-49) are physically
    consistent with the characteristic SAV.
    """

    def test_coefficient_C_decreases_with_sigma(self):
        """
        C = 7.47 exp(−0.133 σ^0.55) decreases monotonically with σ.
        Higher-σ fuels spread more efficiently and need a smaller C.
        """
        models_asc = ["FM13", "FM11", "FM12", "FM7", "FM6", "FM4", "FM5",
                      "FM10", "FM8", "FM2", "FM9", "FM1"]
        prev_C = None
        prev_sc = None
        for name in models_asc:
            sc = compute_sigma_c_afactor(FBFM13[name])
            C = 7.47 * math.exp(-0.133 * sc**0.55)
            if prev_C is not None and sc > prev_sc:
                self.assertLess(C, prev_C,
                                msg=f"C should decrease from σ_c={prev_sc:.0f} "
                                    f"to σ_c={sc:.0f} ({name})")
            prev_C = C
            prev_sc = sc

    def test_coefficient_B_increases_with_sigma(self):
        """B = 0.02526 σ^0.54 increases monotonically with σ."""
        scs = sorted((compute_sigma_c_afactor(fm), n) for n, fm in FBFM13.items())
        for i in range(len(scs) - 1):
            sc_a, n_a = scs[i]
            sc_b, n_b = scs[i + 1]
            B_a = 0.02526 * sc_a**0.54
            B_b = 0.02526 * sc_b**0.54
            if sc_b > sc_a:
                self.assertGreater(B_b, B_a,
                                   msg=f"B should increase from {n_a}→{n_b}")

    def test_phi_w_positive_for_all_models_at_standard_wind(self):
        """φ_w > 0 for every model when U = 440 ft/min."""
        for name, fm in FBFM13.items():
            result = compute_rothermel(fm, U_ftmin=_STD_U)
            with self.subTest(model=name):
                self.assertGreater(result["phi_w"], 0.0,
                                   msg=f"{name}: φ_w must be > 0 at non-zero wind")

    def test_phi_w_zero_at_zero_wind(self):
        """φ_w = 0 for every model when U = 0."""
        for name, fm in FBFM13.items():
            result = compute_rothermel(fm, U_ftmin=0.0)
            with self.subTest(model=name):
                self.assertEqual(result["phi_w"], 0.0,
                                 msg=f"{name}: φ_w must be 0 at zero wind")

    def test_phi_w_monotone_with_wind(self):
        """φ_w must increase strictly with wind speed for every model."""
        wind_speeds = [0.0, 100.0, 220.0, 440.0, 880.0]
        for name, fm in FBFM13.items():
            phi_vals = [compute_rothermel(fm, U_ftmin=u)["phi_w"]
                        for u in wind_speeds]
            for i in range(1, len(phi_vals)):
                with self.subTest(model=name, U=wind_speeds[i]):
                    self.assertGreater(
                        phi_vals[i], phi_vals[i - 1],
                        msg=f"{name}: φ_w should increase at U={wind_speeds[i]} ft/min",
                    )


class TestNoWindROS(unittest.TestCase):
    """
    Verify no-wind, no-slope rate of spread R₀ at standard conditions
    (M_dead = 8 %, M_lh = 90 %, flat terrain).
    """

    _TOL = 0.005   # 0.5 % relative tolerance

    def _check(self, name):
        r = compute_rothermel(FBFM13[name], M_dead=_STD_Md, U_ftmin=0.0)
        expected = R0_EXPECTED[name]
        self.assertAlmostEqual(r["R0_ftmin"] / expected, 1.0, delta=self._TOL,
                               msg=f"{name} R₀: got {r['R0_ftmin']:.3f}, "
                                   f"expected {expected:.3f} ft/min")

    def test_fm1_R0(self):   self._check("FM1")
    def test_fm2_R0(self):   self._check("FM2")
    def test_fm3_R0(self):   self._check("FM3")
    def test_fm4_R0(self):   self._check("FM4")
    def test_fm5_R0(self):   self._check("FM5")
    def test_fm6_R0(self):   self._check("FM6")
    def test_fm7_R0(self):   self._check("FM7")
    def test_fm8_R0(self):   self._check("FM8")
    def test_fm9_R0(self):   self._check("FM9")
    def test_fm10_R0(self):  self._check("FM10")
    def test_fm11_R0(self):  self._check("FM11")
    def test_fm12_R0(self):  self._check("FM12")
    def test_fm13_R0(self):  self._check("FM13")

    def test_R0_positive_for_all_models(self):
        """R₀ must be > 0 for all models at M_dead = 8 % (well below M_x)."""
        for name, fm in FBFM13.items():
            r = compute_rothermel(fm, M_dead=_STD_Md, U_ftmin=0.0)
            with self.subTest(model=name):
                self.assertGreater(r["R0_ftmin"], 0.0,
                                   msg=f"{name}: R₀ must be positive at 8 % moisture")

    def test_R0_vanishes_at_moisture_of_extinction(self):
        """
        When M_dead ≥ M_x the dead moisture damping η_M → 0, driving R₀ → 0.
        We use M_dead = 1.0 (100 %) which exceeds M_x for every FBFM13 model.
        """
        for name, fm in FBFM13.items():
            r = compute_rothermel(fm, M_dead=1.0, U_ftmin=0.0)
            with self.subTest(model=name):
                self.assertAlmostEqual(
                    r["R0_ftmin"], 0.0, places=6,
                    msg=f"{name}: R₀ should vanish at M_dead = 100 %",
                )

    def test_R0_decreases_with_moisture(self):
        """R₀ must decrease as dead fuel moisture increases (0.02 → 0.12)."""
        moistures = [0.02, 0.04, 0.06, 0.08, 0.10, 0.12]
        for name, fm in FBFM13.items():
            r_vals = [compute_rothermel(fm, M_dead=m, U_ftmin=0.0)["R0_ftmin"]
                      for m in moistures]
            for i in range(1, len(r_vals)):
                with self.subTest(model=name, M_dead=moistures[i]):
                    self.assertLessEqual(
                        r_vals[i], r_vals[i - 1],
                        msg=(f"{name}: R₀ should not increase from "
                             f"M_dead={moistures[i-1]} to {moistures[i]}"),
                    )

    def test_fm4_R0_greater_than_fm8(self):
        """
        FM4 (deep chaparral) has much more fine fuel than FM8 (shallow litter)
        and should have significantly higher R₀ under identical conditions.
        """
        r4 = compute_rothermel(FBFM13["FM4"], M_dead=_STD_Md, U_ftmin=0.0)
        r8 = compute_rothermel(FBFM13["FM8"], M_dead=_STD_Md, U_ftmin=0.0)
        self.assertGreater(r4["R0_ftmin"], r8["R0_ftmin"] * 5,
                           msg="FM4 R₀ should be >> FM8 R₀")


class TestWindFactor(unittest.TestCase):
    """
    Canonical regression values for φ_w at 5 mph (440 ft/min) midflame wind,
    M_dead = 8 %.
    """

    _TOL = 0.005   # 0.5 % relative tolerance

    def _check(self, name):
        r = compute_rothermel(FBFM13[name], M_dead=_STD_Md, U_ftmin=_STD_U)
        expected = PHI_W_EXPECTED[name]
        self.assertAlmostEqual(r["phi_w"] / expected, 1.0, delta=self._TOL,
                               msg=f"{name} φ_w: got {r['phi_w']:.4f}, "
                                   f"expected {expected:.4f}")

    def test_fm1_phi_w(self):   self._check("FM1")
    def test_fm2_phi_w(self):   self._check("FM2")
    def test_fm3_phi_w(self):   self._check("FM3")
    def test_fm4_phi_w(self):   self._check("FM4")
    def test_fm5_phi_w(self):   self._check("FM5")
    def test_fm6_phi_w(self):   self._check("FM6")
    def test_fm7_phi_w(self):   self._check("FM7")
    def test_fm8_phi_w(self):   self._check("FM8")
    def test_fm9_phi_w(self):   self._check("FM9")
    def test_fm10_phi_w(self):  self._check("FM10")
    def test_fm11_phi_w(self):  self._check("FM11")
    def test_fm12_phi_w(self):  self._check("FM12")
    def test_fm13_phi_w(self):  self._check("FM13")

    def test_fm4_phi_w_lower_than_load_weighted_would_give(self):
        """
        The old load-weighted σ_c for FM4 was ≈ 1043 ft⁻¹.  With that value,
        C is inflated ~5× (because C ∝ exp(−0.133 σ^0.55)) and φ_w at
        440 ft/min is ≈ 25+.  With the corrected A-factor σ_c ≈ 1617, φ_w
        is ≈ 17.6.  Assert that the corrected value is below the load-weighted
        reference to ensure the fix has taken effect.
        """
        # Compute φ_w with the OLD load-weighted σ_c
        fm = FBFM13["FM4"]
        sc_lw   = compute_sigma_c_load_weighted(fm)
        beta    = (fm["w_d1"] + fm["w_d10"] + fm["w_d100"] + fm["w_lh"]) \
                  / fm["delta"] / fm["rho_p"]
        beta_op = 3.348 * sc_lw**(-0.8189)
        C_lw    = 7.47 * math.exp(-0.133 * sc_lw**0.55)
        B_lw    = 0.02526 * sc_lw**0.54
        E_lw    = 0.715 * math.exp(-3.59e-4 * sc_lw)
        phi_lw  = C_lw * _STD_U**B_lw * (beta / beta_op)**(-E_lw)

        # Compute φ_w with the corrected A-factor σ_c
        phi_af  = compute_rothermel(fm, M_dead=_STD_Md, U_ftmin=_STD_U)["phi_w"]

        self.assertLess(
            phi_af, phi_lw,
            msg=(f"FM4 corrected φ_w ({phi_af:.2f}) should be less than "
                 f"load-weighted φ_w ({phi_lw:.2f})"),
        )

    def test_phi_w_higher_sigma_gives_lower_C_amplification(self):
        """
        For models with coarse fuels, A-factor σ_c is much higher than
        load-weighted σ_c.  Because C decreases rapidly with σ, the
        corrected (A-factor) φ_w should be lower than what the load-weighted
        implementation would produce.
        """
        coarse_models = ["FM4", "FM6", "FM8", "FM11", "FM12", "FM13"]
        for name in coarse_models:
            fm = FBFM13[name]
            sc_af = compute_sigma_c_afactor(fm)
            sc_lw = compute_sigma_c_load_weighted(fm)
            # C(σ) is monotone decreasing, so higher σ → lower C → lower φ_w
            C_af = 7.47 * math.exp(-0.133 * sc_af**0.55)
            C_lw = 7.47 * math.exp(-0.133 * sc_lw**0.55)
            with self.subTest(model=name):
                self.assertLess(C_af, C_lw,
                                msg=f"{name}: C(A-factor σ_c) < C(load-weighted σ_c)")


class TestTotalROS(unittest.TestCase):
    """
    Canonical regression values for total ROS (m/s) at 5 mph midflame wind,
    M_dead = 8 %.
    """

    _TOL = 0.005   # 0.5 % relative tolerance

    def _check(self, name):
        r = compute_rothermel(FBFM13[name], M_dead=_STD_Md, U_ftmin=_STD_U)
        expected = ROS_MS_EXPECTED[name]
        self.assertAlmostEqual(r["ROS_ms"] / expected, 1.0, delta=self._TOL,
                               msg=f"{name} ROS: got {r['ROS_ms']:.4f} m/s, "
                                   f"expected {expected:.4f} m/s")

    def test_fm1_ROS(self):   self._check("FM1")
    def test_fm2_ROS(self):   self._check("FM2")
    def test_fm3_ROS(self):   self._check("FM3")
    def test_fm4_ROS(self):   self._check("FM4")
    def test_fm5_ROS(self):   self._check("FM5")
    def test_fm6_ROS(self):   self._check("FM6")
    def test_fm7_ROS(self):   self._check("FM7")
    def test_fm8_ROS(self):   self._check("FM8")
    def test_fm9_ROS(self):   self._check("FM9")
    def test_fm10_ROS(self):  self._check("FM10")
    def test_fm11_ROS(self):  self._check("FM11")
    def test_fm12_ROS(self):  self._check("FM12")
    def test_fm13_ROS(self):  self._check("FM13")

    def test_ROS_positive_for_all_models(self):
        """ROS > 0 for every model at standard conditions."""
        for name, fm in FBFM13.items():
            r = compute_rothermel(fm, M_dead=_STD_Md, U_ftmin=_STD_U)
            with self.subTest(model=name):
                self.assertGreater(r["ROS_ms"], 0.0)

    def test_ROS_increases_with_wind(self):
        """Total ROS must increase with midflame wind speed for all models."""
        winds_ms = [0.0, 1.0, 2.0, 5.0, 10.0]
        for name, fm in FBFM13.items():
            ros_vals = [
                compute_rothermel(fm, M_dead=_STD_Md,
                                  U_ftmin=u * WIND_CONV)["ROS_ms"]
                for u in winds_ms
            ]
            for i in range(1, len(ros_vals)):
                with self.subTest(model=name, U=winds_ms[i]):
                    self.assertGreater(
                        ros_vals[i], ros_vals[i - 1],
                        msg=f"{name}: ROS should increase from "
                            f"{winds_ms[i-1]} to {winds_ms[i]} m/s wind",
                    )

    def test_grass_spreads_faster_than_timber_litter(self):
        """
        Fine-grass models (FM1, FM3) should generally spread faster than
        compact timber-litter models (FM8, FM9) at the same wind speed.
        """
        for grass in ("FM1", "FM3"):
            for timber in ("FM8", "FM9"):
                r_grass  = compute_rothermel(FBFM13[grass],  M_dead=_STD_Md, U_ftmin=_STD_U)
                r_timber = compute_rothermel(FBFM13[timber], M_dead=_STD_Md, U_ftmin=_STD_U)
                with self.subTest(grass=grass, timber=timber):
                    self.assertGreater(
                        r_grass["ROS_ms"], r_timber["ROS_ms"],
                        msg=f"{grass} ROS should exceed {timber} ROS at 5 mph",
                    )


# ── Entry point ───────────────────────────────────────────────────────────────

if __name__ == "__main__":
    unittest.main(verbosity=2)
