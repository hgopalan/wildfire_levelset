// ============================================================================
// Wildfire-AMR: Global Physical and Numerical Constants
// ============================================================================
// All constants are declared as inline constexpr inside the WildfireConst
// namespace.  Using inline constexpr (C++17) guarantees:
//   1. A single definition across all translation units (ODR-safe).
//   2. The value is a compile-time constant that can be used directly inside
//      AMREX_GPU_DEVICE lambdas without any capture or device-storage issues.
// ============================================================================

#ifndef WILDFIRE_CONSTANTS_H
#define WILDFIRE_CONSTANTS_H

#include <AMReX.H>

namespace WildfireConst {

// ============================================================================
// Mathematical constants
// ============================================================================

/// π (pi)
inline constexpr amrex::Real PI     = amrex::Real(3.14159265358979323846);
/// 2π
inline constexpr amrex::Real TWO_PI = amrex::Real(6.283185307179586);

// ============================================================================
// Physical constants
// ============================================================================

/// Gravitational acceleration [m/s²]
inline constexpr amrex::Real GRAVITY = amrex::Real(9.81);

/// Seconds per day [s]
inline constexpr double DAY_S = 86400.0;

// ============================================================================
// Wind speed thresholds
// ============================================================================

/// Maximum plausible wind speed used to clamp outlier cells [m/s]
inline constexpr amrex::Real MAX_WIND_SPEED_MS = amrex::Real(25.0);

/// Square of MAX_WIND_SPEED_MS [m²/s²], used for fast magnitude checks
inline constexpr amrex::Real MAX_WIND_SPEED_SQ =
    MAX_WIND_SPEED_MS * MAX_WIND_SPEED_MS;

// ============================================================================
// Unit conversions: Imperial → SI
// ============================================================================

/// feet → metres
inline constexpr amrex::Real FT_TO_M         = amrex::Real(0.3048);
/// ft⁻¹ → m⁻¹ (surface-area-to-volume ratio)
inline constexpr amrex::Real FT_INV_TO_M_INV = amrex::Real(3.28084);
/// lb/ft² → kg/m² (fuel loading)
inline constexpr amrex::Real LB_FT2_TO_KG_M2 = amrex::Real(4.8824);
/// lb/ft³ → kg/m³ (bulk density)
inline constexpr amrex::Real LB_FT3_TO_KG_M3 = amrex::Real(16.0185);
/// BTU/lb → J/kg (heat of combustion)
inline constexpr amrex::Real BTU_LB_TO_J_KG  = amrex::Real(2326.0);
/// BTU/lb → kJ/kg (heat of combustion, for kW outputs)
inline constexpr amrex::Real BTU_LB_TO_KJ_KG = amrex::Real(2.326);
/// BTU/(ft²·min) → kW/m² (reaction intensity)
inline constexpr amrex::Real BTU_FT2MIN_TO_KW_M2 = amrex::Real(0.18941);
/// BTU/(ft·s) → kW/m (fireline intensity, Albini / Byram)
inline constexpr amrex::Real BTU_FTS_TO_KW_M = amrex::Real(3.459);
/// BTU/(ft·s) → kW/m (fireline intensity, Van Wagner / crown models)
inline constexpr amrex::Real BTU_FTS_TO_KW_M_VW = amrex::Real(3.4613);

// ============================================================================
// Rothermel (1972) fixed surface-area-to-volume ratios [ft⁻¹]
// ============================================================================

/// SAV ratio for 10-hr dead fuel size class
inline constexpr amrex::Real SIGMA_D10  = amrex::Real(109.0);
/// SAV ratio for 100-hr dead fuel size class
inline constexpr amrex::Real SIGMA_D100 = amrex::Real(30.0);

/// Slope factor coefficient (Rothermel 1972 Eq. 51): 5.275
inline constexpr amrex::Real PHI_S_COEFF    = amrex::Real(5.275);
/// Packing-ratio exponent in slope factor (Rothermel 1972 Eq. 51): -0.3
inline constexpr amrex::Real PHI_S_BETA_EXP = amrex::Real(-0.3);

// ============================================================================
// Balbi (2009) model thresholds
// ============================================================================

/// Minimum buoyancy velocity [m/s] — prevents division-by-zero in tilt kernel
inline constexpr amrex::Real BALBI_V_B_MIN    = amrex::Real(1.0e-3);
/// Minimum dimensionless ignition-energy ratio B* — non-burnable floor
inline constexpr amrex::Real BALBI_B_STAR_MIN = amrex::Real(1.0e-6);

// ============================================================================
// Viegas (2004) model thresholds
// ============================================================================

/// Denominator floor for ROS excess calculation
inline constexpr amrex::Real VIEGAS_EPS    = amrex::Real(1.0e-10);
/// Minimum buoyancy velocity [m/s]
inline constexpr amrex::Real VIEGAS_VB_MIN = amrex::Real(1.0e-3);

// ============================================================================
// Weise & Biging (1996) fire whirl model constants
// ============================================================================

/// Empirical flame-tilt coefficient A: tan(θ_wind) = A × Fr^B
inline constexpr amrex::Real WB_A_COEFF     = amrex::Real(0.382);
/// Empirical flame-tilt exponent B
inline constexpr amrex::Real WB_B_COEFF     = amrex::Real(0.432);
/// Minimum whirl radius [m] — prevents division-by-zero
inline constexpr amrex::Real WB_WHIRL_R_MIN = amrex::Real(0.01);

// ============================================================================
// Cruz, Alexander & Wakimoto (2005) crown fire model constants
// ============================================================================

// ============================================================================
// Per-class exponential fuel burnout model (Albini 1976)
// ============================================================================

/// Particle burnout time-constant coefficient [s·ft⁻¹/(lb/ft³)].
/// tau_i [s] = BURNOUT_ALPHA * rho_p [lb/ft³] / sigma_i [ft⁻¹]
/// Calibrated so that 1-hr fine fuel (sigma ≈ 1739 ft⁻¹, rho_p = 32 lb/ft³)
/// has tau ≈ 66 s, and 100-hr fuel (sigma = 30 ft⁻¹) has tau ≈ 3840 s.
inline constexpr amrex::Real BURNOUT_ALPHA = amrex::Real(3600.0);

/// Live-fuel burnout time-constant multiplier relative to 10-hr dead fuel.
/// Live herbaceous burns roughly as fast as 10-hr dead under fire conditions;
/// live woody burns roughly as fast as 100-hr dead.
inline constexpr amrex::Real BURNOUT_LIVE_HERB_MULT = amrex::Real(1.0);   // herb ≈ 10-hr
inline constexpr amrex::Real BURNOUT_LIVE_WOOD_MULT = amrex::Real(10.0);  // wood ≈ 100-hr

// ============================================================================
// Cruz, Alexander & Wakimoto (2005) crown fire model constants
// ============================================================================

/// Lead coefficient in the Cruz et al. (2005) crown ROS formula
inline constexpr amrex::Real CRUZ_A_LEAD   = amrex::Real(11.02);
/// Canopy bulk density exponent
inline constexpr amrex::Real CRUZ_CBD_EXP  = amrex::Real(0.19);
/// Fuel moisture exponent coefficient
inline constexpr amrex::Real CRUZ_MC_COEFF = amrex::Real(0.17);
/// Wind speed exponent
inline constexpr amrex::Real CRUZ_WIND_EXP = amrex::Real(0.90);

// ============================================================================
// Numerical scheme constants
// ============================================================================

/// WENO3 smoothness-indicator regularization parameter
inline constexpr amrex::Real WENO_EPS  = amrex::Real(1.0e-6);
/// Level-set field upper bound used to detect NaN/Inf in WENO stencils
inline constexpr amrex::Real PHI_BOUND = amrex::Real(1.0e6);
/// Artificial viscosity coefficient for level-set advection (Eq. 2)
inline constexpr amrex::Real ADVECTION_VISC = amrex::Real(0.4);

// ============================================================================
// Wind–terrain model numerical thresholds
// ============================================================================

/// Minimum slope magnitude used as a denominator floor in uphill unit-vector
inline constexpr amrex::Real WT_SLOPE_EPS = amrex::Real(1.0e-8);

// ============================================================================
// Canadian FBP System constants (Forestry Canada 1992)
// ============================================================================

/// FBP ISI wind correction exponent coefficient
inline constexpr amrex::Real FBP_ISI_W_COEFF = amrex::Real(0.05039);
/// FBP ISI moisture pre-factor
inline constexpr amrex::Real FBP_ISI_F_PREFAC = amrex::Real(91.9);
/// FBP ISI moisture exponential coefficient
inline constexpr amrex::Real FBP_ISI_F_EXP = amrex::Real(-0.1386);
/// FBP ISI moisture polynomial power
inline constexpr amrex::Real FBP_ISI_F_POW = amrex::Real(5.31);
/// FBP ISI moisture polynomial divisor
inline constexpr amrex::Real FBP_ISI_F_DIV = amrex::Real(4.93e7);
/// FBP ISI overall multiplier
inline constexpr amrex::Real FBP_ISI_MULT = amrex::Real(0.208);

/// FBP O1a lead coefficient [m/s]: 190 m/min / 60
inline constexpr amrex::Real FBP_O1A_C = amrex::Real(190.0 / 60.0);
/// FBP O1a ISI exponent
inline constexpr amrex::Real FBP_O1A_K = amrex::Real(0.031);
/// FBP O1b lead coefficient [m/s]
inline constexpr amrex::Real FBP_O1B_C = amrex::Real(250.0 / 60.0);
/// FBP O1b ISI exponent
inline constexpr amrex::Real FBP_O1B_K = amrex::Real(0.035);

/// FBP S1 lead coefficient [m/s]
inline constexpr amrex::Real FBP_S1_C = amrex::Real(75.0  / 60.0);
inline constexpr amrex::Real FBP_S1_K = amrex::Real(0.110);
/// FBP S2 lead coefficient [m/s]
inline constexpr amrex::Real FBP_S2_C = amrex::Real(200.0 / 60.0);
inline constexpr amrex::Real FBP_S2_K = amrex::Real(0.062);
/// FBP S3 lead coefficient [m/s]
inline constexpr amrex::Real FBP_S3_C = amrex::Real(320.0 / 60.0);
inline constexpr amrex::Real FBP_S3_K = amrex::Real(0.010);

/// FBP curing factor exponent (O1a/O1b)
inline constexpr amrex::Real FBP_O1_CURING_K = amrex::Real(2.1);

// ============================================================================
// Lautenberger (2013) physics-based spread model constants
// ============================================================================

/// Lautenberger (2013) moisture sensitivity coefficient B_L
inline constexpr amrex::Real LAUT_B_L = amrex::Real(2.5);
/// Lautenberger (2013) wind correction coefficient C_L [(m/s)^-1]
inline constexpr amrex::Real LAUT_C_L = amrex::Real(0.40);
/// Lautenberger (2013) slope correction coefficient D_L [-]
inline constexpr amrex::Real LAUT_D_L = amrex::Real(0.50);
/// Lautenberger (2013) pre-factor A_L [m^2/s] (sigma_m x A_L gives [m/s])
inline constexpr amrex::Real LAUT_A_L = amrex::Real(1.05e-5);

// ============================================================================
// Rothermel (1991) crown fire constant
// ============================================================================

/// Rothermel (1991) active crown fire ROS multiplier: R_crown = 3.34 x R_surface
inline constexpr amrex::Real ROTHERMEL1991_CROWN_MULT = amrex::Real(3.34);

// ============================================================================
// Thomas (1963) / Byram (1959) flame length <-> intensity
// ============================================================================

/// Thomas / Byram (1959) SI coefficient: L_f [m] = BYRAM_A * I_B [kW/m]^BYRAM_B
inline constexpr amrex::Real BYRAM_A    = amrex::Real(0.0775);
/// Thomas / Byram (1959) SI exponent
inline constexpr amrex::Real BYRAM_B    = amrex::Real(0.46);
/// Inverse exponent: 1 / BYRAM_B approximately 2.17391
inline constexpr amrex::Real BYRAM_INVB = amrex::Real(1.0 / 0.46);

} // namespace WildfireConst

// ---------------------------------------------------------------------------
// Portability shim: define M_PI if the platform does not supply it.
// New code should prefer WildfireConst::PI; this shim keeps legacy callers
// compiling without changes.
// ---------------------------------------------------------------------------
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#endif // WILDFIRE_CONSTANTS_H
