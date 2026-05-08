// ============================================================================
// Lautenberger (2013) Physics-Based Cellular Fire Spread Model
// ============================================================================
// Implements the semi-physical rate-of-spread formula from:
//
//   Lautenberger, C. (2013). "Wildland fire modeling with an Eulerian level
//   set method and automated calibration." Fire Safety Journal 62:289-298.
//
// Physical Background:
// --------------------
// The model uses a semi-physical ROS formula that accounts for fuel properties,
// fuel moisture, wind speed, and terrain slope:
//
//   R0 = A_L * sigma_m * exp(-B_L * M_f)      [m/s, base ROS at zero wind/slope]
//   phi_W = C_L * U_ms                          [linear wind correction]
//   phi_S = D_L * tan_slope                     [linear slope correction]
//   R = max(R0 * (1 + phi_W + phi_S), 0)        [m/s]
//
// where:
//   sigma_m = SAV ratio in m^-1 (converted from ft^-1: sigma_ft * 3.28084)
//   M_f     = fuel moisture fraction
//   A_L     = calibration pre-factor [m^2/s]   (default: 1.05e-5)
//   B_L     = moisture sensitivity coefficient (default: 2.5)
//   C_L     = wind correction coefficient [(m/s)^-1] (default: 0.40)
//   D_L     = slope correction coefficient [-]  (default: 0.50)
//
// Parameters:
// -----------
//   LautenbergerParams (from parse_inputs.H):
//     A_L – pre-factor [m^2/s]         (default: 1.05e-5)
//     B_L – moisture sensitivity [-]   (default: 2.5)
//     C_L – wind coefficient [(m/s)^-1] (default: 0.40)
//     D_L – slope coefficient [-]      (default: 0.50)
//
//   RothermelParams (from parse_inputs.H):
//     sigma – SAV ratio [ft^-1] (converted internally to m^-1)
//     M_f   – fuel moisture content [fraction]
// ============================================================================

#ifndef LAUTENBERGER_MODEL_H
#define LAUTENBERGER_MODEL_H

#include "constants.H"
#include "parse_inputs.H"
#include "rothermel_model.H"
#include <AMReX.H>
#include <cmath>

// ============================================================================
// Pre-computed Lautenberger parameters for GPU kernel
// ============================================================================
struct LautenbergerComputed {
    amrex::Real R0;    // base ROS [m/s] at zero wind and slope
    amrex::Real C_L;   // wind correction coefficient [(m/s)^-1]
    amrex::Real D_L;   // slope correction coefficient [-]
};

// ============================================================================
// compute_lautenberger_params
//
// Computes LautenbergerComputed from RothermelParams and LautenbergerParams.
// Called on host; result is captured by GPU kernels.
// ============================================================================
inline LautenbergerComputed compute_lautenberger_params(
    const RothermelParams& rp,
    const InputParameters::LautenbergerParams& lp)
{
    LautenbergerComputed lc;

    // Convert SAV ratio from ft^-1 to m^-1
    const amrex::Real sigma_m = rp.sigma * WildfireConst::FT_INV_TO_M_INV;

    // Base ROS: R0 = A_L * sigma_m * exp(-B_L * M_f)
    // Guard against near-zero sigma (non-burnable)
    if (sigma_m > amrex::Real(1.0e-6)) {
        lc.R0 = amrex::Real(lp.A_L) * sigma_m
                * std::exp(-amrex::Real(lp.B_L) * amrex::Real(rp.M_f));
    } else {
        lc.R0 = amrex::Real(0.0);
    }

    lc.C_L = amrex::Real(lp.C_L);
    lc.D_L = amrex::Real(lp.D_L);

    return lc;
}

// ============================================================================
// lautenberger_ros
//
// GPU-device inline: returns ROS [m/s] given wind speed [m/s] and slope.
//
//   phi_W = C_L * U_ms
//   phi_S = D_L * tan_slope
//   R = max(R0 * (1 + phi_W + phi_S), 0)
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real lautenberger_ros(const LautenbergerComputed& lc,
                              amrex::Real wind_speed_ms,
                              amrex::Real tan_slope = amrex::Real(0.0)) noexcept
{
    const amrex::Real phi_W = lc.C_L * wind_speed_ms;
    const amrex::Real phi_S = lc.D_L * tan_slope;
    return amrex::max(lc.R0 * (amrex::Real(1.0) + phi_W + phi_S),
                      amrex::Real(0.0));
}

#endif // LAUTENBERGER_MODEL_H
