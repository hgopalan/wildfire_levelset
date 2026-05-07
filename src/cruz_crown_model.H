// ============================================================================
// Cruz, Alexander & Wakimoto (2005) Crown Fire Spread Model
// ============================================================================
// Implements the simple algebraic crown fire rate-of-spread model from:
//
//   Cruz, M.G., Alexander, M.E., and Wakimoto, R.H. (2005). "Development and
//   testing of models for predicting crown fire rate of spread in conifer
//   forest stands." Canadian Journal of Forest Research, 35(7), 1626–1639.
//   https://doi.org/10.1139/x05-068
//
// This model was calibrated on active crown fire observations from western
// North American conifer forests, including data from the Sierra Nevada,
// Rocky Mountains, and Pacific Northwest of the continental United States,
// as well as western Canadian forests.
//
// Crown fire rate of spread (R, m/min):
//
//   R = 11.02 × U₁₀^0.90 × CBD^0.19 × exp(−0.17 × MC₁₀)
//
// where:
//   U₁₀  – 10-m open wind speed [km/h]
//   CBD   – canopy bulk density [kg/m³]
//   MC₁₀  – 10-h timelag fine fuel moisture content [%]
//
// The formula is fully algebraic (one closed-form expression) and simple
// to evaluate in each grid cell from local wind speed and two canopy/fuel
// state scalars.
//
// Applicable conditions (from Cruz et al. 2005):
//   U₁₀  : 5–60 km/h
//   CBD   : 0.05–0.40 kg/m³ (typical western US conifer canopies)
//   MC₁₀  : 4–25%
//
// Output is converted to m/s for consistency with the simulation's
// internal length/time units.
// ============================================================================

#ifndef CRUZ_CROWN_MODEL_H
#define CRUZ_CROWN_MODEL_H

#include "constants.H"
#include "parse_inputs.H"
#include <AMReX.H>
#include <AMReX_GpuContainers.H>
#include <cmath>

// ============================================================================
// Pre-computed Cruz–Crown parameters for GPU kernel
// ============================================================================
// Fuel/canopy scalars are folded into pre-computed coefficients so the
// GPU kernel needs only the per-cell wind speed as a dynamic input.
// ============================================================================
struct CruzCrownComputed {
    // Combined constant pre-multiplied by CBD and moisture terms
    // For wind-speed evaluation: R_ms = coeff_base * U_kmh^wind_exp
    // where U_kmh is 10-m wind in km/h and R_ms is ROS in m/s.
    amrex::Real coeff_base;  // 11.02 * CBD^0.19 * exp(-0.17 * MC10) / 60.0  [m/s per (km/h)^0.90]
    amrex::Real wind_exp;    // 0.90 (fixed empirical exponent)
};

// ============================================================================
// Compute Cruz–Crown model coefficients from input parameters
//
// Returns CruzCrownComputed ready for the GPU kernel.
// ============================================================================
inline CruzCrownComputed compute_cruz_crown_params(
    const InputParameters::CruzCrownParams& ccp)
{
    CruzCrownComputed ccc;

    // CBD [kg/m³]: clamp to physically valid range for the model
    const amrex::Real CBD = amrex::max(amrex::Real(0.01), ccp.CBD);

    // MC10 [%]: clamp to valid range
    const amrex::Real MC10 = amrex::max(amrex::Real(0.0), ccp.MC10);

    // Cruz et al. (2005) empirical constants (from WildfireConst)
    constexpr amrex::Real A_const   = WildfireConst::CRUZ_A_LEAD;
    constexpr amrex::Real cbd_exp   = WildfireConst::CRUZ_CBD_EXP;
    constexpr amrex::Real mc_coeff  = WildfireConst::CRUZ_MC_COEFF;
    constexpr amrex::Real wind_exp  = WildfireConst::CRUZ_WIND_EXP;
    constexpr amrex::Real min_per_s = amrex::Real(60.0);  // m/min → m/s conversion

    // Pre-compute the canopy/moisture factor:
    //   (11.02 × CBD^0.19 × exp(-0.17 × MC10)) / 60   [m/s per (km/h)^0.90]
    ccc.coeff_base = A_const
                    * std::pow(CBD,  cbd_exp)
                    * std::exp(-mc_coeff * MC10)
                    / min_per_s;

    ccc.wind_exp = wind_exp;

    return ccc;
}

// ============================================================================
// GPU-device inline: evaluate Cruz–Crown ROS [m/s] from wind speed [m/s]
//
// Called once per cell in the GPU kernel.
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real cruz_crown_ros(const CruzCrownComputed& ccc,
                            amrex::Real wind_speed_ms) noexcept
{
    // Convert simulation wind speed (m/s) to 10-m open wind (km/h)
    const amrex::Real U_kmh = wind_speed_ms * amrex::Real(3.6);

    // Cruz et al. (2005): R [m/min] = 11.02 × U₁₀^0.90 × CBD^0.19 × exp(-0.17 × MC₁₀)
    // In pre-computed form (m/s):
    //   R_ms = coeff_base × U_kmh^wind_exp
    // Wind clamp: minimum 1 km/h to ensure non-zero fire spread even in calm
    // conditions (wind is the dominant driver in this model).
    const amrex::Real U_eff = amrex::max(U_kmh, amrex::Real(1.0));
    const amrex::Real R_ms  = ccc.coeff_base * std::pow(U_eff, ccc.wind_exp);

    return amrex::max(R_ms, amrex::Real(0.0));
}

#endif // CRUZ_CROWN_MODEL_H
