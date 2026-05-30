// ============================================================================
// Effective Wind Speed (Wind + Slope Combined Vector)
// ============================================================================
// Computes a single diagnostic that combines wind and slope effects for easier
// interpretation and comparison with field observations.
//
// The effective wind speed represents the combined forcing from both
// meteorological wind and topographic slope, expressed as an equivalent
// wind speed that would produce the same fire spread rate on level ground.
//
// References:
//   Rothermel, R.C. (1983). "How to predict the spread and intensity of
//     forest and range fires." General Technical Report INT-143,
//     USDA Forest Service.
//   Sharples, J.J. (2008). "Review of formal fire spread models for
//     operational fire management." International Journal of Wildland Fire,
//     17(1):1-17.
//
// Implementation:
//   The slope effect is converted to an equivalent wind speed based on
//   the Rothermel slope factor φ_s, then vector-combined with the actual
//   wind speed:
//
//     U_eff = sqrt(U_wind² + U_slope_equiv²)
//
//   where U_slope_equiv is derived from φ_s such that it produces the
//   same ROS enhancement as the slope.
// ============================================================================

#ifndef EFFECTIVE_WIND_SPEED_H
#define EFFECTIVE_WIND_SPEED_H

#include "constants.H"
#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_GpuLaunch.H>
#include <cmath>

using namespace amrex;

// ============================================================================
// Convert slope factor to equivalent wind speed
//
// The Rothermel slope factor φ_s produces a fractional increase in ROS.
// This function finds the equivalent wind speed that would produce the
// same fractional increase if applied on level ground.
//
// For Rothermel model: φ_s = 5.275 * β^(-0.3) * tan²(θ)
//
// The equivalent wind is approximated as:
//   U_slope_equiv ≈ sqrt(φ_s) * U_ref
//
// where U_ref is a reference wind speed (typically 5 m/s for normalization).
//
// Parameters:
//   phi_slope  - Rothermel slope factor (dimensionless)
//   U_ref      - reference wind speed for scaling [m/s] (default: 5.0)
//
// Returns:
//   U_slope_equiv - equivalent wind speed from slope [m/s]
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real slope_factor_to_equivalent_wind(Real phi_slope, Real U_ref = Real(5.0))
{
    // Ensure non-negative slope factor
    phi_slope = amrex::max(phi_slope, Real(0.0));
    
    // Convert slope factor to equivalent wind speed
    // Empirical relationship: U_equiv ≈ U_ref * sqrt(φ_s)
    Real U_equiv = U_ref * std::sqrt(phi_slope);
    
    return U_equiv;
}

// ============================================================================
// Compute effective wind speed from wind and slope
//
// Combines meteorological wind and slope-induced draft into a single
// effective wind speed diagnostic.
//
// Parameters:
//   U_wind     - wind speed [m/s]
//   phi_slope  - Rothermel slope factor (dimensionless)
//   U_ref      - reference wind speed for slope conversion [m/s] (default: 5.0)
//
// Returns:
//   U_eff - effective wind speed [m/s]
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real compute_effective_wind_speed(
    Real U_wind,
    Real phi_slope,
    Real U_ref = Real(5.0))
{
    // Ensure non-negative wind
    U_wind = amrex::max(U_wind, Real(0.0));
    
    // Convert slope to equivalent wind
    Real U_slope = slope_factor_to_equivalent_wind(phi_slope, U_ref);
    
    // Vector combination (assuming perpendicular worst case)
    Real U_eff = std::sqrt(U_wind * U_wind + U_slope * U_slope);
    
    return U_eff;
}

// ============================================================================
// Compute effective wind speed for entire MultiFab
//
// Parameters:
//   eff_wind_mf   - output: effective wind speed [m/s]
//   vel_mf        - input:  velocity field (component 0 = u, 1 = v) [m/s]
//   slope_mf      - input:  slope magnitude (tan(θ)) or slope factor φ_s
//   use_phi_slope - if true, slope_mf contains φ_s; if false, compute from tan(θ)
//   beta          - packing ratio (needed only if use_phi_slope = false)
//   U_ref         - reference wind speed for slope scaling [m/s]
// ============================================================================
inline void compute_effective_wind_speed_field(
    MultiFab&       eff_wind_mf,
    const MultiFab& vel_mf,
    const MultiFab& slope_mf,
    bool            use_phi_slope = true,
    Real            beta = Real(0.0),
    Real            U_ref = Real(5.0))
{
    for (MFIter mfi(eff_wind_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto eff_wind = eff_wind_mf.array(mfi);
        auto vel      = vel_mf.const_array(mfi);
        auto slope    = slope_mf.const_array(mfi);
        
        const Real beta_val = beta;
        const bool use_phi = use_phi_slope;
        const Real U_ref_val = U_ref;

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            // Compute wind magnitude
            Real u = vel(i, j, k, 0);
            Real v = vel(i, j, k, 1);
            Real U_wind = std::sqrt(u * u + v * v);
            
            // Get slope factor
            Real phi_s;
            if (use_phi) {
                // slope_mf contains φ_s directly
                phi_s = slope(i, j, k);
            } else {
                // slope_mf contains tan(θ); compute φ_s
                // φ_s = 5.275 * β^(-0.3) * tan²(θ)
                Real tan_theta = slope(i, j, k);
                Real beta_safe = amrex::max(beta_val, Real(0.001));
                phi_s = Real(5.275) * std::pow(beta_safe, Real(-0.3))
                      * tan_theta * tan_theta;
            }
            
            // Compute effective wind speed
            eff_wind(i, j, k) = compute_effective_wind_speed(U_wind, phi_s, U_ref_val);
        });
    }
}

#endif // EFFECTIVE_WIND_SPEED_H
