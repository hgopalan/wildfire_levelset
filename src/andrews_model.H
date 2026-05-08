// ============================================================================
// Andrews (2018) Surface Fire Spread - Wind Adjustments
// ============================================================================
// Implements two key refinements from:
//
//   Andrews, P.L. (2018). The Rothermel Surface Fire Spread Model and
//   Associated Developments: A Comprehensive Explanation.
//   Gen. Tech. Rep. RMRS-GTR-371. USDA Forest Service.
//   https://doi.org/10.2737/RMRS-GTR-371
//
// 1. Wind Adjustment Factor (WAF): converts 20-ft (6.1 m) open wind speed
//    to midflame height wind speed for use in the Rothermel model.
//    Enabled via `rothermel.use_waf = 1`.
//    Two sub-modes are available:
//      a. Unsheltered (open fuel beds, default):
//           WAF = 1.83 / ln((20 + 0.36*h) / (0.13*h))   where h = fuel depth [ft]
//      b. Sheltered (closed-canopy fuel beds, FARSITE-style):
//           WAF_sheltered = 0.555 / (sqrt(f_c × h_c) × ln((20 + 0.36*h_c)/(0.13*h_c)))
//           where f_c = canopy cover fraction [0-1], h_c = canopy height [ft]
//           Enabled automatically when canopy_cover > waf_shelter_threshold (default 0.5)
//           and canopy height > 0.
//
// 2. Maximum Effective Wind Speed (MEWS): the wind factor φ_w is bounded
//    above by 0.9 × I_R (Rothermel 1972); the corresponding maximum
//    midflame wind speed U_max caps the wind input to the ROS formula.
//    Enabled via `rothermel.use_wind_limit = 1`.
//
// These two options can be enabled independently of each other and are
// compatible with all fire spread paths (level-set and FARSITE) and all
// fuel models (FBFM13, FBFM40, and user-specified values).
// ============================================================================

#ifndef ANDREWS_MODEL_H
#define ANDREWS_MODEL_H

#include <AMReX.H>
#include <cmath>

// ============================================================================
// Wind Adjustment Factor (WAF) for unsheltered (open) fuel beds
// ============================================================================
// Andrews (2018) Eq. 2.1 (after Albini & Baughman 1979):
//
//   WAF = 1.83 / ln((20 + 0.36 * h) / (0.13 * h))
//
// where h = fuel bed depth [ft].
//
// The WAF converts 20-ft (6.1 m) open wind speed to midflame height wind
// speed. Typical values range from ~0.3 (deep fuel beds, h ≈ 6 ft) to
// ~0.5 (shallow beds, h ≈ 0.5 ft). For the default FM4 fuel (h = 2 ft)
// the WAF is approximately 0.42.
//
// Parameters:
//   delta_ft: fuel bed depth [ft]
//
// Returns:
//   WAF (dimensionless), clamped to [0.1, 1.0]
// ============================================================================
static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_waf_unsheltered(amrex::Real delta_ft)
{
    using namespace amrex;

    if (delta_ft < Real(0.01)) return Real(0.4); // floor for degenerate fuel

    const Real log_arg = (Real(20.0) + Real(0.36) * delta_ft)
                       / (Real(0.13) * delta_ft);
    if (log_arg <= Real(1.0)) return Real(0.4); // safety guard

    const Real waf = Real(1.83) / std::log(log_arg);
    return amrex::max(Real(0.1), amrex::min(waf, Real(1.0)));
}

// ============================================================================
// Maximum Effective Wind Speed (MEWS) for the Rothermel model
// ============================================================================
// Andrews (2018) / Rothermel (1972): the wind factor φ_w is bounded above by
//
//   φ_w_max = 0.9 * I_R   [BTU/ft²/min]
//
// The corresponding maximum midflame wind speed [ft/min] is found by
// inverting the wind factor equation (Rothermel 1972, Eq. 47):
//
//   φ_w = C * U^B * beta_ratio^{-E}
//   => U_max = (φ_w_max / (C * beta_ratio_E))^{1/B}
//
// Applying this cap prevents unrealistically large ROS at very high wind
// speeds where Rothermel's empirical model is no longer calibrated.
//
// Parameters:
//   I_R:          reaction intensity [BTU/ft²/min]
//   C, B:         Rothermel wind factor coefficients (Eqs. 47–48)
//   beta_ratio_E: (beta/beta_op)^{-E}, precomputed packing-ratio factor
//
// Returns:
//   Maximum effective wind speed [ft/min]; 0.0 if not computable.
// ============================================================================
static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_max_ews_ftmin(amrex::Real I_R,
                                   amrex::Real C,
                                   amrex::Real B,
                                   amrex::Real beta_ratio_E)
{
    using namespace amrex;

    if (C <= Real(0.0) || B <= Real(0.0) || beta_ratio_E <= Real(0.0)
        || I_R <= Real(0.0)) {
        return Real(0.0);
    }

    const Real phi_w_max = Real(0.9) * I_R;
    const Real U_max = std::pow(phi_w_max / (C * beta_ratio_E), Real(1.0) / B);
    return amrex::max(U_max, Real(0.0));
}

// ============================================================================
// Wind Adjustment Factor (WAF) for sheltered (closed-canopy) fuel beds
// ============================================================================
// FARSITE / Albini & Baughman (1979) sheltered-canopy WAF formula:
//
//   WAF_sheltered = 0.555 / (sqrt(f_c × h_c) × ln((20 + 0.36*h_c) / (0.13*h_c)))
//
// where:
//   f_c = canopy cover fraction [0, 1]  (0 = open, 1 = fully closed)
//   h_c = canopy height [ft]
//
// This formula differs from the unsheltered version in that the numerator is
// 0.555 (not 1.83) and the denominator includes the square root of the
// fractional canopy cover, which dampens the midflame wind more strongly under
// a closed canopy.  The formula degenerates gracefully:
//   – When f_c → 0 the denominator → 0 → WAF → unsheltered (guarded below).
//   – When h_c → delta (fuel depth) the two formulas converge.
//
// Returns the unsheltered WAF (based on delta_ft) when:
//   • f_c ≤ 0  (no canopy),  OR
//   • h_c ≤ delta_ft (canopy no taller than fuel bed),  OR
//   • the log argument ≤ 1 (degenerate geometry).
//
// Parameters:
//   f_c:       canopy cover fraction [0, 1]
//   h_c_ft:    canopy height [ft]
//   delta_ft:  fuel bed depth [ft] (used as fallback)
//
// Returns:
//   WAF (dimensionless), clamped to [0.05, 1.0]
// ============================================================================
static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_waf_sheltered(amrex::Real f_c,
                                   amrex::Real h_c_ft,
                                   amrex::Real delta_ft)
{
    using namespace amrex;

    // Fall back to unsheltered when canopy is absent or shorter than fuel bed
    if (f_c <= Real(0.0) || h_c_ft <= delta_ft) {
        return compute_waf_unsheltered(delta_ft);
    }

    const Real log_arg = (Real(20.0) + Real(0.36) * h_c_ft)
                       / (Real(0.13) * h_c_ft);
    if (log_arg <= Real(1.0)) return compute_waf_unsheltered(delta_ft);

    const Real denom = std::sqrt(f_c * h_c_ft) * std::log(log_arg);
    if (denom <= Real(1.0e-6)) return compute_waf_unsheltered(delta_ft);

    const Real waf = Real(0.555) / denom;
    // Clamp: never below 0.05 (extremely dense canopy) nor above 1.0
    return amrex::max(Real(0.05), amrex::min(waf, Real(1.0)));
}

// ============================================================================
// Select the appropriate WAF (sheltered vs. unsheltered) per cell
// ============================================================================
// Chooses between compute_waf_sheltered and compute_waf_unsheltered based on
// whether the local canopy cover fraction exceeds a threshold.
//
// Parameters:
//   f_c:               canopy cover fraction [0, 1]
//   h_c_ft:            canopy height [ft]
//   delta_ft:          fuel bed depth [ft]
//   shelter_threshold: minimum canopy cover to invoke sheltered WAF (default: 0.5)
//
// Returns:
//   WAF (dimensionless)
// ============================================================================
static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_waf_per_cell(amrex::Real f_c,
                                  amrex::Real h_c_ft,
                                  amrex::Real delta_ft,
                                  amrex::Real shelter_threshold = amrex::Real(0.5))
{
    if (f_c >= shelter_threshold && h_c_ft > delta_ft) {
        return compute_waf_sheltered(f_c, h_c_ft, delta_ft);
    }
    return compute_waf_unsheltered(delta_ft);
}

#endif // ANDREWS_MODEL_H
