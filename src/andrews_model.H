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
//    The formula is selected via `rothermel.waf_formula`:
//
//      "andrews" (default) — logarithmic Albini & Baughman (1979) formulas:
//        a. Unsheltered (open/shrub fuel beds):
//             WAF = 1.83 / ln((20 + 0.36×h) / (0.13×h))   [h in ft]
//        b. Sheltered (closed-canopy, FARSITE-style):
//             WAF = 0.555 / (sqrt(f_c × h_c) × ln((20 + 0.36×h_c)/(0.13×h_c)))
//             Enabled automatically when canopy_cover > waf_shelter_threshold
//             (default 0.5) and canopy height > 0.
//
//      "behaviorplus" — BehavePlus linear / exponential formulas:
//        a. Unsheltered (open/shrub fuel beds):
//             WAF = 0.36 + 0.004 × h_in   [h_in = h in inches = h_ft × 12]
//        b. Sheltered (forest/closed-canopy fuel beds):
//             WAF = WAF_open(h_c) × exp(−α_c × f_c)
//             where WAF_open is the BehavePlus linear WAF at canopy height h_c,
//             α_c is the canopy attenuation coefficient (default 1.5), and
//             f_c is the canopy cover fraction.
//             Enabled automatically when canopy_cover > waf_shelter_threshold
//             and canopy height > fuel depth.
//             α_c is set via `rothermel.waf_canopy_alpha` (default: 1.5).
//
// 2. Maximum Effective Wind Speed (MEWS): the wind factor φ_w is bounded
//    above by 0.9 × I_R (Rothermel 1972); the corresponding maximum
//    midflame wind speed U_max caps the wind input to the ROS formula.
//    Enabled via `rothermel.use_wind_limit = 1`.
//
// These options can be enabled independently of each other and are
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
//
// See also: compute_waf_behaviorplus for the BehavePlus linear alternative.
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
// Wind Adjustment Factor (WAF) for open/shrub fuels — BehavePlus linear formula
// ============================================================================
// BehavePlus / Andrews (2012) simplified linear WAF for open and shrub fuel
// models:
//
//   WAF = 0.36 + 0.004 × h_in
//
// where h_in = fuel bed depth in inches (= delta_ft × 12).
//
// This is a practical approximation used in BehavePlus that gives a nearly
// flat but slightly increasing WAF with fuel depth.  It is most appropriate
// for open grassland, shrubland, and short-shrub fuel models (FBFM13 1–3,
// FBFM40 GR/GS/SH groups).  Typical values:
//   h = 0.5 ft  ( 6 in): WAF ≈ 0.384
//   h = 2.0 ft  (24 in): WAF ≈ 0.456
//   h = 6.0 ft  (72 in): WAF ≈ 0.648
//   h = 8.0 ft  (96 in): WAF ≈ 0.744
//
// Enabled when rothermel.waf_formula = "behaviorplus".
//
// Parameters:
//   delta_ft: fuel bed depth [ft]
//
// Returns:
//   WAF (dimensionless), clamped to [0.1, 1.0]
// ============================================================================
static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_waf_behaviorplus(amrex::Real delta_ft)
{
    using namespace amrex;
    if (delta_ft < Real(0.01)) return Real(0.36);
    const Real h_in = delta_ft * Real(12.0);  // feet → inches
    const Real waf  = Real(0.36) + Real(0.004) * h_in;
    return amrex::max(Real(0.1), amrex::min(waf, Real(1.0)));
}

// ============================================================================
// Exponential canopy wind-sheltering WAF for forest/closed-canopy fuels
// ============================================================================
// For forest and closed-canopy fuel models the wind at midflame height is
// strongly reduced by the canopy.  This exponential decay model applies a
// Beer–Lambert-style attenuation to the open-fuel WAF:
//
//   WAF_canopy = WAF_open × exp(−α_c × f_c)
//
// where:
//   WAF_open – open-fuel WAF (from compute_waf_behaviorplus or
//              compute_waf_unsheltered) evaluated at the canopy height h_c
//   f_c      – canopy cover fraction [0, 1]
//   α_c      – canopy attenuation coefficient (default 1.5; larger values give
//               stronger sheltering); set via rothermel.waf_canopy_alpha
//
// Physical basis: the exponential profile u(z)/u_h = exp(−α(1 − z/h_c))
// (Massman 1987; Shaw & Pereira 1982) integrated over the lower canopy gives
// a mean wind reduction of approximately exp(−α_c × f_c) relative to the
// above-canopy wind.  This model requires only canopy cover, which is
// available from standard LCP landscape files.
//
// The formula degenerates gracefully:
//   f_c → 0: WAF_canopy → WAF_open (no sheltering)
//   f_c = 1, α_c = 1.5: WAF_canopy ≈ 0.22 × WAF_open (dense canopy)
//
// Enabled when rothermel.waf_formula = "behaviorplus" and canopy cover > 0.
//
// Parameters:
//   waf_open  – open-fuel WAF at reference height (dimensionless)
//   f_c       – canopy cover fraction [0, 1]
//   alpha_c   – canopy attenuation coefficient (default: 1.5)
//
// Returns:
//   WAF (dimensionless), clamped to [0.05, 1.0]
// ============================================================================
static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_waf_canopy_exp(amrex::Real waf_open,
                                    amrex::Real f_c,
                                    amrex::Real alpha_c)
{
    using namespace amrex;
    if (f_c <= Real(0.0) || alpha_c <= Real(0.0)) return waf_open;
    const Real waf = waf_open * std::exp(-alpha_c * f_c);
    return amrex::max(Real(0.05), amrex::min(waf, Real(1.0)));
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
// Chooses the WAF formula based on waf_formula_id and whether the local
// canopy cover fraction exceeds a threshold:
//
//   waf_formula_id = 0  (Andrews, default):
//     Unsheltered: compute_waf_unsheltered    — Andrews (2018) logarithmic
//     Sheltered:   compute_waf_sheltered      — Albini & Baughman (1979)
//
//   waf_formula_id = 1  (BehavePlus):
//     Unsheltered: compute_waf_behaviorplus   — linear  0.36 + 0.004×h_in
//     Sheltered:   compute_waf_canopy_exp     — exponential Beer–Lambert decay
//
// A cell is "sheltered" when f_c >= shelter_threshold AND h_c_ft > delta_ft.
//
// Parameters:
//   f_c:               canopy cover fraction [0, 1]
//   h_c_ft:            canopy height [ft]
//   delta_ft:          fuel bed depth [ft]
//   shelter_threshold: minimum canopy cover to invoke sheltered WAF (default: 0.5)
//   waf_formula_id:    0 = Andrews (default), 1 = BehavePlus
//   waf_canopy_alpha:  canopy attenuation coefficient for BehavePlus path (default: 1.5)
//
// Returns:
//   WAF (dimensionless)
// ============================================================================
static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_waf_per_cell(amrex::Real f_c,
                                  amrex::Real h_c_ft,
                                  amrex::Real delta_ft,
                                  amrex::Real shelter_threshold = amrex::Real(0.5),
                                  int         waf_formula_id    = 0,
                                  amrex::Real waf_canopy_alpha  = amrex::Real(1.5))
{
    const bool sheltered = (f_c >= shelter_threshold && h_c_ft > delta_ft);

    if (waf_formula_id == 1) {
        // BehavePlus path
        if (sheltered) {
            // Exponential canopy attenuation; use canopy height as reference
            // so that the midflame WAF reflects wind reduction through the canopy
            // rather than through a shallow surface fuel bed.
            const amrex::Real waf_open = compute_waf_behaviorplus(h_c_ft);
            return compute_waf_canopy_exp(waf_open, f_c, waf_canopy_alpha);
        }
        return compute_waf_behaviorplus(delta_ft);
    }

    // Andrews / Albini & Baughman path (default)
    if (sheltered) {
        return compute_waf_sheltered(f_c, h_c_ft, delta_ft);
    }
    return compute_waf_unsheltered(delta_ft);
}

#endif // ANDREWS_MODEL_H
