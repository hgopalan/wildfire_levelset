// ============================================================================
// Packing Ratio Diagnostic
// ============================================================================
// Computes the Rothermel (1972) fuel bed packing ratio (β) and the relative
// packing ratio (β / β_opt), which are key parameters in fire spread behavior.
//
// The packing ratio is defined as:
//   β = ρ_b / ρ_p
//
// where:
//   ρ_b - oven-dry bulk density [lb/ft³] = w_0 / δ
//   ρ_p - oven-dry particle density [lb/ft³]
//   w_0 - oven-dry fuel loading [lb/ft²]
//   δ   - fuel bed depth [ft]
//
// The optimum packing ratio is:
//   β_opt = 3.348 × σ^(-0.8189)
//
// where σ is the surface-area-to-volume ratio [ft⁻¹].
//
// The relative packing ratio β / β_opt affects the reaction velocity and
// overall fire behavior. Values close to 1.0 indicate near-optimal fuel
// bed compaction for fire spread.
//
// References:
//   Rothermel, R.C. (1972). "A mathematical model for predicting fire spread
//     in wildland fuels." USDA Forest Service Research Paper INT-115.
//   Andrews, P.L. (2018). "The Rothermel surface fire spread model and
//     associated developments." USDA Forest Service Gen. Tech. Rep. RMRS-GTR-371.
// ============================================================================

#ifndef PACKING_RATIO_H
#define PACKING_RATIO_H

#include "constants.H"
#include "parse_inputs.H"
#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_GpuLaunch.H>

using namespace amrex;

// ============================================================================
// Compute packing ratio and relative packing ratio
//
// Parameters:
//   packing_ratio_mf     – output: packing ratio β [-]
//   relative_packing_mf  – output: relative packing ratio β / β_opt [-]
//   rp                   – Rothermel fuel parameters
// ============================================================================
inline void compute_packing_ratio(MultiFab& packing_ratio_mf,
                                  MultiFab& relative_packing_mf,
                                  const RothermelParams& rp)
{
    // Use multi-class path if per-class loads are specified
    const Real total_per_class = rp.w_d1 + rp.w_d10 + rp.w_d100
                                + rp.w_lh + rp.w_lw;
    const bool use_multiclass = (total_per_class > 0.0);

    // Effective fuel load [lb/ft²]
    const Real w_eff = use_multiclass ? total_per_class : rp.w0;

    // Characteristic SAV [ft⁻¹]
    Real sigma_c;
    if (use_multiclass) {
        // Load-weighted mean SAV (A-factor weighted)
        const Real SIGMA_D10  = WildfireConst::SIGMA_D10;   // 109 ft⁻¹
        const Real SIGMA_D100 = WildfireConst::SIGMA_D100;  // 30 ft⁻¹

        const Real A_d1   = rp.w_d1  * rp.sigma_d1 * std::exp(-Real(138.0) / rp.sigma_d1);
        const Real A_d10  = rp.w_d10 * SIGMA_D10   * std::exp(-Real(138.0) / SIGMA_D10);
        const Real A_d100 = rp.w_d100 * SIGMA_D100 * std::exp(-Real(138.0) / SIGMA_D100);
        const Real A_lh   = rp.w_lh  * rp.sigma_lh * std::exp(-Real(138.0) / rp.sigma_lh);
        const Real A_lw   = rp.w_lw  * rp.sigma_lw * std::exp(-Real(138.0) / rp.sigma_lw);

        const Real A_total = A_d1 + A_d10 + A_d100 + A_lh + A_lw;

        const Real sigma_sum = A_d1 * rp.sigma_d1 / A_total
                             + A_d10 * SIGMA_D10 / A_total
                             + A_d100 * SIGMA_D100 / A_total
                             + A_lh * rp.sigma_lh / A_total
                             + A_lw * rp.sigma_lw / A_total;
        sigma_c = (A_total > 0.0) ? sigma_sum : rp.sigma;
    } else {
        sigma_c = rp.sigma;
    }

    // Bulk density: ρ_b = w_0 / δ [lb/ft³]
    const Real rho_b = w_eff / amrex::max(rp.delta, Real(0.01));

    // Packing ratio: β = ρ_b / ρ_p [-]
    const Real beta = rho_b / amrex::max(rp.rho_p, Real(1.0));

    // Optimum packing ratio: β_opt = 3.348 × σ^(-0.8189)
    const Real beta_opt = Real(3.348) * std::pow(sigma_c, Real(-0.8189));

    // Relative packing ratio: β / β_opt [-]
    const Real rel_beta = beta / amrex::max(beta_opt, Real(1.0e-6));

    // Fill MultiFabs with constant values (uniform fuel parameters)
    for (MFIter mfi(packing_ratio_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto beta_arr = packing_ratio_mf.array(mfi);
        auto rel_arr  = relative_packing_mf.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            beta_arr(i, j, k) = beta;
            rel_arr(i, j, k) = rel_beta;
        });
    }
}

#endif // PACKING_RATIO_H
