// ============================================================================
// Flame Front Depth Diagnostic
// ============================================================================
// Computes the flame front depth (thickness) as the product of the rate of
// spread and the flame residence time.
//
// Flame front depth:
//   D_f = R × τ_res
//
// where:
//   R       - rate of spread [m/s]
//   τ_res   - flame residence time [s]
//
// The flame front depth represents the spatial extent of the actively burning
// zone and is used in heat transfer calculations, backing fire behavior, and
// fire line construction planning.
//
// For uniform fuels, τ_res is typically in the range 20–200 seconds for
// surface fires. The residence time depends on fuel particle size, fuel load,
// and fire intensity (Rothermel 1983):
//
//   τ_res ≈ (α × ρ_p) / σ
//
// where:
//   α    - empirical constant [s·ft/(lb/ft³)] ≈ 3600
//   ρ_p  - particle density [lb/ft³]
//   σ    - surface-area-to-volume ratio [ft⁻¹]
//
// References:
//   Rothermel, R.C. (1983). "How to predict the spread and intensity of
//     forest and range fires." USDA Forest Service Gen. Tech. Rep. INT-143.
//   Andrews, P.L. & Rothermel, R.C. (1982). "Charts for interpreting wildland
//     fire behavior characteristics." USDA Forest Service Gen. Tech. Rep. INT-131.
//   BehavePlus fire modeling system documentation.
// ============================================================================

#ifndef FLAME_FRONT_DEPTH_H
#define FLAME_FRONT_DEPTH_H

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_GpuLaunch.H>

using namespace amrex;

// ============================================================================
// Compute flame front depth from rate of spread and residence time
//
// Parameters:
//   flame_depth_mf       – output: flame front depth [m]
//   R_mf                 – input:  rate of spread [m/s]
//   residence_time_mf    – input:  flame residence time [s]
// ============================================================================
inline void compute_flame_front_depth(MultiFab& flame_depth_mf,
                                      const MultiFab& R_mf,
                                      const MultiFab& residence_time_mf)
{
    for (MFIter mfi(flame_depth_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const R   = R_mf.const_array(mfi);
        auto const tau = residence_time_mf.const_array(mfi);
        auto depth     = flame_depth_mf.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            // Flame front depth [m] = R [m/s] × τ_res [s]
            depth(i, j, k) = amrex::max(R(i, j, k), Real(0.0)) * tau(i, j, k);
        });
    }
}

#endif // FLAME_FRONT_DEPTH_H
