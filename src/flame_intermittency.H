// ============================================================================
// Flame Intermittency Factor Model
// ============================================================================
// Accounts for pulsating/intermittent flames rather than steady burning.
// Real flames flicker and pulse, which affects heat transfer efficiency and
// is important for spotting (firebrands released in pulses).
//
// Physics:
//   - Real flames are not steady; they flicker and pulsate
//   - Intermittency reduces effective heat transfer compared to steady flames
//   - Higher intensity fires have more continuous flames (less intermittent)
//   - Low intensity fires have highly intermittent, flickering flames
//
// Model (from flame stability research):
//   γ = 1 - exp(-k × I_byram)
//
// where:
//   γ         = intermittency factor [0, 1]
//              0 = fully intermittent (no continuous flame)
//              1 = fully continuous flame
//   I_byram   = Byram fireline intensity [kW/m]
//   k         = empirical coefficient ≈ 0.00325 [m/kW]
//
// Application:
//   Q_rad_eff = γ × Q_rad_steady
//
// This modifies the effective radiative heat flux from flames, reducing it
// for low-intensity fires with flickering flames.
//
// References:
//   Finney, M.A., et al. (2015). Role of buoyant flame dynamics in wildfire
//     spread. PNAS 112(32):9833-9838.
//   Frankman, D., et al. (2013). Measurements of convective and radiative
//     heating in wildland fires. International Journal of Wildland Fire
//     22(2):157-167.
// ============================================================================

#ifndef FLAME_INTERMITTENCY_H
#define FLAME_INTERMITTENCY_H

#include "constants.H"
#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_GpuLaunch.H>
#include <cmath>

using namespace amrex;

// Default empirical coefficient for intermittency model [m/kW]
constexpr Real INTERMITTENCY_COEFF_DEFAULT = Real(0.00325);

// ============================================================================
// Compute flame intermittency factor from fireline intensity
//
// γ = 1 - exp(-k × I_byram)
//
// Parameters:
//   I_byram  – Byram fireline intensity [kW/m]
//   k        – empirical coefficient [m/kW] (default: 0.00325)
//
// Returns: intermittency factor γ ∈ [0, 1]
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real compute_flame_intermittency(Real I_byram,
                                 Real k = INTERMITTENCY_COEFF_DEFAULT)
{
    // Intermittency factor: γ = 1 - exp(-k × I_byram)
    // For low intensity: γ → 0 (highly intermittent)
    // For high intensity: γ → 1 (continuous flame)
    
    if (I_byram <= Real(0.0)) {
        return Real(0.0);
    }
    
    return Real(1.0) - std::exp(-k * I_byram);
}

// ============================================================================
// Compute flame intermittency field
//
// Parameters:
//   intermittency_mf   – output: intermittency factor γ [0, 1]
//   fireline_intensity_mf – input: Byram intensity [kW/m]
//   k_coeff            – empirical coefficient [m/kW]
// ============================================================================
inline void compute_flame_intermittency_field(MultiFab& intermittency_mf,
                                              const MultiFab& fireline_intensity_mf,
                                              Real k_coeff = INTERMITTENCY_COEFF_DEFAULT)
{
    for (MFIter mfi(intermittency_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto gamma     = intermittency_mf.array(mfi);
        auto const I_B = fireline_intensity_mf.const_array(mfi);
        const Real k   = k_coeff;

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            gamma(i, j, k) = compute_flame_intermittency(I_B(i, j, k), k);
        });
    }
}

// ============================================================================
// Apply intermittency factor to radiative heat flux
//
// Q_rad_eff = γ × Q_rad_steady
//
// Parameters:
//   Q_rad_eff_mf       – output: effective radiative flux [W/m²]
//   Q_rad_steady_mf    – input: steady-state radiative flux [W/m²]
//   intermittency_mf   – input: intermittency factor γ [0, 1]
// ============================================================================
inline void apply_flame_intermittency(MultiFab& Q_rad_eff_mf,
                                      const MultiFab& Q_rad_steady_mf,
                                      const MultiFab& intermittency_mf)
{
    for (MFIter mfi(Q_rad_eff_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto Q_eff     = Q_rad_eff_mf.array(mfi);
        auto const Q_steady = Q_rad_steady_mf.const_array(mfi);
        auto const gamma    = intermittency_mf.const_array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Q_eff(i, j, k) = gamma(i, j, k) * Q_steady(i, j, k);
        });
    }
}

// ============================================================================
// Combined: compute intermittency and apply to radiative flux
//
// This is a convenience function that combines the two steps above.
//
// Parameters:
//   Q_rad_eff_mf          – output: effective radiative flux [W/m²]
//   Q_rad_steady_mf       – input: steady-state radiative flux [W/m²]
//   fireline_intensity_mf – input: Byram intensity [kW/m]
//   k_coeff               – empirical coefficient [m/kW]
// ============================================================================
inline void apply_intermittency_to_radiation(MultiFab& Q_rad_eff_mf,
                                             const MultiFab& Q_rad_steady_mf,
                                             const MultiFab& fireline_intensity_mf,
                                             Real k_coeff = INTERMITTENCY_COEFF_DEFAULT)
{
    for (MFIter mfi(Q_rad_eff_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto Q_eff     = Q_rad_eff_mf.array(mfi);
        auto const Q_steady = Q_rad_steady_mf.const_array(mfi);
        auto const I_B      = fireline_intensity_mf.const_array(mfi);
        const Real k        = k_coeff;

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const Real gamma = compute_flame_intermittency(I_B(i, j, k), k);
            Q_eff(i, j, k) = gamma * Q_steady(i, j, k);
        });
    }
}

// ============================================================================
// Classify flame stability based on intermittency
//
// Returns a classification:
//   0: Highly intermittent (γ < 0.3) - flickering flames
//   1: Moderate intermittency (0.3 <= γ < 0.7) - pulsating flames
//   2: Continuous (γ >= 0.7) - stable flames
//
// This can be useful for fire behavior analysis and spotting predictions.
//
// Parameters:
//   stability_mf      – output: flame stability class (0, 1, or 2)
//   intermittency_mf  – input: intermittency factor γ [0, 1]
// ============================================================================
inline void classify_flame_stability(MultiFab& stability_mf,
                                     const MultiFab& intermittency_mf)
{
    for (MFIter mfi(stability_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto stab      = stability_mf.array(mfi);
        auto const gamma = intermittency_mf.const_array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const Real g = gamma(i, j, k);
            
            if (g < Real(0.3)) {
                stab(i, j, k) = Real(0.0);  // Highly intermittent
            } else if (g < Real(0.7)) {
                stab(i, j, k) = Real(1.0);  // Moderate
            } else {
                stab(i, j, k) = Real(2.0);  // Continuous
            }
        });
    }
}

#endif // FLAME_INTERMITTENCY_H
