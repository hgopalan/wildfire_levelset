// ============================================================================
// Spot Fire Ignition Delay Model (FARSITE-compatible)
// ============================================================================
// Models the delay between firebrand landing and actual ignition based on
// fuel moisture content. FARSITE includes this delay to better represent
// the reality that firebrands don't ignite instantly upon landing - there's
// a lag time that depends on fuel conditions.
//
// Physical basis:
//   - Wetter fuels take longer to ignite after a firebrand lands
//   - The delay is inversely proportional to the fuel moisture deficit
//   - At extinction moisture (M_x), delay approaches infinity (no ignition)
//   - At very dry conditions (M_f → 0), delay approaches the base value
//
// Model:
//   τ_delay [s] = τ_base × (M_x / (M_x - M_f))
//
// where:
//   τ_base  = base delay time at zero moisture [s] (typically 60-300 s)
//   M_x     = extinction moisture content [fraction]
//   M_f     = current fuel moisture content [fraction]
//
// References:
//   Albini, F.A. (1979). Spot fire distance from burning trees.
//     USDA Forest Service INT-56.
//   Finney, M.A. (1998). FARSITE: Fire Area Simulator.
//     USDA Forest Service RMRS-RP-4.
//   Anderson, H.E. (1970). Forest fuel ignitibility.
//     Fire Technology 6(4):312-319.
// ============================================================================

#ifndef SPOT_IGNITION_DELAY_H
#define SPOT_IGNITION_DELAY_H

#include "constants.H"
#include "rothermel_model.H"
#include <AMReX.H>
#include <cmath>

using namespace amrex;

// Tolerance for fuel moisture comparisons to avoid numerical issues
constexpr Real MOISTURE_TOLERANCE = Real(1e-6);

// ============================================================================
// Compute ignition delay time [s] based on fuel moisture
//
// Returns the delay time between firebrand landing and ignition.
// When M_f approaches M_x (extinction moisture), delay → ∞.
// When M_f → 0 (dry fuel), delay → τ_base.
//
// Parameters:
//   τ_base – base delay time at zero moisture [s]
//   M_x    – extinction moisture content [fraction]
//   M_f    – current fuel moisture content [fraction]
//
// Returns: delay time [s], or -1.0 if ignition is impossible (M_f >= M_x)
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real compute_spot_ignition_delay(Real tau_base, Real M_x, Real M_f)
{
    // Safety check: if fuel moisture at or above extinction, no ignition possible
    if (M_f >= M_x - MOISTURE_TOLERANCE) {
        return Real(-1.0); // Negative indicates no ignition will occur
    }

    // Prevent division by very small numbers
    Real moisture_deficit = M_x - M_f;
    if (moisture_deficit < MOISTURE_TOLERANCE) {
        return Real(1e8); // Effectively infinite delay
    }

    // FARSITE-style delay formula
    Real delay = tau_base * (M_x / moisture_deficit);

    // Cap at a reasonable maximum (1 hour)
    return amrex::min(delay, Real(3600.0));
}

// ============================================================================
// Compute ignition delay using Rothermel parameters
//
// Convenience wrapper that extracts M_x and M_f from RothermelParams.
// Uses 1-hr dead fuel moisture (M_d) as the limiting moisture class.
//
// Parameters:
//   τ_base – base delay time at zero moisture [s]
//   rp     – Rothermel parameters (contains M_x and M_f)
//
// Returns: delay time [s], or -1.0 if ignition is impossible
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real compute_spot_ignition_delay(Real tau_base, const RothermelParams& rp)
{
    return compute_spot_ignition_delay(tau_base, rp.M_x, rp.M_f);
}

// ============================================================================
// Compute spatially-varying ignition delay from spatial moisture field
//
// When spatial moisture is available (from .fms file or diurnal moisture),
// compute per-cell ignition delays based on local 1-hr dead fuel moisture.
//
// Parameters:
//   spot_ignition_delay_mf – output: delay time [s] per cell
//   τ_base                 – base delay time at zero moisture [s]
//   M_x                    – extinction moisture content [fraction]
//   spatial_moisture_mf    – input: per-cell moisture MultiFab (FMS_MOISTURE_NCOMP)
//                            Component 0 = 1-hr dead moisture
//
// Output: spot_ignition_delay_mf filled with delay times
//         -1.0 for cells where ignition is impossible (M_f >= M_x)
// ============================================================================
inline void compute_spatial_spot_ignition_delay(
    MultiFab&       spot_ignition_delay_mf,
    Real            tau_base,
    Real            M_x,
    const MultiFab& spatial_moisture_mf)
{
    for (MFIter mfi(spot_ignition_delay_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto delay_arr = spot_ignition_delay_mf.array(mfi);
        auto moisture_arr = spatial_moisture_mf.const_array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            Real M_f = moisture_arr(i, j, k, 0); // Component 0 = 1-hr dead moisture
            delay_arr(i, j, k) = compute_spot_ignition_delay(tau_base, M_x, M_f);
        });
    }
}

#endif // SPOT_IGNITION_DELAY_H
