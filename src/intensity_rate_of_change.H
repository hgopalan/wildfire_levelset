// ============================================================================
// Fire Line Intensity Rate of Change (dI/dt) Model
// ============================================================================
// Computes the temporal derivative of fireline intensity at each cell, providing
// an indicator of fire acceleration/deceleration. This is useful for:
//   - Identifying dangerous "blow-up" conditions (rapid intensification)
//   - Suppression resource allocation
//   - Fire behavior forecasting
//
// Model:
//   dI/dt [kW/m/s] = (I_current - I_previous) / Δt
//
// where:
//   I_current   = current fireline intensity [kW/m]
//   I_previous  = previous fireline intensity [kW/m]
//   Δt          = time step [s]
//
// Interpretation:
//   dI/dt > 0  : Fire is intensifying (accelerating)
//   dI/dt = 0  : Fire is steady-state
//   dI/dt < 0  : Fire is weakening (decelerating)
//
// Large positive values indicate dangerous rapid fire growth conditions.
//
// References:
//   Rothermel, R.C. (1972). A mathematical model for predicting fire spread
//     in wildland fuels. USDA Forest Service Research Paper INT-115.
//   Byram, G.M. (1959). Combustion of forest fuels. In: Davis, K.P. (ed.)
//     Forest Fire: Control and Use. McGraw-Hill, New York.
// ============================================================================

#ifndef INTENSITY_RATE_OF_CHANGE_H
#define INTENSITY_RATE_OF_CHANGE_H

#include "constants.H"
#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_GpuLaunch.H>

using namespace amrex;

// ============================================================================
// Compute rate of change of fireline intensity
//
// Parameters:
//   dI_dt_mf         – output: temporal derivative of intensity [kW/m/s]
//   I_current_mf     – input:  current fireline intensity [kW/m]
//   I_previous_mf    – input:  previous fireline intensity [kW/m]
//   dt               – time step [s]
//
// Notes:
//   - The previous intensity field should be stored from the last time step
//   - At the first time step, I_previous should be initialized to zero
//   - For stability, consider applying a smoothing filter to dI/dt
// ============================================================================
inline void compute_intensity_rate_of_change(MultiFab& dI_dt_mf,
                                             const MultiFab& I_current_mf,
                                             const MultiFab& I_previous_mf,
                                             Real dt)
{
    // Avoid division by zero
    if (dt <= Real(0.0)) {
        dI_dt_mf.setVal(Real(0.0));
        return;
    }

    const Real inv_dt = Real(1.0) / dt;

    for (MFIter mfi(dI_dt_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const I_curr = I_current_mf.const_array(mfi);
        auto const I_prev = I_previous_mf.const_array(mfi);
        auto dI_dt        = dI_dt_mf.array(mfi);
        const Real inv_dt_local = inv_dt;

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            // Temporal derivative: (I_current - I_previous) / dt
            dI_dt(i, j, k) = (I_curr(i, j, k) - I_prev(i, j, k)) * inv_dt_local;
        });
    }
}

// ============================================================================
// Update previous intensity field (copy current to previous)
//
// This should be called after each time step to prepare for the next
// dI/dt computation.
//
// Parameters:
//   I_previous_mf  – output: storage for previous intensity [kW/m]
//   I_current_mf   – input:  current intensity to be stored [kW/m]
// ============================================================================
inline void update_previous_intensity(MultiFab& I_previous_mf,
                                      const MultiFab& I_current_mf)
{
    MultiFab::Copy(I_previous_mf, I_current_mf, 0, 0, 1, 0);
}

// ============================================================================
// Classify fire behavior based on dI/dt
//
// Returns a classification of the fire behavior:
//   -2: Rapid weakening (dI/dt < -threshold_high)
//   -1: Weakening       (-threshold_high <= dI/dt < -threshold_low)
//    0: Steady-state    (-threshold_low <= dI/dt <= threshold_low)
//    1: Intensifying    (threshold_low < dI/dt <= threshold_high)
//    2: Rapid intensification (dI/dt > threshold_high)
//
// Parameters:
//   fire_behavior_mf  – output: classification (-2 to +2)
//   dI_dt_mf          – input:  rate of change of intensity [kW/m/s]
//   threshold_low     – low threshold [kW/m/s] (default: 10.0)
//   threshold_high    – high threshold [kW/m/s] (default: 50.0)
// ============================================================================
inline void classify_fire_behavior(MultiFab& fire_behavior_mf,
                                   const MultiFab& dI_dt_mf,
                                   Real threshold_low  = Real(10.0),
                                   Real threshold_high = Real(50.0))
{
    for (MFIter mfi(fire_behavior_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const dI_dt = dI_dt_mf.const_array(mfi);
        auto fb          = fire_behavior_mf.array(mfi);
        const Real t_low = threshold_low;
        const Real t_high = threshold_high;

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const Real rate = dI_dt(i, j, k);
            
            if (rate > t_high) {
                fb(i, j, k) = Real(2.0);  // Rapid intensification
            } else if (rate > t_low) {
                fb(i, j, k) = Real(1.0);  // Intensifying
            } else if (rate >= -t_low) {
                fb(i, j, k) = Real(0.0);  // Steady-state
            } else if (rate >= -t_high) {
                fb(i, j, k) = Real(-1.0); // Weakening
            } else {
                fb(i, j, k) = Real(-2.0); // Rapid weakening
            }
        });
    }
}

#endif // INTENSITY_RATE_OF_CHANGE_H
