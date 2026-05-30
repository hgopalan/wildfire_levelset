// ============================================================================
// Simard Exponential Time-Lag Moisture Model
// ============================================================================
// Simple exponential approximation for 10-hr and 100-hr fuel moisture
// equilibration that's easier to parameterize than full differential equations.
//
// Formula: M(t+dt) = M_eq + (M(t) - M_eq) * exp(-dt/tau)
//
// where:
//   M(t)   - current moisture content [fraction]
//   M_eq   - equilibrium moisture content [fraction]
//   tau    - time-lag constant [hours] (10 hr for 10-hr fuels, 100 hr for 100-hr fuels)
//   dt     - time step [hours]
//
// References:
//   Simard, A.J. (1968). "The moisture content of forest fuels." Forestry
//     Branch Info. Rep. FF-X-14, Dept. of Forestry, Canada.
//   Fosberg, M.A. (1975). "Heat and water vapor flux in conifer forest litter
//     and duff: a theoretical model." USDA Forest Service Research Paper RM-152.
//
// Implementation:
//   This provides a simpler alternative to the full time-lag differential
//   equations in fuel_moisture_timelag_de.H for applications where detailed
//   physics-based moisture dynamics are not required.
// ============================================================================

#ifndef SIMARD_MOISTURE_H
#define SIMARD_MOISTURE_H

#include "constants.H"
#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_GpuLaunch.H>
#include <cmath>

using namespace amrex;

// ============================================================================
// Simard exponential time-lag model for single fuel class
//
// Updates moisture content toward equilibrium using exponential relaxation.
//
// Parameters:
//   M_current    - current moisture content [fraction]
//   M_equilibrium - equilibrium moisture content [fraction]
//   dt_hours     - time step [hours]
//   tau_hours    - time-lag constant [hours]
//
// Returns:
//   M_new - updated moisture content [fraction]
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real simard_moisture_update(
    amrex::Real M_current,
    amrex::Real M_equilibrium,
    amrex::Real dt_hours,
    amrex::Real tau_hours) noexcept
{
    using namespace amrex;
    
    // Ensure positive time-lag
    tau_hours = amrex::max(tau_hours, Real(0.1));
    
    // Exponential approach to equilibrium
    // M(t+dt) = M_eq + (M(t) - M_eq) * exp(-dt/tau)
    Real decay_factor = std::exp(-dt_hours / tau_hours);
    Real M_new = M_equilibrium + (M_current - M_equilibrium) * decay_factor;
    
    // Clamp to physically valid range [0, 2] (200% moisture is reasonable max)
    M_new = amrex::max(Real(0.0), amrex::min(Real(2.0), M_new));
    
    return M_new;
}

// ============================================================================
// Simard model with standard time-lag constants
//
// Provides simplified interface with standard time-lag values:
//   - 1-hr  fuel: tau = 1 hour
//   - 10-hr fuel: tau = 10 hours
//   - 100-hr fuel: tau = 100 hours
//
// Parameters:
//   M_current     - current moisture content [fraction]
//   M_equilibrium - equilibrium moisture content from weather [fraction]
//   dt_seconds    - time step [seconds]
//   fuel_class    - 1 (1-hr), 10 (10-hr), or 100 (100-hr)
//
// Returns:
//   M_new - updated moisture content [fraction]
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real simard_moisture_standard(
    amrex::Real M_current,
    amrex::Real M_equilibrium,
    amrex::Real dt_seconds,
    int fuel_class) noexcept
{
    using namespace amrex;
    
    // Convert dt to hours
    Real dt_hours = dt_seconds / Real(3600.0);
    
    // Standard time-lag constants [hours]
    Real tau_hours;
    switch (fuel_class) {
        case 1:   tau_hours = Real(1.0);    break;  // 1-hr fuel
        case 10:  tau_hours = Real(10.0);   break;  // 10-hr fuel
        case 100: tau_hours = Real(100.0);  break;  // 100-hr fuel
        default:  tau_hours = Real(10.0);   break;  // Default to 10-hr
    }
    
    return simard_moisture_update(M_current, M_equilibrium, dt_hours, tau_hours);
}

// ============================================================================
// Compute equilibrium moisture content from relative humidity
//
// Simplified equilibrium moisture model based on relative humidity.
// More sophisticated models would include temperature effects.
//
// Formula (empirical from Simard 1968):
//   M_eq ≈ 0.03 + 0.25 * RH  for RH < 0.5
//   M_eq ≈ 0.14 + 0.44 * RH  for RH >= 0.5
//
// Parameters:
//   RH - relative humidity [0-1]
//
// Returns:
//   M_eq - equilibrium moisture content [fraction]
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_equilibrium_moisture_from_rh(amrex::Real RH) noexcept
{
    using namespace amrex;
    
    // Clamp RH to [0, 1]
    RH = amrex::max(Real(0.0), amrex::min(Real(1.0), RH));
    
    Real M_eq;
    if (RH < Real(0.5)) {
        // Low humidity regime
        M_eq = Real(0.03) + Real(0.25) * RH;
    } else {
        // High humidity regime
        M_eq = Real(0.14) + Real(0.44) * RH;
    }
    
    return M_eq;
}

// ============================================================================
// Update moisture fields using Simard model
//
// Updates all moisture components (1-hr, 10-hr, 100-hr) toward equilibrium
// using exponential relaxation.
//
// Parameters:
//   moisture_mf    - moisture MultiFab (3+ components: M_d1, M_d10, M_d100)
//   M_eq           - equilibrium moisture [fraction]
//   dt_seconds     - time step [seconds]
// ============================================================================
inline void update_moisture_simard(
    MultiFab&   moisture_mf,
    amrex::Real M_eq,
    amrex::Real dt_seconds)
{
    const Real dt_hours = dt_seconds / Real(3600.0);
    const int ncomp = moisture_mf.nComp();
    
    for (MFIter mfi(moisture_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto moisture = moisture_mf.array(mfi);
        const Real M_equilibrium = M_eq;
        const Real dt_h = dt_hours;
        const int nc = ncomp;

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            // Update 1-hr moisture (component 0)
            moisture(i, j, k, 0) = simard_moisture_update(
                moisture(i, j, k, 0), M_equilibrium, dt_h, Real(1.0));
            
            // Update 10-hr moisture (component 1)
            moisture(i, j, k, 1) = simard_moisture_update(
                moisture(i, j, k, 1), M_equilibrium, dt_h, Real(10.0));
            
            // Update 100-hr moisture (component 2)
            if (nc > 2) {
                moisture(i, j, k, 2) = simard_moisture_update(
                    moisture(i, j, k, 2), M_equilibrium, dt_h, Real(100.0));
            }
        });
    }
}

#endif // SIMARD_MOISTURE_H
