// ============================================================================
// Critical Heat Flux for Ignition Model
// ============================================================================
// Replaces binary ignition (φ < 0) with threshold-based ignition from incident
// heat flux. This is more physically realistic and depends on fuel moisture
// (wet fuels need more heat to ignite).
//
// Physics:
//   - Ignition requires a minimum incident heat flux to raise fuel to ignition
//   - Critical heat flux depends on fuel moisture content
//   - Wetter fuels require higher heat flux (more energy to evaporate water)
//   - Used in Balbi, Lautenberger models; enhances Rothermel approach
//
// Model:
//   q_crit [kW/m²] = q_base × (1 + k_moisture × M_fuel)
//
// where:
//   q_base      = base critical heat flux for dry fuel [kW/m²] (typically 10-15)
//   k_moisture  = moisture sensitivity coefficient (typically 2-4)
//   M_fuel      = fuel moisture content [fraction]
//
// Ignition criterion:
//   Ignite when: q_incident > q_crit
//
// Typical values:
//   Dry fuel (M=0.05):     q_crit ≈ 11 kW/m²
//   Moderate (M=0.15):     q_crit ≈ 16 kW/m²
//   Wet fuel (M=0.30):     q_crit ≈ 25 kW/m²
//
// References:
//   Drysdale, D. (2011). An Introduction to Fire Dynamics, 3rd ed. Wiley.
//   Quintiere, J.G. (2006). Fundamentals of Fire Phenomena. Wiley.
//   Balbi, J.H., et al. (2009). A physical model for wildland fires.
//     Combustion and Flame 156:2217-2230.
//   Lautenberger, C. (2013). Wildland fire modeling with an Eulerian level
//     set method. Fire Safety Journal 62:289-298.
// ============================================================================

#ifndef CRITICAL_HEAT_FLUX_H
#define CRITICAL_HEAT_FLUX_H

#include "constants.H"
#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_GpuLaunch.H>

using namespace amrex;

// Default parameters for critical heat flux model
namespace CriticalHeatFluxConst {
    constexpr Real Q_BASE_DEFAULT      = 12.0;  // Base critical heat flux [kW/m²]
    constexpr Real K_MOISTURE_DEFAULT  = 3.0;   // Moisture sensitivity coefficient
    constexpr Real Q_MIN               = 5.0;   // Minimum possible q_crit [kW/m²]
    constexpr Real Q_MAX               = 50.0;  // Maximum possible q_crit [kW/m²]
}

// ============================================================================
// Compute critical heat flux for ignition
//
// q_crit [kW/m²] = q_base × (1 + k_moisture × M_fuel)
//
// Parameters:
//   M_fuel      – fuel moisture content [fraction]
//   q_base      – base critical heat flux [kW/m²]
//   k_moisture  – moisture sensitivity coefficient
//
// Returns: critical heat flux [kW/m²]
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real compute_critical_heat_flux(Real M_fuel,
                                Real q_base    = CriticalHeatFluxConst::Q_BASE_DEFAULT,
                                Real k_moisture = CriticalHeatFluxConst::K_MOISTURE_DEFAULT)
{
    // q_crit = q_base × (1 + k × M_fuel)
    Real q_crit = q_base * (Real(1.0) + k_moisture * amrex::max(M_fuel, Real(0.0)));
    
    // Apply reasonable bounds
    q_crit = amrex::max(q_crit, CriticalHeatFluxConst::Q_MIN);
    q_crit = amrex::min(q_crit, CriticalHeatFluxConst::Q_MAX);
    
    return q_crit;
}

// ============================================================================
// Compute critical heat flux field
//
// Parameters:
//   q_crit_mf    – output: critical heat flux [kW/m²]
//   M_fuel_mf    – input: fuel moisture content [fraction]
//   q_base       – base critical heat flux [kW/m²]
//   k_moisture   – moisture sensitivity coefficient
// ============================================================================
inline void compute_critical_heat_flux_field(MultiFab& q_crit_mf,
                                             const MultiFab& M_fuel_mf,
                                             Real q_base    = CriticalHeatFluxConst::Q_BASE_DEFAULT,
                                             Real k_moisture = CriticalHeatFluxConst::K_MOISTURE_DEFAULT)
{
    for (MFIter mfi(q_crit_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto q_c       = q_crit_mf.array(mfi);
        auto const M_f = M_fuel_mf.const_array(mfi);
        const Real q_b = q_base;
        const Real k_m = k_moisture;

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            q_c(i, j, k) = compute_critical_heat_flux(M_f(i, j, k), q_b, k_m);
        });
    }
}

// ============================================================================
// Check ignition based on incident heat flux threshold
//
// Sets ignition_flag = 1.0 where q_incident > q_crit, 0.0 otherwise.
//
// Parameters:
//   ignition_flag_mf  – output: ignition flag (0 or 1)
//   q_incident_mf     – input: incident heat flux [kW/m²]
//   q_crit_mf         – input: critical heat flux threshold [kW/m²]
// ============================================================================
inline void check_heat_flux_ignition(MultiFab& ignition_flag_mf,
                                     const MultiFab& q_incident_mf,
                                     const MultiFab& q_crit_mf)
{
    for (MFIter mfi(ignition_flag_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto ign       = ignition_flag_mf.array(mfi);
        auto const q_inc = q_incident_mf.const_array(mfi);
        auto const q_c   = q_crit_mf.const_array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            // Ignite when incident flux exceeds critical threshold
            ign(i, j, k) = (q_inc(i, j, k) > q_c(i, j, k)) ? Real(1.0) : Real(0.0);
        });
    }
}

// ============================================================================
// Compute ignition probability from heat flux ratio
//
// For more gradual ignition near the threshold, use a probability model:
//   P_ignition = tanh(α × (q_incident/q_crit - 1))^+
//
// where (...)^+ means max(..., 0), and α controls the sharpness of transition.
//
// Parameters:
//   P_ignition_mf   – output: ignition probability [0, 1]
//   q_incident_mf   – input: incident heat flux [kW/m²]
//   q_crit_mf       – input: critical heat flux [kW/m²]
//   alpha           – sharpness parameter (default: 5.0)
//                     Higher α → sharper transition at q_crit
// ============================================================================
inline void compute_ignition_probability(MultiFab& P_ignition_mf,
                                         const MultiFab& q_incident_mf,
                                         const MultiFab& q_crit_mf,
                                         Real alpha = Real(5.0))
{
    for (MFIter mfi(P_ignition_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto P_ign     = P_ignition_mf.array(mfi);
        auto const q_inc = q_incident_mf.const_array(mfi);
        auto const q_c   = q_crit_mf.const_array(mfi);
        const Real a     = alpha;

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const Real q_ratio = (q_c(i, j, k) > Real(0.0))
                               ? q_inc(i, j, k) / q_c(i, j, k)
                               : Real(0.0);
            
            // Probability: tanh(α × (ratio - 1))^+
            const Real arg = a * (q_ratio - Real(1.0));
            Real prob = std::tanh(amrex::max(arg, Real(0.0)));
            
            // Ensure [0, 1] bounds
            prob = amrex::max(Real(0.0), amrex::min(prob, Real(1.0)));
            
            P_ign(i, j, k) = prob;
        });
    }
}

// ============================================================================
// Combined: compute critical flux and check ignition
//
// This is a convenience function combining the two steps:
//   1. Compute q_crit from moisture
//   2. Check ignition criterion
//
// Parameters:
//   ignition_flag_mf  – output: ignition flag (0 or 1)
//   q_incident_mf     – input: incident heat flux [kW/m²]
//   M_fuel_mf         – input: fuel moisture content [fraction]
//   q_base            – base critical heat flux [kW/m²]
//   k_moisture        – moisture sensitivity coefficient
// ============================================================================
inline void check_moisture_dependent_ignition(MultiFab& ignition_flag_mf,
                                              const MultiFab& q_incident_mf,
                                              const MultiFab& M_fuel_mf,
                                              Real q_base    = CriticalHeatFluxConst::Q_BASE_DEFAULT,
                                              Real k_moisture = CriticalHeatFluxConst::K_MOISTURE_DEFAULT)
{
    for (MFIter mfi(ignition_flag_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto ign       = ignition_flag_mf.array(mfi);
        auto const q_inc = q_incident_mf.const_array(mfi);
        auto const M_f   = M_fuel_mf.const_array(mfi);
        const Real q_b   = q_base;
        const Real k_m   = k_moisture;

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const Real q_crit = compute_critical_heat_flux(M_f(i, j, k), q_b, k_m);
            ign(i, j, k) = (q_inc(i, j, k) > q_crit) ? Real(1.0) : Real(0.0);
        });
    }
}

#endif // CRITICAL_HEAT_FLUX_H
