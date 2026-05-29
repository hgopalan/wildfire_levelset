// ============================================================================
// Duff Moisture and Smoldering Combustion Model
// ============================================================================
// Tracks deep organic layer (duff) moisture and implements smoldering
// combustion in these layers following fire passage.
//
// Physics:
//   - Duff layers (0.1–1+ m thick) have long moisture time constants (days-weeks)
//   - Smoldering combustion can persist for hours-days after flaming front passes
//   - Contributes to long-duration heat release, spotting, and smoke production
//   - FARSITE/FLAMMAP track duff moisture (DMC) and smoldering residual time
//
// Duff Moisture:
//   - Uses Canadian Forest Fire Danger Rating System Duff Moisture Code (DMC)
//   - Similar to existing fuel moisture tracking but on 10-20 day timelag
//   - Affects ignition probability of deep organic layers
//
// Smoldering Combustion:
//   - Exponential decay model: Q_smolder(t) = Q_0 × exp(-t / tau_smolder)
//   - Tau_smolder depends on duff moisture, depth, and packing
//   - Outputs: residual heat release, smolder time, smolder probability
//
// References:
//   Van Wagner, C.E. (1987). Development and structure of the Canadian Forest
//     Fire Weather Index System. Canadian Forestry Service Technical Report 35.
//   Albini, F.A. & Reinhardt, E.D. (1997). "Improved calibration of a large
//     fuel burnout model." International Journal of Wildland Fire, 7(1), 21-28.
//   Poulter, B. et al. (2015). "Contribution of semi-arid ecosystems to
//     interannual variability of the global carbon cycle." Nature, 509, 600-603.
//   Hille, M.G. & Stephens, S.L. (2005). "Mixed conifer forest duff consumption
//     during prescribed fires." Canadian Journal of Forest Research, 35, 2386-2395.
// ============================================================================

#ifndef DUFF_MOISTURE_SMOLDERING_H
#define DUFF_MOISTURE_SMOLDERING_H

#include "constants.H"
#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_GpuLaunch.H>
#include <cmath>

using namespace amrex;

// Default parameters for duff moisture and smoldering combustion
namespace DuffSmolderConst {
    // Duff layer properties
    constexpr Real DUFF_DEPTH_DEFAULT       = 0.05;  // Duff depth [m] (typical: 5-15 cm)
    constexpr Real DUFF_BULK_DENSITY        = 100.0; // Duff bulk density [kg/m³]
    constexpr Real DUFF_HEAT_CONTENT        = 18000.0; // Heat content [kJ/kg]
    
    // Smoldering combustion time constants
    constexpr Real TAU_SMOLDER_DRY          = 3600.0;  // Dry duff smolder time [s] (1 hour)
    constexpr Real TAU_SMOLDER_WET          = 14400.0; // Wet duff smolder time [s] (4 hours)
    constexpr Real MOISTURE_TRANSITION      = 0.30;    // Moisture transition point
    
    // Smoldering combustion fractions
    constexpr Real FRACTION_CONSUMED_DRY    = 0.95;    // Fraction consumed when dry
    constexpr Real FRACTION_CONSUMED_WET    = 0.30;    // Fraction consumed when wet
    
    // Heat release rate scaling
    constexpr Real Q_SMOLDER_SCALE          = 0.1;     // Smolder heat / flaming heat ratio
}

// ============================================================================
// compute_duff_smolder_fraction
//
// Computes the fraction of duff consumed during smoldering based on moisture.
//
// Model:
//   f_c = f_dry + (f_wet - f_dry) × (M_duff / M_transition)^n
//
// where:
//   M_duff       - duff moisture content [fraction]
//   M_transition - moisture transition point (default: 0.30)
//   f_dry        - fraction consumed when dry (default: 0.95)
//   f_wet        - fraction consumed when wet (default: 0.30)
//   n            - exponent (default: 1.5)
//
// Parameters:
//   M_duff - duff moisture content [fraction, 0-1]
//
// Returns: fraction consumed [0-1]
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real compute_duff_smolder_fraction(Real M_duff)
{
    const Real f_dry = DuffSmolderConst::FRACTION_CONSUMED_DRY;
    const Real f_wet = DuffSmolderConst::FRACTION_CONSUMED_WET;
    const Real M_trans = DuffSmolderConst::MOISTURE_TRANSITION;
    const Real exponent = Real(1.5);
    
    if (M_duff <= Real(0.0)) {
        return f_dry;
    } else if (M_duff >= M_trans) {
        return f_wet;
    } else {
        const Real ratio = M_duff / M_trans;
        return f_dry + (f_wet - f_dry) * std::pow(ratio, exponent);
    }
}

// ============================================================================
// compute_duff_smolder_time_constant
//
// Computes the smoldering combustion time constant based on duff moisture.
//
// Model:
//   tau = tau_dry + (tau_wet - tau_dry) × (M_duff / M_transition)
//
// where:
//   tau_dry      - dry duff smolder time constant [s] (default: 3600 s = 1 hour)
//   tau_wet      - wet duff smolder time constant [s] (default: 14400 s = 4 hours)
//   M_duff       - duff moisture content [fraction]
//   M_transition - moisture transition point (default: 0.30)
//
// Parameters:
//   M_duff - duff moisture content [fraction, 0-1]
//
// Returns: time constant [s]
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real compute_duff_smolder_time_constant(Real M_duff)
{
    const Real tau_dry = DuffSmolderConst::TAU_SMOLDER_DRY;
    const Real tau_wet = DuffSmolderConst::TAU_SMOLDER_WET;
    const Real M_trans = DuffSmolderConst::MOISTURE_TRANSITION;
    
    if (M_duff <= Real(0.0)) {
        return tau_dry;
    } else if (M_duff >= M_trans) {
        return tau_wet;
    } else {
        const Real ratio = M_duff / M_trans;
        return tau_dry + (tau_wet - tau_dry) * ratio;
    }
}

// ============================================================================
// compute_duff_smolder_heat_release
//
// Computes the instantaneous heat release rate from smoldering combustion.
//
// Model:
//   Q(t) = Q_0 × exp(-t / tau_smolder)
//
// where:
//   Q_0          - initial heat release rate [kW/m²]
//   t            - elapsed time since ignition [s]
//   tau_smolder  - smolder time constant [s]
//
// Parameters:
//   Q_initial     - initial heat release rate [kW/m²]
//   elapsed_time  - time since ignition [s]
//   tau_smolder   - smoldering time constant [s]
//
// Returns: heat release rate [kW/m²]
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real compute_duff_smolder_heat_release(Real Q_initial, Real elapsed_time, Real tau_smolder)
{
    const Real tau_safe = amrex::max(tau_smolder, Real(1.0));
    return Q_initial * std::exp(-elapsed_time / tau_safe);
}

// ============================================================================
// compute_duff_initial_heat_release
//
// Computes the initial heat release rate from duff smoldering based on
// duff properties and moisture.
//
// Model:
//   Q_0 = (rho_duff × delta_duff × h_duff × f_c) / tau_smolder
//
// where:
//   rho_duff     - duff bulk density [kg/m³]
//   delta_duff   - duff depth [m]
//   h_duff       - duff heat content [kJ/kg]
//   f_c          - fraction consumed
//   tau_smolder  - smolder time constant [s]
//
// Parameters:
//   delta_duff   - duff depth [m]
//   M_duff       - duff moisture content [fraction]
//
// Returns: initial heat release rate [kW/m²]
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real compute_duff_initial_heat_release(Real delta_duff, Real M_duff)
{
    const Real rho_duff = DuffSmolderConst::DUFF_BULK_DENSITY;
    const Real h_duff   = DuffSmolderConst::DUFF_HEAT_CONTENT;
    
    const Real f_c = compute_duff_smolder_fraction(M_duff);
    const Real tau_smolder = compute_duff_smolder_time_constant(M_duff);
    
    // Total energy available [kJ/m²]
    const Real E_total = rho_duff * delta_duff * h_duff * f_c;
    
    // Initial heat release rate [kJ/m²/s = kW/m²]
    const Real Q_0 = E_total / tau_smolder;
    
    return Q_0;
}

// ============================================================================
// compute_duff_smolder_probability
//
// Computes the probability of smoldering ignition based on duff moisture
// and fireline intensity.
//
// Model (empirical):
//   P_smolder = P_max × (1 - M_duff / M_x) × (I_B / I_threshold)^n
//
// where:
//   P_max        - maximum probability (default: 0.95)
//   M_duff       - duff moisture content [fraction]
//   M_x          - moisture of extinction for duff (default: 0.40)
//   I_B          - Byram fireline intensity [kW/m]
//   I_threshold  - threshold intensity for smoldering (default: 100 kW/m)
//   n            - exponent (default: 0.5)
//
// Parameters:
//   M_duff       - duff moisture content [fraction]
//   I_B          - Byram fireline intensity [kW/m]
//
// Returns: probability of smoldering ignition [0-1]
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real compute_duff_smolder_probability(Real M_duff, Real I_B)
{
    constexpr Real P_max = Real(0.95);
    constexpr Real M_x = Real(0.40);
    constexpr Real I_threshold = Real(100.0);
    constexpr Real exponent = Real(0.5);
    
    // Moisture factor (0 when M_duff >= M_x)
    const Real moisture_factor = (M_duff < M_x)
                               ? (Real(1.0) - M_duff / M_x)
                               : Real(0.0);
    
    // Intensity factor (0 when I_B < I_threshold)
    const Real intensity_factor = (I_B > I_threshold)
                                ? std::pow(I_B / I_threshold, exponent)
                                : Real(0.0);
    
    // Combined probability
    Real P_smolder = P_max * moisture_factor * intensity_factor;
    
    // Clamp to [0, 1]
    P_smolder = amrex::max(Real(0.0), amrex::min(Real(1.0), P_smolder));
    
    return P_smolder;
}

// ============================================================================
// update_duff_smolder_state
//
// Updates duff smoldering state for one timestep.
//
// MultiFab layout (duff_smolder_mf):
//   Component 0 - duff moisture content [fraction]
//   Component 1 - smoldering heat release rate [kW/m²]
//   Component 2 - elapsed smolder time [s]
//   Component 3 - smolder probability [-]
//
// Parameters:
//   duff_smolder_mf       - duff smoldering state MultiFab (4 components)
//   fireline_intensity_mf - Byram fireline intensity [kW/m]
//   phi                   - level-set field (phi < 0 = burned)
//   dt                    - timestep [s]
//   duff_depth_mf         - duff depth field [m] (optional; uses default if nullptr)
//   duff_moisture_mf      - duff moisture field [fraction] (optional; uses DMC if nullptr)
// ============================================================================
inline void update_duff_smolder_state(
    MultiFab&       duff_smolder_mf,
    const MultiFab& fireline_intensity_mf,
    const MultiFab& phi,
    Real            dt,
    const MultiFab* duff_depth_mf = nullptr,
    const MultiFab* duff_moisture_mf = nullptr)
{
    for (MFIter mfi(duff_smolder_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto smolder = duff_smolder_mf.array(mfi);
        auto I_B     = fireline_intensity_mf.const_array(mfi);
        auto p       = phi.const_array(mfi);
        
        // Optional per-cell duff depth and moisture
        const bool has_depth = (duff_depth_mf != nullptr);
        const bool has_moisture = (duff_moisture_mf != nullptr);
        
        auto delta_d = has_depth ? duff_depth_mf->const_array(mfi)
                                 : fireline_intensity_mf.const_array(mfi);  // Dummy
        auto M_d     = has_moisture ? duff_moisture_mf->const_array(mfi)
                                    : fireline_intensity_mf.const_array(mfi);  // Dummy
        
        const Real dt_loc = dt;
        const Real default_duff_depth = DuffSmolderConst::DUFF_DEPTH_DEFAULT;
        const Real default_duff_moisture = Real(0.20);  // Default 20% moisture
        
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            // Skip unburned cells
            if (p(i, j, k) >= Real(0.0)) {
                smolder(i, j, k, 0) = default_duff_moisture;  // duff moisture
                smolder(i, j, k, 1) = Real(0.0);              // heat release
                smolder(i, j, k, 2) = Real(0.0);              // elapsed time
                smolder(i, j, k, 3) = Real(0.0);              // probability
                return;
            }
            
            // Get duff properties
            const Real delta_duff = has_depth ? delta_d(i, j, k) : default_duff_depth;
            const Real M_duff     = has_moisture ? M_d(i, j, k) : default_duff_moisture;
            
            // Get current state
            Real elapsed_time = smolder(i, j, k, 2);
            Real Q_current    = smolder(i, j, k, 1);
            
            // Check if smoldering has just started (first timestep after burning)
            const bool just_ignited = (elapsed_time < Real(0.1));
            
            if (just_ignited) {
                // Compute smolder probability
                const Real P_smolder = compute_duff_smolder_probability(M_duff, I_B(i, j, k));
                
                // Initialize smoldering if probability > threshold
                // (In a stochastic version, this would use a random draw)
                if (P_smolder > Real(0.5)) {
                    // Compute initial heat release rate
                    Q_current = compute_duff_initial_heat_release(delta_duff, M_duff);
                    elapsed_time = Real(0.0);
                    
                    smolder(i, j, k, 3) = P_smolder;
                } else {
                    // No smoldering ignition
                    smolder(i, j, k, 1) = Real(0.0);
                    smolder(i, j, k, 2) = Real(0.0);
                    smolder(i, j, k, 3) = P_smolder;
                    return;
                }
            }
            
            // Update elapsed time
            elapsed_time += dt_loc;
            
            // Compute current heat release using exponential decay
            const Real tau_smolder = compute_duff_smolder_time_constant(M_duff);
            Q_current = compute_duff_smolder_heat_release(
                compute_duff_initial_heat_release(delta_duff, M_duff),
                elapsed_time,
                tau_smolder
            );
            
            // Update state
            smolder(i, j, k, 0) = M_duff;
            smolder(i, j, k, 1) = Q_current;
            smolder(i, j, k, 2) = elapsed_time;
            // Component 3 (probability) is already set or unchanged
        });
    }
}

#endif // DUFF_MOISTURE_SMOLDERING_H
