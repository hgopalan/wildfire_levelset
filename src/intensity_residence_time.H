// ============================================================================
// Intensity-Dependent Flame Residence Time Model
// ============================================================================
// Enhances the basic burnout_time.H model with intensity-dependent residence
// time. Higher intensity fires burn more vigorously and complete combustion
// faster than low-intensity fires.
//
// Physics:
//   - Current models use constant τ_residence
//   - Real fires: higher intensity → faster burnout (more vigorous combustion)
//   - Low intensity smoldering: longer residence time
//   - High intensity flaming: shorter residence time
//
// Model (Byram 1959):
//   τ_res [s] = τ_base × (I_base / I_byram)^α
//
// where:
//   τ_base   = baseline residence time [s] at reference intensity
//   I_base   = reference intensity [kW/m] (typically 100-500 kW/m)
//   I_byram  = current fireline intensity [kW/m]
//   α        = intensity exponent (typically 0.3-0.5)
//
// For I_byram > I_base: τ_res < τ_base (faster burnout)
// For I_byram < I_base: τ_res > τ_base (slower burnout)
//
// This refines the existing burnout model by making residence time dynamic
// rather than constant.
//
// References:
//   Byram, G.M. (1959). Combustion of forest fuels. In: Davis, K.P. (ed.)
//     Forest Fire: Control and Use. McGraw-Hill, New York.
//   Rothermel, R.C. (1972). A mathematical model for predicting fire spread
//     in wildland fuels. USDA Forest Service Research Paper INT-115.
//   Anderson, H.E. (1969). Heat transfer and fire spread.
//     USDA Forest Service Research Paper INT-69.
// ============================================================================

#ifndef INTENSITY_RESIDENCE_TIME_H
#define INTENSITY_RESIDENCE_TIME_H

#include "constants.H"
#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_GpuLaunch.H>
#include <cmath>

using namespace amrex;

// Default parameters for intensity-dependent residence time model
namespace IntensityResidenceConst {
    constexpr Real TAU_BASE_DEFAULT  = 60.0;   // Baseline residence time [s]
    constexpr Real I_BASE_DEFAULT    = 300.0;  // Reference intensity [kW/m]
    constexpr Real ALPHA_DEFAULT     = 0.4;    // Intensity exponent
    constexpr Real TAU_MIN           = 10.0;   // Minimum residence time [s]
    constexpr Real TAU_MAX           = 600.0;  // Maximum residence time [s]
}

// ============================================================================
// Compute intensity-dependent residence time
//
// τ_res [s] = τ_base × (I_base / I_byram)^α
//
// Parameters:
//   I_byram   – Byram fireline intensity [kW/m]
//   tau_base  – baseline residence time [s]
//   I_base    – reference intensity [kW/m]
//   alpha     – intensity exponent
//
// Returns: residence time [s]
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real compute_intensity_residence_time(Real I_byram,
                                      Real tau_base = IntensityResidenceConst::TAU_BASE_DEFAULT,
                                      Real I_base   = IntensityResidenceConst::I_BASE_DEFAULT,
                                      Real alpha    = IntensityResidenceConst::ALPHA_DEFAULT)
{
    // Avoid division by zero or negative intensities
    if (I_byram <= Real(0.0)) {
        return IntensityResidenceConst::TAU_MAX;  // Very long residence for no fire
    }
    
    // τ_res = τ_base × (I_base / I_byram)^α
    const Real ratio = I_base / I_byram;
    Real tau_res = tau_base * std::pow(ratio, alpha);
    
    // Apply reasonable bounds
    tau_res = amrex::max(tau_res, IntensityResidenceConst::TAU_MIN);
    tau_res = amrex::min(tau_res, IntensityResidenceConst::TAU_MAX);
    
    return tau_res;
}

// ============================================================================
// Compute intensity-dependent residence time field
//
// Parameters:
//   tau_residence_mf      – output: residence time [s]
//   fireline_intensity_mf – input: Byram intensity [kW/m]
//   tau_base              – baseline residence time [s]
//   I_base                – reference intensity [kW/m]
//   alpha                 – intensity exponent
// ============================================================================
inline void compute_intensity_residence_field(MultiFab& tau_residence_mf,
                                              const MultiFab& fireline_intensity_mf,
                                              Real tau_base = IntensityResidenceConst::TAU_BASE_DEFAULT,
                                              Real I_base   = IntensityResidenceConst::I_BASE_DEFAULT,
                                              Real alpha    = IntensityResidenceConst::ALPHA_DEFAULT)
{
    for (MFIter mfi(tau_residence_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto tau       = tau_residence_mf.array(mfi);
        auto const I_B = fireline_intensity_mf.const_array(mfi);
        const Real t_b = tau_base;
        const Real I_b = I_base;
        const Real a   = alpha;

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            tau(i, j, k) = compute_intensity_residence_time(I_B(i, j, k), t_b, I_b, a);
        });
    }
}

// ============================================================================
// Compute burnout time with intensity-dependent residence
//
// This is an enhanced version of the standard compute_burnout_time that uses
// intensity-dependent residence time instead of a constant value.
//
// Parameters:
//   burnout_time_mf       – output: burnout time [s]; -1.0 = not yet burned
//   arrival_time_mf       – input: simulation time of ignition [s]
//   fireline_intensity_mf – input: Byram intensity [kW/m]
//   tau_base              – baseline residence time [s]
//   I_base                – reference intensity [kW/m]
//   alpha                 – intensity exponent
// ============================================================================
inline void compute_intensity_burnout_time(MultiFab& burnout_time_mf,
                                           const MultiFab& arrival_time_mf,
                                           const MultiFab& fireline_intensity_mf,
                                           Real tau_base = IntensityResidenceConst::TAU_BASE_DEFAULT,
                                           Real I_base   = IntensityResidenceConst::I_BASE_DEFAULT,
                                           Real alpha    = IntensityResidenceConst::ALPHA_DEFAULT)
{
    for (MFIter mfi(burnout_time_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto bot       = burnout_time_mf.array(mfi);
        auto const at  = arrival_time_mf.const_array(mfi);
        auto const I_B = fireline_intensity_mf.const_array(mfi);
        const Real t_b = tau_base;
        const Real I_b = I_base;
        const Real a   = alpha;

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const Real t_arrive = at(i, j, k);
            
            // Not yet ignited
            if (t_arrive < Real(0.0)) {
                bot(i, j, k) = Real(-1.0);
                return;
            }
            
            // Compute intensity-dependent residence time
            const Real tau_res = compute_intensity_residence_time(
                I_B(i, j, k), t_b, I_b, a
            );
            
            // Burnout time = arrival time + residence time
            bot(i, j, k) = t_arrive + tau_res;
        });
    }
}

// ============================================================================
// Classify combustion rate based on residence time
//
// Returns a classification:
//   0: Smoldering (τ > 300 s) - very slow combustion
//   1: Low intensity flaming (120 < τ <= 300 s)
//   2: Moderate flaming (60 < τ <= 120 s)
//   3: High intensity flaming (30 < τ <= 60 s)
//   4: Very intense flaming (τ <= 30 s) - rapid consumption
//
// Parameters:
//   combustion_class_mf  – output: combustion classification (0-4)
//   tau_residence_mf     – input: residence time [s]
// ============================================================================
inline void classify_combustion_rate(MultiFab& combustion_class_mf,
                                     const MultiFab& tau_residence_mf)
{
    for (MFIter mfi(combustion_class_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto cc        = combustion_class_mf.array(mfi);
        auto const tau = tau_residence_mf.const_array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const Real t = tau(i, j, k);
            
            if (t > Real(300.0)) {
                cc(i, j, k) = Real(0.0);  // Smoldering
            } else if (t > Real(120.0)) {
                cc(i, j, k) = Real(1.0);  // Low intensity
            } else if (t > Real(60.0)) {
                cc(i, j, k) = Real(2.0);  // Moderate
            } else if (t > Real(30.0)) {
                cc(i, j, k) = Real(3.0);  // High intensity
            } else {
                cc(i, j, k) = Real(4.0);  // Very intense
            }
        });
    }
}

#endif // INTENSITY_RESIDENCE_TIME_H
