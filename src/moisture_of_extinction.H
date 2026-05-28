// ============================================================================
// Fuel Moisture of Extinction Gradient Model
// ============================================================================
// Makes moisture of extinction (M_x) depend on fuel properties rather than
// using a constant 0.30 value. M_x varies by fuel type (grass vs. timber),
// affecting marginal burning conditions.
//
// Physics:
//   - Moisture of extinction is not universal; it depends on fuel characteristics
//   - Fine fuels (high SAV) have lower M_x (easier to extinguish)
//   - Coarse fuels (low SAV) have higher M_x (harder to extinguish)
//   - Current Rothermel uses constant M_x = 0.30
//   - More accurate fire/no-fire transitions with variable M_x
//
// Model (Rothermel 1972):
//   M_x = M_x_min + k_sigma × (σ_ref / σ_fuel)^n
//
// Alternative simplified form:
//   M_x = 0.12 + 0.28 × (σ_fuel / σ_ref)^(-0.3)
//
// where:
//   σ_fuel  = fuel surface-area-to-volume ratio [ft⁻¹]
//   σ_ref   = reference SAV (typically 1739 ft⁻¹ for 1-hr fuels)
//   M_x_min = minimum moisture of extinction (typically 0.12-0.15)
//   k_sigma = SAV sensitivity coefficient
//   n       = exponent (typically 0.2-0.4)
//
// Typical values:
//   Grass (σ ~ 2000 ft⁻¹):  M_x ≈ 0.25
//   1-hr fuels (σ ~ 1739):   M_x ≈ 0.30
//   10-hr fuels (σ ~ 109):   M_x ≈ 0.40
//
// References:
//   Rothermel, R.C. (1972). A mathematical model for predicting fire spread
//     in wildland fuels. USDA Forest Service Research Paper INT-115.
//   Anderson, H.E. (1970). Forest fuel ignitibility. Fire Technology 6(4):312-319.
//   Albini, F.A. (1976). Estimating wildfire behavior and effects.
//     USDA Forest Service INT-30.
// ============================================================================

#ifndef MOISTURE_OF_EXTINCTION_H
#define MOISTURE_OF_EXTINCTION_H

#include "constants.H"
#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_GpuLaunch.H>
#include <cmath>

using namespace amrex;

// Default parameters for M_x model
namespace MoistureExtinctionConst {
    constexpr Real M_X_MIN_DEFAULT    = 0.12;   // Minimum M_x
    constexpr Real M_X_MAX_DEFAULT    = 0.45;   // Maximum M_x
    constexpr Real SIGMA_REF_DEFAULT  = 1739.0; // Reference SAV [ft⁻¹] (1-hr fuel)
    constexpr Real EXPONENT_DEFAULT   = -0.3;   // SAV exponent
}

// ============================================================================
// Compute fuel-dependent moisture of extinction
//
// M_x = 0.12 + 0.28 × (σ_fuel / σ_ref)^n
//
// Parameters:
//   sigma_fuel  – fuel SAV [ft⁻¹]
//   sigma_ref   – reference SAV [ft⁻¹]
//   M_x_min     – minimum M_x
//   exponent    – SAV exponent
//
// Returns: moisture of extinction [fraction]
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real compute_moisture_of_extinction(Real sigma_fuel,
                                    Real sigma_ref = MoistureExtinctionConst::SIGMA_REF_DEFAULT,
                                    Real M_x_min   = MoistureExtinctionConst::M_X_MIN_DEFAULT,
                                    Real exponent  = MoistureExtinctionConst::EXPONENT_DEFAULT)
{
    // Avoid division by zero
    if (sigma_ref <= Real(0.0) || sigma_fuel <= Real(0.0)) {
        return Real(0.30);  // Default value
    }
    
    // Sigma ratio
    Real sigma_ratio = sigma_fuel / sigma_ref;
    
    // M_x = M_x_min + (M_x_max - M_x_min) × (σ_fuel / σ_ref)^n
    // Simplified: M_x = 0.12 + 0.28 × (σ_fuel / σ_ref)^(-0.3)
    Real M_x = M_x_min + Real(0.28) * std::pow(sigma_ratio, exponent);
    
    // Apply reasonable bounds
    M_x = amrex::max(M_x, MoistureExtinctionConst::M_X_MIN_DEFAULT);
    M_x = amrex::min(M_x, MoistureExtinctionConst::M_X_MAX_DEFAULT);
    
    return M_x;
}

// ============================================================================
// Alternative: Rothermel (1972) original formulation
//
// M_x = M_x_base × f(mineral_content, SAV)
//
// This is a more complete formulation that includes mineral content effects.
//
// Parameters:
//   sigma_fuel     – fuel SAV [ft⁻¹]
//   S_T            – total mineral content [fraction]
//   sigma_ref      – reference SAV [ft⁻¹]
//
// Returns: moisture of extinction [fraction]
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real compute_M_x_rothermel(Real sigma_fuel,
                          Real S_T         = Real(0.0555),  // Default mineral content
                          Real sigma_ref   = MoistureExtinctionConst::SIGMA_REF_DEFAULT)
{
    // Rothermel (1972) original empirical relationship
    // M_x ≈ 0.30 + 0.10 × ln(σ/1000)  [approximate]
    
    if (sigma_fuel <= Real(0.0)) {
        return Real(0.30);
    }
    
    // SAV effect (logarithmic)
    Real sigma_factor = std::log(sigma_fuel / Real(1000.0));
    Real M_x = Real(0.30) + Real(0.05) * sigma_factor;
    
    // Mineral content correction (higher minerals → lower M_x)
    Real mineral_factor = Real(1.0) - Real(0.5) * S_T;
    M_x *= mineral_factor;
    
    // Apply bounds
    M_x = amrex::max(M_x, MoistureExtinctionConst::M_X_MIN_DEFAULT);
    M_x = amrex::min(M_x, MoistureExtinctionConst::M_X_MAX_DEFAULT);
    
    return M_x;
}

// ============================================================================
// Compute M_x field from fuel SAV field
//
// Parameters:
//   M_x_mf      – output: moisture of extinction [fraction]
//   sigma_mf    – input: fuel SAV [ft⁻¹]
//   sigma_ref   – reference SAV [ft⁻¹]
//   M_x_min     – minimum M_x
//   exponent    – SAV exponent
// ============================================================================
inline void compute_M_x_field(MultiFab& M_x_mf,
                              const MultiFab& sigma_mf,
                              Real sigma_ref = MoistureExtinctionConst::SIGMA_REF_DEFAULT,
                              Real M_x_min   = MoistureExtinctionConst::M_X_MIN_DEFAULT,
                              Real exponent  = MoistureExtinctionConst::EXPONENT_DEFAULT)
{
    for (MFIter mfi(M_x_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto M_x       = M_x_mf.array(mfi);
        auto const sig = sigma_mf.const_array(mfi);
        const Real s_r = sigma_ref;
        const Real M_m = M_x_min;
        const Real exp = exponent;

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            M_x(i, j, k) = compute_moisture_of_extinction(
                sig(i, j, k), s_r, M_m, exp
            );
        });
    }
}

// ============================================================================
// Compute M_x field using Rothermel formulation with mineral content
//
// Parameters:
//   M_x_mf      – output: moisture of extinction [fraction]
//   sigma_mf    – input: fuel SAV [ft⁻¹]
//   S_T_mf      – input: total mineral content [fraction]
//   sigma_ref   – reference SAV [ft⁻¹]
// ============================================================================
inline void compute_M_x_field_rothermel(MultiFab& M_x_mf,
                                        const MultiFab& sigma_mf,
                                        const MultiFab& S_T_mf,
                                        Real sigma_ref = MoistureExtinctionConst::SIGMA_REF_DEFAULT)
{
    for (MFIter mfi(M_x_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto M_x       = M_x_mf.array(mfi);
        auto const sig = sigma_mf.const_array(mfi);
        auto const S_T = S_T_mf.const_array(mfi);
        const Real s_r = sigma_ref;

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            M_x(i, j, k) = compute_M_x_rothermel(
                sig(i, j, k), S_T(i, j, k), s_r
            );
        });
    }
}

// ============================================================================
// Classify fuel ignitability based on M_x
//
// Returns a classification:
//   0: Very difficult to ignite (M_x > 0.40)
//   1: Difficult (0.35 < M_x <= 0.40)
//   2: Moderate (0.25 < M_x <= 0.35)
//   3: Easy to ignite (M_x <= 0.25)
//
// Parameters:
//   ignitability_mf  – output: ignitability class (0-3)
//   M_x_mf           – input: moisture of extinction [fraction]
// ============================================================================
inline void classify_fuel_ignitability(MultiFab& ignitability_mf,
                                       const MultiFab& M_x_mf)
{
    for (MFIter mfi(ignitability_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto ign       = ignitability_mf.array(mfi);
        auto const M_x = M_x_mf.const_array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const Real mx = M_x(i, j, k);
            
            if (mx > Real(0.40)) {
                ign(i, j, k) = Real(0.0);  // Very difficult
            } else if (mx > Real(0.35)) {
                ign(i, j, k) = Real(1.0);  // Difficult
            } else if (mx > Real(0.25)) {
                ign(i, j, k) = Real(2.0);  // Moderate
            } else {
                ign(i, j, k) = Real(3.0);  // Easy
            }
        });
    }
}

#endif // MOISTURE_OF_EXTINCTION_H
