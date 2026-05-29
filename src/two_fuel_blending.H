// ============================================================================
// Two-Fuel Model Blending
// ============================================================================
// Blends rate of spread from two different fuel models based on fractional
// coverage. Common scenarios include:
//   - Grass/timber transitions (Rothermel + FBP or Balbi)
//   - Surface/crown fire mixing (already handled by crown_ros_override.H)
//   - Fuel type transitions in heterogeneous landscapes
//
// This extends the existing passive crown fire blending (Van Wagner 1977)
// to general two-fuel scenarios following Finney (1998) / Andrews (2018)
// blend formulas.
//
// Blending methods:
//   1. Linear blend: R = f × R_1 + (1-f) × R_2
//   2. Harmonic mean: 1/R = f/R_1 + (1-f)/R_2  (conservative, fire travels through slowest fuel)
//   3. Maximum: R = max(R_1, R_2)  (optimistic, fire finds fastest path)
//   4. Finney weighted: R = sqrt(f × R_1² + (1-f) × R_2²)  (FARSITE-style vector combination)
//
// References:
//   Finney, M.A. (1998). FARSITE: Fire Area Simulator. USDA Forest Service RMRS-RP-4.
//   Andrews, P.L. (2018). The Rothermel surface fire spread model and associated
//     developments: A comprehensive explanation. USDA Forest Service RMRS-GTR-371.
//   Stocks, B.J. et al. (1989). The Canadian Forest Fire Danger Rating System.
//     Canadian Forestry Service.
// ============================================================================

#ifndef TWO_FUEL_BLENDING_H
#define TWO_FUEL_BLENDING_H

#include "constants.H"
#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_GpuLaunch.H>
#include <cmath>

using namespace amrex;

// Blending method enumeration
enum class BlendMethod {
    LINEAR,      // Simple weighted average
    HARMONIC,    // Harmonic mean (conservative)
    MAXIMUM,     // Take maximum ROS (optimistic)
    FINNEY       // Finney (1998) FARSITE-style weighted quadratic
};

// ============================================================================
// blend_two_fuel_ros_linear
//
// Linear blending: R = f × R_1 + (1-f) × R_2
//
// Simple weighted average. Appropriate when the two fuel types are well-mixed
// at subgrid scale.
//
// Parameters:
//   R_1     - ROS from fuel model 1 [m/s]
//   R_2     - ROS from fuel model 2 [m/s]
//   f       - fractional coverage of fuel model 1 [0-1]
//
// Returns: blended ROS [m/s]
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real blend_two_fuel_ros_linear(Real R_1, Real R_2, Real f)
{
    const Real f_clamped = amrex::max(Real(0.0), amrex::min(Real(1.0), f));
    return f_clamped * R_1 + (Real(1.0) - f_clamped) * R_2;
}

// ============================================================================
// blend_two_fuel_ros_harmonic
//
// Harmonic mean blending: 1/R = f/R_1 + (1-f)/R_2
//
// Conservative estimate. Fire must traverse both fuel types, so the slower
// fuel dominates. Appropriate for patchy fuel distributions where fire must
// pass through both types.
//
// Parameters:
//   R_1     - ROS from fuel model 1 [m/s]
//   R_2     - ROS from fuel model 2 [m/s]
//   f       - fractional coverage of fuel model 1 [0-1]
//
// Returns: blended ROS [m/s]
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real blend_two_fuel_ros_harmonic(Real R_1, Real R_2, Real f)
{
    const Real f_clamped = amrex::max(Real(0.0), amrex::min(Real(1.0), f));
    const Real R_1_safe = amrex::max(R_1, Real(1.0e-8));
    const Real R_2_safe = amrex::max(R_2, Real(1.0e-8));
    
    const Real inv_R = f_clamped / R_1_safe + (Real(1.0) - f_clamped) / R_2_safe;
    return Real(1.0) / inv_R;
}

// ============================================================================
// blend_two_fuel_ros_maximum
//
// Maximum ROS blending: R = max(R_1, R_2)
//
// Optimistic estimate. Fire finds the fastest path through the landscape.
// Appropriate when fuel types are arranged such that fire can preferentially
// spread through the faster-burning fuel.
//
// Note: This ignores the fractional coverage f. It's included for API consistency.
//
// Parameters:
//   R_1     - ROS from fuel model 1 [m/s]
//   R_2     - ROS from fuel model 2 [m/s]
//   f       - fractional coverage of fuel model 1 [0-1] (not used)
//
// Returns: blended ROS [m/s]
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real blend_two_fuel_ros_maximum(Real R_1, Real R_2, Real /*f*/)
{
    return amrex::max(R_1, R_2);
}

// ============================================================================
// blend_two_fuel_ros_finney
//
// Finney (1998) FARSITE-style blending: R = sqrt(f × R_1² + (1-f) × R_2²)
//
// Combines ROS values as vectors (similar to how FARSITE combines slope and
// wind effects). Provides a smooth transition between the two fuels with
// intermediate values.
//
// This is the method used in FARSITE for blending surface and crown fire
// spread rates.
//
// Parameters:
//   R_1     - ROS from fuel model 1 [m/s]
//   R_2     - ROS from fuel model 2 [m/s]
//   f       - fractional coverage of fuel model 1 [0-1]
//
// Returns: blended ROS [m/s]
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real blend_two_fuel_ros_finney(Real R_1, Real R_2, Real f)
{
    const Real f_clamped = amrex::max(Real(0.0), amrex::min(Real(1.0), f));
    const Real R_sq = f_clamped * R_1 * R_1 + (Real(1.0) - f_clamped) * R_2 * R_2;
    return std::sqrt(R_sq);
}

// ============================================================================
// blend_two_fuel_ros
//
// General two-fuel ROS blending function with selectable method.
//
// Parameters:
//   R_1     - ROS from fuel model 1 [m/s]
//   R_2     - ROS from fuel model 2 [m/s]
//   f       - fractional coverage of fuel model 1 [0-1]
//   method  - blending method (LINEAR, HARMONIC, MAXIMUM, FINNEY)
//
// Returns: blended ROS [m/s]
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real blend_two_fuel_ros(Real R_1, Real R_2, Real f, BlendMethod method)
{
    switch (method) {
        case BlendMethod::LINEAR:
            return blend_two_fuel_ros_linear(R_1, R_2, f);
        case BlendMethod::HARMONIC:
            return blend_two_fuel_ros_harmonic(R_1, R_2, f);
        case BlendMethod::MAXIMUM:
            return blend_two_fuel_ros_maximum(R_1, R_2, f);
        case BlendMethod::FINNEY:
            return blend_two_fuel_ros_finney(R_1, R_2, f);
        default:
            return blend_two_fuel_ros_linear(R_1, R_2, f);  // Default to linear
    }
}

// ============================================================================
// apply_two_fuel_blending
//
// Applies two-fuel blending to a MultiFab of ROS values.
//
// This function modifies R_mf in place by blending with R_secondary_mf
// according to the fractional coverage in blend_fraction_mf.
//
// Usage example:
//   // Blend Rothermel surface fire (R_mf) with Balbi model (R_balbi_mf)
//   // using grass coverage fraction (grass_fraction_mf)
//   apply_two_fuel_blending(R_mf, R_balbi_mf, grass_fraction_mf, BlendMethod::LINEAR);
//
// Parameters:
//   R_mf                - primary ROS field [m/s] (modified in place)
//   R_secondary_mf      - secondary ROS field [m/s]
//   blend_fraction_mf   - fractional coverage of primary fuel [0-1]
//   method              - blending method
// ============================================================================
inline void apply_two_fuel_blending(
    MultiFab&       R_mf,
    const MultiFab& R_secondary_mf,
    const MultiFab& blend_fraction_mf,
    BlendMethod     method = BlendMethod::LINEAR)
{
    for (MFIter mfi(R_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto R_1 = R_mf.array(mfi);
        auto R_2 = R_secondary_mf.const_array(mfi);
        auto f   = blend_fraction_mf.const_array(mfi);
        
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const Real R_primary   = R_1(i, j, k);
            const Real R_secondary = R_2(i, j, k);
            const Real fraction    = f(i, j, k);
            
            R_1(i, j, k) = blend_two_fuel_ros(R_primary, R_secondary, fraction, method);
        });
    }
}

#endif // TWO_FUEL_BLENDING_H
