// ============================================================================
// Wind-Fuel Interaction Feedback Model
// ============================================================================
// Adjusts effective wind speed based on fuel structure (sheltering, penetration
// depth). Dense fuels shelter the surface from wind, reducing effective wind
// speed at the fuel bed level.
//
// Physics:
//   - Dense fuels create a sheltering effect
//   - Wind speed decreases exponentially through canopy
//   - Leaf area index (LAI) determines attenuation strength
//   - Current WAF models are height-based only; this adds fuel structure feedback
//
// Model (Wilson 1985, Albini & Baughman 1979):
//   U_eff [m/s] = U_ref × exp(-a × LAI × (1 - z/h))
//
// where:
//   U_ref  = reference wind speed (typically at 20 ft / 6.1 m) [m/s]
//   a      = extinction coefficient (typically 0.5-2.0)
//   LAI    = leaf area index [dimensionless]
//   z      = height of interest within canopy [m] (typically midflame height)
//   h      = canopy height [m]
//
// LAI estimation from fuel bed properties:
//   LAI ≈ k_LAI × δ × β
//
// where:
//   k_LAI  = empirical coefficient (typically 3-6 for wildland fuels)
//   δ      = fuel bed depth [m]
//   β      = packing ratio [dimensionless]
//
// References:
//   Wilson, N.R. (1985). A mathematical model for estimating surface-layer
//     wind speeds in a vegetative canopy. Agricultural and Forest Meteorology
//     34(2-3):191-208.
//   Albini, F.A., & Baughman, R.G. (1979). Estimating windspeeds for predicting
//     wildland fire behavior. USDA Forest Service Research Paper INT-221.
//   Finney, M.A. (2006). WindNinja. USDA Forest Service RMRS.
// ============================================================================

#ifndef WIND_FUEL_INTERACTION_H
#define WIND_FUEL_INTERACTION_H

#include "constants.H"
#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_GpuLaunch.H>
#include <cmath>

using namespace amrex;

// Default parameters for wind-fuel interaction model
namespace WindFuelConst {
    constexpr Real EXTINCTION_COEFF_DEFAULT = 1.0;    // Extinction coefficient
    constexpr Real K_LAI_DEFAULT            = 4.0;    // LAI estimation coefficient
    constexpr Real Z_FRACTION_DEFAULT       = 0.5;    // Midflame height as fraction of canopy
    constexpr Real LAI_MIN                  = 0.1;    // Minimum LAI
    constexpr Real LAI_MAX                  = 10.0;   // Maximum LAI
}

// ============================================================================
// Estimate leaf area index from fuel bed properties
//
// LAI ≈ k_LAI × δ × β
//
// Parameters:
//   fuel_depth     – fuel bed depth [m]
//   packing_ratio  – packing ratio [dimensionless]
//   k_LAI          – empirical coefficient
//
// Returns: leaf area index [dimensionless]
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real estimate_LAI(Real fuel_depth,
                 Real packing_ratio,
                 Real k_LAI = WindFuelConst::K_LAI_DEFAULT)
{
    Real LAI = k_LAI * fuel_depth * packing_ratio;
    
    // Apply reasonable bounds
    LAI = amrex::max(LAI, WindFuelConst::LAI_MIN);
    LAI = amrex::min(LAI, WindFuelConst::LAI_MAX);
    
    return LAI;
}

// ============================================================================
// Compute wind reduction factor through fuel canopy
//
// f_wind = exp(-a × LAI × (1 - z/h))
//
// Parameters:
//   LAI      – leaf area index [dimensionless]
//   z_h      – height fraction z/h (typically 0.5 for midflame)
//   a        – extinction coefficient
//
// Returns: wind reduction factor [0, 1]
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real compute_wind_reduction_factor(Real LAI,
                                   Real z_h = WindFuelConst::Z_FRACTION_DEFAULT,
                                   Real a   = WindFuelConst::EXTINCTION_COEFF_DEFAULT)
{
    // f = exp(-a × LAI × (1 - z/h))
    Real exponent = -a * LAI * (Real(1.0) - z_h);
    Real factor = std::exp(exponent);
    
    // Ensure [0, 1] range
    factor = amrex::max(Real(0.0), amrex::min(factor, Real(1.0)));
    
    return factor;
}

// ============================================================================
// Compute effective wind speed with fuel sheltering
//
// U_eff = U_ref × exp(-a × LAI × (1 - z/h))
//
// Parameters:
//   U_ref        – reference wind speed [m/s]
//   fuel_depth   – fuel bed depth [m]
//   packing_ratio – packing ratio [dimensionless]
//   k_LAI        – LAI estimation coefficient
//   z_h          – height fraction z/h
//   a            – extinction coefficient
//
// Returns: effective wind speed [m/s]
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real compute_sheltered_wind_speed(Real U_ref,
                                  Real fuel_depth,
                                  Real packing_ratio,
                                  Real k_LAI = WindFuelConst::K_LAI_DEFAULT,
                                  Real z_h   = WindFuelConst::Z_FRACTION_DEFAULT,
                                  Real a     = WindFuelConst::EXTINCTION_COEFF_DEFAULT)
{
    // Estimate LAI from fuel properties
    Real LAI = estimate_LAI(fuel_depth, packing_ratio, k_LAI);
    
    // Compute wind reduction factor
    Real f_wind = compute_wind_reduction_factor(LAI, z_h, a);
    
    // Effective wind speed
    return U_ref * f_wind;
}

// ============================================================================
// Compute effective wind speed field with fuel sheltering
//
// Parameters:
//   U_eff_mf       – output: effective wind speed [m/s]
//   U_ref_mf       – input: reference wind speed [m/s]
//   fuel_depth_mf  – input: fuel bed depth [m]
//   packing_mf     – input: packing ratio [dimensionless]
//   k_LAI          – LAI estimation coefficient
//   z_h            – height fraction z/h
//   a              – extinction coefficient
// ============================================================================
inline void compute_sheltered_wind_field(MultiFab& U_eff_mf,
                                        const MultiFab& U_ref_mf,
                                        const MultiFab& fuel_depth_mf,
                                        const MultiFab& packing_mf,
                                        Real k_LAI = WindFuelConst::K_LAI_DEFAULT,
                                        Real z_h   = WindFuelConst::Z_FRACTION_DEFAULT,
                                        Real a     = WindFuelConst::EXTINCTION_COEFF_DEFAULT)
{
    for (MFIter mfi(U_eff_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto U_eff     = U_eff_mf.array(mfi);
        auto const U_ref   = U_ref_mf.const_array(mfi);
        auto const depth   = fuel_depth_mf.const_array(mfi);
        auto const pack    = packing_mf.const_array(mfi);
        const Real k_L     = k_LAI;
        const Real z_h_loc = z_h;
        const Real a_loc   = a;

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            U_eff(i, j, k) = compute_sheltered_wind_speed(
                U_ref(i, j, k),
                depth(i, j, k),
                pack(i, j, k),
                k_L, z_h_loc, a_loc
            );
        });
    }
}

// ============================================================================
// Compute LAI field from fuel properties
//
// Parameters:
//   LAI_mf        – output: leaf area index [dimensionless]
//   fuel_depth_mf – input: fuel bed depth [m]
//   packing_mf    – input: packing ratio [dimensionless]
//   k_LAI         – LAI estimation coefficient
// ============================================================================
inline void compute_LAI_field(MultiFab& LAI_mf,
                              const MultiFab& fuel_depth_mf,
                              const MultiFab& packing_mf,
                              Real k_LAI = WindFuelConst::K_LAI_DEFAULT)
{
    for (MFIter mfi(LAI_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto LAI       = LAI_mf.array(mfi);
        auto const depth = fuel_depth_mf.const_array(mfi);
        auto const pack  = packing_mf.const_array(mfi);
        const Real k_L   = k_LAI;

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            LAI(i, j, k) = estimate_LAI(depth(i, j, k), pack(i, j, k), k_L);
        });
    }
}

// ============================================================================
// Apply wind reduction to existing wind field (in-place)
//
// Modifies the wind field by applying the fuel sheltering effect.
//
// Parameters:
//   U_mf          – input/output: wind speed field [m/s]
//   fuel_depth_mf – input: fuel bed depth [m]
//   packing_mf    – input: packing ratio [dimensionless]
//   k_LAI         – LAI estimation coefficient
//   z_h           – height fraction z/h
//   a             – extinction coefficient
// ============================================================================
inline void apply_wind_fuel_reduction(MultiFab& U_mf,
                                      const MultiFab& fuel_depth_mf,
                                      const MultiFab& packing_mf,
                                      Real k_LAI = WindFuelConst::K_LAI_DEFAULT,
                                      Real z_h   = WindFuelConst::Z_FRACTION_DEFAULT,
                                      Real a     = WindFuelConst::EXTINCTION_COEFF_DEFAULT)
{
    for (MFIter mfi(U_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto U         = U_mf.array(mfi);
        auto const depth = fuel_depth_mf.const_array(mfi);
        auto const pack  = packing_mf.const_array(mfi);
        const Real k_L   = k_LAI;
        const Real z_h_loc = z_h;
        const Real a_loc   = a;

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const Real LAI = estimate_LAI(depth(i, j, k), pack(i, j, k), k_L);
            const Real f_wind = compute_wind_reduction_factor(LAI, z_h_loc, a_loc);
            U(i, j, k) *= f_wind;
        });
    }
}

#endif // WIND_FUEL_INTERACTION_H
