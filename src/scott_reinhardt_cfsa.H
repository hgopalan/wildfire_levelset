// ============================================================================
// Scott & Reinhardt (2001) Crown Fire Surface Area (CFSA)
// ============================================================================
// Computes the active crown fire surface area per unit ground area, which
// represents the effective burning area within the crown layer during active
// crown fire. This is a key parameter for estimating crown fire heat release
// and emissions.
//
// Physics:
//   - Active crown fires consume canopy fuels over a vertical zone
//   - The CFSA accounts for the 3-D structure of the canopy layer
//   - Higher CFSA → more surface area for combustion → higher heat release
//   - CFSA depends on canopy bulk density (CBD), canopy base height (CBH),
//     and canopy height
//
// Model (Scott & Reinhardt 2001):
//   CFSA = min(CBD × (CH - CBH) × k_sa, max_cfsa)
//
// where:
//   CBD    - canopy bulk density [kg/m³]
//   CH     - canopy height (top of live crown) [m]
//   CBH    - canopy base height (bottom of live crown) [m]
//   k_sa   - surface area coefficient [m²/kg] (empirical, ~4-8 m²/kg)
//   max_cfsa - maximum CFSA (typically 2-4 for dense canopies)
//
// The surface area coefficient k_sa relates the fuel loading (CBD × height)
// to the effective burning surface area. It accounts for foliage clumping,
// branch geometry, and crown architecture.
//
// Applications:
//   - Crown fire heat release rate calculations
//   - Emissions modeling (particulate matter from crown fuels)
//   - Radiation heat flux from active crown fires
//   - Plume rise and smoke production estimates
//
// References:
//   Scott, J.H. & Reinhardt, E.D. (2001). "Assessing Crown Fire Potential by
//     Linking Models of Surface and Crown Fire Behavior." USDA Forest Service
//     Research Paper RMRS-RP-29.
//   Reinhardt, E.D., Crookston, N.L. (tech. eds.) (2003). "The Fire and Fuels
//     Extension to the Forest Vegetation Simulator." USDA Forest Service
//     General Technical Report RMRS-GTR-116.
//   Cruz, M.G., Alexander, M.E., & Wakimoto, R.H. (2005). "Development and
//     testing of models for predicting crown fire rate of spread in conifer
//     forest stands." Canadian Journal of Forest Research, 35(7), 1626-1639.
// ============================================================================

#ifndef SCOTT_REINHARDT_CFSA_H
#define SCOTT_REINHARDT_CFSA_H

#include "constants.H"
#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_GpuLaunch.H>
#include <cmath>

using namespace amrex;

// Default parameters for crown fire surface area
namespace CFSAConst {
    // Surface area coefficient [m²/kg]
    // Typical range: 4-8 m²/kg for conifer foliage
    // Lower values: dense, clumped crowns (Douglas-fir, true firs)
    // Higher values: open, dispersed crowns (ponderosa pine, lodgepole pine)
    constexpr Real K_SA_DEFAULT      = 6.0;   // Default coefficient
    constexpr Real K_SA_MIN          = 2.0;   // Minimum (dense clumping)
    constexpr Real K_SA_MAX          = 10.0;  // Maximum (very open crowns)
    
    // Maximum crown fire surface area [dimensionless]
    // Represents the upper limit of burning surface area ratio
    // Typical range: 2-4 for dense conifer canopies
    constexpr Real MAX_CFSA_DEFAULT  = 3.0;   // Default maximum
    constexpr Real MAX_CFSA_RANGE_LOW  = 1.5;   // Lower bound (sparse canopy)
    constexpr Real MAX_CFSA_RANGE_HIGH = 5.0;   // Upper bound (very dense)
    
    // Minimum canopy depth for valid CFSA calculation [m]
    constexpr Real MIN_CANOPY_DEPTH  = 0.5;   // Avoid numerical issues
}

// ============================================================================
// compute_crown_fire_surface_area
//
// Computes the active crown fire surface area per unit ground area.
//
// Model:
//   CFSA = min(CBD × (CH - CBH) × k_sa, max_cfsa)
//
// where:
//   CBD      - canopy bulk density [kg/m³]
//   CH       - canopy height [m]
//   CBH      - canopy base height [m]
//   k_sa     - surface area coefficient [m²/kg]
//   max_cfsa - maximum CFSA [dimensionless]
//
// Physical interpretation:
//   CFSA represents the ratio of burning surface area to ground area.
//   Values typically range from 0 (no crown fire) to 3-4 (dense active crown).
//
// Parameters:
//   CBD      - canopy bulk density [kg/m³]
//   CH       - canopy height [m]
//   CBH      - canopy base height [m]
//   k_sa     - surface area coefficient [m²/kg] (default: 6.0)
//   max_cfsa - maximum CFSA [dimensionless] (default: 3.0)
//
// Returns: crown fire surface area [dimensionless, 0-max_cfsa]
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real compute_crown_fire_surface_area(
    Real CBD,
    Real CH,
    Real CBH,
    Real k_sa = CFSAConst::K_SA_DEFAULT,
    Real max_cfsa = CFSAConst::MAX_CFSA_DEFAULT)
{
    using namespace amrex;
    
    // Ensure positive values
    CBD = amrex::max(CBD, Real(0.0));
    CH  = amrex::max(CH,  Real(0.0));
    CBH = amrex::max(CBH, Real(0.0));
    
    // Canopy depth (vertical extent of live crown)
    const Real canopy_depth = CH - CBH;
    
    // Check for valid canopy depth
    if (canopy_depth < CFSAConst::MIN_CANOPY_DEPTH) {
        return Real(0.0);  // No significant crown layer
    }
    
    // Ensure k_sa is in valid range
    k_sa = amrex::min(amrex::max(k_sa, CFSAConst::K_SA_MIN), CFSAConst::K_SA_MAX);
    
    // Compute CFSA using Scott & Reinhardt (2001) model
    // CFSA = CBD [kg/m³] × depth [m] × k_sa [m²/kg] = [m²/m²] = dimensionless
    Real cfsa = CBD * canopy_depth * k_sa;
    
    // Apply maximum limit
    cfsa = amrex::min(cfsa, max_cfsa);
    
    // Ensure non-negative
    cfsa = amrex::max(cfsa, Real(0.0));
    
    return cfsa;
}

// ============================================================================
// compute_crown_heat_release_with_cfsa
//
// Computes the crown fire heat release rate using CFSA.
//
// Model:
//   Q_crown = CFSA × w_crown × h_crown / tau_crown
//
// where:
//   CFSA     - crown fire surface area [dimensionless]
//   w_crown  - crown fuel loading [kg/m²]
//   h_crown  - heat content of crown fuels [kJ/kg]
//   tau_crown - residence time of crown combustion [s]
//
// This is an alternative to simpler crown heat release models that
// don't account for canopy architecture.
//
// Parameters:
//   cfsa         - crown fire surface area [dimensionless]
//   w_crown      - crown fuel loading [kg/m²]
//   h_crown      - heat content [kJ/kg] (default: 18000 kJ/kg)
//   tau_crown    - residence time [s] (default: 30 s)
//
// Returns: heat release rate [kW/m²]
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real compute_crown_heat_release_with_cfsa(
    Real cfsa,
    Real w_crown,
    Real h_crown = Real(18000.0),    // kJ/kg (typical for conifer foliage)
    Real tau_crown = Real(30.0))     // s (typical crown combustion duration)
{
    using namespace amrex;
    
    // Ensure positive values
    cfsa     = amrex::max(cfsa,     Real(0.0));
    w_crown  = amrex::max(w_crown,  Real(0.0));
    h_crown  = amrex::max(h_crown,  Real(1.0));
    tau_crown = amrex::max(tau_crown, Real(1.0));
    
    // Heat release rate [kW/m²] = CFSA × loading [kg/m²] × heat [kJ/kg] / time [s]
    const Real Q_crown = cfsa * w_crown * h_crown / tau_crown;
    
    return Q_crown;
}

// ============================================================================
// compute_cfsa_from_cbd_depth
//
// Convenience function to compute CFSA from CBD and canopy depth directly.
//
// Parameters:
//   CBD          - canopy bulk density [kg/m³]
//   canopy_depth - vertical extent of live crown [m]
//   k_sa         - surface area coefficient [m²/kg] (default: 6.0)
//   max_cfsa     - maximum CFSA [dimensionless] (default: 3.0)
//
// Returns: crown fire surface area [dimensionless]
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real compute_cfsa_from_cbd_depth(
    Real CBD,
    Real canopy_depth,
    Real k_sa = CFSAConst::K_SA_DEFAULT,
    Real max_cfsa = CFSAConst::MAX_CFSA_DEFAULT)
{
    using namespace amrex;
    
    // Ensure positive values
    CBD = amrex::max(CBD, Real(0.0));
    canopy_depth = amrex::max(canopy_depth, Real(0.0));
    
    // Check for valid canopy depth
    if (canopy_depth < CFSAConst::MIN_CANOPY_DEPTH) {
        return Real(0.0);
    }
    
    // Ensure k_sa is in valid range
    k_sa = amrex::min(amrex::max(k_sa, CFSAConst::K_SA_MIN), CFSAConst::K_SA_MAX);
    
    // Compute CFSA
    Real cfsa = CBD * canopy_depth * k_sa;
    
    // Apply maximum limit
    cfsa = amrex::min(cfsa, max_cfsa);
    
    return amrex::max(cfsa, Real(0.0));
}

// ============================================================================
// estimate_crown_fuel_loading_from_cbd
//
// Estimates crown fuel loading [kg/m²] from canopy bulk density and depth.
//
// Model:
//   w_crown = CBD × canopy_depth
//
// This converts the 3-D canopy fuel density to a 2-D ground-projected loading.
//
// Parameters:
//   CBD          - canopy bulk density [kg/m³]
//   canopy_depth - vertical extent of live crown [m]
//
// Returns: crown fuel loading [kg/m²]
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real estimate_crown_fuel_loading_from_cbd(
    Real CBD,
    Real canopy_depth)
{
    using namespace amrex;
    
    CBD = amrex::max(CBD, Real(0.0));
    canopy_depth = amrex::max(canopy_depth, Real(0.0));
    
    return CBD * canopy_depth;
}

#endif // SCOTT_REINHARDT_CFSA_H
