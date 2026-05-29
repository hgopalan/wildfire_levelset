// ============================================================================
// Multi-Layer Canopy Wind Profile
// ============================================================================
// Implements vertical wind profile through a multi-layer canopy following
// exponential attenuation models. This provides a more physically realistic
// representation of wind speed variation with height in forest environments
// compared to single-layer models.
//
// Physics:
//   - Wind speed decreases exponentially from canopy top to surface
//   - Attenuation depends on canopy density (leaf area index, LAI)
//   - Different layers (understory, mid-canopy, crown) have different profiles
//   - Above canopy: logarithmic boundary layer profile
//   - Within canopy: exponential decay (Massman 1997, Inoue 1963)
//
// Model (Massman 1997 / Inoue 1963):
//   u(z) / u_h = exp(-α × (1 - z/h))    for z < h  (within canopy)
//   u(z) / u_h = (ln((z-d)/z₀) / ln((h-d)/z₀))  for z > h  (above canopy)
//
// where:
//   u(z)  - wind speed at height z [m/s]
//   u_h   - wind speed at canopy top h [m/s]
//   h     - canopy height [m]
//   α     - attenuation coefficient (typically 2-4, depends on LAI)
//   d     - displacement height [m] (typically 0.7h)
//   z₀    - roughness length [m] (typically 0.1h)
//
// The attenuation coefficient α is related to leaf area index (LAI):
//   α ≈ 2.5 × LAI^0.5   (Massman 1997)
//
// Multi-layer extension:
//   - Divide canopy into N vertical layers
//   - Compute wind speed at each layer midpoint
//   - Use layer-specific wind for fire spread calculations
//   - Improves accuracy for surface vs. crown fire transitions
//
// References:
//   Inoue, E. (1963). "On the turbulent structure of airflow within crop
//     canopies." Journal of the Meteorological Society of Japan, 41, 317-326.
//   Massman, W.J. (1997). "An analytical one-dimensional model of momentum
//     transfer by vegetation of arbitrary structure." Boundary-Layer
//     Meteorology, 83(3), 407-421.
//   Cionco, R.M. (1965). "A mathematical model for air flow in a vegetative
//     canopy." Journal of Applied Meteorology, 4(4), 517-522.
//   Yi, C. (2008). "Momentum transfer within canopies." Journal of Applied
//     Meteorology and Climatology, 47, 262-275.
// ============================================================================

#ifndef MULTI_LAYER_CANOPY_WIND_H
#define MULTI_LAYER_CANOPY_WIND_H

#include "constants.H"
#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_GpuLaunch.H>
#include <cmath>

using namespace amrex;

// Default parameters for multi-layer canopy wind profile
namespace CanopyWindConst {
    // Attenuation coefficient [dimensionless]
    // Typical range: 1.5-4.0 depending on canopy density
    constexpr Real ALPHA_DEFAULT       = 2.5;   // Moderate canopy density
    constexpr Real ALPHA_SPARSE        = 1.5;   // Sparse canopy (LAI < 2)
    constexpr Real ALPHA_DENSE         = 4.0;   // Dense canopy (LAI > 6)
    
    // Displacement height ratio [dimensionless]
    // d = d_ratio × h, typically 0.6-0.8 for forests
    constexpr Real D_RATIO_DEFAULT     = 0.7;
    
    // Roughness length ratio [dimensionless]
    // z₀ = z0_ratio × h, typically 0.05-0.15 for forests
    constexpr Real Z0_RATIO_DEFAULT    = 0.1;
    
    // Minimum wind speed [m/s] to avoid division by zero
    constexpr Real MIN_WIND_SPEED      = 0.01;
    
    // Default number of vertical layers
    constexpr int N_LAYERS_DEFAULT     = 5;
    
    // Minimum canopy height for valid calculation [m]
    constexpr Real MIN_CANOPY_HEIGHT   = 0.5;
}

// ============================================================================
// compute_canopy_wind_attenuation_coefficient
//
// Computes the attenuation coefficient α from leaf area index (LAI).
//
// Model (Massman 1997):
//   α = β × LAI^γ
//
// where:
//   β = 2.5 (empirical constant)
//   γ = 0.5 (typical for many forest types)
//
// Parameters:
//   LAI - leaf area index [m²/m²] (typically 0-10)
//
// Returns: attenuation coefficient α [dimensionless]
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real compute_canopy_wind_attenuation_coefficient(Real LAI)
{
    using namespace amrex;
    
    // Ensure positive LAI
    LAI = amrex::max(LAI, Real(0.0));
    
    // Massman (1997) empirical relationship
    constexpr Real beta  = Real(2.5);
    constexpr Real gamma = Real(0.5);
    
    const Real alpha = beta * std::sqrt(LAI);  // Optimized: gamma = 0.5
    
    // Clamp to physically reasonable range
    return amrex::min(amrex::max(alpha, CanopyWindConst::ALPHA_SPARSE),
                      CanopyWindConst::ALPHA_DENSE);
}

// ============================================================================
// compute_wind_speed_at_height_within_canopy
//
// Computes wind speed at height z within the canopy using exponential profile.
//
// Model (Inoue 1963 / Massman 1997):
//   u(z) = u_h × exp(-α × (1 - z/h))
//
// where:
//   u(z)  - wind speed at height z [m/s]
//   u_h   - wind speed at canopy top h [m/s]
//   α     - attenuation coefficient [dimensionless]
//   z     - height above ground [m]
//   h     - canopy height [m]
//
// Parameters:
//   u_h   - wind speed at canopy top [m/s]
//   z     - height above ground [m]
//   h     - canopy height [m]
//   alpha - attenuation coefficient [dimensionless]
//
// Returns: wind speed at height z [m/s]
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real compute_wind_speed_at_height_within_canopy(
    Real u_h,
    Real z,
    Real h,
    Real alpha = CanopyWindConst::ALPHA_DEFAULT)
{
    using namespace amrex;
    
    // Ensure positive values
    u_h   = amrex::max(u_h,   CanopyWindConst::MIN_WIND_SPEED);
    z     = amrex::max(z,     Real(0.0));
    h     = amrex::max(h,     CanopyWindConst::MIN_CANOPY_HEIGHT);
    alpha = amrex::max(alpha, Real(0.1));
    
    // Clamp z to be within canopy (0 ≤ z ≤ h)
    z = amrex::min(z, h);
    
    // Normalized height
    const Real z_norm = z / h;
    
    // Exponential wind profile (Inoue 1963)
    const Real u_z = u_h * std::exp(-alpha * (Real(1.0) - z_norm));
    
    return amrex::max(u_z, CanopyWindConst::MIN_WIND_SPEED);
}

// ============================================================================
// compute_wind_speed_at_height_above_canopy
//
// Computes wind speed at height z above the canopy using logarithmic profile.
//
// Model (standard atmospheric boundary layer):
//   u(z) = u_h × [ln((z-d)/z₀) / ln((h-d)/z₀)]
//
// where:
//   u(z)  - wind speed at height z [m/s]
//   u_h   - wind speed at canopy top h [m/s]
//   d     - displacement height [m] (typically 0.7h)
//   z₀    - roughness length [m] (typically 0.1h)
//
// Parameters:
//   u_h      - wind speed at canopy top [m/s]
//   z        - height above ground [m]
//   h        - canopy height [m]
//   d_ratio  - displacement height ratio d/h [dimensionless]
//   z0_ratio - roughness length ratio z₀/h [dimensionless]
//
// Returns: wind speed at height z [m/s]
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real compute_wind_speed_at_height_above_canopy(
    Real u_h,
    Real z,
    Real h,
    Real d_ratio  = CanopyWindConst::D_RATIO_DEFAULT,
    Real z0_ratio = CanopyWindConst::Z0_RATIO_DEFAULT)
{
    using namespace amrex;
    
    // Ensure positive values
    u_h   = amrex::max(u_h, CanopyWindConst::MIN_WIND_SPEED);
    z     = amrex::max(z,   Real(0.0));
    h     = amrex::max(h,   CanopyWindConst::MIN_CANOPY_HEIGHT);
    
    // Compute displacement height and roughness length
    const Real d  = d_ratio  * h;
    const Real z0 = z0_ratio * h;
    
    // Ensure z > d for valid logarithmic profile
    z = amrex::max(z, d + z0);
    
    // Logarithmic wind profile
    const Real ln_z  = std::log((z - d) / z0);
    const Real ln_h  = std::log((h - d) / z0);
    const Real ratio = (ln_h > Real(1.0e-6)) ? (ln_z / ln_h) : Real(1.0);
    
    const Real u_z = u_h * ratio;
    
    return amrex::max(u_z, CanopyWindConst::MIN_WIND_SPEED);
}

// ============================================================================
// compute_multi_layer_wind_profile
//
// Computes wind speed at multiple heights within and above the canopy.
//
// This function divides the domain into vertical layers and computes
// the wind speed at each layer midpoint. Layers below canopy top use
// exponential profile; layers above use logarithmic profile.
//
// Parameters:
//   u_ref     - reference wind speed (at canopy top or above) [m/s]
//   h         - canopy height [m]
//   z_ref     - reference height for u_ref [m]
//   n_layers  - number of vertical layers
//   layer_heights - output: heights of layer midpoints [m] (length n_layers)
//   layer_winds   - output: wind speeds at layer midpoints [m/s] (length n_layers)
//   alpha     - attenuation coefficient [dimensionless]
//   d_ratio   - displacement height ratio [dimensionless]
//   z0_ratio  - roughness length ratio [dimensionless]
//
// Note: layer_heights and layer_winds must be pre-allocated arrays of length n_layers
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void compute_multi_layer_wind_profile(
    Real u_ref,
    Real h,
    Real z_ref,
    int  n_layers,
    Real* layer_heights,
    Real* layer_winds,
    Real alpha    = CanopyWindConst::ALPHA_DEFAULT,
    Real d_ratio  = CanopyWindConst::D_RATIO_DEFAULT,
    Real z0_ratio = CanopyWindConst::Z0_RATIO_DEFAULT)
{
    using namespace amrex;
    
    // Ensure positive values
    u_ref = amrex::max(u_ref, CanopyWindConst::MIN_WIND_SPEED);
    h     = amrex::max(h,     CanopyWindConst::MIN_CANOPY_HEIGHT);
    z_ref = amrex::max(z_ref, Real(0.0));
    
    // Compute wind speed at canopy top if reference is at different height
    Real u_h = u_ref;
    if (z_ref > h + Real(0.01)) {
        // Reference is above canopy, adjust to canopy top using log profile
        const Real d  = d_ratio  * h;
        const Real z0 = z0_ratio * h;
        const Real ln_zref = std::log((z_ref - d) / z0);
        const Real ln_h    = std::log((h - d) / z0);
        u_h = u_ref * (ln_h / amrex::max(ln_zref, Real(1.0e-6)));
    } else if (z_ref < h - Real(0.01)) {
        // Reference is within canopy, adjust to canopy top using exp profile
        const Real z_norm = z_ref / h;
        u_h = u_ref / std::exp(-alpha * (Real(1.0) - z_norm));
    }
    
    // Compute layer heights and wind speeds
    const Real z_max = h * Real(2.0);  // Extend to 2× canopy height
    const Real dz = z_max / Real(n_layers);
    
    for (int i = 0; i < n_layers; ++i) {
        // Layer midpoint height
        const Real z_mid = (Real(i) + Real(0.5)) * dz;
        layer_heights[i] = z_mid;
        
        // Compute wind speed at this height
        if (z_mid <= h) {
            // Within canopy: exponential profile
            layer_winds[i] = compute_wind_speed_at_height_within_canopy(
                u_h, z_mid, h, alpha);
        } else {
            // Above canopy: logarithmic profile
            layer_winds[i] = compute_wind_speed_at_height_above_canopy(
                u_h, z_mid, h, d_ratio, z0_ratio);
        }
    }
}

// ============================================================================
// compute_effective_wind_speed_for_layer
//
// Computes effective wind speed for fire spread at a given vertical layer.
//
// This function accounts for:
//   - Surface fires: use wind at fuel bed midheight
//   - Crown fires: use wind at crown midpoint (between CBH and CH)
//
// Parameters:
//   u_h       - wind speed at canopy top [m/s]
//   h         - canopy height [m]
//   CBH       - canopy base height [m]
//   layer_type - 0 = surface layer, 1 = crown layer
//   fuel_depth - fuel bed depth for surface fires [m]
//   alpha     - attenuation coefficient [dimensionless]
//
// Returns: effective wind speed for fire spread [m/s]
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real compute_effective_wind_speed_for_layer(
    Real u_h,
    Real h,
    Real CBH,
    int  layer_type,
    Real fuel_depth = Real(0.6),  // Default fuel depth [m] (2 ft)
    Real alpha      = CanopyWindConst::ALPHA_DEFAULT)
{
    using namespace amrex;
    
    Real z_effective;
    
    if (layer_type == 0) {
        // Surface fire: wind at fuel bed midheight
        z_effective = fuel_depth * Real(0.5);
    } else {
        // Crown fire: wind at crown layer midpoint
        z_effective = CBH + (h - CBH) * Real(0.5);
    }
    
    // Compute wind speed at effective height
    return compute_wind_speed_at_height_within_canopy(u_h, z_effective, h, alpha);
}

#endif // MULTI_LAYER_CANOPY_WIND_H
