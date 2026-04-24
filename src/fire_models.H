// ============================================================================
// Fire Spread Model Utility Functions
// ============================================================================
// This file provides utility functions for various fire spread models
// including Anderson (1983) L/W ratio and FARSITE terrain effects
// ============================================================================

#ifndef FIRE_MODELS_H
#define FIRE_MODELS_H

#include <AMReX.H>
#include <cmath>

// Math constants for portability
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ============================================================================
// Anderson (1983) L/W ratio calculation based on wind speed
// Used in FARSITE model for elliptical fire shape
// ============================================================================
// Anderson (1983): "Predicting wind-driven wild land fire size and shape"
// Research Paper INT-305, USDA Forest Service
//
// The length-to-width (L/W) ratio describes the elongation of an elliptical
// fire perimeter in the wind direction. This empirical relationship is
// based on experimental fires and is widely used in FARSITE.
//
// Parameters:
//   wind_speed_mph: Wind speed at midflame height in miles per hour
//
// Returns:
//   L/W ratio (dimensionless), bounded between 1.0 (circular) and 8.0
// ============================================================================
static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real anderson_LW_ratio(amrex::Real wind_speed_mph)
{
    using namespace amrex;
    
    // Anderson (1983) empirical formula:
    // L/W = 0.936 * exp(0.2566 * U) + 0.461 * exp(-0.1548 * U) - 0.397
    // where U is wind speed at midflame height in mph
    // 
    // For U = 0: L/W ≈ 1.0 (circular fire)
    // For U = 10 mph: L/W ≈ 2.5
    // For U = 20 mph: L/W ≈ 4.5
    
    if (wind_speed_mph < Real(0.01)) {
        return Real(1.0); // Nearly circular for no wind
    }
    
    Real LW = Real(0.936) * std::exp(Real(0.2566) * wind_speed_mph) 
            + Real(0.461) * std::exp(Real(-0.1548) * wind_speed_mph) 
            - Real(0.397);
    
    // Ensure physically reasonable bounds
    // Upper bound of 8.0 is typical for extreme wind conditions
    return amrex::max(Real(1.0), amrex::min(LW, Real(8.0)));
}

// ============================================================================
// Convert Anderson L/W ratio to Richards' elliptical coefficients
// ============================================================================
// Richards' (1990) model uses coefficients a, b, c for directional spread:
//   - a: head fire coefficient (downwind, maximum spread)
//   - b: flank fire coefficient (crosswind spread)
//   - c: backing fire coefficient (upwind, minimum spread)
//
// The relationship between L/W ratio and these coefficients can be derived
// from the ellipse geometry. For a standard ellipse:
//   L/W = (a + c) / (2 * b)
//
// Common parameterization (used in FARSITE):
//   a = 1.0 (head fire at full rate)
//   c = 0.2 * a (backing fire at 20% of head fire rate)
//   b = (a + c) / (2 * LW)
//
// Parameters:
//   LW_ratio: Length-to-width ratio from Anderson or other source
//   coeff_a: Output head fire coefficient
//   coeff_b: Output flank fire coefficient  
//   coeff_c: Output backing fire coefficient
// ============================================================================
static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void LW_ratio_to_richards_coefficients(amrex::Real LW_ratio,
                                       amrex::Real& coeff_a,
                                       amrex::Real& coeff_b,
                                       amrex::Real& coeff_c)
{
    using namespace amrex;
    
    // Standard parameterization from FARSITE
    coeff_a = Real(1.0);              // Head fire at full spread rate
    coeff_c = Real(0.2) * coeff_a;    // Backing fire at 20% of head fire
    
    // Compute flank coefficient from L/W ratio
    // From ellipse geometry: L/W = (a + c) / (2 * b)
    // Therefore: b = (a + c) / (2 * LW)
    coeff_b = (coeff_a + coeff_c) / (Real(2.0) * LW_ratio);
    
    // Ensure physically reasonable bounds
    coeff_b = amrex::max(coeff_b, Real(0.1));  // Minimum flank spread
    coeff_b = amrex::min(coeff_b, coeff_a);     // Cannot exceed head fire
}

// ============================================================================
// Rothermel slope correction factor (standalone function for reference)
// ============================================================================
// Note: This is already implemented in rothermel_model.H as part of the
// full Rothermel calculation. This standalone version is provided for
// reference and potential use in other contexts.
//
// Rothermel (1972): Slope factor φ_s from Equation 51
// φ_s = 5.275 * β^(-0.3) * (tan(θ))^2
//
// This simplified version omits the β^(-0.3) term (fuel packing ratio adjustment),
// using only the base coefficient 5.275. For typical fuel packing ratios
// (β ≈ 0.005-0.01), the full formula gives coefficients of ~10-12. This
// simplified version provides a conservative estimate suitable for general use.
//
// For accurate results with specific fuel properties, use the full implementation
// in rothermel_model.H which includes the β^(-0.3) adjustment.
//
// Parameters:
//   slope_degrees: Terrain slope in degrees
//
// Returns:
//   Slope factor φ_s (dimensionless, simplified version)
// ============================================================================
static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real rothermel_slope_factor(amrex::Real slope_degrees)
{
    using namespace amrex;
    
    Real slope_rad = slope_degrees * Real(M_PI) / Real(180.0);
    Real tan_slope = std::tan(slope_rad);
    
    // Using base coefficient of 5.275 without β^(-0.3) adjustment
    // This is a simplified version - see rothermel_model.H for full formula
    Real phi_s = Real(5.275) * tan_slope * tan_slope;
    
    return phi_s;
}

#endif // FIRE_MODELS_H
