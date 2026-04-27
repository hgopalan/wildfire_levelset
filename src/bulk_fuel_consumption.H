// ============================================================================
// Bulk Fuel Consumption Fraction Model
// ============================================================================
// This implements a bulk fuel consumption fraction model for FARSITE-based
// fire spread simulations. The model computes the fraction of available fuel
// that is consumed by the fire based on fire intensity and residence time.
//
// The model is based on:
// - Albini (1976) - "Estimating wildfire behavior and effects"
// - Byram (1959) - Fire intensity and residence time relationships
// - FARSITE documentation - Fuel consumption modeling
//
// Physical Background:
// --------------------
// As a fire spreads through fuel, not all fuel is consumed. The fraction
// consumed depends on:
// 1. Fire intensity (higher intensity → more complete combustion)
// 2. Residence time (longer residence → more complete combustion)
// 3. Fuel moisture (lower moisture → more complete combustion)
// 4. Fuel type and packing
//
// The bulk consumption fraction (f_c) ranges from:
// - f_c_min (typically 0.5): Low-intensity, fast-moving fires
// - f_c_max (typically 0.9): High-intensity, slow-moving fires
//
// This model uses a simple parametric relationship based on the
// characteristic residence time and fire intensity.
// ============================================================================

#ifndef BULK_FUEL_CONSUMPTION_H
#define BULK_FUEL_CONSUMPTION_H

#include <AMReX.H>
#include <cmath>
#include "parse_inputs.H"

// ============================================================================
// Compute bulk fuel consumption fraction based on fire intensity and
// residence time
// ============================================================================
// The model uses a sigmoid-like function to transition between minimum and
// maximum consumption fractions based on the normalized fire intensity.
//
// Parameters:
//   I_R: Reaction intensity [BTU/ft²/min] from Rothermel model
//   tau_residence: Residence time [seconds] - characteristic time fire
//                  remains at a given location
//   f_c_min: Minimum fuel consumption fraction (fast/low-intensity fires)
//   f_c_max: Maximum fuel consumption fraction (slow/high-intensity fires)
//
// Returns:
//   f_c: Bulk fuel consumption fraction (dimensionless, 0.0-1.0)
//
// Physical Interpretation:
//   - Low intensity (I_R < 500 BTU/ft²/min): f_c ≈ f_c_min
//   - Moderate intensity (500-2000 BTU/ft²/min): f_c transitions
//   - High intensity (I_R > 2000 BTU/ft²/min): f_c ≈ f_c_max
//
// The transition is modulated by residence time:
//   - Short residence time: favors lower consumption
//   - Long residence time: favors higher consumption
// ============================================================================
static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_fuel_consumption_fraction(
    amrex::Real I_R,
    amrex::Real tau_residence,
    amrex::Real f_c_min,
    amrex::Real f_c_max)
{
    using namespace amrex;
    
    // Reference intensity for normalization (BTU/ft²/min)
    // This represents a characteristic intensity where consumption
    // is halfway between minimum and maximum
    const Real I_ref = Real(1000.0);
    
    // Reference residence time (seconds)
    // Typical residence times range from 30s (fast grass fires)
    // to 120s (slow-moving crown fires)
    const Real tau_ref = Real(60.0);
    
    // Ensure positive values
    I_R = amrex::max(I_R, Real(0.0));
    tau_residence = amrex::max(tau_residence, Real(1.0));
    
    // Normalized intensity with residence time weighting
    // Higher residence time increases effective intensity
    Real I_normalized = (I_R / I_ref) * std::sqrt(tau_residence / tau_ref);
    
    // Sigmoid-like transition function
    // Uses tanh for smooth transition between f_c_min and f_c_max
    // The factor 0.5 controls the steepness of the transition
    Real transition = Real(0.5) * (Real(1.0) + std::tanh(I_normalized - Real(1.0)));
    
    // Compute consumption fraction
    Real f_c = f_c_min + (f_c_max - f_c_min) * transition;
    
    // Ensure physically reasonable bounds
    f_c = amrex::max(f_c, f_c_min);
    f_c = amrex::min(f_c, f_c_max);
    
    return f_c;
}

// ============================================================================
// Compute residual fuel fraction (complement of consumption fraction)
// ============================================================================
// This is simply 1 - f_c, representing the fraction of fuel that remains
// unconsumed after fire passage.
//
// Parameters:
//   f_c: Fuel consumption fraction (0.0-1.0)
//
// Returns:
//   f_r: Residual fuel fraction (0.0-1.0)
//
// Applications:
//   - Used to update fuel loading for subsequent fire passages
//   - Important for computing post-fire fuel conditions
//   - Affects fire intensity on subsequent burns
// ============================================================================
static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_residual_fuel_fraction(amrex::Real f_c)
{
    using namespace amrex;
    
    // Ensure f_c is in valid range
    f_c = amrex::max(Real(0.0), amrex::min(f_c, Real(1.0)));
    
    return Real(1.0) - f_c;
}

// ============================================================================
// Compute modified fire intensity accounting for fuel consumption
// ============================================================================
// When fuel is not completely consumed, the effective intensity is reduced.
// This can be important for coupled models where intensity affects other
// processes (e.g., firebrand generation, plume rise).
//
// Parameters:
//   I_R: Base reaction intensity [BTU/ft²/min]
//   f_c: Fuel consumption fraction (0.0-1.0)
//
// Returns:
//   I_effective: Effective fire intensity [BTU/ft²/min]
//
// Note: This is a simple linear scaling. More sophisticated models
// might use non-linear relationships.
// ============================================================================
static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_effective_intensity(amrex::Real I_R, amrex::Real f_c)
{
    using namespace amrex;
    
    // Ensure valid inputs
    I_R = amrex::max(I_R, Real(0.0));
    f_c = amrex::max(Real(0.0), amrex::min(f_c, Real(1.0)));
    
    // Linear scaling: effective intensity proportional to consumed fuel
    return I_R * f_c;
}

#endif // BULK_FUEL_CONSUMPTION_H
