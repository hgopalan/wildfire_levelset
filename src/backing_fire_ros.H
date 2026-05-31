// ============================================================================
// Backing Fire Rate of Spread Adjustment
// ============================================================================
// Provides empirical backing fire rate of spread (ROS) ratios for use in
// elliptical fire spread models. The backing fire ROS is significantly lower
// than the heading fire ROS, typically 10-40% of the heading rate.
//
// This module documents the backing fire coefficients used in the FARSITE
// ellipse model and provides helper functions for computing backing fire
// ratios from environmental conditions.
//
// The backing fire is the portion of the fire perimeter spreading against
// the wind direction. It is characterized by:
//   - Lower ROS (typically 10-25% of heading ROS for surface fires)
//   - Lower intensity and shorter flames
//   - More complete fuel consumption
//   - Reduced spotting potential
//
// Backing fire coefficients in Richards (1990) ellipse model:
//   - a: head fire coefficient (maximum spread, downwind)
//   - b: flank fire coefficient (crosswind spread)
//   - c: backing fire coefficient (minimum spread, upwind)
//
// Typical ratios (Anderson 1983, Rothermel 1972):
//   c/a ≈ 0.10-0.40 depending on wind speed and fuel type
//
// References:
//   Rothermel, R.C. (1972). "A mathematical model for predicting fire spread
//     in wildland fuels." USDA Forest Service Research Paper INT-115.
//   Anderson, H.E. (1983). "Predicting wind-driven wildland fire size and shape."
//     USDA Forest Service Research Paper INT-305.
//   Richards, G.D. (1990). "An elliptical growth model of forest fire fronts
//     and its numerical solution." International Journal of Numerical Methods
//     in Engineering, 30(6):1163-1179.
//   Finney, M.A. (1998). "FARSITE: Fire Area Simulator." USDA Forest Service
//     Research Paper RMRS-RP-4.
// ============================================================================

#ifndef BACKING_FIRE_ROS_H
#define BACKING_FIRE_ROS_H

#include <AMReX.H>
#include <cmath>

using namespace amrex;

// ============================================================================
// Empirical backing fire ratios (c/a) from literature
// ============================================================================
namespace BackingFireConst {
    // Default backing fire ratio for moderate wind conditions
    constexpr Real DEFAULT_BACKING_RATIO = 0.20;  // 20% of heading ROS
    
    // Conservative (lower) backing fire ratio for high winds
    constexpr Real LOW_BACKING_RATIO = 0.10;      // 10% of heading ROS
    
    // Liberal (higher) backing fire ratio for light winds
    constexpr Real HIGH_BACKING_RATIO = 0.40;     // 40% of heading ROS
    
    // Wind speed thresholds for backing fire ratio selection [mph]
    constexpr Real WIND_THRESHOLD_LOW = 5.0;      // < 5 mph: light wind
    constexpr Real WIND_THRESHOLD_HIGH = 15.0;    // > 15 mph: strong wind
}

// ============================================================================
// compute_backing_fire_ratio
//
// Computes the backing fire ROS ratio (c/a) based on wind speed.
// Uses empirical relationships from Anderson (1983) and Rothermel (1972).
//
// The backing fire ratio decreases with increasing wind speed:
//   - Light winds (<5 mph):   c/a ≈ 0.30-0.40 (backing fire ~30-40% of heading)
//   - Moderate winds (5-15 mph): c/a ≈ 0.15-0.25 (backing fire ~15-25% of heading)
//   - High winds (>15 mph):   c/a ≈ 0.10-0.15 (backing fire ~10-15% of heading)
//
// Parameters:
//   wind_speed_mph - wind speed [mph]
//   use_conservative - if true, use lower backing ratios (default: false)
//
// Returns: backing fire ratio c/a (dimensionless, 0.1-0.4)
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real compute_backing_fire_ratio(Real wind_speed_mph, bool use_conservative = false)
{
    using namespace BackingFireConst;
    
    // Base backing ratio interpolated from wind speed
    Real backing_ratio;
    
    if (wind_speed_mph < WIND_THRESHOLD_LOW) {
        // Light wind: higher backing fire ratio
        backing_ratio = HIGH_BACKING_RATIO;
    } else if (wind_speed_mph > WIND_THRESHOLD_HIGH) {
        // Strong wind: lower backing fire ratio
        backing_ratio = LOW_BACKING_RATIO;
    } else {
        // Moderate wind: interpolate between high and low
        Real wind_factor = (wind_speed_mph - WIND_THRESHOLD_LOW) / 
                          (WIND_THRESHOLD_HIGH - WIND_THRESHOLD_LOW);
        backing_ratio = HIGH_BACKING_RATIO - wind_factor * 
                       (HIGH_BACKING_RATIO - LOW_BACKING_RATIO);
    }
    
    // Apply conservative adjustment if requested
    if (use_conservative) {
        backing_ratio *= Real(0.75);  // Reduce by 25% for safety
    }
    
    // Clamp to valid range [0.05, 0.50]
    backing_ratio = amrex::max(Real(0.05), amrex::min(backing_ratio, Real(0.50)));
    
    return backing_ratio;
}

// ============================================================================
// compute_backing_ros_from_heading
//
// Computes the backing fire ROS directly from heading ROS and wind speed.
// Convenience function that combines heading ROS with empirical backing ratio.
//
// Parameters:
//   R_heading - heading fire ROS [m/s or any consistent units]
//   wind_speed_mph - wind speed [mph]
//   use_conservative - if true, use lower backing ratios (default: false)
//
// Returns: backing fire ROS in same units as R_heading
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real compute_backing_ros_from_heading(Real R_heading, Real wind_speed_mph, 
                                      bool use_conservative = false)
{
    Real backing_ratio = compute_backing_fire_ratio(wind_speed_mph, use_conservative);
    return R_heading * backing_ratio;
}

// ============================================================================
// Ellipse coefficient adjustment for backing fire
//
// Adjusts Richards (1990) ellipse coefficient c based on wind conditions
// to match empirical backing fire observations.
//
// The FARSITE ellipse model computes spread as:
//   - Heading direction (cos θ > 0): R = (a cos θ + b sin θ) × R_head
//   - Backing direction (cos θ < 0): R = (b sin θ - c cos θ) × R_head
//
// This function ensures c is consistent with empirical backing fire ratios.
//
// Parameters:
//   coeff_a - head fire coefficient (from L/W ratio or user input)
//   wind_speed_mph - wind speed [mph]
//   use_conservative - if true, use lower backing ratios (default: false)
//
// Returns: adjusted backing fire coefficient c
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real adjust_backing_coefficient(Real coeff_a, Real wind_speed_mph,
                                bool use_conservative = false)
{
    Real backing_ratio = compute_backing_fire_ratio(wind_speed_mph, use_conservative);
    
    // Backing coefficient c should produce backing_ratio × R_head when θ = 180°
    // At θ = 180°: cos θ = -1, sin θ = 0
    //   R_back = (-c × (-1)) × R_head = c × R_head
    // Therefore: c = backing_ratio × a
    Real coeff_c = backing_ratio * coeff_a;
    
    return coeff_c;
}

// ============================================================================
// Document empirical backing fire observations
//
// This namespace collects empirical backing fire ROS observations from
// field studies and experimental fires for reference and validation.
// ============================================================================
namespace BackingFireObservations {
    // Anderson (1983) Table 4 - Backing fire ROS as fraction of heading ROS
    // Wind 5-10 mph, various fuel models:
    //   Fuel Model 1 (grass): backing/heading ≈ 0.25
    //   Fuel Model 2 (grass-shrub): backing/heading ≈ 0.20
    //   Fuel Model 4 (chaparral): backing/heading ≈ 0.15
    //   Fuel Model 10 (timber-litter): backing/heading ≈ 0.20
    
    // Rothermel (1972) - Backing fire observations:
    //   Light wind (< 5 mph): backing ≈ 30-40% of heading
    //   Moderate wind (5-15 mph): backing ≈ 15-25% of heading
    //   Strong wind (> 15 mph): backing ≈ 10-15% of heading
    
    // FARSITE default (no wind): backing/heading ≈ 1.0 (circular spread)
    // FARSITE with wind: backing coefficient c computed from L/W ratio
    //   L/W = 1 (no wind): c/a = 1.0
    //   L/W = 2: c/a ≈ 0.33
    //   L/W = 3: c/a ≈ 0.20
    //   L/W = 5: c/a ≈ 0.12
    //   L/W = 10: c/a ≈ 0.06
}

#endif // BACKING_FIRE_ROS_H
