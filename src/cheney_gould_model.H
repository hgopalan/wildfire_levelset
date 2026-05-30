// ============================================================================
// Cheney & Gould (1995) / Cheney, Gould & Catchpole (1998) Grassland Fire
// Spread Model
// ============================================================================
// Implements the empirical rate-of-spread model for grassland fuels introduced
// in:
//
//   Cheney, N.P. & Gould, J.S. (1995). "Fire growth in grassland fuels."
//   Int. J. Wildland Fire, 5(4), 237–247.
//
// with the equilibrium spread rate equations from:
//
//   Cheney, N.P., Gould, J.S. & Catchpole, W.R. (1998). "Prediction of fire
//   spread in grasslands." Int. J. Wildland Fire, 8(1), 1–13.
//
// The model was calibrated specifically for open grassland fuels in Australia
// and produces better agreement with observed grassland fires than Rothermel in
// this fuel type.
//
// Head-fire rate of spread (R, km/h) as a function of:
//   U₁₀  – 10-m open wind speed [km/h]
//   MC   – dead fine fuel moisture content [%]
//   CF   – degree of curing [0–1; 1 = fully cured]
//
// Piecewise-linear wind speed dependence (Cheney et al. 1998, Eq. base):
//   For U₁₀ ≤ 5 km/h:  R = (0.165 + 0.534 × U₁₀) × f_MC × CF
//   For U₁₀ > 5 km/h:  R = (-0.020 + 0.640 × U₁₀) × f_MC × CF
//
// Moisture correction factor (Cheney et al. 1998):
//   f_MC = exp(-0.108 × MC)
//
// Conversion to simulation units (m/s):
//   R [m/s] = R [km/h] / 3.6
//
// Parameters:
//   mc       – dead fine fuel moisture content [%]         default: 10.0
//   curing   – degree of curing (0–1; 1 = fully cured)    default: 1.0
// ============================================================================

#ifndef CHENEY_GOULD_MODEL_H
#define CHENEY_GOULD_MODEL_H

#include "parse_inputs.H"
#include <AMReX.H>
#include <AMReX_GpuContainers.H>
#include <cmath>

// ============================================================================
// Pre-computed Cheney–Gould parameters for GPU kernel
// ============================================================================
// These scalars encode all fuel and curing state so the GPU kernel only
// needs wind speed as input.
// ============================================================================
struct CheneyGouldComputed {
    // Coefficient groups for piecewise wind formula [km/h → km/h]:
    //   For U ≤ 5 km/h:  R_kmh = (a_lo + b_lo * U) * f_mc * curing
    //   For U >  5 km/h: R_kmh = (a_hi + b_hi * U) * f_mc * curing
    amrex::Real a_lo;     // intercept for low-wind regime (km/h)
    amrex::Real b_lo;     // slope    for low-wind regime  (-)
    amrex::Real a_hi;     // intercept for high-wind regime (km/h)
    amrex::Real b_hi;     // slope    for high-wind regime  (-)
    amrex::Real scale;    // combined f_mc * curing / 3.6  [converts km/h to m/s after moisture/curing correction]
};

// ============================================================================
// Compute Cheney–Gould model coefficients from input parameters
//
// Returns CheneyGouldComputed ready for the GPU kernel.
// ============================================================================
inline CheneyGouldComputed compute_cheney_gould_params(
    const InputParameters::CheneyGouldParams& cgp)
{
    CheneyGouldComputed cgc;

    // Piecewise empirical constants (Cheney et al. 1998)
    cgc.a_lo = 0.165;   // [km/h] intercept, U ≤ 5 km/h
    cgc.b_lo = 0.534;   // [-]    slope,     U ≤ 5 km/h
    cgc.a_hi = -0.020;  // [km/h] intercept, U > 5 km/h
    cgc.b_hi = 0.640;   // [-]    slope,     U > 5 km/h

    // Moisture correction factor: f_MC = exp(-0.108 * MC)
    const amrex::Real f_mc = std::exp(-0.108 * cgp.moisture);

    // Curing factor [0–1]: clamp to valid range
    const amrex::Real cf = std::max(amrex::Real(0.0),
                                    std::min(amrex::Real(1.0), cgp.curing));

    // Combined scaling: f_mc * CF * (km/h → m/s conversion, dividing by 3.6)
    cgc.scale = f_mc * cf / amrex::Real(3.6);

    return cgc;
}

// ============================================================================
// GPU-device inline: evaluate Cheney–Gould ROS [m/s] from wind speed [m/s]
//
// Called once per cell in the GPU kernel.
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real cheney_gould_ros(const CheneyGouldComputed& cgc,
                              amrex::Real wind_speed_ms) noexcept
{
    // Convert simulation wind speed (m/s) to 10-m open wind (km/h)
    const amrex::Real U = wind_speed_ms * amrex::Real(3.6);

    // Piecewise empirical spread formula (Cheney et al. 1998)
    amrex::Real R_kmh;
    if (U <= amrex::Real(5.0)) {
        R_kmh = cgc.a_lo + cgc.b_lo * U;
    } else {
        R_kmh = cgc.a_hi + cgc.b_hi * U;
    }

    // Apply moisture and curing correction; convert to m/s
    // Negative base ROS is physically meaningless → clamp to 0
    return amrex::max(R_kmh * cgc.scale, amrex::Real(0.0));
}

// ============================================================================
// CSIRO Grassfire Acceleration Model (Cheney & Gould 1995; Gould et al. 2007)
// ============================================================================
// Models the non-equilibrium fire growth in grasslands, where small fires
// do not immediately achieve their quasi-steady-state ROS.
//
// Formula: a(t) = 1 - exp(-t/t_accel)
//
// where:
//   a(t)     - acceleration factor [0-1]
//   t        - time since ignition [min]
//   t_accel  - acceleration time constant [min] (typically 5-15 min for grass)
//
// References:
//   Cheney, N.P. & Gould, J.S. (1995). "Fire growth in grassland fuels."
//     International Journal of Wildland Fire, 5(4):237-247.
//   Gould, J.S., et al. (2007). "Project Vesta: Fire in Dry Eucalypt Forest."
//     Ensis-CSIRO, Canberra, Australia.
//
// Parameters:
//   time_since_ignition  - elapsed time since ignition [s]
//   t_accel_seconds     - acceleration time constant [s] (default: 600 = 10 min)
//
// Returns:
//   acceleration factor [0-1] to multiply the equilibrium ROS
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real csiro_grassfire_acceleration(
    amrex::Real time_since_ignition,
    amrex::Real t_accel_seconds = amrex::Real(600.0)) noexcept
{
    using namespace amrex;
    
    // Ensure non-negative time
    time_since_ignition = amrex::max(time_since_ignition, Real(0.0));
    
    // Ensure positive acceleration time constant
    t_accel_seconds = amrex::max(t_accel_seconds, Real(1.0));
    
    // Exponential approach to equilibrium
    Real accel_factor = Real(1.0) - std::exp(-time_since_ignition / t_accel_seconds);
    
    // Clamp to [0, 1]
    accel_factor = amrex::max(Real(0.0), amrex::min(Real(1.0), accel_factor));
    
    return accel_factor;
}

// ============================================================================
// Cheney-Gould ROS with CSIRO acceleration
//
// Combines the equilibrium Cheney-Gould ROS with the CSIRO acceleration
// model for non-equilibrium fire growth.
//
// Parameters:
//   cgc                  - Cheney-Gould computed parameters
//   wind_speed_ms        - wind speed [m/s]
//   time_since_ignition  - time since ignition [s]
//   t_accel_seconds     - acceleration time constant [s]
//
// Returns:
//   ROS [m/s] with acceleration applied
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real cheney_gould_ros_with_acceleration(
    const CheneyGouldComputed& cgc,
    amrex::Real wind_speed_ms,
    amrex::Real time_since_ignition,
    amrex::Real t_accel_seconds = amrex::Real(600.0)) noexcept
{
    // Get equilibrium ROS
    amrex::Real R_eq = cheney_gould_ros(cgc, wind_speed_ms);
    
    // Apply acceleration factor
    amrex::Real accel = csiro_grassfire_acceleration(time_since_ignition, t_accel_seconds);
    
    return R_eq * accel;
}

#endif // CHENEY_GOULD_MODEL_H
