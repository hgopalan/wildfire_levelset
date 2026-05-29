// ============================================================================
// Flame Tilt Angle Diagnostic
// ============================================================================
// Computes the flame tilt angle from the horizontal due to wind and buoyancy.
// This diagnostic is useful for understanding flame geometry, radiation heat
// transfer, and fire behavior classification.
//
// Flame tilt angle (degrees from vertical):
//   θ_f = arctan(U / v_buoyancy)
//
// where:
//   U           - horizontal wind speed [m/s]
//   v_buoyancy  - buoyancy-driven vertical velocity [m/s]
//
// The buoyancy velocity is estimated from the Byram fireline intensity:
//   v_b = sqrt((g × I_B) / (ρ_air × C_p × T_a))
//
// where:
//   g       - gravitational acceleration [m/s²]
//   I_B     - Byram fireline intensity [kW/m]
//   ρ_air   - air density [kg/m³]
//   C_p     - specific heat of air [J/(kg·K)]
//   T_a     - ambient temperature [K]
//
// Output:
//   Flame tilt angle in degrees from vertical (0° = vertical, 90° = horizontal)
//
// References:
//   Byram, G.M. (1959). "Combustion of forest fuels."
//   Anderson, H.E. (1968). "Sundance fire: An analysis of fire phenomena."
//     USDA Forest Service Research Paper INT-56.
//   Weise, D.R. & Biging, G.S. (1996). "Effects of wind velocity and slope
//     on flame properties." Canadian Journal of Forest Research, 26(10):1849–1858.
// ============================================================================

#ifndef FLAME_TILT_ANGLE_H
#define FLAME_TILT_ANGLE_H

#include "constants.H"
#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_GpuLaunch.H>

using namespace amrex;

// ============================================================================
// Compute flame tilt angle from fireline intensity and wind velocity
//
// Parameters:
//   flame_tilt_mf         – output: flame tilt angle [degrees from vertical]
//   fireline_intensity_mf – input:  Byram fireline intensity [kW/m]
//   vel                   – input:  wind velocity field (3 components) [m/s]
//   rho_air               – air density [kg/m³] (default: 1.2)
//   C_p                   – specific heat of air [J/(kg·K)] (default: 1005)
//   T_a                   – ambient temperature [K] (default: 300)
// ============================================================================
inline void compute_flame_tilt_angle(
    MultiFab& flame_tilt_mf,
    const MultiFab& fireline_intensity_mf,
    const MultiFab& vel,
    Real rho_air = Real(1.2),
    Real C_p     = Real(1005.0),
    Real T_a     = Real(300.0))
{
    const Real g = WildfireConst::GRAVITY;  // 9.81 m/s²

    // Convert kW to W: 1 kW = 1000 W = 1000 J/s
    const Real kW_to_W = Real(1000.0);

    // Buoyancy velocity denominator: ρ_air × C_p × T_a
    const Real buoy_denom = rho_air * C_p * T_a;

    // Radians to degrees conversion
    const Real rad_to_deg = Real(180.0) / WildfireConst::PI;

    for (MFIter mfi(flame_tilt_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const I_B = fireline_intensity_mf.const_array(mfi);
        auto const u   = vel.const_array(mfi);
        auto tilt      = flame_tilt_mf.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            // Fireline intensity [kW/m] → [W/m]
            const Real I_B_W = I_B(i, j, k) * kW_to_W;

            // Wind speed magnitude [m/s]
            const Real u_x = u(i, j, k, 0);
            const Real u_y = u(i, j, k, 1);
            const Real U_mag = std::sqrt(u_x * u_x + u_y * u_y);

            // Buoyancy velocity [m/s]: v_b = sqrt((g × I_B) / (ρ_air × C_p × T_a))
            // Use a minimum buoyancy velocity to prevent division by zero
            const Real v_b_sq = (I_B_W > Real(0.0))
                                ? (g * I_B_W) / buoy_denom
                                : Real(0.0);
            const Real v_b = std::sqrt(v_b_sq);
            const Real v_b_safe = amrex::max(v_b, Real(0.1));  // minimum 0.1 m/s

            // Flame tilt angle: θ = arctan(U / v_b)
            const Real tilt_rad = std::atan(U_mag / v_b_safe);
            tilt(i, j, k) = tilt_rad * rad_to_deg;  // convert to degrees
        });
    }
}

#endif // FLAME_TILT_ANGLE_H
