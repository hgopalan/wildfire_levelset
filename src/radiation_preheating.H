// ============================================================================
// Radiation-Driven Preheating Distance Model
// ============================================================================
// Calculates the distance ahead of the fire front where fuel is preheated by
// radiant energy from flames. This affects ignition timing and spread rate,
// particularly in non-uniform fuels.
//
// Physics:
//   - Radiant energy from flames preheats fuel ahead of the fire front
//   - The preheating distance depends on flame length, flame tilt, and view factor
//   - Longer, more tilted flames preheat fuel at greater distances
//
// Model:
//   d_preheat [m] = L_f × cos(θ_tilt) × F_v
//
// where:
//   L_f      = Byram flame length [m]
//   θ_tilt   = flame tilt angle from vertical [radians]
//   F_v      = view factor (typically 0.3-0.5)
//
// The flame tilt angle is computed from slope and wind effects:
//   tan(θ_tilt) = U_eff / v_buoyancy
//
// where:
//   U_eff        = effective wind speed [m/s]
//   v_buoyancy   = buoyant rise velocity ~ sqrt(g × L_f) [m/s]
//
// References:
//   Butler, B.W., et al. (2004). A method for evaluating the effectiveness of
//     firefighter escape routes. In: Fire Safety Science—Proceedings of the
//     Seventh International Symposium, pp. 1435-1446.
//   Weber, R.O. (1991). Toward a comprehensive wildfire spread model.
//     International Journal of Wildland Fire 1(4):245-248.
//   Balbi, J.H., et al. (2009). A physical model for wildland fires.
//     Combustion and Flame 156:2217-2230.
// ============================================================================

#ifndef RADIATION_PREHEATING_H
#define RADIATION_PREHEATING_H

#include "constants.H"
#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_GpuLaunch.H>
#include <cmath>

using namespace amrex;

// ============================================================================
// Compute preheating distance field
//
// Parameters:
//   preheating_dist_mf – output: preheating distance [m]
//   flame_length_mf    – input:  Byram flame length [m]
//   vel_eff_x_mf       – input:  effective wind velocity x-component [m/s]
//   vel_eff_y_mf       – input:  effective wind velocity y-component [m/s]
//   view_factor        – radiative view factor (default: 0.4, range: 0.3-0.5)
//
// The effective wind (vel_eff) should already include wind adjustment factor,
// slope effects, and any fire-induced wind components.
// ============================================================================
inline void compute_preheating_distance(MultiFab& preheating_dist_mf,
                                        const MultiFab& flame_length_mf,
                                        const MultiFab& vel_eff_x_mf,
                                        const MultiFab& vel_eff_y_mf,
                                        Real view_factor = Real(0.4))
{
    // Gravitational acceleration [m/s²]
    constexpr Real g = WildfireConst::GRAVITY;

    for (MFIter mfi(preheating_dist_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const Lf  = flame_length_mf.const_array(mfi);
        auto const Ux  = vel_eff_x_mf.const_array(mfi);
        auto const Uy  = vel_eff_y_mf.const_array(mfi);
        auto d_pre     = preheating_dist_mf.array(mfi);
        const Real Fv  = view_factor;

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const Real L_flame = Lf(i, j, k);
            
            // Skip if no flame
            if (L_flame <= Real(0.0)) {
                d_pre(i, j, k) = Real(0.0);
                return;
            }

            // Effective wind speed magnitude [m/s]
            const Real U_mag = std::sqrt(Ux(i, j, k) * Ux(i, j, k) +
                                         Uy(i, j, k) * Uy(i, j, k));

            // Buoyant rise velocity [m/s] from Byram flame length
            // v_b ~ sqrt(g × L_f) is a characteristic buoyancy velocity scale
            const Real v_buoy = std::sqrt(g * L_flame);

            // Flame tilt angle from vertical [radians]
            // tan(θ) = U_eff / v_buoy
            // For small angles: cos(θ) ≈ 1 / sqrt(1 + tan²(θ))
            const Real tan_tilt = (v_buoy > Real(1e-6)) 
                                ? U_mag / v_buoy 
                                : Real(0.0);
            const Real cos_tilt = Real(1.0) / std::sqrt(Real(1.0) + 
                                                        tan_tilt * tan_tilt);

            // Preheating distance [m]
            // d_preheat = L_f × cos(θ_tilt) × F_v
            d_pre(i, j, k) = L_flame * cos_tilt * Fv;
        });
    }
}

// ============================================================================
// Variant: compute preheating distance with scalar wind speed
//
// This version is simpler when only the wind speed magnitude is available
// rather than separate x/y components.
//
// Parameters:
//   preheating_dist_mf – output: preheating distance [m]
//   flame_length_mf    – input:  Byram flame length [m]
//   wind_speed_mf      – input:  effective wind speed magnitude [m/s]
//   view_factor        – radiative view factor (default: 0.4)
// ============================================================================
inline void compute_preheating_distance_scalar(MultiFab& preheating_dist_mf,
                                               const MultiFab& flame_length_mf,
                                               const MultiFab& wind_speed_mf,
                                               Real view_factor = Real(0.4))
{
    constexpr Real g = WildfireConst::GRAVITY;

    for (MFIter mfi(preheating_dist_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const Lf  = flame_length_mf.const_array(mfi);
        auto const U   = wind_speed_mf.const_array(mfi);
        auto d_pre     = preheating_dist_mf.array(mfi);
        const Real Fv  = view_factor;

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const Real L_flame = Lf(i, j, k);
            
            if (L_flame <= Real(0.0)) {
                d_pre(i, j, k) = Real(0.0);
                return;
            }

            const Real U_mag   = U(i, j, k);
            const Real v_buoy  = std::sqrt(g * L_flame);
            const Real tan_tilt = (v_buoy > Real(1e-6)) 
                                 ? U_mag / v_buoy 
                                 : Real(0.0);
            const Real cos_tilt = Real(1.0) / std::sqrt(Real(1.0) + 
                                                        tan_tilt * tan_tilt);
            
            d_pre(i, j, k) = L_flame * cos_tilt * Fv;
        });
    }
}

// ============================================================================
// compute_preheating_distance_with_slope
//
// Enhanced version that includes explicit slope effects on flame tilt.
// Upslope fires have flames tilted forward by slope angle, increasing
// preheating distance.
//
// Flame tilt angle combines wind and slope effects:
//   tan(θ_tilt) = (U_eff + v_slope) / v_buoyancy
//   v_slope = k_slope × tan(slope_angle) × sqrt(g × L_f)
//
// where k_slope is an empirical coefficient (default: 0.5).
//
// Parameters:
//   preheating_dist_mf – output: preheating distance [m]
//   flame_length_mf    – input:  Byram flame length [m]
//   wind_speed_mf      – input:  effective wind speed magnitude [m/s]
//   slope_mf           – input:  slope angle [radians] (positive = upslope)
//   view_factor        – radiative view factor (default: 0.4)
//   k_slope            – slope tilt coefficient (default: 0.5)
//
// References:
//   Anderson, H.E. (1968). Sundance fire: An analysis of fire phenomena.
//     USDA Forest Service Research Paper INT-56.
//   Butler, B.W., et al. (2004). A method for evaluating the effectiveness of
//     firefighter escape routes.
// ============================================================================
inline void compute_preheating_distance_with_slope(
    MultiFab&       preheating_dist_mf,
    const MultiFab& flame_length_mf,
    const MultiFab& wind_speed_mf,
    const MultiFab& slope_mf,
    Real            view_factor = Real(0.4),
    Real            k_slope     = Real(0.5))
{
    constexpr Real g = WildfireConst::GRAVITY;

    for (MFIter mfi(preheating_dist_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const Lf     = flame_length_mf.const_array(mfi);
        auto const U      = wind_speed_mf.const_array(mfi);
        auto const slope  = slope_mf.const_array(mfi);
        auto d_pre        = preheating_dist_mf.array(mfi);
        const Real Fv     = view_factor;
        const Real k_s    = k_slope;

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const Real L_flame = Lf(i, j, k);
            
            if (L_flame <= Real(0.0)) {
                d_pre(i, j, k) = Real(0.0);
                return;
            }

            const Real U_mag      = U(i, j, k);
            const Real slope_rad  = slope(i, j, k);
            const Real v_buoy     = std::sqrt(g * L_flame);
            
            // Slope contribution to forward tilt (Anderson 1968 / Butler 2004)
            // Upslope fires: slope tilts flame forward, increasing preheating
            const Real v_slope = k_s * std::tan(slope_rad) * v_buoy;
            
            // Combined wind + slope tilt
            const Real tan_tilt = (v_buoy > Real(1e-6))
                                 ? (U_mag + v_slope) / v_buoy
                                 : Real(0.0);
            const Real cos_tilt = Real(1.0) / std::sqrt(Real(1.0) + 
                                                        tan_tilt * tan_tilt);
            
            // Preheating distance with slope enhancement
            d_pre(i, j, k) = L_flame * cos_tilt * Fv;
        });
    }
}

#endif // RADIATION_PREHEATING_H
