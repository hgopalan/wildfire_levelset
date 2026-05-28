// ============================================================================
// FARSITE Wind Stream Simulation Model
// ============================================================================
// Implements FARSITE's internal simplified wind flow model around terrain
// using empirical wind deflection rules based on terrain features.
//
// FARSITE (Fire Area Simulator) uses a simplified empirical approach to
// adjust wind speed and direction based on local terrain characteristics
// without requiring full CFD computation. This model implements the key
// empirical relationships used in FARSITE for wind-terrain interactions.
//
// Key Features:
//  1. Ridge acceleration: Wind speeds up over ridge crests
//  2. Lee-side sheltering: Wind slows in protected areas downwind of ridges
//  3. Wind direction deflection: Wind follows terrain contours
//  4. Valley channeling: Wind is funneled through valley corridors
//  5. Slope-aspect interaction: Combined effect of slope orientation and wind
//
// The model uses terrain curvature, slope, and aspect to compute local wind
// modifications. Unlike full CFD models (e.g., WindNinja mass-consistent),
// FARSITE's approach uses simple algebraic relationships calibrated from
// field observations and wind tunnel experiments.
//
// Mathematical Model:
// -------------------
//
// 1. Terrain Analysis:
//    - Compute terrain curvature (2nd derivatives of elevation)
//    - Identify ridge/valley features from curvature sign and magnitude
//    - Calculate terrain aspect (direction of steepest ascent)
//
// 2. Wind Speed Adjustment:
//    Ridge amplification (positive curvature):
//      f_ridge = 1 + k_ridge_farsite * curvature * (1 + alignment)
//    
//    Lee-side sheltering (negative curvature, downwind):
//      f_shelter = 1 - k_shelter * |curvature| * (1 - alignment)
//    
//    Valley channeling (negative curvature):
//      f_valley = 1 + k_valley * |curvature| * channel_factor
//
// 3. Wind Direction Deflection:
//    The wind direction is deflected to follow terrain contours:
//      delta_theta = k_deflection * sin(aspect - wind_dir) * tan(slope)
//
// 4. Combined Effect:
//    U_eff = U * f_speed
//    theta_eff = theta + delta_theta
//
// Parameters (wind_terrain.farsite_* namespace):
//   k_ridge_farsite  – ridge speed-up coefficient (default: 1.5)
//   k_shelter        – lee-side sheltering coefficient (default: 0.6)
//   k_valley         – valley channeling coefficient (default: 0.8)
//   k_deflection     – wind direction deflection coefficient (default: 0.3)
//   curvature_smooth – smoothing radius for curvature calculation [m] (default: 50.0)
//   min_curvature    – minimum curvature threshold to apply corrections (default: 0.0001)
//
// References:
//   Finney, M.A. (1998). FARSITE: Fire Area Simulator—Model Development
//     and Evaluation. USDA Forest Service Research Paper RMRS-RP-4.
//   Finney, M.A. et al. (2011). A method for ensemble wildland fire
//     simulation. Environmental Modeling and Assessment, 16(2), 153-167.
//   Butler, B.W. et al. (1998). Fire behavior associated with the 1994
//     South Canyon Fire. USDA Forest Service Research Paper RMRS-RP-9.
//   Rothermel, R.C. & Rinehart, G.C. (1983). Field procedures for
//     verification and adjustment of fire behavior predictions.
//     USDA Forest Service General Technical Report INT-142.
// ============================================================================

#ifndef FARSITE_WIND_MODEL_H
#define FARSITE_WIND_MODEL_H

#include "constants.H"
#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_GpuLaunch.H>
#include <cmath>

using namespace amrex;

// ============================================================================
// compute_terrain_curvature
//
// Computes terrain curvature (Laplacian of elevation) from slope field.
// Positive curvature indicates ridges (convex), negative indicates valleys
// (concave). Uses centered finite differences on the slope field.
//
// Parameters:
//   curvature_mf   – output: terrain curvature field (1 component)
//   terrain_slopes – input: 2-component MultiFab (slope_x = dz/dx, slope_y = dz/dy)
//   geom           – grid geometry for cell spacing
// ============================================================================
inline void compute_terrain_curvature(
    MultiFab&       curvature_mf,
    const MultiFab& terrain_slopes,
    const Geometry& geom)
{
    const Real* dx_arr = geom.CellSize();
    const Real dx = dx_arr[0];
    const Real dy = dx_arr[1];
    
    for (MFIter mfi(curvature_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto       curv = curvature_mf.array(mfi);
        auto const sl   = terrain_slopes.const_array(mfi);
        
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            // Compute second derivatives of elevation using centered differences
            // on the slope field (slope_x = dz/dx, slope_y = dz/dy)
            
            // d²z/dx² = d(slope_x)/dx
            Real d2z_dx2 = Real(0.0);
            if (i > bx.smallEnd(0) && i < bx.bigEnd(0)) {
                d2z_dx2 = (sl(i+1, j, k, 0) - sl(i-1, j, k, 0)) / (Real(2.0) * dx);
            }
            
            // d²z/dy² = d(slope_y)/dy
            Real d2z_dy2 = Real(0.0);
            if (j > bx.smallEnd(1) && j < bx.bigEnd(1)) {
                d2z_dy2 = (sl(i, j+1, k, 1) - sl(i, j-1, k, 1)) / (Real(2.0) * dy);
            }
            
            // Laplacian (mean curvature indicator)
            curv(i, j, k) = d2z_dx2 + d2z_dy2;
        });
    }
}

// ============================================================================
// apply_farsite_wind_terrain
//
// Applies FARSITE's empirical wind-terrain corrections to the velocity field.
// Modifies both wind speed and direction based on terrain curvature, slope,
// and aspect.
//
// Parameters:
//   vel_eff        – output: effective wind field [m/s] (3 components)
//   vel            – input: ambient wind field [m/s]
//   terrain_slopes – input: 2-component MultiFab (slope_x, slope_y)
//   geom           – grid geometry
//   k_ridge_farsite– ridge speed-up coefficient
//   k_shelter      – lee-side sheltering coefficient
//   k_valley       – valley channeling coefficient
//   k_deflection   – wind direction deflection coefficient
//   min_curvature  – minimum curvature threshold
// ============================================================================
inline void apply_farsite_wind_terrain(
    MultiFab&       vel_eff,
    const MultiFab& vel,
    const MultiFab& terrain_slopes,
    const Geometry& geom,
    Real            k_ridge_farsite,
    Real            k_shelter,
    Real            k_valley,
    Real            k_deflection,
    Real            min_curvature)
{
    // Compute terrain curvature
    MultiFab curvature_mf(vel_eff.boxArray(), vel_eff.DistributionMap(), 1, 0);
    compute_terrain_curvature(curvature_mf, terrain_slopes, geom);
    
    // Domain-mean wind for extreme-cell clamping
    const Long ncells_vel = vel.boxArray().numPts();
    const Real ux_mean = (ncells_vel > 0) ? vel.sum(0) / Real(ncells_vel) : Real(0.0);
    const Real uy_mean = (ncells_vel > 0) ? vel.sum(1) / Real(ncells_vel) : Real(0.0);
    
    const Real FARSITE_WIND_EPS = WildfireConst::WT_SLOPE_EPS;
    
    for (MFIter mfi(vel_eff); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const v    = vel.const_array(mfi);
        auto const sl   = terrain_slopes.const_array(mfi);
        auto const curv = curvature_mf.const_array(mfi);
        auto       ve   = vel_eff.array(mfi);
        
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            // ----------------------------------------------------------------
            // Ambient wind, clamping extreme values
            // ----------------------------------------------------------------
            Real ux = v(i, j, k, 0);
            Real uy = v(i, j, k, 1);
            constexpr Real MAX_WIND_SQ = WildfireConst::MAX_WIND_SPEED_SQ;
            if (ux * ux + uy * uy > MAX_WIND_SQ) {
                ux = ux_mean;
                uy = uy_mean;
            }
            
            const Real wind_spd = std::sqrt(ux * ux + uy * uy);
            
            // ----------------------------------------------------------------
            // Terrain slope and aspect
            // ----------------------------------------------------------------
            const Real sx = sl(i, j, k, 0);  // dz/dx
            const Real sy = sl(i, j, k, 1);  // dz/dy
            const Real tan_phi = std::sqrt(sx * sx + sy * sy);
            
            // Terrain aspect (direction of steepest ascent)
            Real aspect = Real(0.0);
            if (tan_phi > FARSITE_WIND_EPS) {
                aspect = std::atan2(sy, sx);  // radians, -π to π
            }
            
            // ----------------------------------------------------------------
            // Wind direction
            // ----------------------------------------------------------------
            Real wind_dir = Real(0.0);
            if (wind_spd > FARSITE_WIND_EPS) {
                wind_dir = std::atan2(uy, ux);  // radians, -π to π
            }
            
            // ----------------------------------------------------------------
            // Terrain curvature and feature classification
            // ----------------------------------------------------------------
            const Real curvature = curv(i, j, k);
            const Real abs_curvature = std::abs(curvature);
            
            // ----------------------------------------------------------------
            // Wind-terrain alignment
            // Alignment = cos(angle between wind and upslope direction)
            //  +1: wind blowing directly upslope
            //   0: wind perpendicular to slope
            //  -1: wind blowing directly downslope
            // ----------------------------------------------------------------
            Real alignment = Real(0.0);
            if (wind_spd > FARSITE_WIND_EPS && tan_phi > FARSITE_WIND_EPS) {
                // Upslope unit vector
                const Real inv_slope = Real(1.0) / tan_phi;
                const Real uphill_x = sx * inv_slope;
                const Real uphill_y = sy * inv_slope;
                // Alignment: dot product of wind and upslope unit vectors
                alignment = (ux * uphill_x + uy * uphill_y) / wind_spd;
            }
            
            // ----------------------------------------------------------------
            // Speed modification factor (default: no change)
            // ----------------------------------------------------------------
            Real f_speed = Real(1.0);
            
            if (abs_curvature > min_curvature) {
                if (curvature > Real(0.0)) {
                    // -------- Ridge (positive curvature) --------
                    // Speed increases on ridges, especially when wind climbs
                    // f_ridge = 1 + k_ridge * curv * (1 + alignment)
                    // Maximum effect when alignment = +1 (upslope wind)
                    f_speed = Real(1.0) + k_ridge_farsite * curvature 
                              * (Real(1.0) + alignment);
                    
                } else {
                    // -------- Valley (negative curvature) --------
                    // Two competing effects:
                    //  1. Lee-side sheltering (downwind): reduces speed
                    //  2. Valley channeling (cross-wind): increases speed
                    
                    if (alignment < Real(-0.3)) {
                        // Wind blowing downslope into valley → sheltering
                        // f_shelter = 1 - k_shelter * |curv| * (1 - alignment)
                        // Maximum sheltering when alignment = -1 (downslope)
                        f_speed = Real(1.0) - k_shelter * abs_curvature 
                                  * (Real(1.0) - alignment);
                        // Prevent negative or zero speed
                        f_speed = std::max(f_speed, Real(0.2));
                        
                    } else {
                        // Wind across or up valley → channeling
                        // Valley channeling intensifies when wind is perpendicular
                        // to valley axis (alignment near 0) or upslope
                        Real channel_factor = Real(1.0) - alignment * alignment;
                        f_speed = Real(1.0) + k_valley * abs_curvature 
                                  * channel_factor;
                    }
                }
                
                // Limit speed factor to reasonable range
                f_speed = std::max(Real(0.2), std::min(f_speed, Real(3.0)));
            }
            
            // ----------------------------------------------------------------
            // Direction deflection (terrain-following flow)
            // ----------------------------------------------------------------
            Real delta_theta = Real(0.0);
            if (wind_spd > FARSITE_WIND_EPS && tan_phi > FARSITE_WIND_EPS) {
                // Wind deflects to follow terrain contours
                // delta_theta = k_deflection * sin(aspect - wind_dir) * tan(slope)
                // Maximum deflection when wind crosses contours at 90°
                Real angle_diff = aspect - wind_dir;
                // Normalize to [-π, π]
                while (angle_diff > Real(M_PI))  angle_diff -= Real(2.0 * M_PI);
                while (angle_diff < Real(-M_PI)) angle_diff += Real(2.0 * M_PI);
                
                delta_theta = k_deflection * std::sin(angle_diff) * tan_phi;
                // Limit deflection to ±45° per cell
                delta_theta = std::max(Real(-0.7854), std::min(delta_theta, Real(0.7854)));
            }
            
            // ----------------------------------------------------------------
            // Apply speed and direction modifications
            // ----------------------------------------------------------------
            const Real new_wind_dir = wind_dir + delta_theta;
            const Real new_wind_spd = wind_spd * f_speed;
            
            // Convert back to Cartesian components
            Real ux_eff = new_wind_spd * std::cos(new_wind_dir);
            Real uy_eff = new_wind_spd * std::sin(new_wind_dir);
            
            // ----------------------------------------------------------------
            // Store effective velocity (z-component unchanged)
            // ----------------------------------------------------------------
            ve(i, j, k, 0) = ux_eff;
            ve(i, j, k, 1) = uy_eff;
            ve(i, j, k, 2) = v(i, j, k, 2);
        });
    }
}

#endif // FARSITE_WIND_MODEL_H
