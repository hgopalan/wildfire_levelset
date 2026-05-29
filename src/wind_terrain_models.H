// ============================================================================
// Wind–Terrain Feedback Models for Rothermel ROS Enhancement
// ============================================================================
// Provides eight selectable options for how terrain-induced or fire-feedback
// winds modify the Rothermel wind field before (or after) ROS computation.
//
//  Option 1 – "none"                  : default Rothermel, no modification
//  Option 2 – "viegas_ros"            : replace R_mf with max(R_rothermel, R_viegas) in
//                                       eruptive cells (Viegas 2004); vel unchanged
//  Option 3 – "viegas_wind"           : add Viegas buoyancy-induced upslope wind only
//                                       in eruptive cells (tan_phi > tan_phi_c)
//  Option 4 – "canyon_wind"           : Rothermel (1983) canyon wind amplification
//                                       U_eff = U * (1 + k_canyon * tan_phi)
//  Option 5 – "viegas_neto"           : Viegas & Neto (1994) buoyancy-driven upslope
//                                       wind applied at every cell
//  Option 6 – "pimont"                : Pimont et al. (2009) exponential slope correction
//                                       U_eff = U * exp(k_pimont * tan_phi)
//  Option 7 – "windninja_ridge_canyon": WindNinja empirical ridge speed-up and canyon
//                                       channeling based on wind-slope alignment:
//                                         alignment = dot(wind_xy, upslope) / |wind|
//                                         ridge  (alignment > 0): f = 1 + k_ridge * tan_phi * alignment
//                                         canyon (alignment < 0): f = 1 + k_canyon_wn * tan_phi * |alignment|
//  Option 8 – "farsite_wind"          : FARSITE empirical wind stream simulation using
//                                       terrain curvature to identify ridges/valleys and
//                                       modify both wind speed and direction
//
// References:
//   Rothermel, R.C. (1983). How to Predict the Spread and Intensity of Forest
//     and Range Fires. USDA Forest Service General Technical Report INT-143.
//   Viegas, D.X. & Neto, L.P.S. (1994). Wind tunnel study of the convective
//     air flow of slope fires. COMOESTAS Annual Report.
//   Viegas, D.X. (2004). Slope and wind effects on fire propagation.
//     Int. J. Wildland Fire, 13(2), 143-156.
//   Pimont, F., Dupuy, J.L., Linn, R.R., Dupont, S. (2009). Validation of
//     FIRETEC wind-flows over a canopy and fuel-break. Int. J. Wildland Fire,
//     18(7), 775-790.
//   WindNinja: Forthofer, J.M. (2007). Modeling Wind in Complex Terrain for Use
//     in Fire Spread Prediction. Colorado State University MS thesis.
//   Finney, M.A. (1998). FARSITE: Fire Area Simulator—Model Development
//     and Evaluation. USDA Forest Service Research Paper RMRS-RP-4.
// ============================================================================

#ifndef WIND_TERRAIN_MODELS_H
#define WIND_TERRAIN_MODELS_H

#include "constants.H"
#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_GpuLaunch.H>
#include "parse_inputs.H"
#include "farsite_wind_model.H"

using namespace amrex;

// ============================================================================
// wind_terrain_modifies_velocity
//
// Returns true when the selected wind-terrain model (Options 3-8) produces a
// terrain-modified velocity field in vel_effective that differs from the
// ambient vel.  For "none" (Option 1) and "viegas_ros" (Option 2) the
// ambient vel is used directly and this returns false.
// ============================================================================
inline bool wind_terrain_modifies_velocity(const InputParameters& inputs)
{
    const std::string& m = inputs.wind_terrain.model;
    return (m == "viegas_wind" ||
            m == "canyon_wind" ||
            m == "viegas_neto" ||
            m == "pimont"      ||
            m == "windninja_ridge_canyon" ||
            m == "farsite_wind");
}

// ============================================================================
// apply_wind_terrain_velocity
//
// Fills vel_eff with a terrain-modified copy of the ambient wind field vel.
// Called for Options 3 (viegas_wind), 4 (canyon_wind), 5 (viegas_neto),
// 6 (pimont), 7 (windninja_ridge_canyon), and 8 (farsite_wind).
// For "none" and "viegas_ros" this function should not be called.
//
// Parameters:
//   vel_eff        – output: effective wind for Rothermel (3 components, same
//                    layout as vel)
//   vel            – input:  ambient wind field [m/s]
//   terrain_slopes – optional: 2-component MultiFab (slope_x = dz/dx,
//                    slope_y = dz/dy); nullptr = flat terrain (no modification)
//   geom           – grid geometry (required for farsite_wind curvature calc)
//   inputs         – full input parameter set (uses wind_terrain, viegas, and
//                    rothermel sub-structs)
// ============================================================================
inline void apply_wind_terrain_velocity(
    MultiFab&              vel_eff,
    const MultiFab&        vel,
    const MultiFab*        terrain_slopes,
    const Geometry&        geom,
    const InputParameters& inputs)
{
    const std::string& model = inputs.wind_terrain.model;

    // -----------------------------------------------------------------------
    // Option 8 - FARSITE wind stream simulation
    // -----------------------------------------------------------------------
    if (model == "farsite_wind") {
        if (terrain_slopes == nullptr) {
            // No terrain → copy ambient wind unchanged
            MultiFab::Copy(vel_eff, vel, 0, 0, 3, 0);
            return;
        }
        apply_farsite_wind_terrain(
            vel_eff, vel, *terrain_slopes, geom,
            inputs.wind_terrain.k_ridge_farsite,
            inputs.wind_terrain.k_shelter,
            inputs.wind_terrain.k_valley,
            inputs.wind_terrain.k_deflection,
            inputs.wind_terrain.min_curvature);
        return;
    }

    // -----------------------------------------------------------------------
    // Parameters captured for GPU lambda
    // -----------------------------------------------------------------------
    const Real k_canyon  = inputs.wind_terrain.k_canyon;
    const Real k_pimont  = inputs.wind_terrain.k_pimont;
    const Real tan_phi_c = inputs.viegas.tan_phi_c;

    // Buoyancy velocity scale used by viegas_wind and viegas_neto:
    //   v_b = sqrt(g * delta_m * (T_f - T_a) / T_a)
    // delta is stored in feet; multiply by 0.3048 to convert to metres for SI buoyancy formula.
    constexpr Real g_SI   = WildfireConst::GRAVITY;
    const Real delta_m    = inputs.rothermel.delta * WildfireConst::FT_TO_M;  // feet to metres
    const Real T_a        = inputs.viegas.T_a;
    const Real T_f        = inputs.viegas.T_f;
    const Real T_ratio    = (T_f - T_a) / std::max(T_a, Real(1.0));
    const Real v_b        = std::sqrt(std::max(g_SI * delta_m * T_ratio, Real(0.0)));

    const bool has_slope  = (terrain_slopes != nullptr);

    // Enumerate model so GPU lambda captures an int (no std::string on device)
    // 3 = viegas_wind, 4 = canyon_wind, 5 = viegas_neto, 6 = pimont, 
    // 7 = windninja_ridge_canyon, 8 = farsite_wind
    int model_id = 0;
    if      (model == "viegas_wind")            model_id = 3;
    else if (model == "canyon_wind")            model_id = 4;
    else if (model == "viegas_neto")            model_id = 5;
    else if (model == "pimont")                 model_id = 6;
    else if (model == "windninja_ridge_canyon") model_id = 7;
    // Note: farsite_wind (model_id 8) handled separately above

    // WindNinja ridge/canyon coefficients (Option 7)
    const Real k_ridge     = inputs.wind_terrain.k_ridge;
    const Real k_canyon_wn = inputs.wind_terrain.k_canyon_wn;

    // -----------------------------------------------------------------------
    // Domain-mean wind for extreme-cell clamping (same logic as Rothermel)
    // -----------------------------------------------------------------------
    const Long ncells_vel  = vel.boxArray().numPts();
    const Real ux_mean = (ncells_vel > 0) ? vel.sum(0) / Real(ncells_vel) : Real(0.0);
    const Real uy_mean = (ncells_vel > 0) ? vel.sum(1) / Real(ncells_vel) : Real(0.0);
    // Minimum slope magnitude to avoid division-by-zero in uphill unit-vector.
    const Real WT_SLOPE_EPS = WildfireConst::WT_SLOPE_EPS;
    for (MFIter mfi(vel_eff); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const v  = vel.const_array(mfi);
        auto       ve = vel_eff.array(mfi);

        Array4<const Real> sl;
        if (has_slope) sl = terrain_slopes->const_array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            // ----------------------------------------------------------------
            // Ambient wind, clamping extreme values to domain mean
            // ----------------------------------------------------------------
            Real ux = v(i, j, k, 0);
            Real uy = v(i, j, k, 1);
            constexpr Real MAX_WIND_SQ = WildfireConst::MAX_WIND_SPEED_SQ;
            if (ux * ux + uy * uy > MAX_WIND_SQ) {
                ux = ux_mean;
                uy = uy_mean;
            }

            // ----------------------------------------------------------------
            // Terrain slope magnitude  tan(phi) = |grad z|
            // ----------------------------------------------------------------
            Real sx = Real(0.0), sy = Real(0.0), tan_phi = Real(0.0);
            if (has_slope) {
                sx      = sl(i, j, k, 0);
                sy      = sl(i, j, k, 1);
                tan_phi = std::sqrt(sx * sx + sy * sy);
            }

            // ----------------------------------------------------------------
            // Apply selected model
            // ----------------------------------------------------------------
            Real ux_eff = ux;
            Real uy_eff = uy;

            if (model_id == 3) {
                // Option 3 – viegas_wind: Viegas buoyancy-induced upslope wind
                // only in eruptive cells (tan_phi > tan_phi_c).
                //   U_ind = v_b * tan_phi, directed upslope
                if (tan_phi > tan_phi_c) {
                    Real inv_slope = Real(1.0) / std::max(tan_phi, WT_SLOPE_EPS);
                    Real uphill_x  = sx * inv_slope;  // unit vector upslope
                    Real uphill_y  = sy * inv_slope;
                    Real U_ind     = v_b * tan_phi;
                    ux_eff = ux + U_ind * uphill_x;
                    uy_eff = uy + U_ind * uphill_y;
                }
            } else if (model_id == 4) {
                // Option 4 – canyon_wind: Rothermel (1983) terrain channeling.
                //   |U_eff| = |U| * (1 + k_canyon * tan_phi)
                // Both components are scaled to preserve wind direction while
                // amplifying speed proportionally to slope steepness.
                Real factor = Real(1.0) + k_canyon * tan_phi;
                if (factor < Real(1.0)) factor = Real(1.0);  // never suppress wind
                ux_eff = ux * factor;
                uy_eff = uy * factor;
            } else if (model_id == 5) {
                // Option 5 – viegas_neto: Viegas & Neto (1994) buoyancy-driven
                // upslope wind applied at every cell (no eruptive threshold).
                //   U_ind = v_b * tan_phi, directed upslope
                if (tan_phi > WT_SLOPE_EPS) {
                    Real inv_slope = Real(1.0) / tan_phi;
                    Real uphill_x  = sx * inv_slope;
                    Real uphill_y  = sy * inv_slope;
                    Real U_ind     = v_b * tan_phi;
                    ux_eff = ux + U_ind * uphill_x;
                    uy_eff = uy + U_ind * uphill_y;
                }
            } else if (model_id == 6) {
                // Option 6 – pimont: Pimont et al. (2009) exponential slope
                // correction (calibrated from CFD data).
                //   |U_eff| = |U| * exp(k_pimont * tan_phi)
                Real factor = std::exp(k_pimont * tan_phi);
                ux_eff = ux * factor;
                uy_eff = uy * factor;
            } else if (model_id == 7) {
                // Option 7 – windninja_ridge_canyon: WindNinja empirical ridge
                // speed-up and canyon channeling based on wind-slope alignment.
                //
                //   alignment = dot(wind_xy, upslope_unit) / |wind|
                //     > 0: wind is blowing upslope  → ridge speed-up
                //     < 0: wind is blowing downslope → canyon channeling
                //
                //   Ridge (alignment > 0):
                //     f = 1 + k_ridge * tan_phi * alignment
                //     Speed increase proportional to how directly wind climbs
                //     the slope (Forthofer 2007, terrain-driven acceleration).
                //
                //   Canyon (alignment < 0):
                //     f = 1 + k_canyon_wn * tan_phi * |alignment|
                //     Channel the wind downslope (drainage flow amplification).
                //
                // The factor is always >= 1 to avoid suppressing wind speed.
                if (tan_phi > WT_SLOPE_EPS) {
                    const Real wind_spd = std::sqrt(ux * ux + uy * uy);
                    if (wind_spd > WT_SLOPE_EPS) {
                        // Upslope unit vector (points in direction of steepest ascent)
                        const Real inv_slope = Real(1.0) / tan_phi;
                        const Real uphill_x  = sx * inv_slope;
                        const Real uphill_y  = sy * inv_slope;
                        // Alignment: cos(angle between wind and upslope direction)
                        const Real alignment = (ux * uphill_x + uy * uphill_y) / wind_spd;
                        Real factor = Real(1.0);
                        if (alignment > Real(0.0)) {
                            // Wind blowing upslope: ridge speed-up
                            factor = Real(1.0) + k_ridge * tan_phi * alignment;
                        } else if (alignment < Real(0.0)) {
                            // Wind blowing downslope: canyon channeling
                            factor = Real(1.0) + k_canyon_wn * tan_phi * (-alignment);
                        }
                        // Enforce factor >= 1 (never suppress wind)
                        if (factor < Real(1.0)) factor = Real(1.0);
                        ux_eff = ux * factor;
                        uy_eff = uy * factor;
                    }
                }
            }

            // ----------------------------------------------------------------
            // Store effective velocity (z-component unchanged)
            // ----------------------------------------------------------------
            ve(i, j, k, 0) = ux_eff;
            ve(i, j, k, 1) = uy_eff;
            ve(i, j, k, 2) = v(i, j, k, 2);
        });
    }
}

// ============================================================================
// apply_wind_terrain_effective_velocity
//
// Convenience wrapper: when the selected wind-terrain model modifies velocity
// (Options 3-8), calls apply_wind_terrain_velocity to fill vel_eff with the
// terrain-corrected wind.  Otherwise, copies the ambient vel into vel_eff so
// that downstream code can always read from vel_eff uniformly.
// ============================================================================
inline void apply_wind_terrain_effective_velocity(
    MultiFab&              vel_eff,
    const MultiFab&        vel,
    const MultiFab*        terrain_slopes,
    const Geometry&        geom,
    const InputParameters& inputs)
{
    if (wind_terrain_modifies_velocity(inputs)) {
        apply_wind_terrain_velocity(vel_eff, vel, terrain_slopes, geom, inputs);
    } else {
        MultiFab::Copy(vel_eff, vel, 0, 0, 3, 0);
    }
}


inline void apply_upslope_convection_velocity(
    MultiFab&              vel_eff,
    const MultiFab*        terrain_slopes,
    const InputParameters& inputs)
{
    if (inputs.upslope_convection.enable != 1 || terrain_slopes == nullptr) return;

    const Real coeff = inputs.upslope_convection.coefficient;
    const Real tan_min = std::tan(inputs.upslope_convection.min_slope_deg * WildfireConst::PI / Real(180.0));
    const Real max_added = inputs.upslope_convection.max_added_speed;

    for (MFIter mfi(vel_eff); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto ve = vel_eff.array(mfi);
        auto const slopes = terrain_slopes->const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const Real sx = slopes(i, j, k, 0);
            const Real sy = slopes(i, j, k, 1);
            const Real tan_phi = std::sqrt(sx*sx + sy*sy);
            if (tan_phi < tan_min || tan_phi <= Real(1.0e-12)) return;
            const Real added = amrex::min(max_added, amrex::max(Real(0.0), coeff * tan_phi));
            const Real inv = Real(1.0) / tan_phi;
            ve(i, j, k, 0) += added * sx * inv;
            ve(i, j, k, 1) += added * sy * inv;
        });
    }
}

// ============================================================================
// apply_viegas_ros_override
//
// Option 2 ("viegas_ros"): overrides R_mf in eruptive cells so that
// R_final = max(R_rothermel, R_viegas).
//
// Must be called after both compute_rothermel_R and compute_viegas_diagnostics
// have been run, so that R_mf holds the Rothermel ROS and viegas_data holds
// the Viegas ROS (component 0) and eruptive flag (component 1).
//
// Parameters:
//   R_mf        – in/out: Rothermel ROS field [m/s]; overridden in eruptive cells
//   viegas_data – input:  Viegas diagnostic MultiFab (VIEGAS_NCOMP components):
//                   component 0 = viegas_ROS [m/s]
//                   component 1 = viegas_eruptive_flag (1.0 where eruptive)
// ============================================================================
inline void apply_viegas_ros_override(
    MultiFab&       R_mf,
    const MultiFab& viegas_data)
{
    for (MFIter mfi(R_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto       R  = R_mf.array(mfi);
        auto const vd = viegas_data.const_array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            // In eruptive cells take the larger of Rothermel and Viegas ROS.
            // Where Viegas predicts slower spread (typically on gentle slopes),
            // Rothermel is preserved unchanged.
            if (vd(i, j, k, 1) > Real(0.5)) {  // eruptive_flag == 1
                R(i, j, k) = std::max(R(i, j, k), vd(i, j, k, 0));
            }
        });
    }
}

#endif // WIND_TERRAIN_MODELS_H
