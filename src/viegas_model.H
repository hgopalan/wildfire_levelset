// ============================================================================
// Viegas (2004) Eruptive Fire Model - Fire-Terrain Feedback Diagnostics
// ============================================================================
// Implements the eruptive fire spread diagnostics and flame-tilt feedback
// from:
//
//   Viegas, D.X. (2004). "Slope and wind effects on fire propagation."
//   International Journal of Wildland Fire, 13(2), 143-156.
//   https://doi.org/10.1071/WF03046
//
// ----------------------------------------------------------------------------
// Merging rules with Rothermel (1972)
// ----------------------------------------------------------------------------
//  Compatible (implemented here):
//   - Rothermel no-wind no-slope ROS (R_0) is the baseline for Viegas ROS.
//   - Viegas eruptive-regime diagnostics flag cells where Rothermel under-
//     estimates ROS (slope > critical slope phi_c).
//   - Both models receive the same environmental inputs (slope, wind, fuels).
//   - Viegas flame-tilt output is used for hazard assessment only.
//
//  Not compatible (NOT implemented, never used to modify Rothermel):
//   - Viegas induced wind is NOT fed back into Rothermel.
//   - Viegas eruptive acceleration does NOT modify Rothermel's ROS.
//   - Rothermel slope factor phi_s is NOT replaced by Viegas critical slope.
//   - Viegas dynamic ROS is NOT inserted into Rothermel's static formula.
//
// ----------------------------------------------------------------------------
// Viegas (2004) slope enhancement factor
// ----------------------------------------------------------------------------
// The empirical relationship between the dimensionless spread rate r = R/R_0
// and terrain slope angle phi (from Viegas & Neto 1994 and Viegas 2004):
//
//   Phi_s_V = exp(a_V * tan(phi))
//
// where a_V is a fuel-dependent calibration coefficient (Viegas 2004 reports
// a_V ≈ 1.83 for standard laboratory fuels).  The exponential form grows
// significantly faster than Rothermel's quadratic phi_s for large slopes,
// capturing the eruptive acceleration observed in field experiments.
//
// Viegas ROS (parallel, not coupled to Rothermel):
//
//   R_V = R_0 * (1 + phi_w_Rothermel) * Phi_s_V
//
// where R_0 and phi_w_Rothermel are the Rothermel no-slope baseline and wind
// factor respectively.  The wind contribution is shared with Rothermel; only
// the slope factor differs.
//
// ----------------------------------------------------------------------------
// Eruptive regime diagnostics
// ----------------------------------------------------------------------------
// Eruptive fire occurs when the terrain slope exceeds a critical threshold
// phi_c.  Above phi_c the fire does not approach a steady state – it
// accelerates continuously (Viegas 2004, Section 3).  The typical critical
// slope is tan(phi_c) ≈ 0.4 (~22°), user-adjustable via viegas.tan_phi_c.
//
//   eruptive_flag = 1.0  when  tan(phi) > tan_phi_c
//                  0.0  otherwise
//
// The ROS excess ratio indicates where Rothermel under-predicts:
//
//   ROS_excess = (R_V - R_Rothermel) / max(R_Rothermel, epsilon)
//
// A positive value means Viegas predicts faster spread than Rothermel.
//
// ----------------------------------------------------------------------------
// Viegas flame-tilt feedback (hazard assessment only)
// ----------------------------------------------------------------------------
// The flame-tilt angle from vertical (Viegas 2004 and Balbi 2009 framework):
//
//   tan(alpha) = U / v_b + tan(phi)
//   alpha      = atan(tan(alpha))
//
// where:
//   U   = wind speed [m/s]
//   v_b = buoyancy velocity scale = sqrt(g * delta_m * (T_f - T_a) / T_a)
//   phi = terrain slope angle
//
// alpha is output for hazard assessment (flame projection onto unburned fuel).
// It is NOT used to modify either model's ROS computation.
//
// ----------------------------------------------------------------------------
// Output MultiFab viegas_data (VIEGAS_NCOMP = 5 components)
// ----------------------------------------------------------------------------
//   0 – viegas_ROS          : Viegas ROS [m/s]
//   1 – viegas_eruptive_flag: 1.0 when eruptive conditions met, 0.0 otherwise
//   2 – viegas_ROS_excess   : (R_V - R_Rothermel) / R_Rothermel ratio
//   3 – viegas_flame_tilt   : Viegas flame tilt angle from vertical [rad]
//   4 – viegas_slope_factor : Viegas slope enhancement factor Phi_s_V
// ============================================================================

#ifndef VIEGAS_MODEL_H
#define VIEGAS_MODEL_H

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_GpuLaunch.H>
#include "parse_inputs.H"
#include "rothermel_model.H"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace amrex;

// Number of output components in the viegas_data MultiFab.
static constexpr int VIEGAS_NCOMP = 5;

// Tiny positive value used as denominator floor when computing ROS excess.
static constexpr amrex::Real VIEGAS_EPS = Real(1.0e-10);

// Minimum buoyancy velocity [m/s] to prevent division-by-zero in tilt calc.
static constexpr amrex::Real VIEGAS_VB_MIN = Real(1.0e-3);

// ============================================================================
// compute_viegas_diagnostics
//
// Fills viegas_data (VIEGAS_NCOMP components) using the Viegas (2004) model.
//
// This function is a DIAGNOSTIC: it runs in parallel with Rothermel and does
// NOT modify R_mf or any Rothermel-internal quantity.
//
// Parameters:
//   viegas_data    – output MultiFab (VIEGAS_NCOMP components, no ghost cells)
//   R_mf           – input: Rothermel ROS field [m/s] (1 component, already
//                    computed by compute_rothermel_R)
//   vel            – input: wind velocity field [m/s]
//   rp             – Rothermel fuel parameters (provides R_0 baseline, delta,
//                    wind conversion, C/B/beta_ratio_E for phi_w)
//   vp             – Viegas model parameters
//   terrain_slopes – optional: 2-component MultiFab (slope_x, slope_y) in
//                    tan(slope) units; nullptr = flat terrain (all zeros)
//   fuel_model_mf  – optional: per-cell integer fuel codes
//   fuel_table     – optional: GPU array of RothermelComputed indexed by code
//   fuel_table_size– number of entries in fuel_table
// ============================================================================
static void compute_viegas_diagnostics(
    MultiFab&              viegas_data,
    const MultiFab&        R_mf,
    const MultiFab&        vel,
    const RothermelParams& rp,
    const InputParameters::ViegasParams& vp,
    const MultiFab*        terrain_slopes = nullptr,
    const MultiFab*        fuel_model_mf  = nullptr,
    const RothermelComputed* fuel_table   = nullptr,
    int                    fuel_table_size = 0)
{
    // Physical constant
    constexpr Real g = Real(9.81);  // gravitational acceleration [m/s²]

    // -----------------------------------------------------------------------
    // Pre-compute global Rothermel parameters for the default fuel
    // -----------------------------------------------------------------------
    const RothermelComputed rc_global = compute_rothermel_params(rp);

    // Global defaults captured for GPU lambda
    const Real R0_global        = rc_global.R0;
    const Real C_global         = rc_global.C;
    const Real B_global         = rc_global.B;
    const Real bRE_global       = rc_global.beta_ratio_E;
    const Real waf_global       = rc_global.waf;
    const Real U_max_global     = rc_global.U_max_ftmin;

    // Unit conversions and flags from RothermelParams
    const Real wind_conv        = rp.wind_conv;
    const Real ros_conv         = rp.ros_conv;
    const int  use_waf_flag     = rp.use_waf;
    const int  use_wind_limit   = rp.use_wind_limit;

    // Fuel bed depth in metres (for buoyancy velocity)
    const Real delta_m = rp.delta * Real(0.3048);  // ft → m

    // -----------------------------------------------------------------------
    // Viegas parameters (captured by value for GPU lambda)
    // -----------------------------------------------------------------------
    const Real a_V         = vp.a_V;
    const Real tan_phi_c   = vp.tan_phi_c;
    const Real T_a         = vp.T_a;
    const Real T_f         = vp.T_f;

    // Buoyancy velocity [m/s] using global fuel depth.
    // Per-cell delta is not tracked independently so the global value is used;
    // the per-cell wind & slope corrections are the primary source of spatial
    // variation in flame tilt.
    const Real T_ratio = (T_f - T_a) / std::max(T_a, Real(1.0));
    const Real v_b_global = std::sqrt(std::max(g * delta_m * T_ratio, Real(0.0)));

    // -----------------------------------------------------------------------
    // Fuel table availability
    // -----------------------------------------------------------------------
    const bool use_fuel_table = (fuel_model_mf != nullptr &&
                                 fuel_table != nullptr &&
                                 fuel_table_size > 0);

    // -----------------------------------------------------------------------
    // Domain-mean wind for extreme-cell clamping (same logic as Rothermel)
    // -----------------------------------------------------------------------
    const Long ncells_vel   = vel.boxArray().numPts();
    const Real ux_mean = (ncells_vel > 0) ? vel.sum(0) / Real(ncells_vel) : Real(0.0);
    const Real uy_mean = (ncells_vel > 0) ? vel.sum(1) / Real(ncells_vel) : Real(0.0);

    const bool has_slope = (terrain_slopes != nullptr);

    for (MFIter mfi(viegas_data); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();

        auto const Rr  = R_mf.const_array(mfi);     // Rothermel ROS [m/s]
        auto const v   = vel.const_array(mfi);

        // Optional terrain slopes
        Array4<const Real> sl;
        if (has_slope) sl = terrain_slopes->const_array(mfi);

        // Optional per-cell fuel model codes
        Array4<const Real> fm;
        if (use_fuel_table) fm = fuel_model_mf->const_array(mfi);

        auto vd = viegas_data.array(mfi);

        // Capture scalars for lambda
        const Real R0_def    = R0_global;
        const Real C_def     = C_global;
        const Real B_def     = B_global;
        const Real bRE_def   = bRE_global;
        const Real waf_def   = waf_global;
        const Real Umax_def  = U_max_global;
        const Real v_b_def   = v_b_global;

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {

            // ----------------------------------------------------------------
            // Zero all outputs by default
            // ----------------------------------------------------------------
            vd(i, j, k, 0) = Real(0.0);  // viegas_ROS
            vd(i, j, k, 1) = Real(0.0);  // viegas_eruptive_flag
            vd(i, j, k, 2) = Real(0.0);  // viegas_ROS_excess
            vd(i, j, k, 3) = Real(0.0);  // viegas_flame_tilt
            vd(i, j, k, 4) = Real(1.0);  // viegas_slope_factor (flat terrain = 1)

            // ----------------------------------------------------------------
            // Wind components – clamp extreme values to domain mean
            // ----------------------------------------------------------------
            Real ux = v(i, j, k, 0);
            Real uy = v(i, j, k, 1);
            constexpr Real MAX_WIND_SQ = Real(25.0) * Real(25.0);
            if (ux * ux + uy * uy > MAX_WIND_SQ) {
                ux = ux_mean;
                uy = uy_mean;
            }
            const Real U = std::sqrt(ux * ux + uy * uy);  // wind speed [m/s]

            // ----------------------------------------------------------------
            // Per-cell or global Rothermel parameters
            // ----------------------------------------------------------------
            Real cell_R0   = R0_def;
            Real cell_C    = C_def;
            Real cell_B    = B_def;
            Real cell_bRE  = bRE_def;
            Real cell_waf  = waf_def;
            Real cell_Umax = Umax_def;

            if (use_fuel_table) {
                int code = static_cast<int>(fm(i, j, k));
                if (code >= 0 && code < fuel_table_size) {
                    const RothermelComputed& frc = fuel_table[code];
                    cell_R0   = frc.R0;
                    cell_C    = frc.C;
                    cell_B    = frc.B;
                    cell_bRE  = frc.beta_ratio_E;
                    cell_waf  = frc.waf;
                    cell_Umax = frc.U_max_ftmin;
                }
            }

            // ----------------------------------------------------------------
            // Rothermel wind factor phi_w (same computation as compute_rothermel_R)
            // ----------------------------------------------------------------
            Real U_ftmin = U * wind_conv;
            if (use_waf_flag == 1) {
                U_ftmin *= cell_waf;
            }
            if (use_wind_limit == 1 && cell_Umax > Real(0.0)) {
                U_ftmin = std::min(U_ftmin, cell_Umax);
            }
            const Real phi_w = (U_ftmin > Real(0.0))
                               ? cell_C * std::pow(U_ftmin, cell_B) * cell_bRE
                               : Real(0.0);

            // ----------------------------------------------------------------
            // Terrain slope magnitude  tan(phi) = |grad z|
            // ----------------------------------------------------------------
            Real tan_phi = Real(0.0);
            if (has_slope) {
                const Real sx = sl(i, j, k, 0);
                const Real sy = sl(i, j, k, 1);
                tan_phi = std::sqrt(sx * sx + sy * sy);
            }

            // ----------------------------------------------------------------
            // Viegas slope enhancement factor (Viegas 2004, exponential form)
            //   Phi_s_V = exp(a_V * tan(phi))
            // ----------------------------------------------------------------
            const Real Phi_s_V = std::exp(a_V * tan_phi);

            // ----------------------------------------------------------------
            // Viegas ROS [ft/min → m/s]
            //   R_V = R_0 * (1 + phi_w) * Phi_s_V
            // Uses Rothermel's R_0 (baseline) and phi_w (wind factor); only
            // the slope factor is replaced by Viegas' exponential form.
            // ----------------------------------------------------------------
            const Real R_V = std::max(cell_R0 * (Real(1.0) + phi_w) * Phi_s_V, Real(0.0))
                             * ros_conv;

            // ----------------------------------------------------------------
            // Eruptive regime flag (Viegas 2004, Section 3)
            //   Flag = 1 when slope exceeds critical value tan_phi_c
            // ----------------------------------------------------------------
            const Real eruptive_flag = (tan_phi > tan_phi_c) ? Real(1.0) : Real(0.0);

            // ----------------------------------------------------------------
            // ROS excess ratio: positive where Rothermel under-predicts
            //   excess = (R_V - R_Rothermel) / R_Rothermel
            // ----------------------------------------------------------------
            const Real R_roth = Rr(i, j, k);
            const Real ROS_excess = (R_V - R_roth) / std::max(R_roth, VIEGAS_EPS);

            // ----------------------------------------------------------------
            // Viegas flame-tilt angle (hazard assessment, NOT used for ROS)
            //   tan(alpha) = U / v_b + tan(phi)
            //   alpha      = atan(tan(alpha))
            // ----------------------------------------------------------------
            const Real v_b = std::max(v_b_def, VIEGAS_VB_MIN);
            const Real tan_alpha = U / v_b + tan_phi;
            const Real flame_tilt = std::atan(tan_alpha);

            // ----------------------------------------------------------------
            // Store outputs
            // ----------------------------------------------------------------
            vd(i, j, k, 0) = R_V;
            vd(i, j, k, 1) = eruptive_flag;
            vd(i, j, k, 2) = ROS_excess;
            vd(i, j, k, 3) = flame_tilt;
            vd(i, j, k, 4) = Phi_s_V;
        });
    }
}

#endif // VIEGAS_MODEL_H
