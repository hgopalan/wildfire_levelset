// ============================================================================
// Weise & Biging (1996) Fire Whirl Model
// ============================================================================
// Implements the fire whirl characterisation model presented in:
//
//   Weise, D.R. and Biging, G.S. (1996). "Effects of wind velocity and slope
//   on flame properties." Canadian Journal of Forest Research, 26(10),
//   1849–1858. https://doi.org/10.1139/x26-210
//
// The model predicts fire whirl formation and characteristics from local
// fireline intensity, ambient wind speed, and terrain slope.  It extends
// the Byram (1959) flame-length formula with the Weise & Biging empirical
// flame-tilt relationship and then derives whirl geometry and kinematics
// from Rankine-vortex scaling.
//
// ----------------------------------------------------------------------------
// Key intermediate quantities
// ----------------------------------------------------------------------------
// Modified Froude number (Weise & Biging 1996):
//   Fr = U² / (g × L_f)
//
// Flame tilt angle θ from vertical (Weise & Biging 1996, empirical fit):
//   tan(θ) = 0.382 × Fr^0.432   [wind only]
//   tan(θ) = 0.382 × Fr^0.432 + tan(β)   [wind + slope, β = slope angle]
//
// Vertical flame height [m]:
//   H_f = L_f × cos(θ)
//
// ----------------------------------------------------------------------------
// Unique outputs written to weise_data MultiFab (6 components)
// ----------------------------------------------------------------------------
// Component 0 – weise_flame_height     : vertical flame height H_f [m]
// Component 1 – weise_flame_tilt       : flame tilt angle θ from vertical [rad]
// Component 2 – weise_whirl_height     : fire whirl column height H_w [m]
// Component 3 – weise_whirl_radius     : fire whirl core radius r_w [m]
// Component 4 – weise_angular_velocity : angular velocity Ω [rad/s]
// Component 5 – weise_max_tang_vel     : maximum tangential velocity v_θ [m/s]
//
// ----------------------------------------------------------------------------
// Fire whirl geometry (Rankine-vortex scaling, Weise & Biging 1996)
// ----------------------------------------------------------------------------
// Whirl height (matches the tilted-flame vertical extent):
//   H_w = H_f
//
// Whirl core radius (empirical proportionality, user-configurable c_r):
//   r_w = max(c_r × H_w, r_min)   where r_min = 0.01 m
//
// Circulation set by the ambient wind × whirl height:
//   Γ = U × H_w
//
// Angular velocity of the Rankine vortex core:
//   Ω = Γ / (2π × r_w²) = U × H_w / (2π × r_w²)
//
// Maximum tangential velocity at the core edge:
//   v_θ = Ω × r_w = U × H_w / (2π × r_w)
//
// All outputs are set to 0 where fireline intensity I_B ≤ I_B_min to avoid
// spurious whirl predictions in quiescent (unburned) cells.
// ============================================================================

#ifndef WEISE_BIGING_WHIRL_H
#define WEISE_BIGING_WHIRL_H

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuLaunch.H>
#include "parse_inputs.H"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace amrex;

// Number of output components in the weise_data MultiFab.
// Components (see header documentation above):
//   0 – weise_flame_height
//   1 – weise_flame_tilt
//   2 – weise_whirl_height
//   3 – weise_whirl_radius
//   4 – weise_angular_velocity
//   5 – weise_max_tang_vel
static constexpr int WEISE_NCOMP = 6;

// ============================================================================
// compute_weise_biging_whirl
//
// Fills weise_data (WEISE_NCOMP components) using the Weise & Biging (1996)
// fire whirl model.
//
// Parameters:
//   weise_data           – output MultiFab (WEISE_NCOMP components, no ghost cells)
//   fireline_intensity_mf – input: Byram fireline intensity [kW/m] (1 component)
//   flame_length_mf       – input: Byram flame length [m]          (1 component)
//   vel                   – input: wind velocity field [m/s]
//   terrain_slopes        – optional: 2-component MultiFab (slope_x, slope_y)
//                           in tan(slope) units; nullptr = flat terrain
//   wp                    – Weise & Biging model parameters
// ============================================================================
static void compute_weise_biging_whirl(
    MultiFab&       weise_data,
    const MultiFab& fireline_intensity_mf,
    const MultiFab& flame_length_mf,
    const MultiFab& vel,
    const MultiFab* terrain_slopes,
    const InputParameters::WeiseBigingParams& wp)
{
    // Physical constants
    constexpr Real g      = Real(9.81);     // gravitational acceleration [m/s²]
    constexpr Real TWO_PI = Real(2.0 * M_PI);
    constexpr Real r_min  = Real(0.01);     // minimum whirl radius to prevent division by zero [m]

    // Model parameters (captured by value for GPU lambda)
    const Real c_r    = wp.c_r;
    const Real I_min  = wp.I_B_min;

    // Empirical Weise & Biging (1996) flame-tilt coefficients:
    //   tan(θ_wind) = A_wb × Fr^B_wb
    constexpr Real A_wb = Real(0.382);
    constexpr Real B_wb = Real(0.432);

    // Domain-mean wind components (used to clamp extreme per-cell values)
    const Long ncells_vel   = vel.boxArray().numPts();
    const Real ux_mean = (ncells_vel > 0) ? vel.sum(0) / Real(ncells_vel) : Real(0.0);
    const Real uy_mean = (ncells_vel > 0) ? vel.sum(1) / Real(ncells_vel) : Real(0.0);

    const bool has_slope = (terrain_slopes != nullptr);

    for (MFIter mfi(weise_data); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const I_B = fireline_intensity_mf.const_array(mfi);
        auto const L_f_arr = flame_length_mf.const_array(mfi);
        auto const v   = vel.const_array(mfi);

        // Optional slope arrays; provide a dummy reference when unavailable.
        // The Array4 is only accessed when has_slope == true.
        Array4<Real const> sl;
        if (has_slope) sl = terrain_slopes->const_array(mfi);

        auto wd = weise_data.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            // Default: zero all outputs for non-burning cells
            wd(i, j, k, 0) = Real(0.0);
            wd(i, j, k, 1) = Real(0.0);
            wd(i, j, k, 2) = Real(0.0);
            wd(i, j, k, 3) = Real(0.0);
            wd(i, j, k, 4) = Real(0.0);
            wd(i, j, k, 5) = Real(0.0);

            // Skip cells with negligible fireline intensity
            if (I_B(i, j, k) <= I_min) return;

            // Flame length [m] from Byram formula; must be positive to proceed
            const Real L_f = L_f_arr(i, j, k);
            if (L_f <= Real(0.0)) return;

            // ------------------------------------------------------------------
            // Wind speed [m/s] – clamp extreme values (> 25 m/s) to domain mean
            // ------------------------------------------------------------------
            Real ux = v(i, j, k, 0);
            Real uy = v(i, j, k, 1);
            constexpr Real MAX_WIND_SQ = Real(25.0) * Real(25.0);
            if (ux * ux + uy * uy > MAX_WIND_SQ) {
                ux = ux_mean;
                uy = uy_mean;
            }
            const Real U = std::sqrt(ux * ux + uy * uy);

            // ------------------------------------------------------------------
            // Modified Froude number (Weise & Biging 1996)
            //   Fr = U² / (g × L_f)
            // ------------------------------------------------------------------
            const Real Fr = (U * U) / (g * L_f);

            // ------------------------------------------------------------------
            // Flame tilt angle θ from vertical
            //   tan(θ) = A_wb × Fr^B_wb           [wind component]
            //           + tan(β)                   [slope component, if any]
            // ------------------------------------------------------------------
            Real tan_theta = (Fr > Real(0.0))
                             ? A_wb * std::pow(Fr, B_wb)
                             : Real(0.0);

            if (has_slope) {
                const Real sx = sl(i, j, k, 0);  // tan(slope) in x-direction
                const Real sy = sl(i, j, k, 1);  // tan(slope) in y-direction
                // Project slope onto wind direction; fall back to magnitude when
                // wind is negligible to avoid dividing by zero.
                Real tan_beta;
                if (U > Real(1.0e-3)) {
                    tan_beta = (ux * sx + uy * sy) / U;
                } else {
                    tan_beta = std::sqrt(sx * sx + sy * sy);
                }
                tan_theta += tan_beta;
            }

            // Clamp to avoid θ ≥ 90° (numerically non-physical)
            tan_theta = amrex::min(tan_theta, Real(1.0e6));

            const Real theta = std::atan(tan_theta);   // [rad]  θ ∈ [0, π/2)

            // ------------------------------------------------------------------
            // Vertical flame height H_f [m]
            // ------------------------------------------------------------------
            const Real H_f = L_f * std::cos(theta);

            // ------------------------------------------------------------------
            // Fire whirl height H_w [m]
            // Whirl column matches the vertical flame extent (Weise & Biging).
            // ------------------------------------------------------------------
            const Real H_w = H_f;

            // ------------------------------------------------------------------
            // Fire whirl core radius r_w [m]
            //   r_w = max(c_r × H_w, r_min)
            // ------------------------------------------------------------------
            const Real r_w = amrex::max(c_r * H_w, r_min);

            // ------------------------------------------------------------------
            // Circulation (Rankine vortex)
            //   Γ = U × H_w
            // ------------------------------------------------------------------
            const Real Gamma = U * H_w;

            // ------------------------------------------------------------------
            // Angular velocity Ω [rad/s]
            //   Ω = Γ / (2π × r_w²)
            // ------------------------------------------------------------------
            const Real Omega = Gamma / (TWO_PI * r_w * r_w);

            // ------------------------------------------------------------------
            // Maximum tangential velocity v_θ [m/s]
            //   v_θ = Ω × r_w
            // ------------------------------------------------------------------
            const Real v_theta = Omega * r_w;

            // Write outputs
            wd(i, j, k, 0) = H_f;       // weise_flame_height   [m]
            wd(i, j, k, 1) = theta;     // weise_flame_tilt     [rad]
            wd(i, j, k, 2) = H_w;       // weise_whirl_height   [m]
            wd(i, j, k, 3) = r_w;       // weise_whirl_radius   [m]
            wd(i, j, k, 4) = Omega;     // weise_angular_velocity [rad/s]
            wd(i, j, k, 5) = v_theta;   // weise_max_tang_vel   [m/s]
        });
    }
}

#endif // WEISE_BIGING_WHIRL_H
