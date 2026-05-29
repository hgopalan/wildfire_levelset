// crown_ros_override.H
// Applies the active crown fire ROS override to R_mf.
//
// For cells classified as active crown fire (ecology crown_activity ≥ 2),
// replaces the surface ROS in R_mf with the active crown fire ROS from one of:
//   • Rothermel (1991) : R_crown = 3.34 × R_surface
//   • Cruz, Alexander & Wakimoto (2005) : empirical wind/CBD regression
//   • Van Wagner (1977) : R_crown = 3 / CBD × moisture_factor
// Optional passive blending (Scott & Reinhardt 2001) smoothes the transition.
//
// Extracted from main.cpp to avoid repeating the identical GPU kernel twice.
#ifndef CROWN_ROS_OVERRIDE_H
#define CROWN_ROS_OVERRIDE_H

#include "parse_inputs.H"
#include "crown_initiation.H"
#include "cruz_crown_model.H"

using namespace amrex;

// Apply the active crown fire ROS override to R_mf.
//
// Parameters:
//   R_mf                 – ROS field [m/s]; modified in-place for active crown cells
//   ecology_mf           – fire ecology diagnostics (component 3 = crown_activity)
//   vel                  – wind velocity field (for Cruz 2005 wind magnitude)
//   fireline_intensity_mf– Byram fireline intensity [kW/m] (for passive blending)
//   has_spatial_crown    – whether per-cell cbd_mf is valid
//   cbd_mf               – crown bulk density MultiFab (used when has_spatial_crown)
//   crown                – global crown parameters (CBD, FMC, CBH, flags)
//   cruz_crown           – Cruz crown parameters (MC10)
//   ccc_global           – pre-computed Cruz crown coefficients
//
// Returns the global Van Wagner crown ROS estimate [m/s] computed from the
// global CBD/FMC values.  The caller can use this to guard a dt recomputation
// (recompute only when the returned value > 0).
inline Real apply_crown_ros_override(
    MultiFab& R_mf,
    const MultiFab& ecology_mf,
    const MultiFab& vel,
    const MultiFab& fireline_intensity_mf,
    bool has_spatial_crown,
    const MultiFab& cbh_mf,
    const MultiFab& cbd_mf,
    const InputParameters::CrownInitiationParams& crown,
    const InputParameters::CruzCrownParams& cruz_crown,
    const CruzCrownComputed& /*ccc_global*/)
{
    const Real CBD_global = Real(crown.CBD);
    const Real FMC_global = Real(crown.FMC);
    const Real mf_factor  = amrex::max(Real(0.3),
                                amrex::min(Real(1.0),
                                Real(1.0) - (FMC_global - Real(100.0)) / Real(200.0)));
    // Global Van Wagner crown ROS estimate for dt guard
    const Real R_crown_g_ms = (Real(3.0) / amrex::max(CBD_global, Real(0.01)))
                              * mf_factor / Real(60.0);

    const bool use_cruz     = (crown.use_cruz_crown == 1);
    const bool use_roth1991 = (crown.use_rothermel1991_crown == 1);
    const bool use_passive  = (crown.use_passive_blend == 1);
    const Real MC10         = Real(cruz_crown.MC10);
    const Real ladder_height = Real(crown.ladder_fuel_height);
    const Real ladder_coeff  = Real(crown.ladder_fuel_coefficient);
    const Real I_o_global    = compute_critical_surface_intensity(
        Real(crown.CBH), FMC_global, ladder_height, ladder_coeff);

    for (MFIter mfi(R_mf); mfi.isValid(); ++mfi) {
        const Box& bx  = mfi.validbox();
        auto       R   = R_mf.array(mfi);
        auto const eco = ecology_mf.const_array(mfi);
        auto const v   = vel.const_array(mfi);
        auto const fi  = fireline_intensity_mf.const_array(mfi);

        const bool use_sp = has_spatial_crown;
        Array4<const Real> cbh_arr;
        Array4<const Real> cbd_arr;
        if (use_sp) {
            cbh_arr = cbh_mf.const_array(mfi);
            cbd_arr = cbd_mf.const_array(mfi);
        }

        const Real CBD_g         = CBD_global;
        const Real CBH_g         = Real(crown.CBH);
        const Real mf_val        = mf_factor;
        const bool use_roth91_k  = use_roth1991;
        const bool use_cruz_k    = use_cruz;
        const bool use_passive_k = use_passive;
        const Real FMC_k         = FMC_global;
        const Real ladder_h_k    = ladder_height;
        const Real ladder_c_k    = ladder_coeff;
        const Real I_o_k_global  = I_o_global;

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            if (eco(i, j, k, 3) < Real(1.5)) return; // surface or passive crown

            Real R_surface = R(i, j, k);
            Real R_crown_ms;

            if (use_roth91_k) {
                // Rothermel (1991): R_crown = 3.34 × R_surface
                R_crown_ms = compute_rothermel_1991_crown_ros(R_surface);
            } else if (use_cruz_k) {
                // Cruz, Alexander & Wakimoto (2005)
                Real ux = v(i, j, k, 0);
                Real uy = v(i, j, k, 1);
                Real wind_mag = std::sqrt(ux * ux + uy * uy);
                Real CBD_c = use_sp ? cbd_arr(i, j, k) : CBD_g;
                R_crown_ms = compute_crown_fire_spread_rate_cruz(wind_mag, CBD_c, MC10);
            } else {
                // Van Wagner (1977) proxy: R_crown = 3 / CBD × mf_factor
                Real CBD_c = use_sp ? cbd_arr(i, j, k) : CBD_g;
                CBD_c = amrex::max(CBD_c, Real(0.01));
                R_crown_ms = (Real(3.0) / CBD_c) * mf_val / Real(60.0);
            }

            if (use_passive_k) {
                const Real I_B_kwm = fi(i, j, k);
                const Real CBH_c = use_sp ? cbh_arr(i, j, k) : CBH_g;
                const Real I_o_k = use_sp
                    ? compute_critical_surface_intensity(CBH_c, FMC_k, ladder_h_k, ladder_c_k)
                    : I_o_k_global;
                R(i, j, k) = compute_van_wagner_passive_blend(
                    R_surface, R_crown_ms, I_B_kwm, I_o_k);
            } else {
                R(i, j, k) = amrex::max(R_surface, R_crown_ms);
            }
        });
    }

    return R_crown_g_ms;
}

#endif // CROWN_ROS_OVERRIDE_H
