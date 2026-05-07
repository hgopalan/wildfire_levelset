// ============================================================================
// Compute Rothermel wind speed R and store in MultiFab
// ============================================================================
#ifndef COMPUTE_ROTHERMEL_R_H
#define COMPUTE_ROTHERMEL_R_H

#include "constants.H"
#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_Gpu.H>
#include "parse_inputs.H"
#include "rothermel_model.H"

using namespace amrex;

// ============================================================================
// Compute Rothermel rate of spread R and store in MultiFab
//
// Optional parameters for spatially-varying fuel models from a landscape file:
//   fuel_model_mf  – per-cell integer fuel codes (from landscape file)
//   fuel_table     – GPU-accessible array of RothermelComputed indexed by code
//   fuel_table_size – number of entries in fuel_table
//
// When fuel_model_mf and fuel_table are provided, each cell uses the Rothermel
// parameters (R0, C, B, beta_ratio_E, beta) corresponding to its local fuel
// code.  Codes outside [0, fuel_table_size) fall back to the global default.
// ============================================================================
static void compute_rothermel_R(MultiFab& R_mf,
                                const MultiFab& vel,
                                const Geometry& /*geom*/,
                                const RothermelParams& rp,
                                const MultiFab* slopes = nullptr,
                                const MultiFab* fuel_model_mf = nullptr,
                                const RothermelComputed* fuel_table = nullptr,
                                int fuel_table_size = 0,
                                const MultiFab* canopy_cover_mf = nullptr)
{
    // Compute Rothermel (1972) fire spread parameters (global default)
    const RothermelComputed rc = compute_rothermel_params(rp);
    const Real R0 = rc.R0;
    const Real C = rc.C;
    const Real B = rc.B;
    const Real beta_ratio_E = rc.beta_ratio_E;
    const Real phi_s_const = rc.phi_s;  // Constant phi_s for backward compatibility
    const Real beta = rc.beta;  // Packing ratio for terrain-based phi_s calculation
    const Real wind_conv = rc.wind_conv;
    const Real ros_conv = rc.ros_conv;
    // Andrews (2018) wind adjustment defaults (global fuel)
    const Real waf_default      = rc.waf;
    const Real U_max_default    = rc.U_max_ftmin;
    const int  use_waf_flag     = rp.use_waf;
    const int  use_wind_limit_flag = rp.use_wind_limit;

    const bool use_fuel_table = (fuel_model_mf != nullptr &&
                                  fuel_table != nullptr &&
                                  fuel_table_size > 0);

    // Per-cell canopy-cover-weighted WAF (Item 4):
    // When canopy_cover_mf is provided and use_waf = 1, the WAF is modulated
    // by the local canopy cover fraction:
    //   waf_cell = 1.0 - (cc/100) * (1.0 - waf_fuel)
    // cc = 0%  → waf_cell = 1.0 (open, no wind reduction)
    // cc = 100% → waf_cell = waf_fuel (fully sheltered, full reduction)
    const bool use_cc_waf = (canopy_cover_mf != nullptr && use_waf_flag == 1);

    // Compute the domain-mean wind components so that extreme wind-speed cells
    // (> 25 m/s) can be replaced with the global mean of the flow field rather
    // than only the mean of their four immediate neighbours.  The latter fails
    // when multiple adjacent cells are all extreme (e.g. near SRTM tile seams).
    const Long ncells_vel = vel.boxArray().numPts();
    const Real ux_global_mean = (ncells_vel > 0) ? vel.sum(0) / Real(ncells_vel) : Real(0.0);
    const Real uy_global_mean = (ncells_vel > 0) ? vel.sum(1) / Real(ncells_vel) : Real(0.0);

    for (MFIter mfi(R_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const v = vel.const_array(mfi);
        auto const R = R_mf.array(mfi);

        // Get slopes array if terrain file was provided
        const bool use_terrain = (slopes != nullptr);
        Array4<const Real> s;
        if (use_terrain) {
            s = slopes->const_array(mfi);
        }

        // Get per-cell fuel model codes if available
        Array4<const Real> fm;
        if (use_fuel_table) {
            fm = fuel_model_mf->const_array(mfi);
        }

        // Get per-cell canopy cover for WAF weighting (if provided)
        Array4<const Real> cc_arr;
        if (use_cc_waf) {
            cc_arr = canopy_cover_mf->const_array(mfi);
        }

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            // Horizontal wind components
            Real ux = v(i,j,k,0);
            Real uy = v(i,j,k,1);

            // Clamp extreme wind speeds (> 25 m/s) with the mean of the entire
            // flow field.  Using the global mean (rather than just the four
            // immediate neighbours) handles clusters of adjacent extreme cells
            // where a local-neighbour average would still be unreliable.
            constexpr Real MAX_WIND_SPEED_SQ = WildfireConst::MAX_WIND_SPEED_SQ;
            if (ux*ux + uy*uy > MAX_WIND_SPEED_SQ) {
                ux = ux_global_mean;
                uy = uy_global_mean;
            }

            // Per-cell or global Rothermel parameters
            Real cell_R0 = R0;
            Real cell_C  = C;
            Real cell_B  = B;
            Real cell_beta_ratio_E = beta_ratio_E;
            Real cell_beta = beta;
            Real cell_waf      = waf_default;
            Real cell_U_max    = U_max_default;

            if (use_fuel_table) {
                int code = static_cast<int>(fm(i, j, k));
                if (code >= 0 && code < fuel_table_size) {
                    const RothermelComputed& frc = fuel_table[code];
                    cell_R0          = frc.R0;
                    cell_C           = frc.C;
                    cell_B           = frc.B;
                    cell_beta_ratio_E = frc.beta_ratio_E;
                    cell_beta        = frc.beta;
                    cell_waf         = frc.waf;
                    cell_U_max       = frc.U_max_ftmin;
                }
            }

            // Wind speed in ft/min (Rothermel units)
            Real U_ftmin = std::sqrt(ux*ux + uy*uy) * wind_conv;

            // Apply Wind Adjustment Factor if enabled (Andrews 2018):
            // converts 20-ft open wind to midflame height wind.
            // When canopy_cover_mf is supplied, blend per-cell:
            //   waf_cell = 1 - (cc/100) * (1 - waf_fuel)
            if (use_waf_flag == 1) {
                Real waf_cell = cell_waf;
                if (use_cc_waf) {
                    Real cc_frac = amrex::max(Real(0.0),
                                  amrex::min(Real(100.0), cc_arr(i,j,k))) / Real(100.0);
                    waf_cell = Real(1.0) - cc_frac * (Real(1.0) - cell_waf);
                }
                U_ftmin *= waf_cell;
            }

            // Apply maximum effective wind speed limit if enabled (Andrews 2018)
            if (use_wind_limit_flag == 1 && cell_U_max > Real(0.0)) {
                U_ftmin = std::min(U_ftmin, cell_U_max);
            }

            // Wind factor φ_w (Rothermel Eq. 47)
            Real phi_w = (U_ftmin > 0.0) ? cell_C * std::pow(U_ftmin, cell_B) * cell_beta_ratio_E : 0.0;

            // Slope factor φ_s: use terrain-based if available, otherwise constant
            Real phi_s;
            if (use_terrain) {
                // Compute phi_s from terrain slopes at this cell
                // Using full Rothermel Eq. 51: φ_s = 5.275 * β^(-0.3) * tan²(φ)
                constexpr Real PHI_S_COEFFICIENT = WildfireConst::PHI_S_COEFF;
                constexpr Real BETA_EXPONENT     = WildfireConst::PHI_S_BETA_EXP;
                Real slope_x = s(i,j,k,0);
                Real slope_y = s(i,j,k,1);
                Real tan_phi_sq = slope_x * slope_x + slope_y * slope_y;
                phi_s = std::max(PHI_S_COEFFICIENT * std::pow(cell_beta, BETA_EXPONENT) * tan_phi_sq, Real(0.0));
            } else {
                // Use constant phi_s from RothermelComputed
                phi_s = phi_s_const;
            }

            // Rate of spread [ft/min] → simulation units
            R(i,j,k) = std::max(cell_R0 * (1.0 + phi_w + phi_s), 0.0) * ros_conv;
        });
    }
}

#endif // COMPUTE_ROTHERMEL_R_H
