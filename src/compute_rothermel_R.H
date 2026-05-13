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
inline void compute_rothermel_R(MultiFab& R_mf,
                                const MultiFab& vel,
                                const Geometry& /*geom*/,
                                const RothermelParams& rp,
                                const MultiFab* slopes = nullptr,
                                const MultiFab* fuel_model_mf = nullptr,
                                const RothermelComputed* fuel_table = nullptr,
                                int fuel_table_size = 0,
                                const MultiFab* canopy_cover_mf = nullptr,
                                const MultiFab* canopy_height_mf = nullptr)
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
    // Slope/wind cross-term
    const int  use_sw_cross  = rc.use_slope_wind_cross;
    const Real k_sw_cross    = rc.k_slope_wind_cross;
    // FARSITE-style vectorial slope/wind combination
    const int  use_sw_vectors = rc.use_slope_wind_vectors;
    // Andrews (2018) wind adjustment defaults (global fuel)
    const Real waf_default      = rc.waf;
    const Real U_max_default    = rc.U_max_ftmin;
    const int  use_waf_flag     = rp.use_waf;
    const int  use_wind_limit_flag = rp.use_wind_limit;
    // WAF formula selection (0 = Andrews logarithmic, 1 = BehavePlus linear/exponential)
    const int  waf_formula_id   = (rp.waf_formula == "behaviorplus") ? 1 : 0;
    const Real waf_canopy_alpha = rp.waf_canopy_alpha;
    // Global fuel bed depth [ft] used as WAF baseline when no canopy height is available
    const Real delta_ft_global = rp.delta;

    const bool use_fuel_table = (fuel_model_mf != nullptr &&
                                  fuel_table != nullptr &&
                                  fuel_table_size > 0);

    // Per-cell canopy WAF:
    // When canopy_height_mf is provided and use_waf = 1, the effective midflame
    // wind height is blended from the fuel-bed depth (open areas) and the
    // canopy height (covered areas) using the local canopy cover fraction:
    //   h_eff [ft] = (1 - cc) * delta_ft + cc * canopy_height_ft
    //   waf_cell = 1.83 / ln((20 + 0.36*h_eff) / (0.13*h_eff))
    // This replaces the previous two-step approach (global WAF + cover weighting)
    // with a single, continuous per-cell WAF that responds to actual canopy height.
    // When canopy_height_mf is absent but canopy_cover_mf is present the old
    // blending formula is retained as a fallback.
    const bool use_canopy_ht_waf = (canopy_height_mf != nullptr && canopy_cover_mf != nullptr && use_waf_flag == 1);
    const bool use_cc_waf        = (!use_canopy_ht_waf && canopy_cover_mf != nullptr && use_waf_flag == 1);

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

        // Get per-cell canopy cover and height arrays for per-cell WAF
        // Safety invariant: cht_arr (and cc_arr) are accessed in the GPU kernel
        // only when use_canopy_ht_waf (or use_cc_waf) is true.  These booleans
        // are evaluated before the MFIter loop and captured by value in the
        // ParallelFor lambda; when true, the Array4 is always initialized here
        // before the kernel is launched.  Default-constructed Array4 (null-pointer
        // state) is never dereferenced because the same boolean guards both the
        // host assignment and the device access path.
        Array4<const Real> cc_arr;
        if (use_cc_waf || use_canopy_ht_waf) {
            cc_arr = canopy_cover_mf->const_array(mfi);
        }
        Array4<const Real> cht_arr;
        if (use_canopy_ht_waf) {
            cht_arr = canopy_height_mf->const_array(mfi);
        }

        const Real delta_ft = delta_ft_global;

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
                // Non-burnable landscape codes (FBFM40: 91–99; FBFM13: NB equivalent)
                // Zero ROS prevents fire from crossing water, rock, urban areas, etc.
                if (code == 91 || code == 92 || code == 93 ||
                    code == 98 || code == 99) {
                    R(i, j, k) = Real(0.0);
                    return;
                }
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

            // Apply Wind Adjustment Factor if enabled (Andrews 2018)
            if (use_waf_flag == 1) {
                Real waf_cell;
                if (use_canopy_ht_waf) {
                    // Per-cell WAF using the selected formula (Andrews or BehavePlus)
                    // based on local canopy cover fraction and canopy height.
                    //
                    // Andrews path (waf_formula_id == 0):
                    //   Unsheltered: WAF = 1.83 / ln((20+0.36*h)/(0.13*h))
                    //   Sheltered:   WAF = 0.555 / (sqrt(f_c*h_c)*ln((20+0.36*h_c)/(0.13*h_c)))
                    //
                    // BehavePlus path (waf_formula_id == 1):
                    //   Unsheltered: WAF = 0.36 + 0.004 * h_in
                    //   Sheltered:   WAF = (0.36 + 0.004*h_c_in) * exp(-alpha_c * f_c)
                    //
                    // Switch to sheltered when cc >= 0.5 and canopy height > fuel depth.
                    const Real cc_frac = amrex::max(Real(0.0),
                                         amrex::min(Real(100.0), cc_arr(i,j,k))) / Real(100.0);
                    const Real canopy_ht_ft = cht_arr(i,j,k) / Real(0.3048); // m → ft
                    waf_cell = compute_waf_per_cell(cc_frac, canopy_ht_ft, delta_ft,
                                                    Real(0.5) /* shelter_threshold */,
                                                    waf_formula_id, waf_canopy_alpha);
                } else if (use_cc_waf) {
                    // Fallback: only canopy cover fraction is available (no canopy height).
                    // Andrews path: linearly blend global WAF toward 1.0 with cover fraction.
                    // BehavePlus path: apply exponential canopy attenuation to open WAF.
                    const Real cc_frac = amrex::max(Real(0.0),
                                         amrex::min(Real(100.0), cc_arr(i,j,k))) / Real(100.0);
                    if (waf_formula_id == 1) {
                        waf_cell = compute_waf_canopy_exp(cell_waf, cc_frac, waf_canopy_alpha);
                    } else {
                        waf_cell = Real(1.0) - cc_frac * (Real(1.0) - cell_waf);
                    }
                } else {
                    waf_cell = cell_waf;
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
            // Three modes (mutually exclusive):
            //
            //  A. FARSITE vectorial (use_slope_wind_vectors = 1):
            //       φ_w_vec = φ_w × wind_unit_vector
            //       φ_s_vec = φ_s × upslope_unit_vector
            //       φ_combined = |φ_w_vec + φ_s_vec|
            //       R = R₀(1 + φ_combined)
            //     Naturally captures alignment effect: constructive when wind
            //     blows upslope, partially cancelling when it blows downslope.
            //     Requires terrain data and non-zero wind; falls back to
            //     scalar additive when either is absent.
            //
            //  B. Scalar cross-term (use_slope_wind_cross = 1):
            //       R = R₀(1 + φ_w + φ_s + k × φ_w × φ_s)
            //     Amplifies both effects when both are non-zero; always ≥ A.
            //
            //  C. Default scalar additive:
            //       R = R₀(1 + φ_w + φ_s)
            Real phi_combined;
            if (use_sw_vectors == 1 && use_terrain) {
                // FARSITE vectorial mode
                Real slope_x = s(i,j,k,0);
                Real slope_y = s(i,j,k,1);
                Real tan_phi = std::sqrt(slope_x*slope_x + slope_y*slope_y);
                Real U_mag   = std::sqrt(ux*ux + uy*uy);

                if (phi_w > Real(0.0) && tan_phi > Real(0.0)) {
                    // Both wind and slope present: true vectorial sum
                    Real wind_hat_x  = ux / U_mag;
                    Real wind_hat_y  = uy / U_mag;
                    Real slope_hat_x = slope_x / tan_phi;
                    Real slope_hat_y = slope_y / tan_phi;
                    Real comb_x = phi_w * wind_hat_x + phi_s * slope_hat_x;
                    Real comb_y = phi_w * wind_hat_y + phi_s * slope_hat_y;
                    phi_combined = std::sqrt(comb_x*comb_x + comb_y*comb_y);
                } else {
                    // One factor is zero: magnitude equals the non-zero one
                    phi_combined = phi_w + phi_s;
                }
            } else {
                // Scalar additive (default) or cross-term
                Real cross = (use_sw_cross == 1) ? k_sw_cross * phi_w * phi_s : Real(0.0);
                phi_combined = phi_w + phi_s + cross;
            }

            // Apply FARSITE-style stall floor: cells below ROS_MIN are zeroed
            const Real R_raw = std::max(cell_R0 * (Real(1.0) + phi_combined), Real(0.0)) * ros_conv;
            R(i,j,k) = (R_raw >= WildfireConst::ROS_MIN_MS) ? R_raw : Real(0.0);
        });
    }
}

#endif // COMPUTE_ROTHERMEL_R_H
