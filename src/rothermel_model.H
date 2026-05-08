// ============================================================================
// Rothermel (1972) Fire Spread Model - Common Implementation
// ============================================================================
// This file provides the Rothermel fire spread model calculations that are
// shared between the level-set advection and FARSITE ellipse modules.
//
// Two computation paths are available:
//
//  Single-class (fallback):
//    Used when all per-class fuel loads in RothermelParams are zero.
//    Reproduces the original single-aggregate (w0, sigma, M_f, M_x) behaviour
//    exactly, preserving backward compatibility.
//
//  Multi-class (Rothermel 1972 full formulation):
//    Used when at least one per-class fuel load (w_d1, w_d10, w_d100, w_lh,
//    w_lw) is non-zero.  Each size class has its own load, SAV, and moisture
//    content.  Dead-class composite moisture damping and the live-fuel
//    moisture-of-extinction formula (Rothermel 1972, Eq. 88) are computed
//    explicitly.  The characteristic SAV sigma_c (load-weighted mean) drives
//    Gamma', C, B, E, xi, eps_h, and Q_ig as in the single-class path.
//
// Fixed SAV constants (Rothermel 1972):
//   10-hr dead: 109 ft⁻¹   100-hr dead: 30 ft⁻¹
// ============================================================================

#ifndef ROTHERMEL_MODEL_H
#define ROTHERMEL_MODEL_H

#include "constants.H"
#include "parse_inputs.H"
#include "fuel_database.H"
#include "andrews_model.H"
#include <AMReX.H>

// Structure to hold computed Rothermel parameters for GPU kernel
struct RothermelComputed {
    amrex::Real R0;              // No-wind, no-slope rate of spread [ft/min]
    amrex::Real I_R;             // Reaction intensity [BTU/ft²/min]
    amrex::Real beta_ratio_E;    // Packing-ratio part of wind factor
    amrex::Real C;               // Wind factor coefficient C (Eq. 47)
    amrex::Real B;               // Wind factor coefficient B (Eq. 48)
    amrex::Real phi_s;           // Slope factor (Eq. 51)
    amrex::Real beta;            // Packing ratio (for terrain-based phi_s calculation)
    amrex::Real wind_conv;       // Unit conversion for wind
    amrex::Real ros_conv;        // Unit conversion for ROS
    amrex::Real waf;             // Wind Adjustment Factor (Andrews 2018): 20-ft→midflame
    amrex::Real U_max_ftmin;     // Maximum effective wind speed [ft/min] (Andrews 2018)
};

// ============================================================================
// Compute Rothermel (1972) fire spread parameters
// ============================================================================
inline RothermelComputed compute_rothermel_params(const RothermelParams& rp)
{
    using namespace amrex;
    
    RothermelComputed rc;

    // Fixed SAV constants (ft⁻¹) for coarser dead fuel size classes
    constexpr Real SIGMA_D10  = WildfireConst::SIGMA_D10;
    constexpr Real SIGMA_D100 = WildfireConst::SIGMA_D100;

    // -----------------------------------------------------------------------
    // Determine which computation path to use
    // -----------------------------------------------------------------------
    const Real total_per_class = rp.w_d1 + rp.w_d10 + rp.w_d100
                                + rp.w_lh + rp.w_lw;
    const bool use_multiclass = (total_per_class > 0.0);

    // -----------------------------------------------------------------------
    // Characteristic SAV, bulk density, and packing ratio
    // -----------------------------------------------------------------------
    Real sigma_c;   // characteristic (load-weighted) SAV [ft⁻¹]
    Real rho_b;     // oven-dry bulk density [lb/ft³]
    Real beta;      // packing ratio [-]

    if (use_multiclass) {
        // Load-weighted mean SAV across all size classes (Rothermel 1972 Eq. 16)
        const Real w_total = total_per_class;
        sigma_c = (rp.w_d1   * rp.sigma_d1
                 + rp.w_d10  * SIGMA_D10
                 + rp.w_d100 * SIGMA_D100
                 + rp.w_lh   * rp.sigma_lh
                 + rp.w_lw   * rp.sigma_lw) / w_total;
        rho_b = w_total / rp.delta;
    } else {
        sigma_c = rp.sigma;
        rho_b   = rp.w0 / rp.delta;      // Eq. 24 (single-class)
    }
    beta = rho_b / rp.rho_p;

    // -----------------------------------------------------------------------
    // Optimum packing ratio and reaction velocity (Eqs. 36–38)
    // -----------------------------------------------------------------------
    const Real beta_op  = 3.348 * std::pow(sigma_c, -0.8189);          // Eq. 37
    const Real sigma_c15 = std::pow(sigma_c, 1.5);
    const Real Gamma_max = sigma_c15 / (495.0 + 0.0594 * sigma_c15);  // Eq. 36
    const Real A         = 133.0 * std::pow(sigma_c, -0.7913);         // Eq. 38
    const Real beta_ratio = beta / beta_op;
    const Real Gamma_prime = Gamma_max * std::pow(beta_ratio, A)
                             * std::exp(A * (1.0 - beta_ratio));       // Eq. 38

    // -----------------------------------------------------------------------
    // Moisture damping and reaction intensity
    // -----------------------------------------------------------------------
    Real I_R;

    if (use_multiclass) {
        // Per-class weighting factors: A_i = w_i * sigma_i * exp(-138/sigma_i)
        // These weight each class's contribution to the composite moisture ratio.
        const Real A_d1   = (rp.w_d1   > 0.0)
                            ? rp.w_d1   * rp.sigma_d1  * std::exp(-138.0 / rp.sigma_d1)
                            : 0.0;
        const Real A_d10  = (rp.w_d10  > 0.0)
                            ? rp.w_d10  * SIGMA_D10    * std::exp(-138.0 / SIGMA_D10)
                            : 0.0;
        const Real A_d100 = (rp.w_d100 > 0.0)
                            ? rp.w_d100 * SIGMA_D100   * std::exp(-138.0 / SIGMA_D100)
                            : 0.0;
        const Real A_lh   = (rp.w_lh   > 0.0 && rp.sigma_lh > 0.0)
                            ? rp.w_lh   * rp.sigma_lh  * std::exp(-138.0 / rp.sigma_lh)
                            : 0.0;
        const Real A_lw   = (rp.w_lw   > 0.0 && rp.sigma_lw > 0.0)
                            ? rp.w_lw   * rp.sigma_lw  * std::exp(-138.0 / rp.sigma_lw)
                            : 0.0;

        const Real A_dead = A_d1 + A_d10 + A_d100;
        const Real A_live = A_lh + A_lw;
        const Real A_tot  = A_dead + A_live;

        // Dead-class composite moisture ratio and damping coefficient
        Real eta_M_dead = 0.0;
        Real rm_dead_val = 0.0;  // saved for Eq. 88 live M'_x calculation
        if (A_dead > 0.0) {
            Real rm_dead = (A_d1 * rp.M_d1 + A_d10 * rp.M_d10 + A_d100 * rp.M_d100)
                           / (A_dead * rp.M_x);
            rm_dead = std::min(rm_dead, 1.0);
            rm_dead_val = rm_dead;
            eta_M_dead = std::max(1.0 - 2.59*rm_dead
                                      + 5.11*rm_dead*rm_dead
                                      - 3.52*rm_dead*rm_dead*rm_dead, 0.0); // Eq. 29
        }

        // Live-class moisture of extinction: Rothermel (1972) Eq. 88
        //   M'_x = 2.9 * (W'_nd / W'_nl) * (1 - M'_f / M'_x) - 0.226
        // where W'_nd and W'_nl are the dead and live net fuel loads [lb/ft²],
        // and M'_f / M'_x is the dead weighted-mean moisture ratio (= rm_dead_val).
        // Net loads are computed once here and reused below for I_R.
        const Real w_n_dead = (rp.w_d1 + rp.w_d10 + rp.w_d100) * (1.0 - rp.S_T);
        const Real w_n_live = (rp.w_lh + rp.w_lw)               * (1.0 - rp.S_T);

        Real eta_M_live = 0.0;
        if (A_live > 0.0) {
            Real M_x_live = 0.30; // minimum live moisture of extinction
            if (A_dead > 0.0 && w_n_live > 1.0e-10) {
                // Eq. 88: correct dead/live ratio uses net fuel loads, not A-factors
                const Real mx_live_calc = 2.9 * (w_n_dead / w_n_live)
                                        * (1.0 - rm_dead_val) - 0.226;
                M_x_live = std::max(mx_live_calc, 0.30);
            }
            Real rm_live = (A_lh * rp.M_lh + A_lw * rp.M_lw)
                           / (A_live * M_x_live);
            rm_live = std::min(rm_live, 1.0);
            eta_M_live = std::max(1.0 - 2.59*rm_live
                                      + 5.11*rm_live*rm_live
                                      - 3.52*rm_live*rm_live*rm_live, 0.0); // Eq. 29
        }

        // Mineral damping (Eq. 30) — same for all classes
        const Real eta_s = 0.174 * std::pow(rp.S_e, -0.19);

        // Reaction intensity: sum of dead and live contributions (Eq. 27)
        // w_n_dead and w_n_live were already computed above for Eq. 88
        I_R = Gamma_prime * (w_n_dead * rp.h_heat * eta_M_dead
                           + w_n_live * rp.h_heat * eta_M_live) * eta_s;

        // Use dead 1-hr moisture for heat of preignition (fine fuel preignition)
        const Real M_f_preignition = (rp.w_d1 > 0.0) ? rp.M_d1 : rp.M_f;
        const Real Q_ig  = 250.0 + 1116.0 * M_f_preignition;                  // Eq. 12
        const Real eps_h = std::exp(-138.0 / sigma_c);                         // Eq. 14
        const Real xi    = std::exp((0.792 + 0.681 * std::sqrt(sigma_c))
                                    * (beta + 0.1))
                           / (192.0 + 0.2595 * sigma_c);                       // Eq. 42
        rc.R0 = (A_tot > 0.0 && rho_b > 0.0 && eps_h > 0.0 && Q_ig > 0.0)
                ? I_R * xi / (rho_b * eps_h * Q_ig)
                : 0.0;                                                          // Eq. 1

    } else {
        // ----------------------------------------------------------------
        // Single-class path (original formulation, unchanged)
        // ----------------------------------------------------------------
        const Real w0     = rp.w0;
        const Real M_f    = rp.M_f;
        const Real M_x    = rp.M_x;

        // Net fuel load (Eq. 24)
        const Real w_n = w0 * (1.0 - rp.S_T);

        // Moisture damping (Eq. 29)
        const Real rm    = std::min(M_f / M_x, 1.0);
        const Real eta_M = std::max(1.0 - 2.59*rm + 5.11*rm*rm - 3.52*rm*rm*rm, 0.0);

        // Mineral damping (Eq. 30)
        const Real eta_s = 0.174 * std::pow(rp.S_e, -0.19);

        // Reaction intensity (Eq. 27)
        I_R = Gamma_prime * w_n * rp.h_heat * eta_M * eta_s;

        // Propagating flux ratio (Eq. 42)
        const Real xi    = std::exp((0.792 + 0.681 * std::sqrt(sigma_c)) * (beta + 0.1))
                           / (192.0 + 0.2595 * sigma_c);

        // Effective heating number (Eq. 14)
        const Real eps_h = std::exp(-138.0 / sigma_c);

        // Heat of preignition (Eq. 12)
        const Real Q_ig  = 250.0 + 1116.0 * M_f;

        // No-wind, no-slope rate of spread (Eq. 1)
        rc.R0 = I_R * xi / (rho_b * eps_h * Q_ig);
    }

    rc.I_R = I_R;

    // -----------------------------------------------------------------------
    // Wind factor coefficients (Eqs. 47–48) — use sigma_c for both paths
    // -----------------------------------------------------------------------
    rc.C = 7.47  * std::exp(-0.133  * std::pow(sigma_c, 0.55));
    rc.B = 0.02526 * std::pow(sigma_c, 0.54);
    const Real E = 0.715  * std::exp(-3.59e-4 * sigma_c);
    rc.beta_ratio_E = std::pow(beta_ratio, -E);

    // Slope factor (Eq. 51)
    const Real tan_phi_sq = rp.slope_x * rp.slope_x + rp.slope_y * rp.slope_y;
    rc.phi_s = std::max(5.275 * std::pow(beta, -0.3) * tan_phi_sq, 0.0);

    // Store beta for terrain-based phi_s recalculation in GPU kernel
    rc.beta = beta;

    // Unit-conversion scalars for the GPU kernel
    rc.wind_conv = rp.wind_conv;
    rc.ros_conv  = rp.ros_conv;

    // -----------------------------------------------------------------------
    // Andrews (2018) wind adjustments
    // -----------------------------------------------------------------------
    // Wind Adjustment Factor: converts 20-ft open wind to midflame wind
    rc.waf = compute_waf_unsheltered(rp.delta);
    // Maximum effective wind speed: prevents unrealistically large phi_w
    rc.U_max_ftmin = compute_max_ews_ftmin(I_R, rc.C, rc.B, rc.beta_ratio_E);

    return rc;
}

// ============================================================================
// Build a per-code lookup table of precomputed Rothermel parameters
// ============================================================================
// Creates a flat array of RothermelComputed indexed by fuel integer code, so
// GPU kernels can look up per-cell Rothermel parameters by fuel code.
//
// The table size is chosen to accommodate all valid codes:
//   FBFM13: indices 0–13   (size 14)
//   FBFM40: indices 0–204  (size 205)
//
// Entries for non-burnable or unknown codes are set to the default values
// derived from *base_params* (the global rothermel parameter set).
//
// The moisture contents (M_f, M_d1, …, M_lw), wind conversion, and ROS
// conversion are taken from *base_params* and shared across all per-code
// entries; only the structural fuel parameters (w0, sigma, delta, M_x,
// h_heat, S_T, S_e, rho_p, and per-class loads/SAVs) are varied per code.
// ============================================================================

inline std::vector<RothermelComputed> build_fuel_rothermel_table(
    const RothermelParams& base_params,
    const std::string& fuel_type)
{
    const int table_size = (fuel_type == "40") ? 205 : 14;
    std::vector<RothermelComputed> table(table_size);

    // Fill with default (base_params) values first
    RothermelComputed default_rc = compute_rothermel_params(base_params);
    for (auto& entry : table) entry = default_rc;

    // Override entries for each known burnable fuel code
    for (int code = 0; code < table_size; ++code) {
        FuelModel model;
        if (lookup_fuel_model_by_int_code(code, fuel_type, model)) {
            RothermelParams rp = base_params;
            apply_fuel_model(model, rp);
            // Per-class moistures remain as specified in base_params so that
            // the user-supplied moisture values apply to all fuel codes.
            table[code] = compute_rothermel_params(rp);
        }
    }

    return table;
}

#endif // ROTHERMEL_MODEL_H
