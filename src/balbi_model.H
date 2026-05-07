// ============================================================================
// Balbi (2009) Physical Fire Spread Model
// ============================================================================
// Implements the physical rate-of-spread model introduced in:
//
//   Balbi, J.-H., Rossi, J.-L., Marcelli, T., and Santoni, P.-A. (2009).
//   "A physical model for wildland fires."
//   Combustion and Flame, 156(12), 2217–2230.
//   https://doi.org/10.1016/j.combustflame.2009.07.010
//
// The model derives the rate of spread from the radiant energy flux received
// by unburned fuel from a tilted flame front. Key physical mechanisms:
//   1. Radiation fraction χ – fraction of combustion heat radiated forward
//   2. Ignition energy B* – dimensionless energy needed to ignite the fuel
//   3. Buoyancy velocity v_b – vertical flame velocity set by thermal forcing
//   4. Flame tilt angle α – from cross-wind (U) and terrain slope (θ)
//
// Rate-of-spread formula:
//   R = A_coeff · (1 + sin α − cos α)
//
// where
//   A_coeff = χ · σ_m · δ_m / (2 · τ₀ · B*)    [m/s]
//   χ       = r₀₀ · σ_m / (1 + r₀₀ · σ_m)       [−]
//   B*      = (C_pf·(T_i − T_a) + M_f·Λ) / h     [−]
//   tan α   = U / v_b + tan θ
//   v_b     = sqrt(g · δ_m · (T_f − T_a) / T_a)  [m/s]
//
// Parameters in SI units (fuel parameters converted from US Rothermel units):
//   σ_m  [m⁻¹] – surface-area-to-volume ratio  (= σ_ft × 3.28084)
//   δ_m  [m]   – fuel bed depth                (= δ_ft × 0.3048)
//   h    [J/kg] – heat of combustion            (= h_BTU/lb × 2326)
//   M_f  [−]   – fuel moisture content (fraction)
//   T_a  [K]   – ambient temperature            (parmparse: balbi.T_a, default 300)
//   T_f  [K]   – mean flame temperature         (parmparse: balbi.T_f, default 1000)
//   T_i  [K]   – ignition temperature           (parmparse: balbi.T_i, default 600)
//   Λ    [J/kg] – latent heat of vaporisation   (parmparse: balbi.delta_H, default 2.26e6)
//   C_pf [J/(kg·K)] – specific heat of dry fuel (parmparse: balbi.C_pf, default 1800)
//   r₀₀  [m]   – radiation length scale         (parmparse: balbi.r_00, default 2.5e-4)
//   τ₀   [s/m] – residence-time coefficient     (parmparse: balbi.tau_0, default 75591)
//
// Gravity g = 9.81 m/s² is embedded as a physical constant.
// ============================================================================

#ifndef BALBI_MODEL_H
#define BALBI_MODEL_H

#include "constants.H"
#include "parse_inputs.H"
#include "fuel_database.H"
#include <AMReX.H>
#include <cmath>
#include <iomanip>
#include <vector>

// ============================================================================
// Pre-computed Balbi parameters for GPU kernel (one entry per fuel code)
// ============================================================================
struct BalbiComputed {
    amrex::Real A_coeff;  // ROS amplitude coefficient A_coeff [m/s]
    amrex::Real v_b;      // Buoyancy velocity scale v_b [m/s]
};

// ============================================================================
// Unit-conversion constants (US/Rothermel → SI)
// Aliases into WildfireConst for use with 'using namespace BalbiUnits'.
// ============================================================================
namespace BalbiUnits {
    inline constexpr amrex::Real FT_TO_M          = WildfireConst::FT_TO_M;
    inline constexpr amrex::Real FT_INV_TO_M_INV  = WildfireConst::FT_INV_TO_M_INV;
    inline constexpr amrex::Real LB_FT2_TO_KG_M2  = WildfireConst::LB_FT2_TO_KG_M2;
    inline constexpr amrex::Real LB_FT3_TO_KG_M3  = WildfireConst::LB_FT3_TO_KG_M3;
    inline constexpr amrex::Real BTU_LB_TO_J_KG   = WildfireConst::BTU_LB_TO_J_KG;
    inline constexpr amrex::Real G                = WildfireConst::GRAVITY;
}

// Minimum safe buoyancy velocity [m/s] used as a numerical floor.
// Prevents division-by-zero in the GPU tilt-angle kernel while remaining
// orders of magnitude below physically meaningful v_b values (≥ 1 m/s).
inline constexpr amrex::Real BALBI_V_B_MIN    = WildfireConst::BALBI_V_B_MIN;

// Minimum dimensionless ignition-energy ratio B* before the cell is treated
// as non-burnable (A_coeff = 0).  Numerically, 1e-6 is well below any
// physically plausible value for B* (typically 0.01–0.10).
inline constexpr amrex::Real BALBI_B_STAR_MIN = WildfireConst::BALBI_B_STAR_MIN;

// ============================================================================
// Compute Balbi (2009) model coefficients from fuel + thermal parameters
//
// Inputs:
//   rp – Rothermel fuel parameters (provides σ, δ, w₀, ρ_p, h_heat, M_f)
//   bp – Balbi thermal parameters  (provides T_a, T_f, T_i, Λ, C_pf, r₀₀, τ₀)
//
// Returns: BalbiComputed with A_coeff and v_b ready for the GPU kernel.
// ============================================================================
inline BalbiComputed compute_balbi_params(const RothermelParams& rp,
                                          const InputParameters::BalbiParams& bp)
{
    using namespace amrex;
    using namespace BalbiUnits;

    // -----------------------------------------------------------------------
    // Convert fuel parameters from US/Rothermel units to SI
    // -----------------------------------------------------------------------
    const Real sigma_m  = rp.sigma  * FT_INV_TO_M_INV;  // [m⁻¹]
    const Real delta_m  = rp.delta  * FT_TO_M;           // [m]
    const Real h_si     = rp.h_heat * BTU_LB_TO_J_KG;    // [J/kg]

    // -----------------------------------------------------------------------
    // Guard against degenerate fuel beds
    // -----------------------------------------------------------------------
    if (delta_m <= 0.0 || sigma_m <= 0.0 || h_si <= 0.0) {
        BalbiComputed bc;
        bc.A_coeff = 0.0;
        bc.v_b     = 1.0;   // avoid div-by-zero in GPU kernel
        return bc;
    }

    // -----------------------------------------------------------------------
    // Radiation fraction χ
    //   χ₀ = r₀₀ · σ_m        (r₀₀ [m] × σ_m [m⁻¹] → dimensionless)
    //   χ  = χ₀ / (1 + χ₀)
    // -----------------------------------------------------------------------
    const Real chi_0 = bp.r_00 * sigma_m;
    const Real chi   = chi_0 / (1.0 + chi_0);

    // -----------------------------------------------------------------------
    // Dimensionless ignition energy B*
    //   E_ig = C_pf·(T_i − T_a) + M_f·Λ      [J/kg]
    //   B*   = E_ig / h                        [−]
    // -----------------------------------------------------------------------
    const Real E_ig  = bp.C_pf * (bp.T_i - bp.T_a) + rp.M_f * bp.delta_H;
    const Real B_star = E_ig / h_si;

    // -----------------------------------------------------------------------
    // Buoyancy velocity scale
    //   v_b = sqrt(g · δ_m · (T_f − T_a) / T_a)   [m/s]
    // -----------------------------------------------------------------------
    const Real dT_f = bp.T_f - bp.T_a;
    const Real v_b  = (dT_f > 0.0)
                      ? std::sqrt(G * delta_m * dT_f / bp.T_a)
                      : BALBI_V_B_MIN;  // non-positive temperature delta: use floor

    // -----------------------------------------------------------------------
    // ROS amplitude coefficient
    //   A_coeff = χ · σ_m · δ_m / (2 · τ₀ · B*)   [m/s]
    //   (σ_m·δ_m dimensionless; τ₀ in [s/m] gives [m/s] after division)
    // -----------------------------------------------------------------------
    const Real A_coeff = (B_star > BALBI_B_STAR_MIN)
                         ? chi * sigma_m * delta_m / (2.0 * bp.tau_0 * B_star)
                         : 0.0;

    BalbiComputed bc;
    bc.A_coeff = std::max(A_coeff, 0.0);
    bc.v_b     = std::max(v_b, BALBI_V_B_MIN);  // floor to prevent divide-by-zero
    return bc;
}

// ============================================================================
// Build per-fuel-code lookup table of BalbiComputed entries
//
// Mirrors build_fuel_rothermel_table() but uses the Balbi model.
// Table is indexed by integer fuel code; non-burnable codes receive the
// global default entry computed from *base_params*.
// ============================================================================
inline std::vector<BalbiComputed> build_fuel_balbi_table(
    const RothermelParams& base_params,
    const InputParameters::BalbiParams& bp,
    const std::string& fuel_type)
{
    const int table_size = (fuel_type == "40") ? 205 : 14;
    std::vector<BalbiComputed> table(table_size);

    // Fill with global default first
    BalbiComputed default_bc = compute_balbi_params(base_params, bp);
    for (auto& entry : table) entry = default_bc;

    // Override entries for each known burnable fuel code
    for (int code = 0; code < table_size; ++code) {
        FuelModel model;
        if (lookup_fuel_model_by_int_code(code, fuel_type, model)) {
            RothermelParams rp = base_params;
            apply_fuel_model(model, rp);
            table[code] = compute_balbi_params(rp, bp);
        }
    }

    return table;
}

// ============================================================================
// Print Balbi fuel parameter table for all fuel codes in the database
//
// Outputs a human-readable table showing:
//   Fuel | Name | sigma_m | delta_m | chi | B_star | v_b | A_coeff | R@5m/s
//
// This auto-generates the Balbi fuel parameter reference table.
// ============================================================================
inline void print_balbi_fuel_table(const RothermelParams& base_params,
                                   const InputParameters::BalbiParams& bp,
                                   const std::string& fuel_type)
{
    using namespace amrex;
    using namespace BalbiUnits;

    const bool is40 = (fuel_type == "40");
    const int  table_size = is40 ? 205 : 14;

    Print() << "\n";
    Print() << "==========================================================================\n";
    Print() << " Balbi (2009) Fuel Parameter Table"
            << " (" << (is40 ? "FBFM40" : "FBFM13") << ")\n";
    Print() << "==========================================================================\n";
    Print() << " Balbi thermal parameters:\n";
    Print() << "   T_a=" << bp.T_a << " K  T_f=" << bp.T_f
            << " K  T_i=" << bp.T_i << " K\n";
    Print() << "   C_pf=" << bp.C_pf << " J/(kg·K)"
            << "  r_00=" << bp.r_00 << " m"
            << "  tau_0=" << bp.tau_0 << " s/m"
            << "  delta_H=" << bp.delta_H << " J/kg\n";
    Print() << "--------------------------------------------------------------------------\n";
    Print() << std::left
            << std::setw(6)  << "Code"
            << std::setw(8)  << "Name"
            << std::setw(10) << "sig[m⁻¹]"
            << std::setw(9)  << "del[m]"
            << std::setw(8)  << "chi[-]"
            << std::setw(9)  << "B*[-]"
            << std::setw(9)  << "v_b[m/s]"
            << std::setw(11) << "A[m/s]"
            << std::setw(11) << "R@5m/s[m/s]"
            << "\n";
    Print() << "--------------------------------------------------------------------------\n";

    for (int code = 0; code < table_size; ++code) {
        FuelModel model;
        if (!lookup_fuel_model_by_int_code(code, fuel_type, model)) continue;

        RothermelParams rp = base_params;
        apply_fuel_model(model, rp);

        // Derived SI quantities
        const Real sigma_m  = rp.sigma  * FT_INV_TO_M_INV;
        const Real delta_m  = rp.delta  * FT_TO_M;
        const Real h_si     = rp.h_heat * BTU_LB_TO_J_KG;

        const Real chi_0  = bp.r_00 * sigma_m;
        const Real chi    = chi_0 / (1.0 + chi_0);
        const Real E_ig   = bp.C_pf*(bp.T_i - bp.T_a) + rp.M_f * bp.delta_H;
        const Real B_star = (h_si > 0.0) ? E_ig / h_si : 0.0;
        const Real dT_f   = bp.T_f - bp.T_a;
        const Real v_b    = (dT_f > 0.0 && delta_m > 0.0)
                            ? std::sqrt(G * delta_m * dT_f / bp.T_a)
                            : 0.0;
        const Real A_coeff = (B_star > BALBI_B_STAR_MIN && sigma_m > 0.0 && delta_m > 0.0)
                             ? chi * sigma_m * delta_m / (2.0 * bp.tau_0 * B_star)
                             : 0.0;

        // Predicted ROS at U=5 m/s (typical moderate wind), flat terrain
        const Real vb_safe = std::max(v_b, 1e-3);
        const Real tan_a_5 = 5.0 / vb_safe;
        const Real denom_5 = std::sqrt(1.0 + tan_a_5*tan_a_5);
        const Real sin_a   = tan_a_5 / denom_5;
        const Real cos_a   = 1.0 / denom_5;
        const Real R_5ms   = A_coeff * std::max(1.0 + sin_a - cos_a, 0.0);

        Print() << std::left << std::fixed << std::setprecision(3)
                << std::setw(6)  << code
                << std::setw(8)  << model.name.substr(0, 7)
                << std::setw(10) << sigma_m
                << std::setw(9)  << delta_m
                << std::setw(8)  << chi
                << std::setw(9)  << B_star
                << std::setw(9)  << v_b
                << std::setw(11) << A_coeff
                << std::setw(11) << R_5ms
                << "\n";
    }
    Print() << "==========================================================================\n\n";
}

#endif // BALBI_MODEL_H
