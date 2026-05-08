// ============================================================================
// Fire Ecology Diagnostics
// ============================================================================
// Computes six per-cell fire ecology and behavior fields written to every
// plotfile as additional diagnostic outputs:
//
//   Component 0 - scorch_height      [m]   Van Wagner (1973)
//   Component 1 - prob_ignition      [-]   Anderson (1970) / Rothermel (1983)
//   Component 2 - tree_mortality     [-]   Ryan-Reinhardt (1988) style logistic
//   Component 3 - crown_activity     [-]   0=surface, 1=passive, 2=active
//                                          Van Wagner (1977) / Scott & Reinhardt (2001)
//   Component 4 - torching_ratio     [-]   I_B / I_o  (Scott & Reinhardt 2001 TI proxy)
//                                          ≥ 1 → passive / active crown fire possible
//   Component 5 - crowning_ratio     [-]   R / R'_SA  (Scott & Reinhardt 2001 CI proxy)
//                                          ≥ 1 → active crown fire possible
//
// ----------------------------------------------------------------------------
// Scott & Reinhardt (2001) Torching Index / Crowning Index proxies
// ----------------------------------------------------------------------------
// The original TI and CI are defined as the minimum open (20-ft) wind speed
// at which torching or active crowning begins.  Computing them exactly requires
// inverting the Rothermel ROS function with respect to wind speed, which is
// numerically expensive per cell.
//
// Instead this module computes dimensionless ratios that carry the same
// diagnostic information:
//
//   torching_ratio  = I_B / I_o
//     I_B = Byram fireline intensity [kW/m]
//     I_o = Van Wagner (1977) crown initiation threshold [kW/m]
//         = 0.010 * CBH * (460 + 25.9 * FMC)
//   torching_ratio ≥ 1.0  ↔  torching / passive crown fire conditions met
//
//   crowning_ratio  = R / R'_SA
//     R    = current surface/crown ROS [m/s]
//     R'SA = Van Wagner (1977) critical active-crown ROS [m/s]
//          = 3.0 / CBD / 60
//   crowning_ratio ≥ 1.0  ↔  active crown fire conditions met
//
// Reference:
//   Scott, J.H. & Reinhardt, E.D. (2001). "Assessing Crown Fire Potential by
//     Linking Models of Surface and Crown Fire Behavior." USDA Forest Service
//     Research Paper RMRS-RP-29.
//
// ----------------------------------------------------------------------------
// Scorch height (Van Wagner 1973)
// ----------------------------------------------------------------------------
//   H_s [m] = 0.1483 * I_B [kW/m]^(2/3)
//
// Empirical formula for flame-induced scorch height under standard ambient
// conditions (Scott & Reinhardt 2001).  I_B is the Byram (1959) fireline
// intensity already computed in compute_fire_behavior.H.
//
// ----------------------------------------------------------------------------
// Probability of ignition (Anderson 1970 / Rothermel 1983)
// ----------------------------------------------------------------------------
//   P_i = min(1, max(0, 0.000048 * T_fuel_F^1.4 * exp(-0.07 * MC%)))
//
// Probability that a firebrand landing on fine fuel will cause sustained
// ignition.  T_fuel_F is fine fuel temperature in °F (ambient + solar
// increment, default +25°F).  MC% is 1-hr dead fuel moisture content in %.
//
// ----------------------------------------------------------------------------
// Tree mortality fraction (Ryan & Reinhardt 1988)
// ----------------------------------------------------------------------------
//   CRS = (H_s - CBH) / (H_tree - CBH)   clipped to [0, 1]
//   P_kill = 1 / (1 + exp(-4 * (CRS - 0.6)))   when H_s > CBH
//          = 0                                   when H_s <= CBH
//
// Crown Scorch Ratio (CRS) drives a logistic mortality curve calibrated to
// give ~50% mortality when 60% of the live crown is scorched.  Uses the same
// CBH from crown initiation parameters.
//
// ----------------------------------------------------------------------------
// Crown activity class (Van Wagner 1977 / Scott & Reinhardt 2001)
// ----------------------------------------------------------------------------
//   0 - surface fire     : I_B < I_o
//   1 - passive crown    : I_B >= I_o  AND  R < R'_SA
//   2 - active crown     : I_B >= I_o  AND  R >= R'_SA
//
// Crown initiation intensity (Van Wagner 1977, as in crown_initiation.H):
//   I_o [kW/m] = 0.010 * CBH * (460 + 25.9 * FMC)
//
// Critical active-crown ROS (Van Wagner 1977):
//   R'_SA [m/s] = 3.0 / CBD / 60
//
// ----------------------------------------------------------------------------
// References
// ----------------------------------------------------------------------------
//   Van Wagner, C.E. (1973). "Height of crown scorch in forest fires."
//     Canadian Journal of Forest Research, 3(3), 373-378.
//   Anderson, H.E. (1970). "Forest fuel ignitibility."  Fire Technology, 6(4).
//   Rothermel, R.C. (1983).  USDA Forest Service GTR INT-143.
//   Ryan, K.C. & Reinhardt, E.D. (1988). "Predicting postfire mortality of
//     seven western conifers." Canadian Journal of Forest Research, 18, 1291-1297.
//   Van Wagner, C.E. (1977). "Conditions for the start and spread of crown fire."
//     Canadian Journal of Forest Research, 7(1), 23-34.
//   Scott, J.H. & Reinhardt, E.D. (2001). "Assessing crown fire potential by
//     linking models of surface and crown fire behavior." USDA Forest Service
//     Research Paper RMRS-RP-29.
// ============================================================================

#ifndef FIRE_ECOLOGY_MODEL_H
#define FIRE_ECOLOGY_MODEL_H

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_GpuLaunch.H>
#include "parse_inputs.H"
#include "rothermel_model.H"
#include "constants.H"
#include <cmath>

using namespace amrex;

// Number of components in the ecology MultiFab
constexpr int ECOLOGY_NCOMP = 6;

// ============================================================================
// compute_fire_ecology
//
// Parameters:
//   ecology_mf             - output: ECOLOGY_NCOMP-component MultiFab
//   fireline_intensity_mf  - input:  Byram fireline intensity [kW/m]
//   R_mf                   - input:  rate of spread [m/s]
//   rp                     - Rothermel parameters (provides M_f for P_i)
//   fep                    - fire ecology parameters
//   cp                     - crown initiation parameters (CBH, CBD, FMC,
//                            tree_height for mortality; crown.enable not
//                            checked here — ecology is always computed)
// ============================================================================
inline void compute_fire_ecology(
    MultiFab&                                    ecology_mf,
    const MultiFab&                              fireline_intensity_mf,
    const MultiFab&                              R_mf,
    const RothermelParams&                       rp,
    const InputParameters::FireEcologyParams&    fep,
    const InputParameters::CrownInitiationParams& cp)
{
    // -----------------------------------------------------------------------
    // Ambient and fuel temperature (°F) for probability of ignition
    // -----------------------------------------------------------------------
    // T_fuel = ambient [°C] → °F + solar heating increment [°F]
    const Real T_a_F     = Real(fep.T_a_C) * Real(1.8) + Real(32.0);
    const Real T_fuel_F  = T_a_F + Real(fep.solar_heating_F);
    // Anderson (1970) formula uses T_fuel^1.4 - precompute for GPU
    const Real T_fuel_14 = (T_fuel_F > Real(0.0))
                           ? std::pow(T_fuel_F, Real(1.4))
                           : Real(0.0);
    // 1-hr dead fuel moisture [%]
    const Real MC_pct = rp.M_f * Real(100.0);

    // -----------------------------------------------------------------------
    // Crown activity thresholds
    // -----------------------------------------------------------------------
    const Real CBH = Real(cp.CBH);   // canopy base height [m]
    const Real CBD = Real(cp.CBD);   // canopy bulk density [kg/m³]
    const Real FMC = Real(cp.FMC);   // foliar moisture content [%]

    // Crown initiation intensity [kW/m] (Van Wagner 1977, consistent with
    // the existing crown_initiation.H implementation)
    const Real I_o = Real(0.010) * CBH * (Real(460.0) + Real(25.9) * FMC);

    // Critical active crown fire ROS [m/s] (Van Wagner 1977)
    // R'_SA [m/min] = 3.0 / CBD  →  [m/s] = 3.0 / CBD / 60
    const Real R_SA = Real(3.0) / amrex::max(CBD, Real(0.01)) / Real(60.0);

    // -----------------------------------------------------------------------
    // Tree mortality parameters
    // -----------------------------------------------------------------------
    const Real H_tree = amrex::max(Real(fep.tree_height),
                                   CBH + Real(0.1));  // avoid degenerate case

    // -----------------------------------------------------------------------
    // GPU-kernel locals (scalars captured by value)
    // -----------------------------------------------------------------------
    for (MFIter mfi(ecology_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const fi  = fireline_intensity_mf.const_array(mfi);
        auto const R   = R_mf.const_array(mfi);
        auto       eco = ecology_mf.array(mfi);

        const Real Tf14 = T_fuel_14;
        const Real MC   = MC_pct;
        const Real Io   = I_o;
        const Real Rsa  = R_SA;
        const Real Ht   = H_tree;
        const Real Cbh  = CBH;

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const Real IB = fi(i, j, k);
            const Real Ri = R(i, j, k);

            // ----------------------------------------------------------------
            // (0) Scorch height [m]  -  Van Wagner (1973)
            //     H_s = 0.1483 * I_B^(2/3)   for I_B in kW/m
            // ----------------------------------------------------------------
            const Real H_s = (IB > Real(0.0))
                             ? Real(0.1483) * std::pow(IB, Real(2.0 / 3.0))
                             : Real(0.0);
            eco(i, j, k, 0) = H_s;

            // ----------------------------------------------------------------
            // (1) Probability of ignition [-]  -  Anderson (1970)
            //     P_i = 0.000048 * T_fuel_F^1.4 * exp(-0.07 * MC%)
            // ----------------------------------------------------------------
            const Real P_i = amrex::min(Real(1.0),
                             amrex::max(Real(0.0),
                             Real(0.000048) * Tf14
                             * std::exp(Real(-0.07) * MC)));
            eco(i, j, k, 1) = P_i;

            // ----------------------------------------------------------------
            // (2) Tree mortality fraction [-]  -  Ryan-Reinhardt (1988) style
            //     CRS = (H_s - CBH) / (H_tree - CBH)  clipped [0,1]
            //     P_kill = 1 / (1 + exp(-4 * (CRS - 0.6)))
            // ----------------------------------------------------------------
            Real P_kill = Real(0.0);
            if (H_s > Cbh) {
                const Real CRS = amrex::max(Real(0.0),
                                 amrex::min(Real(1.0),
                                 (H_s - Cbh) / amrex::max(Ht - Cbh, Real(0.1))));
                P_kill = Real(1.0) / (Real(1.0)
                         + std::exp(Real(-4.0) * (CRS - Real(0.6))));
            }
            eco(i, j, k, 2) = P_kill;

            // ----------------------------------------------------------------
            // (3) Crown activity class [-]
            //     0 = surface fire     (I_B < I_o)
            //     1 = passive crown    (I_B >= I_o, R < R'_SA)
            //     2 = active crown     (I_B >= I_o, R >= R'_SA)
            // ----------------------------------------------------------------
            Real activity = Real(0.0);
            if (IB >= Io) {
                activity = (Ri >= Rsa) ? Real(2.0) : Real(1.0);
            }
            eco(i, j, k, 3) = activity;

            // ----------------------------------------------------------------
            // (4) Torching ratio [-]  -  Scott & Reinhardt (2001) TI proxy
            //     torching_ratio = I_B / I_o
            //     ≥ 1.0 → passive or active crown fire possible
            //     < 1.0 → surface fire only
            // ----------------------------------------------------------------
            eco(i, j, k, 4) = IB / amrex::max(Io, Real(1.0e-6));

            // ----------------------------------------------------------------
            // (5) Crowning ratio [-]  -  Scott & Reinhardt (2001) CI proxy
            //     crowning_ratio = R / R'_SA
            //     ≥ 1.0 → active crown fire conditions met
            //     < 1.0 → passive crown fire at most
            // ----------------------------------------------------------------
            eco(i, j, k, 5) = Ri / amrex::max(Rsa, Real(1.0e-8));
        });
    }
}

// ============================================================================
// apply_ecology_p_ignition_feedback
//
// Scales the rate-of-spread field R_mf by the local probability of ignition
// (ecology component 1) when fire_ecology.couple_to_ros == 1.
//
// Physical motivation:
//   P_ignition (Anderson 1970) gives the probability that a firebrand or
//   advancing flame tip successfully initiates sustained combustion at a given
//   point.  Multiplying the ROS by this probability captures the observation
//   that fire spreads more slowly through wet or cold fuel where ignition is
//   uncertain.  A floor value (p_ignition_floor) prevents complete quenching
//   of the fire front in marginally-burnable fuels.
//
// Only cells that have not yet burned (phi >= 0) are modified; already-burned
// cells (phi < 0) retain their pre-computed R value unchanged.
//
// Parameters:
//   R_mf        – rate of spread [m/s], modified in-place
//   ecology_mf  – ECOLOGY_NCOMP ecology fields (component 1 = P_ignition)
//   phi         – level-set field (phi < 0 = burned)
//   fep         – fire ecology parameters (couple_to_ros, p_ignition_floor)
// ============================================================================
inline void apply_ecology_p_ignition_feedback(
    MultiFab&                                    R_mf,
    const MultiFab&                              ecology_mf,
    const MultiFab&                              phi,
    const InputParameters::FireEcologyParams&    fep)
{
    if (fep.couple_to_ros != 1) return;

    const Real p_floor = Real(fep.p_ignition_floor);

    for (MFIter mfi(R_mf); mfi.isValid(); ++mfi) {
        const Box& bx  = mfi.validbox();
        auto       R   = R_mf.array(mfi);
        auto const eco = ecology_mf.const_array(mfi);
        auto const p   = phi.const_array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            // Do not scale ROS in already-burned cells
            if (p(i, j, k) < Real(0.0)) return;

            // Component 1 of ecology_mf is P_ignition ∈ [0, 1]
            const Real P_i = eco(i, j, k, 1);
            const Real scale = amrex::max(P_i, p_floor);
            R(i, j, k) *= scale;
        });
    }
}

// ============================================================================
// Scott & Reinhardt (2001) full Torching Index (TI) and Crowning Index (CI)
// ============================================================================
// The Torching Index is the minimum 20-ft open wind speed [km/h] at which
// the surface fire intensity equals the Van Wagner (1977) crown initiation
// threshold I_o.  The Crowning Index is the minimum 20-ft wind speed at
// which the surface ROS equals the Van Wagner (1977) critical active-crown
// ROS R'_SA.
//
// Both are computed by bisecting over U_20 using the Rothermel formula.
//
// Reference:
//   Scott, J.H. & Reinhardt, E.D. (2001). "Assessing Crown Fire Potential by
//   Linking Models of Surface and Crown Fire Behavior." USDA Forest Service
//   Research Paper RMRS-RP-29.
// ============================================================================

// Returns the minimum 20-ft wind speed [km/h] at which I_B = I_o
// Returns -1.0 if I_o is never reached within [0, U_max_kmh]
inline amrex::Real compute_torching_index_kmh(
    const RothermelComputed& rc,
    const RothermelParams& rp,
    amrex::Real I_o_kwm,
    amrex::Real h_kJ_kg,
    amrex::Real w_net_kg_m2,
    amrex::Real U_max_kmh = 160.0,
    int max_iter = 50)
{
    using namespace amrex;
    if (I_o_kwm <= Real(0.0)) return Real(0.0);

    // Convert max wind to ft/min for Rothermel
    // U [km/h] -> [ft/min]: * 3.28084 * 1000 / 60 = * 54.6807
    const Real U_conv = Real(54.6807);
    const Real waf_val = rc.waf;

    // Bisection over U_20 [km/h]
    Real lo = Real(0.0);
    Real hi = U_max_kmh;

    // Evaluate I_B at hi
    auto eval_IB = [&](Real U_kmh) -> Real {
        Real U_ftmin = U_kmh * U_conv * waf_val;
        U_ftmin = amrex::max(U_ftmin, Real(0.0));
        if (rc.U_max_ftmin > Real(0.0))
            U_ftmin = amrex::min(U_ftmin, rc.U_max_ftmin);
        const Real phi_w = rc.C * std::pow(U_ftmin, rc.B) * rc.beta_ratio_E;
        const Real R_ftmin = rc.R0 * (Real(1.0) + phi_w + rc.phi_s);
        const Real R_ms = amrex::max(R_ftmin * rp.ros_conv, Real(0.0));
        return h_kJ_kg * w_net_kg_m2 * R_ms;  // kW/m
    };

    const Real IB_hi = eval_IB(hi);
    if (IB_hi < I_o_kwm) return Real(-1.0);  // never reached

    for (int iter = 0; iter < max_iter; ++iter) {
        const Real mid = Real(0.5) * (lo + hi);
        const Real IB_mid = eval_IB(mid);
        if (IB_mid < I_o_kwm) lo = mid;
        else                   hi = mid;
        if (hi - lo < Real(0.01)) break;
    }
    return Real(0.5) * (lo + hi);
}

// Returns the minimum 20-ft wind speed [km/h] at which R(U) >= R_SA_ms
// Returns -1.0 if never reached
inline amrex::Real compute_crowning_index_kmh(
    const RothermelComputed& rc,
    const RothermelParams& rp,
    amrex::Real R_SA_ms,
    amrex::Real U_max_kmh = 160.0,
    int max_iter = 50)
{
    using namespace amrex;
    if (R_SA_ms <= Real(0.0)) return Real(0.0);

    const Real U_conv  = Real(54.6807);
    const Real waf_val = rc.waf;

    auto eval_R = [&](Real U_kmh) -> Real {
        Real U_ftmin = U_kmh * U_conv * waf_val;
        U_ftmin = amrex::max(U_ftmin, Real(0.0));
        if (rc.U_max_ftmin > Real(0.0))
            U_ftmin = amrex::min(U_ftmin, rc.U_max_ftmin);
        const Real phi_w = rc.C * std::pow(U_ftmin, rc.B) * rc.beta_ratio_E;
        const Real R_ftmin = rc.R0 * (Real(1.0) + phi_w + rc.phi_s);
        return amrex::max(R_ftmin * rp.ros_conv, Real(0.0));
    };

    const Real R_hi = eval_R(U_max_kmh);
    if (R_hi < R_SA_ms) return Real(-1.0);

    Real lo = Real(0.0);
    Real hi = U_max_kmh;
    for (int iter = 0; iter < max_iter; ++iter) {
        const Real mid = Real(0.5) * (lo + hi);
        const Real R_mid = eval_R(mid);
        if (R_mid < R_SA_ms) lo = mid;
        else                  hi = mid;
        if (hi - lo < Real(0.01)) break;
    }
    return Real(0.5) * (lo + hi);
}

// ============================================================================
// compute_full_scott_reinhardt
//
// Fills ti_mf and ci_mf with the full TI/CI values [km/h] for each cell.
// For uniform fuel runs (no landscape file): constant values everywhere.
// For per-cell fuel tables: evaluate per fuel code (not yet implemented;
// fills uniform value using global RothermelParams).
// ============================================================================
inline void compute_full_scott_reinhardt(
    MultiFab& ti_mf,
    MultiFab& ci_mf,
    const RothermelParams& rp,
    const InputParameters::CrownInitiationParams& cp,
    const MultiFab* /*fuel_model_mf*/ = nullptr,
    const void*     /*fuel_table*/     = nullptr,
    int             /*fuel_table_size*/= 0,
    amrex::Real     U_max_kmh          = amrex::Real(200.0))
{
    using namespace amrex;

    // Compute Rothermel params for global fuel
    const RothermelComputed rc = compute_rothermel_params(rp);

    // Crown initiation threshold
    const Real CBH = Real(cp.CBH);
    const Real CBD = Real(cp.CBD);
    const Real FMC = Real(cp.FMC);
    const Real I_o = Real(0.010) * CBH * (Real(460.0) + Real(25.9) * FMC);

    // Critical active crown ROS [m/s]
    const Real R_SA = Real(3.0) / amrex::max(CBD, Real(0.01)) / Real(60.0);

    // Heat content and net fuel load for Byram intensity
    const Real h_kJ_kg    = Real(rp.h_heat) * WildfireConst::BTU_LB_TO_KJ_KG;
    const Real total_pc   = rp.w_d1 + rp.w_d10 + rp.w_d100 + rp.w_lh + rp.w_lw;
    const Real w_eff      = (total_pc > 0.0) ? total_pc : rp.w0;
    const Real w_net_kg_m2 = w_eff * (Real(1.0) - Real(rp.S_T))
                            * WildfireConst::LB_FT2_TO_KG_M2;

    const Real U_max = amrex::max(U_max_kmh, Real(10.0));

    // Compute scalar TI/CI using bisection
    const Real TI_val = compute_torching_index_kmh(rc, rp, I_o, h_kJ_kg, w_net_kg_m2, U_max);
    const Real CI_val = compute_crowning_index_kmh(rc, rp, R_SA, U_max);

    // Fill MultiFabs with constant value
    ti_mf.setVal(TI_val);
    ci_mf.setVal(CI_val);
}

#endif // FIRE_ECOLOGY_MODEL_H
