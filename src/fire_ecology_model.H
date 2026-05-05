// ============================================================================
// Fire Ecology Diagnostics
// ============================================================================
// Computes four per-cell fire ecology and behavior fields written to every
// plotfile as additional diagnostic outputs:
//
//   Component 0 – scorch_height   [m]   Van Wagner (1973)
//   Component 1 – prob_ignition   [-]   Anderson (1970) / Rothermel (1983)
//   Component 2 – tree_mortality  [-]   Ryan-Reinhardt (1988) style logistic
//   Component 3 – crown_activity  [-]   0=surface, 1=passive, 2=active
//                                       Van Wagner (1977) / Scott & Reinhardt (2001)
//
// ----------------------------------------------------------------------------
// Scorch height (Van Wagner 1973)
// ----------------------------------------------------------------------------
//   H_s [m] = 0.1483 × I_B [kW/m]^(2/3)
//
// Empirical formula for flame-induced scorch height under standard ambient
// conditions (Scott & Reinhardt 2001).  I_B is the Byram (1959) fireline
// intensity already computed in compute_fire_behavior.H.
//
// ----------------------------------------------------------------------------
// Probability of ignition (Anderson 1970 / Rothermel 1983)
// ----------------------------------------------------------------------------
//   P_i = min(1, max(0, 0.000048 × T_fuel_F^1.4 × exp(-0.07 × MC%)))
//
// Probability that a firebrand landing on fine fuel will cause sustained
// ignition.  T_fuel_F is fine fuel temperature in °F (ambient + solar
// increment, default +25°F).  MC% is 1-hr dead fuel moisture content in %.
//
// ----------------------------------------------------------------------------
// Tree mortality fraction (Ryan & Reinhardt 1988)
// ----------------------------------------------------------------------------
//   CRS = (H_s - CBH) / (H_tree - CBH)   clipped to [0, 1]
//   P_kill = 1 / (1 + exp(-4 × (CRS - 0.6)))   when H_s > CBH
//          = 0                                   when H_s <= CBH
//
// Crown Scorch Ratio (CRS) drives a logistic mortality curve calibrated to
// give ~50% mortality when 60% of the live crown is scorched.  Uses the same
// CBH from crown initiation parameters.
//
// ----------------------------------------------------------------------------
// Crown activity class (Van Wagner 1977 / Scott & Reinhardt 2001)
// ----------------------------------------------------------------------------
//   0 – surface fire     : I_B < I_o
//   1 – passive crown    : I_B >= I_o  AND  R < R'_SA
//   2 – active crown     : I_B >= I_o  AND  R >= R'_SA
//
// Crown initiation intensity (Van Wagner 1977, as in crown_initiation.H):
//   I_o [kW/m] = 0.010 × CBH × (460 + 25.9 × FMC)
//
// Critical active-crown ROS (Van Wagner 1977):
//   R'_SA [m/s] = 3.0 / CBD / 60
//
// ----------------------------------------------------------------------------
// References
// ----------------------------------------------------------------------------
//   Van Wagner, C.E. (1973). "Height of crown scorch in forest fires."
//     Canadian Journal of Forest Research, 3(3), 373–378.
//   Anderson, H.E. (1970). "Forest fuel ignitibility."  Fire Technology, 6(4).
//   Rothermel, R.C. (1983).  USDA Forest Service GTR INT-143.
//   Ryan, K.C. & Reinhardt, E.D. (1988). "Predicting postfire mortality of
//     seven western conifers." Canadian Journal of Forest Research, 18, 1291–1297.
//   Van Wagner, C.E. (1977). "Conditions for the start and spread of crown fire."
//     Canadian Journal of Forest Research, 7(1), 23–34.
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
#include <cmath>

using namespace amrex;

// Number of components in the ecology MultiFab
constexpr int ECOLOGY_NCOMP = 4;

// ============================================================================
// compute_fire_ecology
//
// Parameters:
//   ecology_mf             – output: ECOLOGY_NCOMP-component MultiFab
//   fireline_intensity_mf  – input:  Byram fireline intensity [kW/m]
//   R_mf                   – input:  rate of spread [m/s]
//   rp                     – Rothermel parameters (provides M_f for P_i)
//   fep                    – fire ecology parameters
//   cp                     – crown initiation parameters (CBH, CBD, FMC,
//                            tree_height for mortality; crown.enable not
//                            checked here — ecology is always computed)
// ============================================================================
static void compute_fire_ecology(
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
    // Anderson (1970) formula uses T_fuel^1.4 – precompute for GPU
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
            // (0) Scorch height [m]  –  Van Wagner (1973)
            //     H_s = 0.1483 × I_B^(2/3)   for I_B in kW/m
            // ----------------------------------------------------------------
            const Real H_s = (IB > Real(0.0))
                             ? Real(0.1483) * std::pow(IB, Real(2.0 / 3.0))
                             : Real(0.0);
            eco(i, j, k, 0) = H_s;

            // ----------------------------------------------------------------
            // (1) Probability of ignition [-]  –  Anderson (1970)
            //     P_i = 0.000048 × T_fuel_F^1.4 × exp(-0.07 × MC%)
            // ----------------------------------------------------------------
            const Real P_i = amrex::min(Real(1.0),
                             amrex::max(Real(0.0),
                             Real(0.000048) * Tf14
                             * std::exp(Real(-0.07) * MC)));
            eco(i, j, k, 1) = P_i;

            // ----------------------------------------------------------------
            // (2) Tree mortality fraction [-]  –  Ryan-Reinhardt (1988) style
            //     CRS = (H_s - CBH) / (H_tree - CBH)  clipped [0,1]
            //     P_kill = 1 / (1 + exp(-4 × (CRS - 0.6)))
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
        });
    }
}

#endif // FIRE_ECOLOGY_MODEL_H
