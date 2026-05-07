// ============================================================================
// Fire Behavior Diagnostics: Fireline Intensity and Flame Length
// ============================================================================
// Computes Byram (1959) fireline intensity (kW/m) and flame length (m) from
// the Rothermel rate of spread field and fuel parameters.
//
// References:
//   Byram, G.M. (1959). "Combustion of forest fuels." In: Davis, K.P. (ed.)
//   Forest Fire: Control and Use. McGraw-Hill, New York. pp. 61–89.
//
// Fireline intensity (Byram 1959):
//   I_B [kW/m] = H [kJ/kg] * w_a [kg/m²] * R [m/s]
//
// where H is the low heat of combustion, w_a is the available (net) fuel load
// per unit area, and R is the rate of spread.
//
// Flame length (Byram 1959, SI form):
//   L_f [m] = 0.0775 * I_B [kW/m]^0.46
// ============================================================================

#ifndef COMPUTE_FIRE_BEHAVIOR_H
#define COMPUTE_FIRE_BEHAVIOR_H

#include "constants.H"
#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_GpuLaunch.H>
#include "parse_inputs.H"

using namespace amrex;

// ============================================================================
// Compute fireline intensity and flame length from the ROS field
//
// Parameters:
//   fireline_intensity_mf  – output: Byram fireline intensity [kW/m]
//   flame_length_mf        – output: Byram flame length [m]
//   R_mf                   – input:  rate of spread [m/s] (simulation units)
//   rp                     – Rothermel fuel parameters (provides fuel load
//                            and heat content)
// ============================================================================
static void compute_fire_behavior(MultiFab& fireline_intensity_mf,
                                  MultiFab& flame_length_mf,
                                  const MultiFab& R_mf,
                                  const RothermelParams& rp)
{
    // Heat content [kJ/kg]
    const Real h_kJ_kg = rp.h_heat * WildfireConst::BTU_LB_TO_KJ_KG;

    // Available (net) fuel load [kg/m²]
    // Use per-class total when multi-class fuel is active; fall back to w0.
    const Real total_per_class = rp.w_d1 + rp.w_d10 + rp.w_d100
                                + rp.w_lh + rp.w_lw;
    const Real w_eff_lb_ft2 = (total_per_class > 0.0) ? total_per_class : rp.w0;
    const Real w_net_kg_m2  = w_eff_lb_ft2 * (1.0 - rp.S_T) * WildfireConst::LB_FT2_TO_KG_M2;

    const Real h_times_w = h_kJ_kg * w_net_kg_m2;  // kJ/m²

    for (MFIter mfi(fireline_intensity_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const R  = R_mf.const_array(mfi);
        auto fi       = fireline_intensity_mf.array(mfi);
        auto fl       = flame_length_mf.array(mfi);
        const Real hw = h_times_w;

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            // Fireline intensity [kW/m] = H[kJ/kg] * w_net[kg/m²] * R[m/s]
            Real I_B = hw * amrex::max(R(i, j, k), Real(0.0));
            fi(i, j, k) = I_B;

            // Flame length [m] (Byram 1959, SI): L_f = 0.0775 * I_B^0.46
            fl(i, j, k) = (I_B > Real(0.0))
                          ? Real(0.0775) * std::pow(I_B, Real(0.46))
                          : Real(0.0);
        });
    }
}

#endif // COMPUTE_FIRE_BEHAVIOR_H
