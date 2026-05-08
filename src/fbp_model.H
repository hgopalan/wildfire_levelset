// ============================================================================
// Canadian FBP System — O1a, O1b Grassland and S1, S2, S3 Logging Slash
// ============================================================================
// Implements the Canadian Forest Fire Behavior Prediction (FBP) System head-fire
// rate of spread for five fuel types applicable to border US states (MT, ND, MN,
// WA, ID) and interior grasslands:
//
//   O1a – matted/continuous grass  (Forestry Canada 1992, Table 2)
//   O1b – standing grass           (Forestry Canada 1992, Table 2)
//   S1  – light logging slash      (Forestry Canada 1992, Table 2)
//   S2  – moderate logging slash   (Forestry Canada 1992, Table 2)
//   S3  – heavy logging slash      (Forestry Canada 1992, Table 2)
//
// Reference:
//   Forestry Canada Fire Danger Group (1992). "Development and Structure of
//   the Canadian Forest Fire Behavior Prediction System." Forestry Canada
//   Information Report ST-X-3.
//
// ISI (Initial Spread Index) is derived directly from wind speed and fuel
// moisture (Van Wagner 1985), without requiring the full FWI system:
//
//   f_W = exp(0.05039 * W_10)             [W_10 in km/h]
//   f_F = 91.9 * exp(-0.1386 * m_0)
//         * (1 + m_0^5.31 / 4.93e7)       [m_0 = moisture %]
//   ISI = 0.208 * f_W * f_F
//
// Head-fire ROS [m/min] (Forestry Canada 1992 Table 2):
//   O1a: 190 * exp(0.031 * ISI) * CF     CF = (1 - exp(-2.1*PC/100))^2
//   O1b: 250 * exp(0.035 * ISI) * CF
//   S1:   75 * exp(0.110 * ISI)
//   S2:  200 * exp(0.062 * ISI)
//   S3:  320 * exp(0.010 * ISI)
//
// Output is in m/s (divide m/min by 60).
// ============================================================================

#ifndef FBP_MODEL_H
#define FBP_MODEL_H

#include "constants.H"
#include "parse_inputs.H"
#include <AMReX.H>
#include <AMReX_GpuContainers.H>
#include <cmath>
#include <string>

// ============================================================================
// Pre-computed FBP parameters for GPU kernel
// ============================================================================
struct FBPComputed {
    amrex::Real f_F;     // fuel moisture effect (ISI moisture factor, pre-computed)
    amrex::Real C_ros;   // lead coefficient for the fuel type [m/s]
    amrex::Real k_isi;   // ISI exponent for the fuel type
    amrex::Real CF;      // curing factor [0,1] (O1a/O1b only; 1.0 for S-types)
};

// ============================================================================
// compute_fbp_params
//
// Computes the FBPComputed struct from the FBPParams input struct.
// Called on host; result is captured by GPU kernels.
// ============================================================================
inline FBPComputed compute_fbp_params(const InputParameters::FBPParams& fp)
{
    using namespace WildfireConst;
    FBPComputed fc;

    // moisture as percent
    const amrex::Real m0 = amrex::max(fp.moisture, amrex::Real(0.0));

    // Fine-fuel moisture effect (Van Wagner 1985)
    // f_F = 91.9 * exp(-0.1386 * m0) * (1 + m0^5.31 / 4.93e7)
    const amrex::Real exp_part = FBP_ISI_F_PREFAC * std::exp(FBP_ISI_F_EXP * m0);
    const amrex::Real poly_part = amrex::Real(1.0)
                                  + std::pow(m0, FBP_ISI_F_POW) / FBP_ISI_F_DIV;
    fc.f_F = amrex::max(exp_part * poly_part, amrex::Real(0.0));

    // Fuel type: set C_ros and k_isi; curing factor for O1 types
    const std::string& ft = fp.fuel_type;
    if (ft == "o1b") {
        fc.C_ros = FBP_O1B_C;
        fc.k_isi = FBP_O1B_K;
    } else if (ft == "s1") {
        fc.C_ros = FBP_S1_C;
        fc.k_isi = FBP_S1_K;
        fc.CF    = amrex::Real(1.0);
    } else if (ft == "s2") {
        fc.C_ros = FBP_S2_C;
        fc.k_isi = FBP_S2_K;
        fc.CF    = amrex::Real(1.0);
    } else if (ft == "s3") {
        fc.C_ros = FBP_S3_C;
        fc.k_isi = FBP_S3_K;
        fc.CF    = amrex::Real(1.0);
    } else {
        // default: o1a
        fc.C_ros = FBP_O1A_C;
        fc.k_isi = FBP_O1A_K;
    }

    // Curing factor for O1a/O1b: CF = (1 - exp(-2.1 * PC/100))^2
    if (ft == "o1a" || ft == "o1b") {
        const amrex::Real pc_frac = amrex::max(amrex::Real(0.0),
                                    amrex::min(amrex::Real(100.0),
                                    fp.curing)) / amrex::Real(100.0);
        const amrex::Real cf_inner = amrex::Real(1.0)
                                     - std::exp(-FBP_O1_CURING_K * pc_frac);
        fc.CF = cf_inner * cf_inner;
    }

    return fc;
}

// ============================================================================
// fbp_ros
//
// GPU-device inline: returns head-fire ROS [m/s] given wind speed [m/s].
//
// f_W = exp(0.05039 * U_kmh),  ISI = 0.208 * f_W * f_F
// R_ms = C_ros * exp(k_isi * ISI) * CF
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real fbp_ros(const FBPComputed& fc, amrex::Real wind_speed_ms) noexcept
{
    using namespace WildfireConst;
    // Convert m/s → km/h
    const amrex::Real U_kmh = wind_speed_ms * amrex::Real(3.6);
    const amrex::Real f_W   = std::exp(FBP_ISI_W_COEFF * U_kmh);
    const amrex::Real ISI   = FBP_ISI_MULT * f_W * fc.f_F;
    return amrex::max(fc.C_ros * std::exp(fc.k_isi * ISI) * fc.CF,
                      amrex::Real(0.0));
}

#endif // FBP_MODEL_H
