// ============================================================================
// Compute Canadian FBP System rate of spread and store in a MultiFab
// ============================================================================
// Mirrors compute_cheney_gould_R.H but uses the FBP O1a/O1b/S1/S2/S3 models.
// Outputs R in m/s.
// ============================================================================

#ifndef COMPUTE_FBP_R_H
#define COMPUTE_FBP_R_H

#include "constants.H"
#include "parse_inputs.H"
#include "fbp_model.H"
#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_GpuLaunch.H>

using namespace amrex;

// ============================================================================
// compute_fbp_R
//
// Fills R_mf [m/s] using the Canadian FBP System fire spread model.
//
// Parameters:
//   R_mf – output MultiFab (1 component, no ghost cells)
//   vel  – velocity field (wind, m/s)
//   fp   – FBP model parameters (fuel_type, moisture, curing)
//   fwi_mf – optional dynamic FWI MultiFab (6 components)
//   fwi_enable – optional dynamic FWI enable flag
// ============================================================================
inline void compute_fbp_R(MultiFab& R_mf,
                           const MultiFab& vel,
                           const InputParameters::FBPParams& fp,
                           const MultiFab* fwi_mf = nullptr,
                           bool fwi_enable = false)
{
    // Pre-compute scalar coefficients on host; captured by GPU lambda
    const FBPComputed fc = compute_fbp_params(fp);

    // Domain-mean wind for clamping extreme cells
    const Long ncells_vel = vel.boxArray().numPts();
    const Real ux_global_mean = (ncells_vel > 0) ? vel.sum(0) / Real(ncells_vel) : Real(0.0);
    const Real uy_global_mean = (ncells_vel > 0) ? vel.sum(1) / Real(ncells_vel) : Real(0.0);

    const bool use_fwi = fwi_enable && (fwi_mf != nullptr);

    for (MFIter mfi(R_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const v  = vel.const_array(mfi);
        auto       R  = R_mf.array(mfi);
        auto const fwi = use_fwi ? fwi_mf->const_array(mfi) : Array4<const Real>();

        // Map fuel type to an integer code for device kernel
        // 0: o1a, 1: o1b, 2: s1, 3: s2, 4: s3
        int fuel_type_code = 0;
        if (fp.fuel_type == "o1b") fuel_type_code = 1;
        else if (fp.fuel_type == "s1") fuel_type_code = 2;
        else if (fp.fuel_type == "s2") fuel_type_code = 3;
        else if (fp.fuel_type == "s3") fuel_type_code = 4;

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real ux = v(i,j,k,0);
            Real uy = v(i,j,k,1);
            // Clamp extreme wind speeds
            constexpr Real MAX_WIND_SQ = WildfireConst::MAX_WIND_SPEED_SQ;
            if (ux*ux + uy*uy > MAX_WIND_SQ) {
                ux = ux_global_mean;
                uy = uy_global_mean;
            }
            const Real wind_speed = std::sqrt(ux*ux + uy*uy);
            
            if (use_fwi) {
                // Get dynamic FFMC and BUI from the FWI MultiFab
                const Real ffmc = fwi(i, j, k, 0);
                const Real bui  = fwi(i, j, k, 4);

                // Convert FFMC to moisture percent (Van Wagner 1985)
                const Real m0 = Real(147.2) * (Real(101.0) - ffmc) / (Real(59.5) + ffmc);

                // Re-evaluate fine-fuel moisture effect f_F dynamically
                const Real exp_part = WildfireConst::FBP_ISI_F_PREFAC * std::exp(WildfireConst::FBP_ISI_F_EXP * m0);
                const Real poly_part = amrex::Real(1.0) + std::pow(m0, WildfireConst::FBP_ISI_F_POW) / WildfireConst::FBP_ISI_F_DIV;
                const Real f_F = amrex::max(exp_part * poly_part, amrex::Real(0.0));

                // Recalculate ISI with dynamic f_F and local wind speed
                const Real U_kmh = wind_speed * amrex::Real(3.6);
                const Real f_W   = std::exp(WildfireConst::FBP_ISI_W_COEFF * U_kmh);
                const Real ISI   = WildfireConst::FBP_ISI_MULT * f_W * f_F;

                // Calculate ROS with dynamic ISI and Buildup Effect
                Real BE = amrex::Real(1.0);
                if (bui > amrex::Real(0.0)) {
                    if (fuel_type_code == 2) { // s1
                        BE = std::exp(-amrex::Real(50.0) * (amrex::Real(1.0) / bui - amrex::Real(1.0) / amrex::Real(38.0)));
                    } else if (fuel_type_code == 3) { // s2
                        BE = std::exp(-amrex::Real(50.0) * (amrex::Real(1.0) / bui - amrex::Real(1.0) / amrex::Real(63.0)));
                    } else if (fuel_type_code == 4) { // s3
                        BE = std::exp(-amrex::Real(50.0) * (amrex::Real(1.0) / bui - amrex::Real(1.0) / amrex::Real(43.0)));
                    }
                }

                // Curing factor and ROS
                const Real CF = fc.CF;
                const Real ROS = fc.C_ros * std::exp(fc.k_isi * ISI) * CF * BE;
                R(i,j,k) = amrex::max(ROS, amrex::Real(0.0));
            } else {
                R(i,j,k) = fbp_ros(fc, wind_speed);
            }
        });
    }
}

#endif // COMPUTE_FBP_R_H
