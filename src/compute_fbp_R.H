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
// ============================================================================
inline void compute_fbp_R(MultiFab& R_mf,
                           const MultiFab& vel,
                           const InputParameters::FBPParams& fp)
{
    // Pre-compute scalar coefficients on host; captured by GPU lambda
    const FBPComputed fc = compute_fbp_params(fp);

    // Domain-mean wind for clamping extreme cells
    const Long ncells_vel = vel.boxArray().numPts();
    const Real ux_global_mean = (ncells_vel > 0) ? vel.sum(0) / Real(ncells_vel) : Real(0.0);
    const Real uy_global_mean = (ncells_vel > 0) ? vel.sum(1) / Real(ncells_vel) : Real(0.0);

    for (MFIter mfi(R_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const v  = vel.const_array(mfi);
        auto       R  = R_mf.array(mfi);

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
            R(i,j,k) = fbp_ros(fc, wind_speed);
        });
    }
}

#endif // COMPUTE_FBP_R_H
