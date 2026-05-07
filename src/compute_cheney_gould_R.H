// ============================================================================
// Compute Cheney & Gould (1995) / Cheney, Gould & Catchpole (1998)
// grassland rate of spread and store in a MultiFab
// ============================================================================
// Mirrors compute_balbi_R.H but uses the Cheney–Gould empirical model.
// Outputs R in simulation length/time units (m/s by default).
//
// The Cheney–Gould model is a pure wind-speed model; terrain slope is not
// accounted for in the original empirical formulation and is intentionally
// omitted here.
// ============================================================================

#ifndef COMPUTE_CHENEY_GOULD_R_H
#define COMPUTE_CHENEY_GOULD_R_H

#include "constants.H"
#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuLaunch.H>
#include "parse_inputs.H"
#include "cheney_gould_model.H"

using namespace amrex;

// ============================================================================
// compute_cheney_gould_R
//
// Fills R_mf [m/s] using the Cheney–Gould grassland spread model.
//
// Parameters:
//   R_mf   – output MultiFab (1 component, no ghost cells)
//   vel    – velocity field (wind, m/s)
//   cgp    – Cheney–Gould model parameters (moisture, curing)
// ============================================================================
static void compute_cheney_gould_R(MultiFab& R_mf,
                                    const MultiFab& vel,
                                    const InputParameters::CheneyGouldParams& cgp)
{
    // Pre-compute scalar coefficients on the host; copy to GPU via capture
    const CheneyGouldComputed cgc = compute_cheney_gould_params(cgp);

    // Domain-mean wind for clamping extreme cells (same strategy as Rothermel
    // and Balbi implementations)
    const Long ncells_vel = vel.boxArray().numPts();
    const Real ux_global_mean = (ncells_vel > 0) ? vel.sum(0) / Real(ncells_vel) : Real(0.0);
    const Real uy_global_mean = (ncells_vel > 0) ? vel.sum(1) / Real(ncells_vel) : Real(0.0);

    for (MFIter mfi(R_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const v  = vel.const_array(mfi);
        auto       R  = R_mf.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            // Wind components – clamp extreme wind speeds (> 25 m/s)
            Real ux = v(i,j,k,0);
            Real uy = v(i,j,k,1);
            constexpr Real MAX_WIND_SQ = WildfireConst::MAX_WIND_SPEED_SQ;
            if (ux*ux + uy*uy > MAX_WIND_SQ) {
                ux = ux_global_mean;
                uy = uy_global_mean;
            }

            // Wind speed magnitude [m/s]
            const Real wind_speed = std::sqrt(ux*ux + uy*uy);

            // Cheney–Gould ROS [m/s]
            R(i,j,k) = cheney_gould_ros(cgc, wind_speed);
        });
    }
}

#endif // COMPUTE_CHENEY_GOULD_R_H
