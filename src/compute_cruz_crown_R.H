// ============================================================================
// Compute Cruz, Alexander & Wakimoto (2005) crown fire rate of spread
// and store in a MultiFab
// ============================================================================
// Mirrors compute_cheney_gould_R.H but uses the Cruz–Crown empirical model.
// Outputs R in simulation length/time units (m/s by default).
//
// The Cruz–Crown model is a pure wind-speed model (with canopy bulk density
// and 10-h fuel moisture as fuel-state scalars). Terrain slope is not
// accounted for in the original formulation and is intentionally omitted.
// ============================================================================

#ifndef COMPUTE_CRUZ_CROWN_R_H
#define COMPUTE_CRUZ_CROWN_R_H

#include "constants.H"
#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuLaunch.H>
#include "parse_inputs.H"
#include "cruz_crown_model.H"

using namespace amrex;

// ============================================================================
// compute_cruz_crown_R
//
// Fills R_mf [m/s] using the Cruz–Crown (2005) algebraic crown fire spread
// model.
//
// Parameters:
//   R_mf   – output MultiFab (1 component, no ghost cells)
//   vel    – velocity field (wind, m/s)
//   ccp    – Cruz–Crown model parameters (CBD, MC10)
// ============================================================================
inline void compute_cruz_crown_R(MultiFab& R_mf,
                                  const MultiFab& vel,
                                  const InputParameters::CruzCrownParams& ccp)
{
    // Pre-compute scalar coefficients on the host; copy to GPU via capture
    const CruzCrownComputed ccc = compute_cruz_crown_params(ccp);

    // Domain-mean wind for clamping extreme cells (same strategy as
    // Cheney–Gould and Balbi implementations)
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

            // Cruz–Crown ROS [m/s]
            R(i,j,k) = cruz_crown_ros(ccc, wind_speed);
        });
    }
}

#endif // COMPUTE_CRUZ_CROWN_R_H
