// ============================================================================
// Compute Lautenberger (2013) physics-based rate of spread and store in a MultiFab
// ============================================================================
// Mirrors compute_cheney_gould_R.H but uses the Lautenberger semi-physical model.
// Optionally accepts a terrain slopes MultiFab for slope correction.
// Outputs R in m/s.
// ============================================================================

#ifndef COMPUTE_LAUTENBERGER_R_H
#define COMPUTE_LAUTENBERGER_R_H

#include "constants.H"
#include "parse_inputs.H"
#include "lautenberger_model.H"
#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_GpuLaunch.H>

using namespace amrex;

// ============================================================================
// compute_lautenberger_R
//
// Fills R_mf [m/s] using the Lautenberger (2013) semi-physical spread model.
//
// Parameters:
//   R_mf   – output MultiFab (1 component, no ghost cells)
//   vel    – velocity field (wind, m/s)
//   rp     – Rothermel parameters (provides sigma, M_f)
//   lp     – Lautenberger model parameters (A_L, B_L, C_L, D_L)
//   slopes – optional MultiFab with terrain slopes (component 0 = combined slope
//            magnitude tan(phi)), or nullptr for flat terrain
// ============================================================================
inline void compute_lautenberger_R(MultiFab& R_mf,
                                    const MultiFab& vel,
                                    const RothermelParams& rp,
                                    const InputParameters::LautenbergerParams& lp,
                                    const MultiFab* slopes = nullptr)
{
    // Pre-compute scalar coefficients on host; captured by GPU lambda
    const LautenbergerComputed lc = compute_lautenberger_params(rp, lp);

    // Domain-mean wind for clamping extreme cells
    const Long ncells_vel = vel.boxArray().numPts();
    const Real ux_global_mean = (ncells_vel > 0) ? vel.sum(0) / Real(ncells_vel) : Real(0.0);
    const Real uy_global_mean = (ncells_vel > 0) ? vel.sum(1) / Real(ncells_vel) : Real(0.0);

    const bool have_slopes = (slopes != nullptr);

    for (MFIter mfi(R_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const v  = vel.const_array(mfi);
        auto       R  = R_mf.array(mfi);

        Array4<const Real> slp_arr;
        if (have_slopes) {
            slp_arr = slopes->const_array(mfi);
        }

        const bool use_slp = have_slopes;

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

            // Terrain slope: use combined slope magnitude if available
            Real tan_slp = Real(0.0);
            if (use_slp) {
                const Real sx = slp_arr(i,j,k,0);
                const Real sy = slp_arr(i,j,k,1);
                tan_slp = std::sqrt(sx*sx + sy*sy);
            }

            R(i,j,k) = lautenberger_ros(lc, wind_speed, tan_slp);
        });
    }
}

#endif // COMPUTE_LAUTENBERGER_R_H
