// ============================================================================
// Compute Rothermel wind speed R and store in MultiFab
// ============================================================================
#ifndef COMPUTE_ROTHERMEL_R_H
#define COMPUTE_ROTHERMEL_R_H

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include "parse_inputs.H"
#include "rothermel_model.H"

using namespace amrex;

// ============================================================================
// Compute Rothermel rate of spread R and store in MultiFab
// ============================================================================
static void compute_rothermel_R(MultiFab& R_mf,
                                const MultiFab& vel,
                                const Geometry& /*geom*/,
                                const RothermelParams& rp,
                                const MultiFab* slopes = nullptr)
{
    // Compute Rothermel (1972) fire spread parameters
    const RothermelComputed rc = compute_rothermel_params(rp);
    const Real R0 = rc.R0;
    const Real C = rc.C;
    const Real B = rc.B;
    const Real beta_ratio_E = rc.beta_ratio_E;
    const Real phi_s_const = rc.phi_s;  // Constant phi_s for backward compatibility
    const Real beta = rc.beta;  // Packing ratio for terrain-based phi_s calculation
    const Real wind_conv = rc.wind_conv;
    const Real ros_conv = rc.ros_conv;

    for (MFIter mfi(R_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const v = vel.const_array(mfi);
        auto const R = R_mf.array(mfi);
        
        // Get slopes array if terrain file was provided
        const bool use_terrain = (slopes != nullptr);
        Array4<const Real> s;
        if (use_terrain) {
            s = slopes->const_array(mfi);
        }

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            // Horizontal wind components
            Real ux = v(i,j,k,0);
            Real uy = v(i,j,k,1);

            // Wind speed in ft/min (Rothermel units)
            Real U_ftmin = std::sqrt(ux*ux + uy*uy) * wind_conv;

            // Wind factor φ_w (Rothermel Eq. 47)
            Real phi_w = (U_ftmin > 0.0) ? C * std::pow(U_ftmin, B) * beta_ratio_E : 0.0;

            // Slope factor φ_s: use terrain-based if available, otherwise constant
            Real phi_s;
            if (use_terrain) {
                // Compute phi_s from terrain slopes at this cell
                // Using full Rothermel Eq. 51: φ_s = 5.275 * β^(-0.3) * tan²(φ)
                constexpr Real PHI_S_COEFFICIENT = 5.275;
                constexpr Real BETA_EXPONENT = -0.3;
                Real slope_x = s(i,j,k,0);
                Real slope_y = s(i,j,k,1);
                Real tan_phi_sq = slope_x * slope_x + slope_y * slope_y;
                phi_s = std::max(PHI_S_COEFFICIENT * std::pow(beta, BETA_EXPONENT) * tan_phi_sq, Real(0.0));
            } else {
                // Use constant phi_s from RothermelComputed
                phi_s = phi_s_const;
            }

            // Rate of spread [ft/min] → simulation units
            R(i,j,k) = std::max(R0 * (1.0 + phi_w + phi_s), 0.0) * ros_conv;
            amrex::Print() << "Cell (" << i << "," << j << "): "<< R(i,j,k)<<std::endl;
        });
    }
}

#endif // COMPUTE_ROTHERMEL_R_H
