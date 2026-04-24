// ============================================================================
// FARSITE Ellipse Model (Richards 1990)
// ============================================================================
// This implements Richard's FARSITE elliptical fire spread model.
// The model finds locations where phi ≈ 0 (fire front) and uses an ellipse
// equation to compute where the spread will be at the next time step.
// ============================================================================

#ifndef FARSITE_ELLIPSE_H
#define FARSITE_ELLIPSE_H

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuLaunch.H>
#include "parse_inputs.H"
#include "rothermel_model.H"

using namespace amrex;

// ============================================================================
// Find locations where phi is close to 0 (fire front) and compute
// ellipse-based spread using Richard's FARSITE model with Rothermel ROS
// ============================================================================
static void compute_farsite_spread(const MultiFab& phi,
                                   const MultiFab& vel,
                                   MultiFab& farsite_spread,
                                   const Geometry& geom,
                                   Real dt,
                                   const RothermelParams& rp,
                                   const InputParameters::FARSITEParams& fp)
{
    if (fp.enable == 0) {
        return;  // FARSITE disabled
    }

    const auto dx = geom.CellSize();
#if (AMREX_SPACEDIM == 3)
    GpuArray<Real,AMREX_SPACEDIM> gdx{dx[0], dx[1], dx[2]};
#else
    GpuArray<Real,AMREX_SPACEDIM> gdx{dx[0], dx[1]};
#endif

    const Real phi_thresh = fp.phi_threshold;
    // const Real LW_ratio = fp.length_to_width_ratio;  // Reserved for future ellipse implementation

    // Compute Rothermel (1972) fire spread parameters
    const RothermelComputed rc = compute_rothermel_params(rp);
    const Real R0 = rc.R0;
    const Real C = rc.C;
    const Real B = rc.B;
    const Real beta_ratio_E = rc.beta_ratio_E;
    const Real phi_s = rc.phi_s;
    const Real wind_conv = rc.wind_conv;
    const Real ros_conv = rc.ros_conv;

    // Compute spread for each cell
    for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const p = phi.const_array(mfi);
        auto const v = vel.const_array(mfi);
        auto spread = farsite_spread.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            // Check if this cell is near the fire front (phi ≈ 0)
            Real phi_val = p(i,j,k);
            
            // Initialize spread to zero
#if (AMREX_SPACEDIM == 3)
            spread(i,j,k,0) = 0.0;
            spread(i,j,k,1) = 0.0;
            spread(i,j,k,2) = 0.0;
#else
            spread(i,j,k,0) = 0.0;
            spread(i,j,k,1) = 0.0;
#endif

            // Only compute spread for cells near the fire front
            if (std::abs(phi_val) > phi_thresh) {
                return;
            }

            // Compute gradient of phi to get normal direction
            Real dphi_dx = (p(i+1,j,k) - p(i-1,j,k)) / (2.0 * gdx[0]);
            Real dphi_dy = (p(i,j+1,k) - p(i,j-1,k)) / (2.0 * gdx[1]);
#if (AMREX_SPACEDIM == 3)
            Real dphi_dz = (p(i,j,k+1) - p(i,j,k-1)) / (2.0 * gdx[2]);
            Real grad_mag = std::sqrt(dphi_dx*dphi_dx + dphi_dy*dphi_dy + dphi_dz*dphi_dz);
#else
            Real grad_mag = std::sqrt(dphi_dx*dphi_dx + dphi_dy*dphi_dy);
#endif

            // Avoid division by zero
            if (grad_mag < 1.0e-10) {
                return;
            }

            // Normal vector (points toward unburned region, phi > 0)
            Real nx = dphi_dx / grad_mag;
            Real ny = dphi_dy / grad_mag;
#if (AMREX_SPACEDIM == 3)
            Real nz = dphi_dz / grad_mag;
#endif

            // Get wind velocity components
            Real ux = v(i,j,k,0);
            Real uy = v(i,j,k,1);
#if (AMREX_SPACEDIM == 3)
            Real uz = v(i,j,k,2);
#endif

            // Wind speed magnitude in simulation units
            Real wind_mag = std::sqrt(ux*ux + uy*uy
#if (AMREX_SPACEDIM == 3)
                + uz*uz
#endif
            );

            // Convert wind speed to Rothermel units (ft/min)
            Real U_ftmin = wind_mag * wind_conv;

            // Wind factor φ_w (Rothermel Eq. 47)
            Real phi_w = (U_ftmin > 0.0) ? C * std::pow(U_ftmin, B) * beta_ratio_E : 0.0;

            // Rate of spread using full Rothermel model [ft/min] → simulation units
            Real R_head = std::max(R0 * (1.0 + phi_w + phi_s), 0.0) * ros_conv;

            // Future work: Use ellipse semi-axes for more accurate directional spread
            // Real a = R_head * dt;  // semi-major axis (downwind direction)
            // Real b = a / LW_ratio;  // semi-minor axis (crosswind direction)

            // Spread distance in the direction of the fire normal
            // This represents where the fire front will move
            Real spread_dist = R_head * dt;

            // Compute new position (displacement from current position)
            spread(i,j,k,0) = spread_dist * nx;
            spread(i,j,k,1) = spread_dist * ny;
#if (AMREX_SPACEDIM == 3)
            spread(i,j,k,2) = spread_dist * nz;
#endif
        });
    }
}

#endif // FARSITE_ELLIPSE_H
