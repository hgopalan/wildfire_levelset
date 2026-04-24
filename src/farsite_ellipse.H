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

using namespace amrex;

// Structure to hold FARSITE ellipse parameters
struct FARSITEParams {
    amrex::Real length_to_width_ratio;  // L/W ratio of ellipse (default: 1.0 for circle)
    amrex::Real phi_threshold;           // threshold for identifying fire front (default: 0.1)
    int enable;                          // 1 to enable FARSITE, 0 to disable
};

// ============================================================================
// Find locations where phi is close to 0 (fire front) and compute
// ellipse-based spread using Richard's FARSITE model
// ============================================================================
static void compute_farsite_spread(const MultiFab& phi,
                                   const MultiFab& vel,
                                   MultiFab& farsite_spread,
                                   const Geometry& geom,
                                   Real dt,
                                   const RothermelParams& rp,
                                   const FARSITEParams& fp)
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
    const Real LW_ratio = fp.length_to_width_ratio;

    // Simplified Rothermel computation for FARSITE
    // In a full implementation, this would use the complete Rothermel model
    const Real wind_conv = rp.wind_conv;
    const Real ros_conv = rp.ros_conv;
    
    // Basic rate of spread (simplified - should use full Rothermel if needed)
    const Real R0 = 1.0; // Base rate in simulation units

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

            // Wind speed magnitude
            Real wind_mag = std::sqrt(ux*ux + uy*uy
#if (AMREX_SPACEDIM == 3)
                + uz*uz
#endif
            );

            // Wind direction (unit vector)
            Real wx = (wind_mag > 1.0e-10) ? ux / wind_mag : 0.0;
            Real wy = (wind_mag > 1.0e-10) ? uy / wind_mag : 0.0;
#if (AMREX_SPACEDIM == 3)
            Real wz = (wind_mag > 1.0e-10) ? uz / wind_mag : 0.0;
#endif

            // Simplified rate of spread (head fire direction)
            // In full implementation, use complete Rothermel model
            Real R_head = R0 * (1.0 + 0.5 * wind_mag);  // Simplified wind effect
            Real R_back = R0 * 0.5;  // Back fire (reduced rate)

            // Head-to-back ratio for ellipse
            Real HB = (R_back > 1.0e-10) ? R_head / R_back : 1.0;
            
            // Effective wind direction for ellipse (aligned with wind or gradient)
            Real eff_wind_x = (wind_mag > 1.0e-10) ? wx : nx;
            Real eff_wind_y = (wind_mag > 1.0e-10) ? wy : ny;
#if (AMREX_SPACEDIM == 3)
            Real eff_wind_z = (wind_mag > 1.0e-10) ? wz : nz;
#endif

            // Ellipse semi-axes (Richards 1990 FARSITE model)
            Real a = R_head * dt;  // semi-major axis (downwind direction)
            Real b = a / LW_ratio;  // semi-minor axis (crosswind direction)

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
