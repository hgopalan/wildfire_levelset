// ============================================================================
// FARSITE Ellipse Model (Richards 1990)
// ============================================================================
// This implements Richards' (1990) elliptical fire spread model with 
// coefficients a, b, and c for directional spread rates.
//
// The model can use either:
// 1. Fixed Richards' coefficients (a, b, c) from input parameters
// 2. Anderson (1983) L/W ratio based on wind speed, converted to coefficients
//
// The model:
// 1. Finds locations where phi ≈ 0 (fire front)
// 2. Uses Rothermel (1972) to compute base rate of spread (ROS)
// 3. Applies Richards' elliptical expansion with coefficients:
//    - a: head fire coefficient (downwind, maximum spread)
//    - b: flank fire coefficient (crosswind spread)
//    - c: backing fire coefficient (upwind, minimum spread)
// 4. Computes directional spread based on angle between fire normal and wind
// ============================================================================

#ifndef FARSITE_ELLIPSE_H
#define FARSITE_ELLIPSE_H

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuLaunch.H>
#include "parse_inputs.H"
#include "rothermel_model.H"
#include "fire_models.H"

using namespace amrex;

// ============================================================================
// Find locations where phi is close to 0 (fire front) and compute
// ellipse-based spread using Richards' (1990) model with coefficients a, b, c
// combined with Rothermel (1972) ROS
// ============================================================================
static void compute_farsite_spread(MultiFab& phi,
                                   const MultiFab& vel,
                                   MultiFab& farsite_spread,
                                   const Geometry& geom,
                                   Real dt,
                                   const RothermelParams& rp,
                                   const InputParameters::FARSITEParams& fp,
                                   const MultiFab* slopes = nullptr)
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
    
    // Dynamic arrays to collect positions where phi < threshold
    std::vector<Real> xarray;
    std::vector<Real> yarray;
#if (AMREX_SPACEDIM == 3)
    std::vector<Real> zarray;
#endif
    
    // Get physical domain information for coordinate conversion
    const Real* prob_lo = geom.ProbLo();
    
    const int use_anderson_LW = fp.use_anderson_LW;
    
    // Richards' (1990) elliptical expansion coefficients
    // These may be overridden by Anderson L/W ratio if enabled
    Real coeff_a = fp.coeff_a;  // head fire coefficient
    Real coeff_b = fp.coeff_b;  // flank fire coefficient
    Real coeff_c = fp.coeff_c;  // backing fire coefficient

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

    // Compute spread for each cell
    for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const p = phi.const_array(mfi);
        auto const v = vel.const_array(mfi);
        auto spread = farsite_spread.array(mfi);
        
        // Get slopes array if terrain file was provided
        const bool use_terrain = (slopes != nullptr);
        Array4<const Real> s;
        if (use_terrain) {
            s = slopes->const_array(mfi);
        }

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

            // Slope factor φ_s: use terrain-based if available, otherwise constant
            Real phi_s;
            if (use_terrain) {
                // Compute phi_s from terrain slopes at this cell
                // Using full Rothermel Eq. 51: φ_s = 5.275 * β^(-0.3) * tan²(φ)
                Real slope_x = s(i,j,k,0);
                Real slope_y = s(i,j,k,1);
                Real tan_phi_sq = slope_x * slope_x + slope_y * slope_y;
                phi_s = std::max(Real(5.275) * std::pow(beta, Real(-0.3)) * tan_phi_sq, Real(0.0));
            } else {
                // Use constant phi_s from RothermelComputed
                phi_s = phi_s_const;
            }

            // Rate of spread using full Rothermel model [ft/min] → simulation units
            Real R_head = std::max(R0 * (1.0 + phi_w + phi_s), 0.0) * ros_conv;
            
            // Determine elliptical coefficients
            // If Anderson L/W mode is enabled, compute coefficients from wind speed
            Real local_coeff_a = coeff_a;
            Real local_coeff_b = coeff_b;
            Real local_coeff_c = coeff_c;
            
            if (use_anderson_LW == 1) {
                // Convert wind speed to mph for Anderson formula
                // wind_conv converts from sim units to ft/min
                // ft/min to mph: multiply by 60/5280 = 0.01136364
                Real U_mph = U_ftmin * Real(0.01136364);
                
                // Compute Anderson (1983) L/W ratio based on wind speed
                Real LW_ratio = anderson_LW_ratio(U_mph);
                
                // Convert L/W ratio to Richards' coefficients
                LW_ratio_to_richards_coefficients(LW_ratio, 
                                                 local_coeff_a, 
                                                 local_coeff_b, 
                                                 local_coeff_c);
            }
            
            // Richards' (1990) elliptical expansion model
            // Compute spread distances for each direction (not geometric ellipse axes)
            Real head_spread = local_coeff_a * R_head * dt;    // head fire (downwind, maximum spread)
            Real flank_spread = local_coeff_b * R_head * dt;   // flank fire (crosswind)
            Real backing_spread = local_coeff_c * R_head * dt; // backing fire (upwind, minimum spread)
            
            // Determine wind direction (unit vector)
            Real wx = 0.0, wy = 0.0;
#if (AMREX_SPACEDIM == 3)
            Real wz = 0.0;
#endif
            if (wind_mag > 1.0e-10) {
                wx = ux / wind_mag;
                wy = uy / wind_mag;
#if (AMREX_SPACEDIM == 3)
                wz = uz / wind_mag;
#endif
            }
            
            // Compute angle between fire normal and wind direction
            // cos(theta) = dot(normal, wind)
            Real cos_theta = nx * wx + ny * wy;
#if (AMREX_SPACEDIM == 3)
            cos_theta += nz * wz;
#endif
            
            // Richards' elliptical spread rate formula
            // Spread rate R(theta) depends on angle from head fire (wind) direction
            // The formula interpolates between directional spread rates based on cos(theta)
            Real spread_dist = 0.0;
            
            if (wind_mag > 1.0e-10) {
                // With wind: use elliptical expansion based on angle to wind
                // The interpolation follows an elliptical pattern:
                //   theta=0°   (cos=1):  full head fire spread
                //   theta=90°  (cos=0):  full flank fire spread  
                //   theta=180° (cos=-1): full backing fire spread
                
                if (cos_theta > 0.0) {
                    // Forward hemisphere (0° to 90°): interpolate head to flank
                    // Formula: R = head_spread * cos(θ) + flank_spread * sin(θ)
                    // This represents the ellipse equation in polar coordinates
                    Real sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);
                    spread_dist = head_spread * cos_theta + flank_spread * sin_theta;
                } else {
                    // Backward hemisphere (90° to 180°): interpolate flank to backing
                    // Formula: R = flank_spread * sin(θ) + backing_spread * |cos(θ)|
                    // Note: we use -cos_theta since cos is negative in this range
                    Real sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);
                    spread_dist = flank_spread * sin_theta - backing_spread * cos_theta;
                }
            } else {
                // No wind: uniform spread in all directions using base rate
                // Without wind, there's no preferred direction, so use uniform expansion
                spread_dist = R_head * dt;
            }

            // Compute new position (displacement from current position)
            spread(i,j,k,0) = spread_dist * nx;
            spread(i,j,k,1) = spread_dist * ny;
#if (AMREX_SPACEDIM == 3)
            spread(i,j,k,2) = spread_dist * nz;
#endif
        });
    }
    
    // Second pass: collect positions where phi < threshold and add spread
    for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const p = phi.const_array(mfi);
        auto spread = farsite_spread.array(mfi);
        
        // Loop over cells (on CPU) to collect positions
        const auto lo = lbound(bx);
        const auto hi = ubound(bx);
        
        for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
        for (int i = lo.x; i <= hi.x; ++i) {
            Real phi_val = p(i,j,k);
            
            // If phi < threshold, add position + spread to arrays
            if (phi_val < phi_thresh) {
                // Convert cell indices to physical coordinates
                Real x_cell = prob_lo[0] + (i + 0.5) * gdx[0];
                Real y_cell = prob_lo[1] + (j + 0.5) * gdx[1];
                
                // Add spread displacement to get new position
                Real x_new = x_cell + spread(i,j,k,0);
                Real y_new = y_cell + spread(i,j,k,1);
                
                xarray.push_back(x_new);
                yarray.push_back(y_new);
                
#if (AMREX_SPACEDIM == 3)
                Real z_cell = prob_lo[2] + (k + 0.5) * gdx[2];
                Real z_new = z_cell + spread(i,j,k,2);
                zarray.push_back(z_new);
#endif
            }
        }
        }
        }
    }
    
    // Set phi to zero everywhere
    phi.setVal(0.0);
    
    // Set phi at collected positions to a value much smaller than threshold
    const Real PHI_RESET_MULTIPLIER = -10.0;  // Much smaller (more negative) than threshold
    const Real phi_new = phi_thresh * PHI_RESET_MULTIPLIER;
    
    for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto p = phi.array(mfi);
        
        // Loop over all collected positions
        for (size_t idx = 0; idx < xarray.size(); ++idx) {
            Real x_pos = xarray[idx];
            Real y_pos = yarray[idx];
            
            // Convert physical coordinates back to cell indices
            // Use rounding to handle the 0.5 cell-center offset correctly
            int i = static_cast<int>(std::round((x_pos - prob_lo[0]) / gdx[0] - 0.5));
            int j = static_cast<int>(std::round((y_pos - prob_lo[1]) / gdx[1] - 0.5));
#if (AMREX_SPACEDIM == 3)
            Real z_pos = zarray[idx];
            int k = static_cast<int>(std::round((z_pos - prob_lo[2]) / gdx[2] - 0.5));
#else
            int k = 0;  // For 2D, k is always 0
#endif
            
            // Check if the cell is in this box
            if (bx.contains(IntVect(AMREX_D_DECL(i, j, k)))) {
                p(i,j,k) = phi_new;
            }
        }
    }
}

#endif // FARSITE_ELLIPSE_H
