// ============================================================================
// Fuel Boundary Smoothing
// ============================================================================
// Eliminates sharp discontinuities at fuel model boundaries by applying
// distance-weighted ROS blending in cells near fuel type transitions.
//
// This extends the existing two_fuel_blending.H logic from temporal blending
// to spatial transitions at fuel boundaries, improving realism by smoothing
// the transition zone (2-3 cells wide).
//
// References:
//   Finney, M.A. (1998). "FARSITE: Fire Area Simulator—Model development
//     and evaluation." USDA Forest Service Research Paper RMRS-RP-4.
//   Finney, M.A. (2006). "An overview of FlamMap fire modeling capabilities."
//     USDA Forest Service Proceedings RMRS-P-41, pp. 213-220.
// ============================================================================

#ifndef FUEL_BOUNDARY_SMOOTHING_H
#define FUEL_BOUNDARY_SMOOTHING_H

#include "constants.H"
#include "two_fuel_blending.H"
#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_GpuLaunch.H>
#include <cmath>

using namespace amrex;

// ============================================================================
// Detect if cell is near a fuel type boundary
//
// A cell is considered "near a boundary" if any of its neighbors has a
// different fuel model code.
//
// Parameters:
//   fuel_model_mf - MultiFab containing fuel model codes
//   i, j, k       - cell indices
//
// Returns:
//   true if cell is within transition zone (has neighbors with different fuel)
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool is_near_fuel_boundary(
    const Array4<const Real>& fuel_arr,
    int i, int j, int k,
    int ncells_x, int ncells_y, int ncells_z) noexcept
{
    // Get center cell fuel code
    Real center_fuel = fuel_arr(i, j, k);
    
    // Check all 8 neighbors in 2D (or 26 in 3D)
    // For simplicity, check immediate neighbors (4 in 2D, 6 in 3D)
    
    // Check x-neighbors
    if (i > 0 && fuel_arr(i-1, j, k) != center_fuel) return true;
    if (i < ncells_x-1 && fuel_arr(i+1, j, k) != center_fuel) return true;
    
    // Check y-neighbors
    if (j > 0 && fuel_arr(i, j-1, k) != center_fuel) return true;
    if (j < ncells_y-1 && fuel_arr(i, j+1, k) != center_fuel) return true;
    
    // Check z-neighbors (for 3D)
    if (ncells_z > 1) {
        if (k > 0 && fuel_arr(i, j, k-1) != center_fuel) return true;
        if (k < ncells_z-1 && fuel_arr(i, j, k+1) != center_fuel) return true;
    }
    
    return false;
}

// ============================================================================
// Compute distance-weighted blend factor for fuel transition
//
// When a cell is near a fuel boundary, computes a blend factor based on
// the fraction of neighbors with each fuel type.
//
// Parameters:
//   fuel_model_mf - MultiFab containing fuel model codes
//   i, j, k       - cell indices
//   target_fuel   - fuel code to compute weight for
//
// Returns:
//   weight - fraction [0-1] representing prevalence of target_fuel in vicinity
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real compute_fuel_blend_weight(
    const Array4<const Real>& fuel_arr,
    int i, int j, int k,
    int ncells_x, int ncells_y, int ncells_z,
    Real target_fuel) noexcept
{
    int count_total = 0;
    int count_target = 0;
    
    // Check all neighbors in a 3x3 (2D) or 3x3x3 (3D) stencil
    for (int ii = amrex::max(i-1, 0); ii <= amrex::min(i+1, ncells_x-1); ++ii) {
        for (int jj = amrex::max(j-1, 0); jj <= amrex::min(j+1, ncells_y-1); ++jj) {
            if (ncells_z > 1) {
                for (int kk = amrex::max(k-1, 0); kk <= amrex::min(k+1, ncells_z-1); ++kk) {
                    ++count_total;
                    if (fuel_arr(ii, jj, kk) == target_fuel) ++count_target;
                }
            } else {
                ++count_total;
                if (fuel_arr(ii, jj, k) == target_fuel) ++count_target;
            }
        }
    }
    
    if (count_total == 0) return Real(0.0);
    
    return static_cast<Real>(count_target) / static_cast<Real>(count_total);
}

// ============================================================================
// Apply fuel boundary smoothing to ROS field
//
// For cells near fuel type boundaries, blends ROS values from adjacent
// fuel types using distance-weighted averaging.
//
// Parameters:
//   R_mf          - ROS field [m/s] (modified in place)
//   fuel_model_mf - fuel model codes (integer)
//   ncells_x/y/z  - domain dimensions
//   blend_method  - blending method from two_fuel_blending.H
// ============================================================================
inline void apply_fuel_boundary_smoothing(
    MultiFab&       R_mf,
    const MultiFab& fuel_model_mf,
    int ncells_x, int ncells_y, int ncells_z,
    BlendMethod blend_method = BlendMethod::LINEAR)
{
    for (MFIter mfi(R_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto R_arr    = R_mf.array(mfi);
        auto fuel_arr = fuel_model_mf.const_array(mfi);
        
        const int nx = ncells_x;
        const int ny = ncells_y;
        const int nz = ncells_z;

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            // Check if cell is near a fuel boundary
            if (!is_near_fuel_boundary(fuel_arr, i, j, k, nx, ny, nz)) {
                // Not near boundary - no smoothing needed
                return;
            }
            
            // Get center cell values
            Real R_center = R_arr(i, j, k);
            Real fuel_center = fuel_arr(i, j, k);
            
            // Compute average ROS from neighbors with different fuel types
            Real R_sum = Real(0.0);
            Real weight_sum = Real(0.0);
            int n_neighbors = 0;
            
            // Check immediate neighbors (4 in 2D, 6 in 3D)
            if (i > 0 && fuel_arr(i-1, j, k) != fuel_center) {
                R_sum += R_arr(i-1, j, k);
                weight_sum += Real(1.0);
                ++n_neighbors;
            }
            if (i < nx-1 && fuel_arr(i+1, j, k) != fuel_center) {
                R_sum += R_arr(i+1, j, k);
                weight_sum += Real(1.0);
                ++n_neighbors;
            }
            if (j > 0 && fuel_arr(i, j-1, k) != fuel_center) {
                R_sum += R_arr(i, j-1, k);
                weight_sum += Real(1.0);
                ++n_neighbors;
            }
            if (j < ny-1 && fuel_arr(i, j+1, k) != fuel_center) {
                R_sum += R_arr(i, j+1, k);
                weight_sum += Real(1.0);
                ++n_neighbors;
            }
            if (nz > 1) {
                if (k > 0 && fuel_arr(i, j, k-1) != fuel_center) {
                    R_sum += R_arr(i, j, k-1);
                    weight_sum += Real(1.0);
                    ++n_neighbors;
                }
                if (k < nz-1 && fuel_arr(i, j, k+1) != fuel_center) {
                    R_sum += R_arr(i, j, k+1);
                    weight_sum += Real(1.0);
                    ++n_neighbors;
                }
            }
            
            if (n_neighbors == 0) return;  // No different neighbors
            
            // Compute average neighbor ROS
            Real R_neighbor_avg = R_sum / weight_sum;
            
            // Blend factor: more neighbors with different fuel → more blending
            Real blend_factor = static_cast<Real>(n_neighbors) / Real(6.0);  // Normalize by max neighbors
            blend_factor = amrex::min(blend_factor, Real(1.0));
            
            // Apply blending based on method
            Real R_blended;
            if (blend_method == BlendMethod::LINEAR) {
                R_blended = blend_two_fuel_ros_linear(R_center, R_neighbor_avg, Real(1.0) - blend_factor);
            } else if (blend_method == BlendMethod::HARMONIC) {
                R_blended = blend_two_fuel_ros_harmonic(R_center, R_neighbor_avg, Real(1.0) - blend_factor);
            } else if (blend_method == BlendMethod::FINNEY) {
                R_blended = blend_two_fuel_ros_finney(R_center, R_neighbor_avg, Real(1.0) - blend_factor);
            } else {
                // MAXIMUM or unknown
                R_blended = amrex::max(R_center, R_neighbor_avg);
            }
            
            // Update ROS with smoothed value
            R_arr(i, j, k) = R_blended;
        });
    }
}

#endif // FUEL_BOUNDARY_SMOOTHING_H
