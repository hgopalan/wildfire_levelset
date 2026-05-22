// ============================================================================
// Cell Size Effects Correction for Fire Spread Rate
// ============================================================================
// Empirical correction for fire spread rates based on cell/grid resolution.
// Accounts for sub-grid fire behavior: finer grids capture more fire-front
// geometry, leading to higher effective spread rates.
//
// Physical Basis:
// Fire spread models like Rothermel assume a representative fuel bed at a
// given scale. When using coarser grids, numerical diffusion and averaging
// effects cause under-prediction of fire spread. Conversely, finer grids
// better resolve fire-front topology and propagate faster (more accurately).
//
// FARSITE empirical correction (Finney 1998):
// Cell size effects are modeled as: R_corrected = R_base × (dx_ref / dx_actual)^α
// where α is typically in [0.0, 0.3], empirically determined from observations.
//
// Default exponent α = 0.1 balances realism and stability:
//   - At 10m vs 30m:   correction factor = (30/10)^0.1 = 1.116 (+11.6%)
//   - At 30m vs 30m:   correction factor = (30/30)^0.1 = 1.0   (no change)
//   - At 100m vs 30m:  correction factor = (30/100)^0.1 = 0.885 (-11.5%)
//
// Reference:
//   Finney, M.A. (1998). FARSITE: Fire Area Simulator. USDA Forest Service
//   RMRS-RP-4.
// ============================================================================

#ifndef CELL_SIZE_CORRECTION_H
#define CELL_SIZE_CORRECTION_H

#include <AMReX.H>
#include <cmath>

using namespace amrex;

// ============================================================================
// Compute cell size correction factor
// ============================================================================
// Parameters:
//   dx_min         – minimum cell size in domain [m] (min of dx, dy [, dz])
//   dx_ref         – reference cell size [m] (typically 30.0)
//   correction_exp – power-law exponent (typically 0.1)
//
// Returns:
//   Multiplicative correction factor for ROS: correction = (dx_ref / dx_min)^exponent
//
// This function assumes all parameters are positive and valid (validation
// is performed in parse_inputs.cpp).
//
// GPU-compatible: can be called from device code
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real compute_cell_size_correction_factor(amrex::Real dx_min,
                                                amrex::Real dx_ref,
                                                amrex::Real correction_exp)
{
    // Guard against division by zero and invalid inputs
    if (dx_min <= 0.0 || dx_ref <= 0.0) {
        return 1.0;  // No correction if invalid inputs
    }
    
    // Compute ratio of reference to actual cell size
    amrex::Real ratio = dx_ref / dx_min;
    
    // Apply power-law correction: factor = ratio^exponent
    amrex::Real factor = std::pow(ratio, correction_exp);
    
    return factor;
}

#endif // CELL_SIZE_CORRECTION_H
