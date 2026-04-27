// ============================================================================
// Firebrand Spotting Model
// ============================================================================
// This implements a probability-based firebrand spotting model that generates
// new ignition points ahead of the fire front based on wind conditions and
// fire intensity. Works within the FARSITE framework.
//
// The model:
// 1. Identifies fire front locations where phi ≈ 0
// 2. Computes spotting probability based on wind speed, fire intensity, and fuel moisture
// 3. Generates new spot fires using stochastic sampling
// 4. Places new ignition points in the phi field
// ============================================================================

#ifndef FIREBRAND_SPOTTING_H
#define FIREBRAND_SPOTTING_H

#include "fire_models.H"
#include "parse_inputs.H"
#include "rothermel_model.H"
#include <AMReX.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_MultiFab.H>
#include <random>
#include <vector>
#include <cmath>

using namespace amrex;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ============================================================================
// Compute spotting probability field based on wind, fire intensity, and moisture
// ============================================================================
static void compute_spotting_probability(
    MultiFab &spotting_data,
    const MultiFab &phi,
    const MultiFab &vel,
    const Geometry &geom,
    const RothermelParams &rp,
    const InputParameters::SpottingParams &sp,
    const MultiFab *slopes = nullptr)
{
    if (sp.enable == 0) {
        return; // Spotting disabled
    }

    const Real phi_thresh = 0.1; // Threshold for identifying fire front
    
    // Compute Rothermel parameters for intensity calculation
    const RothermelComputed rc = compute_rothermel_params(rp);
    const Real R0 = rc.R0;
    const Real C = rc.C;
    const Real B = rc.B;
    const Real beta_ratio_E = rc.beta_ratio_E;
    const Real wind_conv = rc.wind_conv;
    
    // Loop over all boxes
    for (MFIter mfi(spotting_data); mfi.isValid(); ++mfi) {
        const Box &bx = mfi.validbox();
        auto const p = phi.const_array(mfi);
        auto const v = vel.const_array(mfi);
        auto spot = spotting_data.array(mfi);
        
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            // Initialize all components to zero
            spot(i, j, k, 0) = 0.0; // Probability
            spot(i, j, k, 1) = 0.0; // Count
            spot(i, j, k, 2) = 0.0; // Distance
            spot(i, j, k, 3) = 0.0; // Active flag
            
            // Check if this cell is near the fire front
            Real phi_val = p(i, j, k);
            if (phi_val > phi_thresh) {
                return; // Not near fire front
            }
            
            // Compute wind speed
            Real ux = v(i, j, k, 0);
            Real uy = v(i, j, k, 1);
#if (AMREX_SPACEDIM == 3)
            Real uz = v(i, j, k, 2);
            Real wind_mag = std::sqrt(ux * ux + uy * uy + uz * uz);
#else
            Real wind_mag = std::sqrt(ux * ux + uy * uy);
#endif
            
            // Convert wind speed to Rothermel units (ft/min)
            Real U_ftmin = wind_mag * wind_conv;
            
            // Wind factor for probability: f(wind) = 1 - exp(-k_wind * U^2)
            Real f_wind = Real(1.0) - std::exp(-sp.k_wind * U_ftmin * U_ftmin / Real(100.0));
            
            // Moisture factor: f(moisture) = max(0, (M_x - M_f) / M_x)
            Real f_moisture = amrex::max(Real(0.0), (rp.M_x - rp.M_f) / rp.M_x);
            
            // Intensity factor: approximated from rate of spread
            // φ_w (wind factor from Rothermel)
            Real phi_w = (U_ftmin > 0.0) ? C * std::pow(U_ftmin, B) * beta_ratio_E : Real(0.0);
            Real R = R0 * (Real(1.0) + phi_w);
            Real intensity = R * rp.h_heat; // Simplified intensity (BTU/ft²/min)
            Real f_intensity = amrex::min(Real(1.0), intensity / sp.I_critical);
            
            // Combined probability
            Real P_spot = sp.P_base * f_wind * f_intensity * f_moisture;
            P_spot = amrex::min(Real(1.0), amrex::max(Real(0.0), P_spot));
            
            // Store probability in component 0
            spot(i, j, k, 0) = P_spot;
        });
    }
}

// ============================================================================
// Generate firebrand spots using stochastic sampling
// ============================================================================
static void generate_firebrand_spots(
    MultiFab &phi,
    MultiFab &spotting_data,
    const MultiFab &vel,
    const Geometry &geom,
    const InputParameters::SpottingParams &sp,
    int step)
{
    if (sp.enable == 0) {
        return; // Spotting disabled
    }
    
#if defined(AMREX_USE_MPI) || defined(AMREX_USE_OMP) || defined(AMREX_USE_GPU)
    amrex::Abort("Firebrand spotting data update requires serial execution. "
                 "Parallel execution (MPI/OMP/GPU) is not yet supported for spotting_data updates. "
                 "Please rebuild with AMReX_MPI=OFF, AMReX_OMP=OFF, AMReX_GPU_BACKEND=NONE.");
#endif
    
    const Real phi_thresh = 0.1; // Threshold for identifying fire front
    const auto dx = geom.CellSize();
    const Real *prob_lo = geom.ProbLo();
    const Real *prob_hi = geom.ProbHi();
    Real domain_length = prob_hi[0] - prob_lo[0];
    
#if (AMREX_SPACEDIM == 3)
    GpuArray<Real, AMREX_SPACEDIM> gdx{dx[0], dx[1], dx[2]};
#else
    GpuArray<Real, AMREX_SPACEDIM> gdx{dx[0], dx[1]};
#endif
    
    // Setup random number generator
    std::random_device rd;
    std::mt19937 gen(sp.random_seed == 0 ? rd() : sp.random_seed + step);
    std::uniform_real_distribution<Real> uniform(0.0, 1.0);
    std::normal_distribution<Real> normal(0.0, 1.0);
    
    // For lognormal distribution: ln(d) ~ N(μ, σ²)
    // where μ = ln(d_mean) - 0.5*σ² to get mean = d_mean
    Real mu_log = std::log(sp.d_mean) - 0.5 * sp.d_sigma * sp.d_sigma;
    std::lognormal_distribution<Real> lognormal_dist(mu_log, sp.d_sigma);
    std::exponential_distribution<Real> exponential_dist(sp.d_lambda);
    
    // Collect fire front positions and their probabilities
    std::vector<Real> xarray, yarray;
#if (AMREX_SPACEDIM == 3)
    std::vector<Real> zarray;
#endif
    std::vector<Real> probabilities;
    std::vector<Real> wind_x, wind_y;
#if (AMREX_SPACEDIM == 3)
    std::vector<Real> wind_z;
#endif
    
    // First pass: collect fire front cells
    for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
        const Box &bx = mfi.validbox();
        auto const p = phi.const_array(mfi);
        auto const spot = spotting_data.const_array(mfi);
        auto const v = vel.const_array(mfi);
        
        const auto lo = lbound(bx);
        const auto hi = ubound(bx);
        
        for (int k = lo.z; k <= hi.z; ++k) {
            for (int j = lo.y; j <= hi.y; ++j) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    Real phi_val = p(i, j, k);
                    Real prob = spot(i, j, k, 0);
                    
                    if (phi_val <= phi_thresh && prob > 0.0) {
                        Real x = prob_lo[0] + (i + 0.5) * gdx[0];
                        Real y = prob_lo[1] + (j + 0.5) * gdx[1];
                        xarray.push_back(x);
                        yarray.push_back(y);
                        probabilities.push_back(prob);
                        wind_x.push_back(v(i, j, k, 0));
                        wind_y.push_back(v(i, j, k, 1));
#if (AMREX_SPACEDIM == 3)
                        Real z = prob_lo[2] + (k + 0.5) * gdx[2];
                        zarray.push_back(z);
                        wind_z.push_back(v(i, j, k, 2));
#endif
                    }
                }
            }
        }
    }
    
    // Second pass: generate spots stochastically
    std::vector<Real> spot_x, spot_y;
#if (AMREX_SPACEDIM == 3)
    std::vector<Real> spot_z;
#endif
    std::vector<Real> spot_distances; // Store distances for each spot
    std::vector<int> source_i, source_j, source_k; // Store source cell indices
    
    const int npts = static_cast<int>(xarray.size());
    int total_spots = 0;
    
    for (int idx = 0; idx < npts; ++idx) {
        Real prob = probabilities[idx];
        Real rand_val = uniform(gen);
        
        // Check if spotting occurs at this location
        if (rand_val < prob) {
            // Determine spotting distance
            Real distance;
            if (sp.distance_model == "lognormal") {
                distance = lognormal_dist(gen);
            } else {
                distance = exponential_dist(gen);
            }
            
            // Clamp distance to reasonable bounds
            distance = amrex::min(distance, Real(0.5) * domain_length); // Max 50% of domain
            distance = amrex::max(distance, sp.spot_radius * Real(2.0)); // Min 2x spot radius
            
            // Compute wind direction
            Real wx = wind_x[idx];
            Real wy = wind_y[idx];
#if (AMREX_SPACEDIM == 3)
            Real wz = wind_z[idx];
            Real wind_mag = std::sqrt(wx * wx + wy * wy + wz * wz);
#else
            Real wind_mag = std::sqrt(wx * wx + wy * wy);
#endif
            
            if (wind_mag < 1.0e-10) {
                continue; // No wind, no spotting
            }
            
            // Normalized wind direction
            wx /= wind_mag;
            wy /= wind_mag;
#if (AMREX_SPACEDIM == 3)
            wz /= wind_mag;
#endif
            
            // Add lateral dispersion (perpendicular to wind)
            Real lateral_angle_rad = sp.lateral_spread_angle * M_PI / 180.0;
            Real lateral_spread = std::tan(lateral_angle_rad) * distance;
            Real lateral_offset = normal(gen) * lateral_spread;
            
            // In 2D, perpendicular direction is (-wy, wx)
            Real perp_x = -wy;
            Real perp_y = wx;
            
            // Final spot location: downwind + lateral offset
            Real new_x = xarray[idx] + distance * wx + lateral_offset * perp_x;
            Real new_y = yarray[idx] + distance * wy + lateral_offset * perp_y;
            
            // Compute source cell indices
            int src_i = static_cast<int>(std::floor((xarray[idx] - prob_lo[0]) / gdx[0]));
            int src_j = static_cast<int>(std::floor((yarray[idx] - prob_lo[1]) / gdx[1]));
#if (AMREX_SPACEDIM == 3)
            int src_k = static_cast<int>(std::floor((zarray[idx] - prob_lo[2]) / gdx[2]));
#else
            int src_k = 0;
#endif
            
            spot_x.push_back(new_x);
            spot_y.push_back(new_y);
            spot_distances.push_back(distance);
            source_i.push_back(src_i);
            source_j.push_back(src_j);
            source_k.push_back(src_k);
#if (AMREX_SPACEDIM == 3)
            Real new_z = zarray[idx] + distance * wz;
            spot_z.push_back(new_z);
#endif
            total_spots++;
        }
    }
    
    amrex::Print() << "Generated " << total_spots << " firebrand spots at step " << step << "\n";
    
    // Third pass: apply spots to phi field by setting negative signed distance values
    for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
        const Box &bx = mfi.validbox();
        auto p = phi.array(mfi);
        
        // Loop over all generated spots
        for (int idx = 0; idx < total_spots; ++idx) {
            Real x_spot = spot_x[idx];
            Real y_spot = spot_y[idx];
#if (AMREX_SPACEDIM == 3)
            Real z_spot = spot_z[idx];
#endif
            
            // Find grid indices for spot center
            int i_center = static_cast<int>(std::floor((x_spot - prob_lo[0]) / gdx[0]));
            int j_center = static_cast<int>(std::floor((y_spot - prob_lo[1]) / gdx[1]));
#if (AMREX_SPACEDIM == 3)
            int k_center = static_cast<int>(std::floor((z_spot - prob_lo[2]) / gdx[2]));
#else
            int k_center = 0;
#endif
            
            // Determine radius in grid cells
            int radius_cells = static_cast<int>(std::ceil(sp.spot_radius / gdx[0])) + 1;
            
            const auto lo = lbound(bx);
            const auto hi = ubound(bx);
            
            // Create a circular/spherical ignition zone
            for (int k = amrex::max(k_center - radius_cells, lo.z); 
                     k <= amrex::min(k_center + radius_cells, hi.z); ++k) {
                for (int j = amrex::max(j_center - radius_cells, lo.y); 
                         j <= amrex::min(j_center + radius_cells, hi.y); ++j) {
                    for (int i = amrex::max(i_center - radius_cells, lo.x); 
                             i <= amrex::min(i_center + radius_cells, hi.x); ++i) {
                        
                        // Compute distance from spot center
                        Real x = prob_lo[0] + (i + 0.5) * gdx[0];
                        Real y = prob_lo[1] + (j + 0.5) * gdx[1];
                        Real dist = std::sqrt((x - x_spot) * (x - x_spot) + 
                                            (y - y_spot) * (y - y_spot));
#if (AMREX_SPACEDIM == 3)
                        Real z = prob_lo[2] + (k + 0.5) * gdx[2];
                        dist = std::sqrt(dist * dist + (z - z_spot) * (z - z_spot));
#endif
                        
                        // If within spot radius, set phi to negative
                        if (dist < sp.spot_radius) {
                            // Set to signed distance from spot center
                            p(i, j, k) = amrex::min(p(i, j, k), dist - sp.spot_radius);
                        }
                    }
                }
            }
        }
    }
    
    // Fourth pass: update spotting_data with count, distance, and active flag
    // Reset all spotting data components first
    for (MFIter mfi(spotting_data); mfi.isValid(); ++mfi) {
        const Box &bx = mfi.validbox();
        auto spot = spotting_data.array(mfi);
        
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            spot(i, j, k, 1) = 0.0; // Reset count
            spot(i, j, k, 2) = 0.0; // Reset distance
            spot(i, j, k, 3) = 0.0; // Reset active flag
        });
    }
    
    // Update spotting data based on generated spots
    // Iterate over boxes in outer loop for better performance
    // NOTE: This code assumes serial execution (AMReX_MPI=OFF, AMReX_OMP=OFF, AMReX_GPU=OFF).
    // If parallelization is enabled, atomic operations (Gpu::Atomic::Add, Gpu::Atomic::Max)
    // would be required to avoid race conditions when multiple threads update the same cell.
    for (MFIter mfi(spotting_data); mfi.isValid(); ++mfi) {
        const Box &bx = mfi.validbox();
        auto spot = spotting_data.array(mfi);
        
        const auto lo = lbound(bx);
        const auto hi = ubound(bx);
        
        // Process all spots for this box
        for (int idx = 0; idx < total_spots; ++idx) {
            int src_i = source_i[idx];
            int src_j = source_j[idx];
            int src_k = source_k[idx];
            Real distance = spot_distances[idx];
            
            Real x_spot = spot_x[idx];
            Real y_spot = spot_y[idx];
#if (AMREX_SPACEDIM == 3)
            Real z_spot = spot_z[idx];
#endif
            
            // Check if source cell is in this box
            if (src_i >= lo.x && src_i <= hi.x &&
                src_j >= lo.y && src_j <= hi.y &&
                src_k >= lo.z && src_k <= hi.z) {
                
                // Increment count (component 1)
                spot(src_i, src_j, src_k, 1) += 1.0;
                
                // Store maximum distance (component 2) - useful for visualization
                spot(src_i, src_j, src_k, 2) = amrex::max(spot(src_i, src_j, src_k, 2), distance);
            }
            
            // Mark cells where spot fires were placed (component 3)
            int i_center = static_cast<int>(std::floor((x_spot - prob_lo[0]) / gdx[0]));
            int j_center = static_cast<int>(std::floor((y_spot - prob_lo[1]) / gdx[1]));
#if (AMREX_SPACEDIM == 3)
            int k_center = static_cast<int>(std::floor((z_spot - prob_lo[2]) / gdx[2]));
#else
            int k_center = 0;
#endif
            
            // Check if spot center is in this box
            if (i_center >= lo.x && i_center <= hi.x &&
                j_center >= lo.y && j_center <= hi.y &&
                k_center >= lo.z && k_center <= hi.z) {
                
                // Mark as active spot fire location (component 3)
                spot(i_center, j_center, k_center, 3) = 1.0;
            }
        }
    }
}

#endif // FIREBRAND_SPOTTING_H
