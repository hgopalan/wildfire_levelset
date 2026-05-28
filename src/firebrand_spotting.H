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
#include "scott_spotting_table.H"
#include <AMReX.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <random>
#include <vector>
#include <cmath>

using namespace amrex;

// ============================================================================
// Crown Fire Scaling Multiplier
// ============================================================================
// Computes intensity scaling factors for firebrand spotting based on crown
// fire activity state. FARSITE model increases spotting probability and
// distance during active crown fire.
//
// Crown activity values from ecology_mf component 3:
//   0 = surface fire     → multiplier = 1.0 (baseline)
//   1 = passive crown    → multiplier = 1.25 (25% enhancement)
//   2 = active crown     → multiplier = 1.75 (75% enhancement)
//
// This function is called per-cell to get the appropriate scaling factor.
// ============================================================================
static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_crown_fire_spotting_multiplier(amrex::Real crown_activity)
{
    // Round crown_activity to nearest integer: 0, 1, or 2
    // Using 2.5 upper bound ensures rounding works correctly for floating-point values
    // that may slightly exceed 2.0 due to numerical precision
    int ca = static_cast<int>(amrex::max(amrex::Real(0.0), amrex::min(amrex::Real(2.5), crown_activity)) + amrex::Real(0.5));
    ca = amrex::max(0, amrex::min(2, ca));
    
    if (ca == 0) {
        return amrex::Real(1.0);    // Surface fire: baseline spotting
    } else if (ca == 1) {
        return amrex::Real(1.25);   // Passive crown fire: 25% enhancement
    } else {
        return amrex::Real(1.75);   // Active crown fire: 75% enhancement
    }
}

// ============================================================================
// Compute per-cell spot-fire catching (ignition) probability
//
// Returns a spatially-varying raster:
//   P_catch_local(i,j) = P_catch × max(0, (M_x − M_f_local(i,j)) / M_x)
//
// The Anderson (1970) / Rothermel (1972) moisture-damping factor accounts for
// the fact that wet or near-extinction-moisture fuels are less likely to be
// ignited by a landed firebrand.  When spatial_moisture_mf is provided the
// per-cell 1-hr dead moisture M_d1 (component 0) is used; otherwise the
// global M_f from RothermelParams is applied uniformly.
//
// Parameters:
//   spot_catch_prob_mf  – output MultiFab (1 comp); filled in-place
//   rp                  – Rothermel parameters supplying M_x and M_f defaults
//   P_catch             – global max catching probability [0–1]
//   spatial_moisture_mf – optional per-cell moisture MultiFab (FMS_MOISTURE_NCOMP)
// ============================================================================
inline void compute_spot_catch_probability(
    MultiFab&            spot_catch_prob_mf,
    const RothermelParams& rp,
    Real                 P_catch,
    const MultiFab*      spatial_moisture_mf = nullptr)
{
    const Real M_x   = rp.M_x;
    const Real M_f   = rp.M_f;
    const Real P_c   = amrex::max(Real(0.0), amrex::min(Real(1.0), P_catch));
    const bool has_sm = (spatial_moisture_mf != nullptr);

    for (MFIter mfi(spot_catch_prob_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto cp = spot_catch_prob_mf.array(mfi);

        Array4<const Real> sm_arr;
        if (has_sm) sm_arr = spatial_moisture_mf->const_array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            Real M_f_local = has_sm ? sm_arr(i, j, k, 0) : M_f;
            Real f_moisture = (M_x > Real(0.0))
                ? amrex::max(Real(0.0), (M_x - M_f_local) / M_x)
                : Real(0.0);
            cp(i, j, k) = amrex::min(Real(1.0), P_c * f_moisture);
        });
    }
}


inline void compute_spotting_probability(
    MultiFab &spotting_data,
    const MultiFab &phi,
    const MultiFab &vel,
    const Geometry & /*geom*/,
    const RothermelParams &rp,
    const InputParameters::SpottingParams &sp,
    const MultiFab * /*slopes*/ = nullptr,
    const MultiFab* ecology_mf = nullptr)
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
    
    // Extract scalar fields from structs that contain std::string members.
    // SYCL requires that all values captured by a device lambda are trivially
    // copyable ("device copyable").  SpottingParams and RothermelParams each
    // contain std::string fields, so they cannot be captured directly.
    const Real sp_k_wind      = sp.k_wind;
    const Real sp_I_critical  = sp.I_critical;
    const Real sp_P_base      = sp.P_base;
    const Real rp_M_x         = rp.M_x;
    const Real rp_M_f         = rp.M_f;
    const Real rp_h_heat      = rp.h_heat;
    const bool has_ecology    = (ecology_mf != nullptr);

    // Loop over all boxes
    for (MFIter mfi(spotting_data); mfi.isValid(); ++mfi) {
        const Box &bx = mfi.validbox();
        auto const p = phi.const_array(mfi);
        auto const v = vel.const_array(mfi);
        auto spot = spotting_data.array(mfi);
        
        Array4<const Real> eco_arr;
        if (has_ecology) eco_arr = ecology_mf->const_array(mfi);
        
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
            Real f_wind = Real(1.0) - std::exp(-sp_k_wind * U_ftmin * U_ftmin / Real(100.0));
            
            // Moisture factor: f(moisture) = max(0, (M_x - M_f) / M_x)
            Real f_moisture = amrex::max(Real(0.0), (rp_M_x - rp_M_f) / rp_M_x);
            
            // Intensity factor: approximated from rate of spread
            // φ_w (wind factor from Rothermel)
            Real phi_w = (U_ftmin > 0.0) ? C * std::pow(U_ftmin, B) * beta_ratio_E : Real(0.0);
            Real R = R0 * (Real(1.0) + phi_w);
            Real intensity = R * rp_h_heat; // Simplified intensity (BTU/ft²/min)
            Real f_intensity = amrex::min(Real(1.0), intensity / sp_I_critical);
            
            // Crown fire scaling: increase spotting probability during active crown fire
            Real crown_mult = Real(1.0);
            if (has_ecology) {
                Real crown_activity = eco_arr(i, j, k, 3);
                crown_mult = compute_crown_fire_spotting_multiplier(crown_activity);
            }
            
            // Combined probability with crown fire enhancement
            Real P_spot = crown_mult * sp_P_base * f_wind * f_intensity * f_moisture;
            P_spot = amrex::min(Real(1.0), amrex::max(Real(0.0), P_spot));
            
            // Store probability in component 0
            spot(i, j, k, 0) = P_spot;
        });
    }
}

// ============================================================================
// Generate firebrand spots using stochastic sampling
// ============================================================================
inline void generate_firebrand_spots(
    MultiFab &phi,
    MultiFab &spotting_data,
    const MultiFab &vel,
    const Geometry &geom,
    const InputParameters::SpottingParams &sp,
    int step,
    const MultiFab* fuel_model_mf    = nullptr,  // per-cell fuel model code (optional)
    const std::string* fuel_system   = nullptr,  // "13" or "40" (optional)
    const MultiFab* residual_fuel_mf = nullptr,  // per-cell residual fuel fraction (optional)
    Real reentry_fuel_threshold      = Real(0.05), // min fuel fraction to allow ignition
    MultiFab* spotting_lineage_mf    = nullptr,  // lineage output: src (x,y) per landing cell
    const MultiFab* ecology_mf       = nullptr)  // per-cell ecology including crown_activity (optional)
{
    if (sp.enable == 0) {
        return; // Spotting disabled
    }
    
#if defined(AMREX_USE_OMP) || defined(AMREX_USE_GPU)
    amrex::Abort("Firebrand spotting data update requires serial execution. "
                 "Parallel execution (OMP/GPU) is not yet supported for spotting_data updates. "
                 "Please rebuild with AMReX_OMP=OFF, AMReX_GPU_BACKEND=NONE.");
#endif
    
    const Real phi_thresh = 0.1; // Threshold for identifying fire front
    const auto dx = geom.CellSize();
    const Real *prob_lo = geom.ProbLo();
    const Real *prob_hi = geom.ProbHi();
    Real domain_length = prob_hi[0] - prob_lo[0];
    const bool has_ecology = (ecology_mf != nullptr);
    
#if (AMREX_SPACEDIM == 3)
    GpuArray<Real, AMREX_SPACEDIM> gdx{dx[0], dx[1], dx[2]};
#else
    GpuArray<Real, AMREX_SPACEDIM> gdx{dx[0], dx[1]};
#endif
    
    // Setup random number generator with a seed that is identical on all MPI ranks.
    unsigned rng_seed;
    if (sp.random_seed == 0) {
#ifdef AMREX_USE_MPI
        if (ParallelDescriptor::IOProcessor()) {
            rng_seed = std::random_device{}();
        }
        MPI_Bcast(&rng_seed, 1,
                  ParallelDescriptor::Mpi_typemap<unsigned>::type(),
                  ParallelDescriptor::IOProcessorNumber(),
                  ParallelDescriptor::Communicator());
#else
        rng_seed = std::random_device{}();
#endif
    } else {
        rng_seed = static_cast<unsigned>(sp.random_seed + step);
    }
    std::mt19937 gen(rng_seed);
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
    std::vector<int> src_fuel_model;  // fuel model code at each source cell
    std::vector<Real> crown_activity_arr;  // crown activity at each source cell

    const bool has_fm_sp = (fuel_model_mf != nullptr);
    
    // First pass: collect fire front cells
    for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
        const Box &bx = mfi.validbox();
        auto const p = phi.const_array(mfi);
        auto const spot = spotting_data.const_array(mfi);
        auto const v = vel.const_array(mfi);

        Array4<const Real> fm_arr;
        if (has_fm_sp) fm_arr = fuel_model_mf->const_array(mfi);
        
        Array4<const Real> eco_arr;
        if (has_ecology) eco_arr = ecology_mf->const_array(mfi);

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
                        src_fuel_model.push_back(has_fm_sp
                            ? static_cast<int>(fm_arr(i, j, k)) : 0);
                        crown_activity_arr.push_back(has_ecology
                            ? eco_arr(i, j, k, 3) : Real(0.0));
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
    
#ifdef AMREX_USE_MPI
    // Gather fire-front data from all MPI ranks so every rank generates
    // the same global set of firebrand spots (deterministic with shared seed).
    {
        int local_n = static_cast<int>(xarray.size());
        int nprocs = ParallelDescriptor::NProcs();
        std::vector<int> counts(nprocs), displs(nprocs);
        MPI_Allgather(&local_n, 1, MPI_INT,
                      counts.data(), 1, MPI_INT,
                      ParallelDescriptor::Communicator());
        int total = 0;
        for (int r = 0; r < nprocs; ++r) { displs[r] = total; total += counts[r]; }

        auto allgatherv_real = [&](std::vector<Real>& v) {
            std::vector<Real> all(total);
            MPI_Allgatherv(v.data(), local_n,
                           ParallelDescriptor::Mpi_typemap<Real>::type(),
                           all.data(), counts.data(), displs.data(),
                           ParallelDescriptor::Mpi_typemap<Real>::type(),
                           ParallelDescriptor::Communicator());
            v = std::move(all);
        };
        allgatherv_real(xarray);
        allgatherv_real(yarray);
        allgatherv_real(probabilities);
        allgatherv_real(wind_x);
        allgatherv_real(wind_y);
        allgatherv_real(crown_activity_arr);
#if (AMREX_SPACEDIM == 3)
        allgatherv_real(zarray);
        allgatherv_real(wind_z);
#endif

        // allgather fuel model codes (int)
        {
            std::vector<int> all_fm(total);
            MPI_Allgatherv(src_fuel_model.data(), local_n, MPI_INT,
                           all_fm.data(), counts.data(), displs.data(), MPI_INT,
                           ParallelDescriptor::Communicator());
            src_fuel_model = std::move(all_fm);
        }
    }
#endif

    // Second pass: generate spots stochastically
    // Note: All arrays (xarray, yarray, probabilities, wind_x, wind_y, crown_activity_arr, etc.)
    // are populated in a synchronized manner during fire-front collection.
    // Array synchronization is maintained through a single collection loop and MPI AllGather operations.
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

            // Apply crown fire scaling to spotting distance
            // Active crown fire increases spotting distance based on crown_activity at source
            // crown_activity_arr is synchronized with xarray/yarray/probabilities (all populated together)
            Real crown_mult = amrex::Real(1.0);
            if (has_ecology && idx < static_cast<int>(crown_activity_arr.size())) {
                crown_mult = compute_crown_fire_spotting_multiplier(crown_activity_arr[idx]);
            }
            distance *= crown_mult;

            // ---- Scott/Albini (1979) per-fuel-model maximum spotting distance cap ----
            if (has_fm_sp && fuel_system != nullptr) {
                int fm_code = (idx < static_cast<int>(src_fuel_model.size()))
                              ? src_fuel_model[idx] : 0;
                if (fm_code > 0) {
                    float d_max = (*fuel_system == "40")
                                  ? get_max_spot_dist_fbfm40(fm_code)
                                  : get_max_spot_dist_fbfm13(fm_code);
                    if (d_max > 0.0f && distance > static_cast<Real>(d_max)) {
                        distance = static_cast<Real>(d_max);
                    }
                }
            }
            
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

    // ---- Catching (ignition) probability filter ----
    // When sp.P_catch < 1.0, each landed firebrand ignites the receiving fuel
    // only with probability P_catch × f_moisture, where f_moisture mirrors the
    // Rothermel (1972) moisture damping factor:
    //   f_moisture = max(0, (M_x - M_f) / M_x)
    // This is the Anderson (1970) probability-of-ignition model applied to
    // spot-fire catching: wet or near-extinction-moisture fuels are less likely
    // to sustain ignition from a landed brand.
    // P_catch = 1.0 (default) → always catch; original behaviour is preserved.
    // Apply the filter by removing spots that fail the catch test.
    if (sp.P_catch < Real(1.0)) {
        // Compute receiver moisture factor from SpottingParams
        // (M_x and M_f are captured indirectly via the probability already
        //  factored into spotting_data; here we apply only the P_catch gate)
        const Real P_c = std::max(Real(0.0), std::min(Real(1.0), sp.P_catch));
        std::vector<Real> new_spot_x, new_spot_y, new_spot_dist;
        std::vector<int>  new_src_i,  new_src_j,  new_src_k;
#if (AMREX_SPACEDIM == 3)
        std::vector<Real> new_spot_z;
#endif
        std::uniform_real_distribution<Real> catch_roll(0.0, 1.0);
        for (int idx = 0; idx < total_spots; ++idx) {
            if (catch_roll(gen) < P_c) {
                new_spot_x.push_back(spot_x[idx]);
                new_spot_y.push_back(spot_y[idx]);
                new_spot_dist.push_back(spot_distances[idx]);
                new_src_i.push_back(source_i[idx]);
                new_src_j.push_back(source_j[idx]);
                new_src_k.push_back(source_k[idx]);
#if (AMREX_SPACEDIM == 3)
                new_spot_z.push_back(spot_z[idx]);
#endif
            }
        }
        spot_x        = std::move(new_spot_x);
        spot_y        = std::move(new_spot_y);
        spot_distances= std::move(new_spot_dist);
        source_i      = std::move(new_src_i);
        source_j      = std::move(new_src_j);
        source_k      = std::move(new_src_k);
#if (AMREX_SPACEDIM == 3)
        spot_z        = std::move(new_spot_z);
#endif
        const int caught = static_cast<int>(spot_x.size());
        amrex::Print() << "  After P_catch=" << P_c << " filter: "
                       << caught << " of " << total_spots << " spots caught fire\n";
        total_spots = caught;
    }

    // ---- Post-fire re-entry fuel adjustment ----
    // When residual_fuel_mf is provided, spots landing in previously-burned cells
    // with low residual fuel are filtered out or have reduced ignition probability.
    // This prevents spot fires from igniting in areas where fuel has been depleted.
    if (residual_fuel_mf != nullptr && total_spots > 0) {
        std::vector<Real> new_spot_x, new_spot_y, new_spot_dist;
        std::vector<int>  new_src_i,  new_src_j,  new_src_k;
#if (AMREX_SPACEDIM == 3)
        std::vector<Real> new_spot_z;
#endif
        std::uniform_real_distribution<Real> fuel_roll(0.0, 1.0);

        for (int idx = 0; idx < total_spots; ++idx) {
            // Find the landing cell indices
            int li = static_cast<int>(std::floor((spot_x[idx] - prob_lo[0]) / gdx[0]));
            int lj = static_cast<int>(std::floor((spot_y[idx] - prob_lo[1]) / gdx[1]));

            // Look up residual fuel at the landing cell (only on IO processor for serial code)
            Real f_res = Real(1.0);
            if (ParallelDescriptor::IOProcessor()) {
                for (MFIter mfi(*residual_fuel_mf); mfi.isValid(); ++mfi) {
                    const Box& bx = mfi.validbox();
                    if (bx.contains(IntVect(li, lj))) {
                        auto const rf = residual_fuel_mf->array(mfi);
                        f_res = rf(li, lj, 0);
                        break;
                    }
                }
            }
            // Broadcast f_res to all ranks
#ifdef AMREX_USE_MPI
            MPI_Bcast(&f_res, 1,
                      ParallelDescriptor::Mpi_typemap<Real>::type(),
                      ParallelDescriptor::IOProcessorNumber(),
                      ParallelDescriptor::Communicator());
#endif
            // Skip spots in cells below the threshold (no fuel to burn)
            if (f_res < reentry_fuel_threshold) continue;

            // Scale P_catch by residual fuel fraction: more depleted → less likely to ignite
            if (fuel_roll(gen) < f_res) {
                new_spot_x.push_back(spot_x[idx]);
                new_spot_y.push_back(spot_y[idx]);
                new_spot_dist.push_back(spot_distances[idx]);
                new_src_i.push_back(source_i[idx]);
                new_src_j.push_back(source_j[idx]);
                new_src_k.push_back(source_k[idx]);
#if (AMREX_SPACEDIM == 3)
                new_spot_z.push_back(spot_z[idx]);
#endif
            }
        }
        const int kept = static_cast<int>(new_spot_x.size());
        if (kept < total_spots)
            amrex::Print() << "  Re-entry fuel filter: " << kept << " of "
                           << total_spots << " spots in fueled area\n";
        spot_x         = std::move(new_spot_x);
        spot_y         = std::move(new_spot_y);
        spot_distances = std::move(new_spot_dist);
        source_i       = std::move(new_src_i);
        source_j       = std::move(new_src_j);
        source_k       = std::move(new_src_k);
#if (AMREX_SPACEDIM == 3)
        spot_z         = std::move(new_spot_z);
#endif
        total_spots = kept;
    }

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
    // Iterate over boxes in outer loop (better cache locality than iterating over spots first)
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
                
                // Store maximum distance (component 2) - shows farthest firebrand travel from this source cell
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

    // ---- Lineage tracking: record source-cell coordinates at each landing cell ----
    // Components 0,1 of spotting_lineage_mf: (firebrand_src_x, firebrand_src_y)
    // When multiple brands land in the same cell the most recent source wins.
    if (spotting_lineage_mf != nullptr && total_spots > 0) {
        for (MFIter mfi(*spotting_lineage_mf); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            auto lin = spotting_lineage_mf->array(mfi);
            const auto lo = lbound(bx);
            const auto hi = ubound(bx);

            for (int idx = 0; idx < total_spots; ++idx) {
                Real x_spot = spot_x[idx];
                Real y_spot = spot_y[idx];

                int i_land = static_cast<int>(std::floor((x_spot - prob_lo[0]) / gdx[0]));
                int j_land = static_cast<int>(std::floor((y_spot - prob_lo[1]) / gdx[1]));
                int k_land = 0;

                if (i_land >= lo.x && i_land <= hi.x &&
                    j_land >= lo.y && j_land <= hi.y &&
                    k_land >= lo.z && k_land <= hi.z) {
                    // Source cell physical centre
                    Real src_x = prob_lo[0] + (source_i[idx] + Real(0.5)) * gdx[0];
                    Real src_y = prob_lo[1] + (source_j[idx] + Real(0.5)) * gdx[1];
                    lin(i_land, j_land, k_land, 0) = src_x;
                    lin(i_land, j_land, k_land, 1) = src_y;
                }
            }
        }
    }
}

#endif // FIREBRAND_SPOTTING_H
