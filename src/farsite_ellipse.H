// ============================================================================
// FARSITE Ellipse Model (Richards 1990)
// ============================================================================
// This implements Richards' (1990) elliptical fire spread model with
// coefficients a, b, and c for directional spread rates.
//
// The model can use either:
// 1. Fixed Richards' coefficients (a, b, c) from input parameters
// 2. Anderson (1983) L/W ratio based on wind speed, converted to coefficients
// 3. Van Wagner (1977) crown fire initiation to modify spread rates
//
// The model:
// 1. Finds locations where phi ≈ 0 (fire front)
// 2. Uses the pre-computed ROS field (R_mf) from the configured firespread model
// 3. Applies Van Wagner crown initiation if enabled
// 4. Applies Richards' elliptical expansion with coefficients:
//    - a: head fire coefficient (downwind, maximum spread)
//    - b: flank fire coefficient (crosswind spread)
//    - c: backing fire coefficient (upwind, minimum spread)
// 5. Computes directional spread based on angle between fire normal and wind
// 6. Optionally computes bulk fuel consumption fraction
// ============================================================================

#ifndef FARSITE_ELLIPSE_H
#define FARSITE_ELLIPSE_H

#include "fire_models.H"
#include "parse_inputs.H"
#include "rothermel_model.H"
#include "bulk_fuel_consumption.H"
#include "crown_initiation.H"
#include "cruz_crown_model.H"
#include <AMReX.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>

using namespace amrex;

// Constant for phi reset multiplier
// When resetting phi at fire front positions in FARSITE mode, set to -1 (fire)
constexpr Real PHI_RESET_MULTIPLIER = -1.0;

// ============================================================================
// Find locations where phi is close to 0 (fire front) and compute
// ellipse-based spread using Richards' (1990) model with coefficients a, b, c
// combined with the pre-computed ROS field R_mf (which may come from any
// configured fire-spread model: Rothermel, Balbi, FBP, Cheney-Gould, etc.).
// Optionally applies Van Wagner (1977) crown fire initiation with either
//   Van Wagner 3/CBD proxy (use_cruz_crown=0) or Cruz et al. (2005) active
//   crown ROS (use_cruz_crown=1, requires ccc to be pre-computed).
// Optionally computes per-class fuel consumption fraction if enabled.
// Optionally computes per-cell WAF from canopy height when canopy_height_mf
//   is provided (replaces global WAF + cover-weighted blending).
// ============================================================================
inline void compute_farsite_spread(MultiFab &phi, const MultiFab &vel,
                                   MultiFab &farsite_spread,
                                   const Geometry &geom, Real dt,
                                   const RothermelParams &rp,
                                   const InputParameters::FARSITEParams &fp,
                                   const InputParameters::CrownInitiationParams &cp,
                                   const MultiFab &R_mf,
                                   const MultiFab * /*slopes*/ = nullptr,
                                   MultiFab *fuel_consumption = nullptr,
                                   MultiFab *crown_fire_fraction = nullptr,
                                   const MultiFab *canopy_cover_mf = nullptr,
                                   const MultiFab *canopy_height_mf = nullptr,
                                   const CruzCrownComputed *ccc = nullptr) {
  // (FARSITE is always active when propagation_method = farsite; guard removed)


  const auto dx = geom.CellSize();
#if (AMREX_SPACEDIM == 3)
  GpuArray<Real, AMREX_SPACEDIM> gdx{dx[0], dx[1], dx[2]};
#else
  GpuArray<Real, AMREX_SPACEDIM> gdx{dx[0], dx[1]};
#endif

  const Real phi_thresh = fp.phi_threshold;

  // Dynamic arrays to collect positions where phi < threshold
  std::vector<Real> xarray;
  std::vector<Real> yarray;
#if (AMREX_SPACEDIM == 3)
  std::vector<Real> zarray;
#endif

  // Get physical domain information for coordinate conversion
  const Real *prob_lo = geom.ProbLo();

  const int use_anderson_LW = fp.use_anderson_LW;

  // Fire shape model: convert string to int for GPU-safe capture
  // 0 = richards (default), 1 = catchpole_demestre, 2 = wilson, 3 = lemniscate
  const int fire_shape_id = (fp.fire_shape_model == "catchpole_demestre") ? 1
                          : (fp.fire_shape_model == "wilson")             ? 2
                          : (fp.fire_shape_model == "lemniscate")         ? 3
                          :                                                  0;

  // Richards' (1990) elliptical expansion coefficients
  // These may be overridden by Anderson L/W ratio if enabled
  Real coeff_a = fp.coeff_a; // head fire coefficient
  Real coeff_b = fp.coeff_b; // flank fire coefficient
  Real coeff_c = fp.coeff_c; // backing fire coefficient

  // Compute Rothermel (1972) fire spread parameters
  // (used for crown fire intensity, WAF, and wind-unit conversion only;
  //  the base ROS comes from the pre-computed R_mf field)
  const RothermelComputed rc = compute_rothermel_params(rp);
  const Real I_R = rc.I_R; // Reaction intensity for crown fire intensity
  const Real wind_conv = rc.wind_conv;
  const Real ros_conv = rc.ros_conv;
  // Andrews (2018) wind adjustments
  const Real waf_default         = rc.waf;
  const Real U_max_default       = rc.U_max_ftmin;
  const int  use_waf_flag        = rp.use_waf;
  const int  use_wind_limit_flag = rp.use_wind_limit;

  // Per-cell canopy WAF flag
  const bool use_canopy_ht_waf = (canopy_height_mf != nullptr && canopy_cover_mf != nullptr && use_waf_flag == 1);
  const bool use_cc_waf        = (!use_canopy_ht_waf && canopy_cover_mf != nullptr && use_waf_flag == 1);
  // Global fuel-bed depth [ft] for WAF baseline
  const Real delta_ft_global = rp.delta;
  
  // Cruz crown ROS: enabled when ccc pointer is non-null and crown.use_cruz_crown = 1
  const int use_cruz_crown_flag = cp.use_cruz_crown;
  const bool have_cruz = (ccc != nullptr && use_cruz_crown_flag == 1);
  // Copy CruzCrownComputed scalars by value only when Cruz model is active,
  // so the GPU kernel captures a plain struct (no host pointer indirection).
  CruzCrownComputed ccc_local;
  if (have_cruz) { ccc_local = *ccc; }
  
  // Bulk fuel consumption model parameters
  const int use_fuel_consumption = fp.use_bulk_fuel_consumption;
  const Real tau_residence = fp.tau_residence;
  const Real f_c_min = fp.f_consumed_min;
  const Real f_c_max = fp.f_consumed_max;
  // Per-class fuel load and SAV scalars for per-class burnout model
  const Real rp_w_d1    = rp.w_d1;
  const Real rp_sigma_d1= rp.sigma_d1;
  const Real rp_w_d10   = rp.w_d10;
  const Real rp_w_d100  = rp.w_d100;
  const Real rp_w_lh    = rp.w_lh;
  const Real rp_sigma_lh= rp.sigma_lh;
  const Real rp_w_lw    = rp.w_lw;
  const Real rp_sigma_lw= rp.sigma_lw;
  const Real rp_w0      = rp.w0;
  const Real rp_sigma   = rp.sigma;
  const Real rp_rho_p   = rp.rho_p;
  
  // Crown fire initiation model parameters
  const int use_crown_initiation = cp.enable;
  const Real CBH = cp.CBH;
  const Real CBD = cp.CBD;
  const Real FMC = cp.FMC;
  const Real crown_fraction_weight = cp.crown_fraction_weight;
  const int use_metric_units = cp.use_metric_units;

  // Feature 11: critical active crown ROS [m/s] (Van Wagner 1977)
  // R'_SA = 3.0/CBD [m/min] → [m/s] = 3.0/CBD/60
  // Crown fire is ACTIVE only when R_surface >= R_SA_crit;
  // when I_B >= I_o but R < R_SA_crit it is PASSIVE (surface ROS only).
  const Real R_SA_crit = Real(3.0) / amrex::max(Real(CBD), Real(0.01)) / Real(60.0);

  // Feature 12: ellipse scaling for active crown fire (backward compatible; default off)
  const int scale_ellipse_crown_flag = fp.scale_ellipse_with_crown;
  const Real crown_lw_scale_val      = fp.crown_lw_scale;

  // Passive crown fire ellipse coefficients (Feature 9 in the FARSITE audit).
  // When crown fire is detected (crowning fraction CF > 0), blend the ellipse
  // a/b/c coefficients between surface and crown values proportionally to CF.
  // Use cp.a_crown / b_crown / c_crown when scale_ellipse_crown_flag = 1.
  const Real crown_coeff_a = cp.a_crown;
  const Real crown_coeff_b = cp.b_crown;
  const Real crown_coeff_c = cp.c_crown;

  // Compute critical surface intensity for crown initiation (if enabled)
  Real I_o = 0.0;
  if (use_crown_initiation == 1) {
    I_o = compute_critical_surface_intensity(CBH, FMC);
  }

  // -------------------------------------------------------------------------
  // First pass: compute elliptical spread vector at each fire-front cell
  // -------------------------------------------------------------------------
  for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
    const Box &bx = mfi.validbox();
    auto const p = phi.const_array(mfi);
    auto const v = vel.const_array(mfi);
    auto spread = farsite_spread.array(mfi);
    auto const R_arr = R_mf.const_array(mfi);
    
    // Get fuel consumption array if provided
    const bool compute_fuel_consumption = (fuel_consumption != nullptr && use_fuel_consumption == 1);
    Array4<Real> fc;
    if (compute_fuel_consumption) {
      fc = fuel_consumption->array(mfi);
    }
    
    // Get crown fire fraction array if provided
    const bool should_compute_crown_fraction = (crown_fire_fraction != nullptr && use_crown_initiation == 1);
    Array4<Real> cff;
    if (should_compute_crown_fraction) {
      cff = crown_fire_fraction->array(mfi);
    }

    // Get canopy cover and height arrays for per-cell WAF (if provided).
    // Safety invariant: cht_arr is only dereferenced in the GPU kernel when
    // use_canopy_ht_waf is true (a const bool captured by value).  That same
    // bool gates the initialization below, so cht_arr is never read in a
    // default-constructed (null-pointer) state.
    Array4<const Real> cc_arr;
    if (use_cc_waf || use_canopy_ht_waf) {
      cc_arr = canopy_cover_mf->const_array(mfi);
    }
    Array4<const Real> cht_arr;
    if (use_canopy_ht_waf) {
      cht_arr = canopy_height_mf->const_array(mfi);
    }

    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      // Check if this cell is near the fire front (phi ≈ 0)
      Real phi_val = p(i, j, k);

      // Initialize spread to zero
#if (AMREX_SPACEDIM == 3)
      spread(i, j, k, 0) = 0.0;
      spread(i, j, k, 1) = 0.0;
      spread(i, j, k, 2) = 0.0;
#else
            spread(i,j,k,0) = 0.0;
            spread(i,j,k,1) = 0.0;
#endif

      // Only compute spread for cells near the fire front
      if (phi_val > phi_thresh) {
        return;
      }

      // Compute gradient of phi to get normal direction
      Real dphi_dx = (p(i + 1, j, k) - p(i - 1, j, k)) / (2.0 * gdx[0]);
      Real dphi_dy = (p(i, j + 1, k) - p(i, j - 1, k)) / (2.0 * gdx[1]);
#if (AMREX_SPACEDIM == 3)
      Real dphi_dz = (p(i, j, k + 1) - p(i, j, k - 1)) / (2.0 * gdx[2]);
      Real grad_mag =
          std::sqrt(dphi_dx * dphi_dx + dphi_dy * dphi_dy + dphi_dz * dphi_dz);
#else
            Real grad_mag = std::sqrt(dphi_dx*dphi_dx + dphi_dy*dphi_dy);
#endif

      // Avoid division by zero
      if (grad_mag < 1.0e-50) {
        return;
      }

      // Normal vector (points toward unburned region, phi > 0)
      Real nx = dphi_dx / grad_mag;
      Real ny = dphi_dy / grad_mag;
#if (AMREX_SPACEDIM == 3)
      Real nz = dphi_dz / grad_mag;
#endif

      // Get wind velocity components
      Real ux = v(i, j, k, 0);
      Real uy = v(i, j, k, 1);
#if (AMREX_SPACEDIM == 3)
      Real uz = v(i, j, k, 2);
#endif

      // Wind speed magnitude in simulation units
      Real wind_mag = std::sqrt(ux * ux + uy * uy
#if (AMREX_SPACEDIM == 3)
                                + uz * uz
#endif
      );

      // Convert wind speed to Rothermel units (ft/min)
      Real U_ftmin = wind_mag * wind_conv;

      // Apply Wind Adjustment Factor if enabled (Andrews 2018):
      // converts 20-ft open wind to midflame height wind.
      // When canopy_height_mf is supplied: compute per-cell WAF from blended
      //   fuel-bed depth (open) + canopy height (covered areas).
      // When only canopy_cover_mf is supplied: blend global WAF with cover fraction.
      if (use_waf_flag == 1) {
          Real waf_cell;
          if (use_canopy_ht_waf) {
              const Real cc_frac = amrex::max(Real(0.0),
                                   amrex::min(Real(100.0), cc_arr(i,j,k))) / Real(100.0);
              const Real canopy_ht_ft = cht_arr(i,j,k) / Real(0.3048); // m → ft
              const Real h_eff = (Real(1.0) - cc_frac) * delta_ft_global
                                + cc_frac * amrex::max(canopy_ht_ft, delta_ft_global);
              waf_cell = compute_waf_unsheltered(h_eff);
          } else if (use_cc_waf) {
              Real cc_frac = amrex::max(Real(0.0),
                             amrex::min(Real(100.0), cc_arr(i,j,k))) / Real(100.0);
              waf_cell = Real(1.0) - cc_frac * (Real(1.0) - waf_default);
          } else {
              waf_cell = waf_default;
          }
          U_ftmin *= waf_cell;
      }

      // Apply maximum effective wind speed limit if enabled (Andrews 2018)
      if (use_wind_limit_flag == 1 && U_max_default > Real(0.0)) {
          U_ftmin = std::min(U_ftmin, U_max_default);
      }

      // Wind factor and slope factor are already baked into R_mf by the
      // configured firespread model (Rothermel, Balbi, FBP, etc.).
      // Read the pre-computed surface ROS [simulation units, m/s] directly.
      Real R_surface = R_arr(i, j, k);
      // Convert to ft/min for crown fire fireline-intensity computation.
      Real R_surface_ftmin = R_surface / ros_conv;
      
      // Apply Van Wagner crown fire initiation if enabled
      Real R_head = R_surface;  // Default: surface fire rate
      Real local_crown_fraction = 0.0;
      bool crown_is_active = false;  // Feature 11/12: track active crown state

      if (use_crown_initiation == 1) {
        // Compute fireline intensity from Rothermel reaction intensity
        // I_B [BTU/ft/s] = I_R * R_surface_ftmin / 60
        Real I_B_imperial = I_R * R_surface_ftmin / Real(60.0);

        // Convert to metric if needed for Van Wagner model
        Real I_B = use_metric_units ? convert_intensity_to_metric(I_B_imperial) : I_B_imperial;

        // Check if crown fire initiates (I_B >= I_o)
        bool crown_initiated = will_crown_initiate(I_B, I_o);

        if (crown_initiated) {
          // Feature 11: Distinguish passive crown (I_B >= I_o but R < R'_SA)
          // from active crown (I_B >= I_o AND R >= R'_SA).
          // Passive crown fire uses surface ROS — no crown-spread enhancement.
          crown_is_active = (R_surface >= R_SA_crit);

          if (crown_is_active) {
            // Active crown fire: compute enhanced crown ROS [m/s]
            Real R_crown;
            if (have_cruz) {
              // Cruz, Alexander & Wakimoto (2005)
              R_crown = cruz_crown_ros(ccc_local, wind_mag);  // already in m/s
            } else {
              // Van Wagner (1977) simplified proxy
              Real R_crown_mmin = compute_crown_fire_spread_rate(CBD, FMC);
              Real R_crown_ftmin = R_crown_mmin * Real(3.28084);
              R_crown = R_crown_ftmin * ros_conv;
            }

            // Compute combined rate (crown fire dominates)
            R_head = compute_total_spread_rate(R_surface, R_crown, true, crown_fraction_weight);

            // Compute crown fire fraction for output
            local_crown_fraction = compute_crown_fraction(R_crown, R_surface, true);
          }
          // Passive crown: R_head = R_surface, local_crown_fraction = 0
        }
      }

      // Determine elliptical coefficients
      Real local_coeff_a = coeff_a;
      Real local_coeff_b = coeff_b;
      Real local_coeff_c = coeff_c;

      if (use_anderson_LW == 1) {
        // ft/min to mph: multiply by 60/5280 = 0.01136364
        Real U_mph = U_ftmin * Real(0.01136364);
        Real LW_ratio = anderson_LW_ratio(U_mph);
        // Feature 12: scale L/W ratio for active crown fire (backward compat; off by default)
        if (scale_ellipse_crown_flag == 1 && crown_is_active) {
          LW_ratio *= crown_lw_scale_val;
        }
        LW_ratio_to_richards_coefficients(LW_ratio, local_coeff_a,
                                          local_coeff_b, local_coeff_c);
      } else if (scale_ellipse_crown_flag == 1) {
        if (crown_is_active) {
          // Feature 12: scale fixed coefficients for active crown.
          // Use head/backing ratio as a proxy for L/W, then reconstruct a/b/c.
          const Real LW_approx = (local_coeff_c > Real(1.0e-10))
                                  ? local_coeff_a / local_coeff_c : Real(3.0);
          const Real LW_crown = LW_approx * crown_lw_scale_val;
          LW_ratio_to_richards_coefficients(LW_crown, local_coeff_a,
                                            local_coeff_b, local_coeff_c);
        } else if (local_crown_fraction > Real(0.0)) {
          // Feature 9 (passive crown fire ellipse): blend a/b/c between surface
          // and crown coefficients proportionally to the crowning fraction CF.
          const Real CF = amrex::min(local_crown_fraction, Real(1.0));
          local_coeff_a = coeff_a + CF * (crown_coeff_a - coeff_a);
          local_coeff_b = coeff_b + CF * (crown_coeff_b - coeff_b);
          local_coeff_c = coeff_c + CF * (crown_coeff_c - coeff_c);
        }
      }

      // Wind unit vector
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

      // cos(theta) = dot(normal, wind)
      Real cos_theta = nx * wx + ny * wy;
#if (AMREX_SPACEDIM == 3)
      cos_theta += nz * wz;
#endif
      // Clamp to [-1, 1] to guard against floating-point rounding before
      // computing sin_theta = sqrt(1 - cos²θ).  This is cheap and benefits
      // all shape models, including Richards (which also calls sqrt below).
      cos_theta = amrex::max(Real(-1.0), amrex::min(Real(1.0), cos_theta));
      Real sin_theta = std::sqrt(Real(1.0) - cos_theta * cos_theta);

      Real spread_dist = 0.0;
      if (wind_mag > 1.0e-10) {
        // Dispatch to selected fire shape model
        Real scale_factor = 0.0;
        if (fire_shape_id == 1) {
          // Catchpole & de Mestre (1986) double-ellipse
          scale_factor = catchpole_demestre_spread(cos_theta, sin_theta,
                                                   local_coeff_a, local_coeff_b,
                                                   local_coeff_c);
        } else if (fire_shape_id == 2) {
          // Wilson (1988) single ellipse from rear focus
          scale_factor = wilson_spread(cos_theta, local_coeff_a, local_coeff_c);
        } else if (fire_shape_id == 3) {
          // Alexander et al. lemniscate (Limaçon)
          scale_factor = alexander_lemniscate_spread(cos_theta, local_coeff_a,
                                                     local_coeff_b, local_coeff_c);
        } else {
          // Richards (1990) – default linear head/flank/backing combination
          if (cos_theta > 0.0) {
            scale_factor = local_coeff_a * cos_theta + local_coeff_b * sin_theta;
          } else {
            scale_factor = local_coeff_b * sin_theta - local_coeff_c * cos_theta;
          }
        }
        spread_dist = amrex::max(scale_factor, Real(0.0)) * R_head * dt;
      } else {
        spread_dist = R_head * dt;
      }

      spread(i, j, k, 0) = spread_dist * nx;
      spread(i, j, k, 1) = spread_dist * ny;
#if (AMREX_SPACEDIM == 3)
      spread(i, j, k, 2) = spread_dist * nz;
#endif

      // Compute fuel consumption fraction if enabled
      if (compute_fuel_consumption) {
        // Per-class exponential burnout model (Albini 1976):
        //   f_c_i = 1 - exp(-tau_res / tau_i)  where tau_i = BURNOUT_ALPHA * rho_p / sigma_i
        // Falls back to single-class (using aggregate sigma) when per-class loads are zero.
        Real f_c = compute_per_class_fuel_consumption_frac(
            tau_residence,
            rp_w_d1, rp_sigma_d1, rp_w_d10, rp_w_d100,
            rp_w_lh, rp_sigma_lh, rp_w_lw, rp_sigma_lw,
            rp_w0, rp_sigma, rp_rho_p);
        // Clamp to user-specified bounds
        f_c = amrex::max(f_c_min, amrex::min(f_c_max, f_c));
        fc(i, j, k) = f_c;
      }
      
      // Store crown fire fraction if enabled
      if (should_compute_crown_fraction) {
        cff(i, j, k) = local_crown_fraction;
      }
    });
  }
  amrex::Print()<<" Maximum Spread (m) = "<<farsite_spread.norminf(0)<<" "<<farsite_spread.norminf(1)
#if (AMREX_SPACEDIM == 3)
             <<" "<<farsite_spread.norminf(2) 
#endif  
              <<"\n";
  amrex::Print()<<" Spread Rate (m/s) = "<<R_mf.norminf(0)<<std::endl;
  // -------------------------------------------------------------------------
  // Synchronize so the GPU kernel results are visible on the host before
  // we read phi and farsite_spread in the collection loop below.
  // -------------------------------------------------------------------------
  Gpu::streamSynchronize();

  // -------------------------------------------------------------------------
  // Second pass: collect positions where phi < threshold and add spread
  // -------------------------------------------------------------------------
  for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
    const Box &bx = mfi.validbox();
    auto const spread = farsite_spread.const_array(mfi);

    const auto lo = lbound(bx);
    const auto hi = ubound(bx);

    for (int k = lo.z; k <= hi.z; ++k) {
      for (int j = lo.y; j <= hi.y; ++j) {
        for (int i = lo.x; i <= hi.x; ++i) {
          if (std::abs(spread(i, j, k,0)) > 0 ||
              std::abs(spread(i, j, k, 1)) > 0) {
            Real x = prob_lo[0] + (i + 0.5) * gdx[0];
            Real y = prob_lo[1] + (j + 0.5) * gdx[1];
            Real x_new = x + spread(i, j, k, 0);
            Real y_new = y + spread(i, j, k, 1);
            xarray.push_back(x_new);
            yarray.push_back(y_new);
#if (AMREX_SPACEDIM == 3)
            Real z_new = prob_lo[2] + (k + 0.5) * gdx[2] + spread(i, j, k, 2);
            zarray.push_back(z_new);
#endif
        }
      }
    }
  }
  }

#ifdef AMREX_USE_MPI
  // Gather spread positions from all MPI ranks so every rank stamps the
  // complete global fire front into its local boxes.
  {
    int local_n = static_cast<int>(xarray.size());
    int nprocs = ParallelDescriptor::NProcs();
    std::vector<int> counts(nprocs), displs(nprocs);
    MPI_Allgather(&local_n, 1, MPI_INT,
                  counts.data(), 1, MPI_INT,
                  ParallelDescriptor::Communicator());
    int total = 0;
    for (int r = 0; r < nprocs; ++r) { displs[r] = total; total += counts[r]; }

    std::vector<Real> all_x(total), all_y(total);
    MPI_Allgatherv(xarray.data(), local_n,
                   ParallelDescriptor::Mpi_typemap<Real>::type(),
                   all_x.data(), counts.data(), displs.data(),
                   ParallelDescriptor::Mpi_typemap<Real>::type(),
                   ParallelDescriptor::Communicator());
    MPI_Allgatherv(yarray.data(), local_n,
                   ParallelDescriptor::Mpi_typemap<Real>::type(),
                   all_y.data(), counts.data(), displs.data(),
                   ParallelDescriptor::Mpi_typemap<Real>::type(),
                   ParallelDescriptor::Communicator());
    xarray = std::move(all_x);
    yarray = std::move(all_y);
#if (AMREX_SPACEDIM == 3)
    std::vector<Real> all_z(total);
    MPI_Allgatherv(zarray.data(), local_n,
                   ParallelDescriptor::Mpi_typemap<Real>::type(),
                   all_z.data(), counts.data(), displs.data(),
                   ParallelDescriptor::Mpi_typemap<Real>::type(),
                   ParallelDescriptor::Communicator());
    zarray = std::move(all_z);
#endif
  }
#endif

  // -------------------------------------------------------------------------
  // Reset phi everywhere, then stamp new fire-front positions
  // -------------------------------------------------------------------------
  phi.setVal(0.0);
  const Real phi_new = PHI_RESET_MULTIPLIER;
  const int npts = static_cast<int>(xarray.size());
  amrex::Print() << "Updating phi at " << npts
                 << " positions based on FARSITE spread\n";

  // FIX: PODVector (Gpu::DeviceVector) does not support iterator-range
  // construction.  Allocate by size, then async-copy from host to device.
  Gpu::DeviceVector<Real> d_xarray(npts);
  Gpu::DeviceVector<Real> d_yarray(npts);
#if (AMREX_SPACEDIM == 3)
  Gpu::DeviceVector<Real> d_zarray(npts);
#endif

  Gpu::copyAsync(Gpu::hostToDevice, xarray.begin(), xarray.end(),
                 d_xarray.begin());
  Gpu::copyAsync(Gpu::hostToDevice, yarray.begin(), yarray.end(),
                 d_yarray.begin());
#if (AMREX_SPACEDIM == 3)
  Gpu::copyAsync(Gpu::hostToDevice, zarray.begin(), zarray.end(),
                 d_zarray.begin());
#endif
  Gpu::streamSynchronize(); // ensure copies are complete before ParallelFor

  const Real *px = d_xarray.dataPtr();
  const Real *py = d_yarray.dataPtr();
#if (AMREX_SPACEDIM == 3)
  const Real *pz = d_zarray.dataPtr();
#endif

  // Determine whether Gaussian/SDF smoothing is requested.
  // fp.gaussian_sigma < 0  →  original single-cell stamping (default).
  // fp.gaussian_sigma = 0  →  auto sigma = 3 * min(dx, dy [, dz]).
  // fp.gaussian_sigma > 0  →  user-supplied sigma [m].
  const Real raw_gsigma = fp.gaussian_sigma;
  const bool use_gaussian = (raw_gsigma >= Real(0.0));
  Real eff_sigma = raw_gsigma;
  if (use_gaussian && eff_sigma == Real(0.0)) {
#if (AMREX_SPACEDIM == 3)
    eff_sigma = Real(3.0) * amrex::min({gdx[0], gdx[1], gdx[2]});
#else
    eff_sigma = Real(3.0) * amrex::min(gdx[0], gdx[1]);
#endif
    amrex::Print() << "FARSITE Gaussian stamping: auto sigma = " << eff_sigma << " m\n";
  }

  for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
    const Box &bx = mfi.validbox();
    auto p = phi.array(mfi);

    if (!use_gaussian) {
      // ---------------------------------------------------------------
      // Original behaviour: stamp a single cell per propagated position.
      // ---------------------------------------------------------------
      // ParallelFor over collected points — write goes through the correct
      // (device) execution path in GPU builds.
      ParallelFor(npts, [=] AMREX_GPU_DEVICE(int idx) noexcept {
        Real x_pos = px[idx];
        Real y_pos = py[idx];
        int ii =
            static_cast<int>(std::round((x_pos - prob_lo[0]) / gdx[0] - 0.5));
        int jj =
            static_cast<int>(std::round((y_pos - prob_lo[1]) / gdx[1] - 0.5));
#if (AMREX_SPACEDIM == 3)
        Real z_pos = pz[idx];
        int kk =
            static_cast<int>(std::round((z_pos - prob_lo[2]) / gdx[2] - 0.5));
#else
              int kk = 0;
#endif
        if (bx.contains(IntVect(AMREX_D_DECL(ii, jj, kk)))) {
          p(ii, jj, kk) = phi_new;
        }
      });
    } else {
      // ---------------------------------------------------------------
      // Point-driven SDF stamping (replaces the previous O(cells × npts)
      // cell-driven approach).
      //
      // For each propagated fire-front point, update only the cells within
      // a bounding box of radius ceil(sig/dx)+2 around that point.
      // The per-cell SDF contribution is:
      //   phi_sdf = dist(cell_centre, point) - sig
      // and phi is updated as min(phi, phi_sdf).
      //
      // Complexity: O(npts × (sig/dx)²) vs. the previous O(cells × npts),
      // yielding orders-of-magnitude speedup when sigma << domain_size.
      //
      // AMReX ParallelFor kernels within the same GPU stream execute in
      // order (same-stream ordering guarantees), so the sequential
      // per-point min-update to p(i,j,k) is race-free without atomics.
      // ---------------------------------------------------------------
      const Real sig = eff_sigma;
      const int r_i = static_cast<int>(sig / gdx[0]) + 2;
      const int r_j = static_cast<int>(sig / gdx[1]) + 2;
#if (AMREX_SPACEDIM == 3)
      const int r_k = static_cast<int>(sig / gdx[2]) + 2;
#endif
      for (int n = 0; n < npts; ++n) {
        const Real xn = xarray[n];
        const Real yn = yarray[n];
#if (AMREX_SPACEDIM == 3)
        const Real zn = zarray[n];
        const int kk = static_cast<int>((zn - prob_lo[2]) / gdx[2]);
#else
        const int kk = 0;
#endif
        const int ii = static_cast<int>((xn - prob_lo[0]) / gdx[0]);
        const int jj = static_cast<int>((yn - prob_lo[1]) / gdx[1]);

        Box local_bx(IntVect(AMREX_D_DECL(ii - r_i, jj - r_j, kk - r_k)),
                     IntVect(AMREX_D_DECL(ii + r_i, jj + r_j, kk + r_k)));
        local_bx &= bx;
        if (local_bx.isEmpty()) { continue; }

        ParallelFor(local_bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          const Real x = prob_lo[0] + (i + Real(0.5)) * gdx[0];
          const Real y = prob_lo[1] + (j + Real(0.5)) * gdx[1];
          const Real ddx = x - xn;
          const Real ddy = y - yn;
#if (AMREX_SPACEDIM == 3)
          const Real z = prob_lo[2] + (k + Real(0.5)) * gdx[2];
          const Real ddz = z - zn;
          const Real phi_sdf = std::sqrt(ddx*ddx + ddy*ddy + ddz*ddz) - sig;
#else
          const Real phi_sdf = std::sqrt(ddx*ddx + ddy*ddy) - sig;
#endif
          p(i, j, k) = amrex::min(p(i, j, k), phi_sdf);
        });
      }
    }
  }

  }

#endif // FARSITE_ELLIPSE_H
