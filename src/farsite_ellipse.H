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
// 2. Uses Rothermel (1972) to compute base rate of spread (ROS)
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
#include <AMReX.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_MultiFab.H>

using namespace amrex;

// Constant for phi reset multiplier
// When resetting phi at fire front positions in FARSITE mode, set to -1 (fire)
constexpr Real PHI_RESET_MULTIPLIER = -1.0;

// ============================================================================
// Find locations where phi is close to 0 (fire front) and compute
// ellipse-based spread using Richards' (1990) model with coefficients a, b, c
// combined with Rothermel (1972) ROS
// Optionally applies Van Wagner (1977) crown fire initiation
// Optionally computes bulk fuel consumption fraction if enabled
// ============================================================================
static void compute_farsite_spread(MultiFab &phi, const MultiFab &vel,
                                   MultiFab &farsite_spread,
                                   const Geometry &geom, Real dt,
                                   const RothermelParams &rp,
                                   const InputParameters::FARSITEParams &fp,
                                   const InputParameters::CrownInitiationParams &cp,
                                   const MultiFab *slopes = nullptr,
                                   MultiFab *fuel_consumption = nullptr,
                                   MultiFab *crown_fire_fraction = nullptr) {
  if (fp.enable == 0) {
    return; // FARSITE disabled
  }

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

  // Richards' (1990) elliptical expansion coefficients
  // These may be overridden by Anderson L/W ratio if enabled
  Real coeff_a = fp.coeff_a; // head fire coefficient
  Real coeff_b = fp.coeff_b; // flank fire coefficient
  Real coeff_c = fp.coeff_c; // backing fire coefficient

  // Compute Rothermel (1972) fire spread parameters
  const RothermelComputed rc = compute_rothermel_params(rp);
  const Real R0 = rc.R0;
  const Real I_R = rc.I_R; // Reaction intensity for fuel consumption
  const Real C = rc.C;
  const Real B = rc.B;
  const Real beta_ratio_E = rc.beta_ratio_E;
  const Real phi_s_const =
      rc.phi_s; // Constant phi_s for backward compatibility
  const Real beta =
      rc.beta; // Packing ratio for terrain-based phi_s calculation
  const Real wind_conv = rc.wind_conv;
  const Real ros_conv = rc.ros_conv;
  
  // Bulk fuel consumption model parameters
  const int use_fuel_consumption = fp.use_bulk_fuel_consumption;
  const Real tau_residence = fp.tau_residence;
  const Real f_c_min = fp.f_consumed_min;
  const Real f_c_max = fp.f_consumed_max;
  
  // Crown fire initiation model parameters
  const int use_crown_initiation = cp.enable;
  const Real CBH = cp.CBH;
  const Real CBD = cp.CBD;
  const Real FMC = cp.FMC;
  const Real crown_fraction_weight = cp.crown_fraction_weight;
  const int use_metric_units = cp.use_metric_units;
  
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

    // Get slopes array if terrain file was provided
    const bool use_terrain = (slopes != nullptr);
    Array4<const Real> s;
    if (use_terrain) {
      s = slopes->const_array(mfi);
    }
    
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

      // Wind factor φ_w (Rothermel Eq. 47)
      Real phi_w =
          (U_ftmin > 0.0) ? C * std::pow(U_ftmin, B) * beta_ratio_E : 0.0;

      // Slope factor φ_s: use terrain-based if available, otherwise constant
      Real phi_s;
      if (use_terrain) {
        // Using full Rothermel Eq. 51: φ_s = 5.275 * β^(-0.3) * tan²(φ)
        Real slope_x = s(i, j, k, 0);
        Real slope_y = s(i, j, k, 1);
        Real tan_phi_sq = slope_x * slope_x + slope_y * slope_y;
        phi_s = std::max(Real(5.275) * std::pow(beta, Real(-0.3)) * tan_phi_sq,
                         Real(0.0));
      } else {
        phi_s = phi_s_const;
      }

      // Rate of spread using full Rothermel model [ft/min] → simulation units
      Real R_surface_ftmin = std::max(R0 * (1.0 + phi_w + phi_s), 0.0);
      Real R_surface = R_surface_ftmin * ros_conv;
      
      // Apply Van Wagner crown fire initiation if enabled
      Real R_head = R_surface;  // Default: surface fire rate
      Real local_crown_fraction = 0.0;
      
      if (use_crown_initiation == 1) {
        // Compute fireline intensity from Rothermel reaction intensity
        // Assuming w_consumed ≈ w0 * (1 - S_T) from Rothermel params
        // I_B [BTU/ft/s] = h_heat * w_consumed * R_surface_ftmin / 60
        Real I_B_imperial = I_R * R_surface_ftmin / Real(60.0);
        
        // Convert to metric if needed for Van Wagner model
        Real I_B = use_metric_units ? convert_intensity_to_metric(I_B_imperial) : I_B_imperial;
        
        // Check if crown fire initiates
        bool crown_initiated = will_crown_initiate(I_B, I_o);
        
        if (crown_initiated) {
          // Compute crown fire spread rate [m/min]
          Real R_crown_mmin = compute_crown_fire_spread_rate(CBD, FMC);
          
          // Convert to simulation units (assuming ros_conv is for ft/min to sim units)
          // Need to convert m/min to ft/min first: 1 m = 3.28084 ft
          Real R_crown_ftmin = R_crown_mmin * Real(3.28084);
          Real R_crown = R_crown_ftmin * ros_conv;
          
          // Compute combined rate (crown fire dominates)
          R_head = compute_total_spread_rate(R_surface, R_crown, crown_initiated, crown_fraction_weight);
          
          // Compute crown fire fraction for output
          local_crown_fraction = compute_crown_fraction(R_crown, R_surface, crown_initiated);
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
        LW_ratio_to_richards_coefficients(LW_ratio, local_coeff_a,
                                          local_coeff_b, local_coeff_c);
      }

      // Richards' (1990) elliptical expansion model
      Real head_spread = local_coeff_a * R_head * dt;
      Real flank_spread = local_coeff_b * R_head * dt;
      Real backing_spread = local_coeff_c * R_head * dt;

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

      Real spread_dist = 0.0;
      if (wind_mag > 1.0e-10) {
        if (cos_theta > 0.0) {
          Real sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);
          spread_dist = head_spread * cos_theta + flank_spread * sin_theta;
        } else {
          Real sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);
          spread_dist = flank_spread * sin_theta - backing_spread * cos_theta;
        }
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
        Real f_c = compute_fuel_consumption_fraction(I_R, tau_residence, f_c_min, f_c_max);
        fc(i, j, k) = f_c;
      }
      
      // Store crown fire fraction if enabled
      if (should_compute_crown_fraction) {
        cff(i, j, k) = local_crown_fraction;
      }
    });
  }
  amrex::Print()<<" Maximum Spread (cells) = "<<farsite_spread.norminf(0)<<" "<<farsite_spread.norminf(1)
#if (AMREX_SPACEDIM == 3)
             <<" "<<farsite_spread.norminf(2) 
#endif  
              <<"\n";
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

  for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
    const Box &bx = mfi.validbox();
    auto p = phi.array(mfi);

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
  }

  }

#endif // FARSITE_ELLIPSE_H
