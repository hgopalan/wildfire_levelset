// ============================================================================
// Solar Geometry, Terrain Shading, and Shade-Adjusted EMC
// ============================================================================
// Implements the FARSITE solar radiation model:
//
//   1. Solar position (elevation, azimuth) from date/time and geographic
//      coordinates.  Uses the simplified NOAA solar-position algorithm
//      (accurate to ~1° for dates within a century of J2000).
//
//   2. Per-cell terrain shading.  For each cell with slope β [deg] and
//      aspect α [deg, FARSITE convention: 0 = N, CW], the incidence angle
//      between the sun vector and the surface normal is tested:
//
//        cos_i = cos(θ_s)·sin(β)·cos(φ_s − α) + sin(θ_s)·cos(β)
//
//      where θ_s is the solar elevation and φ_s is the solar azimuth (same
//      FARSITE compass convention).
//      shade_terrain = 1 if θ_s ≤ 0 OR cos_i ≤ 0, else 0.
//
//      Optional canopy shading: shade_canopy = canopy_cover [%] / 100.
//      Total shade fraction f = max(shade_terrain, shade_canopy), in [0, 1].
//
//   3. Shade-adjusted EMC.  Solar heating raises fuel temperature above the
//      ambient air temperature.  Shade reduces this increment:
//
//        T_fuel = T_air + ΔT_solar · (1 − f)
//
//      EMC is then recomputed with T_fuel using the Simard (1968) equations
//      (same as in fuel_moisture_scheduler.H – compute_diurnal_emc).
//      Timelag corrections follow Nelson (2000):
//        M_1hr   = EMC(T_fuel, RH)
//        M_10hr  = EMC · 1.10
//        M_100hr = EMC · 1.30
//
//      The adjusted dead-fuel moistures overwrite components 0-2 of the
//      per-cell spatial_moisture_mf, preserving live-fuel moistures (3-4).
//
// Usage:
//   // Once per timestep: compute sun position
//   SolarPosition sun = compute_solar_position(lat, lon, year, month, day,
//                                              hour_local, tz_offset_h);
//   // Per-cell shade fraction
//   compute_shade_fraction_mf(shade_mf, slope_mf, aspect_mf, &cc_mf, sun, 1);
//   // Per-cell shade-adjusted EMC written to spatial_moisture_mf (comps 0-2)
//   apply_solar_emc_to_spatial_moisture(spatial_moisture_mf, shade_mf,
//                                       T_C_air, RH_pct, solar_heating_C);
//
// References:
//   Finney, M.A. (2004). FARSITE: Fire Area Simulator.
//     USDA FS Research Paper RMRS-RP-4 (Revised).
//   Simard, A.J. (1968). The moisture content of forest fuels.
//     Forestry Branch Info. Rep. FF-X-14, Dept. of Forestry, Canada.
//   Nelson, R.M. Jr. (2000). Prediction of diurnal change in 10-h fuel
//     stick moisture content. Canadian J. of Forest Research, 30(7).
//   Blanco-Chao, R. et al. / Spencer (1971): equation of time.
// ============================================================================

#ifndef SOLAR_RADIATION_H
#define SOLAR_RADIATION_H

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_GpuLaunch.H>
#include <cmath>
#include <string>

#include "constants.H"

using namespace amrex;

// ============================================================================
// Solar-date helper functions
// ============================================================================

// Returns true when the given year is a leap year (Gregorian calendar).
inline bool solar_is_leap_year(int year) noexcept
{
    return (year % 4 == 0 && year % 100 != 0) || (year % 400 == 0);
}

// Returns the number of days in each month for the given year.
// Index 0 = January, ..., 11 = December.
inline void solar_days_in_months(int year, int (&dim)[12]) noexcept
{
    const int feb = solar_is_leap_year(year) ? 29 : 28;
    const int d[] = {31, feb, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    for (int i = 0; i < 12; ++i) dim[i] = d[i];
}

// ============================================================================
// SolarPosition
// Holds the sun's elevation angle and azimuth at a given instant.
// Azimuth follows the FARSITE compass convention: 0 = North, increasing CW.
// ============================================================================
struct SolarPosition {
    amrex::Real elevation_rad;  // solar elevation above horizon [rad]; ≤ 0 = below horizon
    amrex::Real azimuth_rad;    // solar azimuth [rad], FARSITE convention: 0=N, CW
};

// ============================================================================
// compute_solar_position
//
// Returns the solar elevation and azimuth for a given date, local time, and
// geographic location using the simplified NOAA/Spencer algorithm.
//
// Parameters:
//   lat_deg       – site latitude [deg N; negative = S]
//   lon_deg       – site longitude [deg E; negative = W]
//   year          – calendar year (e.g. 2024)
//   month         – calendar month (1-12)
//   day           – calendar day (1-31)
//   hour_local    – local civil time [decimal hours from midnight, 0-24)
//   tz_offset_h   – UTC offset [hours] (e.g., −8 for PST, +1 for CET)
// ============================================================================
inline SolarPosition compute_solar_position(
    amrex::Real lat_deg,
    amrex::Real lon_deg,
    int         year,
    int         month,
    int         day,
    amrex::Real hour_local,
    amrex::Real tz_offset_h)
{
    using std::sin;  using std::cos;  using std::asin;  using std::acos;
    using std::fmod; using std::sqrt;

    constexpr amrex::Real DEG2RAD = WildfireConst::PI / amrex::Real(180.0);
    constexpr amrex::Real RAD2DEG = amrex::Real(180.0) / WildfireConst::PI;

    // ---- Day of year ----
    // Follows the formula of Spencer (1971) / ASHRAE via days in each month.
    int doy = day;
    int dim_pos[12];
    solar_days_in_months(year, dim_pos);
    for (int m = 1; m < month; ++m) doy += dim_pos[m - 1];

    // ---- Equation of time (Spencer 1971, minutes) ----
    amrex::Real B = DEG2RAD * amrex::Real(360.0) * (doy - 1) / amrex::Real(365.0);
    amrex::Real E_min = amrex::Real(229.18) *
        (amrex::Real(0.000075)
         + amrex::Real(0.001868) * cos(B)
         - amrex::Real(0.032077) * sin(B)
         - amrex::Real(0.014615) * cos(amrex::Real(2.0) * B)
         - amrex::Real(0.04089)  * sin(amrex::Real(2.0) * B));

    // ---- Solar declination (Cooper 1969) [rad] ----
    amrex::Real decl_rad = DEG2RAD * amrex::Real(23.45) *
        sin(DEG2RAD * amrex::Real(360.0) * (doy + 284) / amrex::Real(365.0));

    // ---- Local Solar Time Meridian & Time Correction [min] ----
    amrex::Real LSTM_deg = amrex::Real(15.0) * tz_offset_h;
    amrex::Real TC_min   = amrex::Real(4.0) * (lon_deg - LSTM_deg) + E_min;

    // ---- Local Solar Time [h] ----
    amrex::Real LST = hour_local + TC_min / amrex::Real(60.0);

    // ---- Hour Angle [rad]: positive in the afternoon ----
    amrex::Real HRA_rad = DEG2RAD * amrex::Real(15.0) * (LST - amrex::Real(12.0));

    // ---- Solar elevation (altitude) [rad] ----
    amrex::Real lat_rad = DEG2RAD * lat_deg;
    amrex::Real sin_elev = sin(lat_rad) * sin(decl_rad)
                         + cos(lat_rad) * cos(decl_rad) * cos(HRA_rad);
    // Clamp for asin numerical safety
    sin_elev = amrex::max(amrex::Real(-1.0), amrex::min(amrex::Real(1.0), sin_elev));
    amrex::Real elev_rad = asin(sin_elev);

    // ---- Solar azimuth [deg, geographic N=0, CW] ----
    // Compute raw azimuth from N via S (standard geographic):
    //   cos_az = (sin(decl) - sin(elev)·sin(lat)) / (cos(elev)·cos(lat))
    amrex::Real azimuth_rad = amrex::Real(0.0);
    amrex::Real cos_elev = cos(elev_rad);
    if (cos_elev > amrex::Real(1.0e-6)) {
        amrex::Real cos_az = (sin(decl_rad) - sin_elev * sin(lat_rad))
                            / (cos_elev * cos(lat_rad));
        cos_az = amrex::max(amrex::Real(-1.0), amrex::min(amrex::Real(1.0), cos_az));
        amrex::Real az_rad = acos(cos_az);
        // HRA > 0 means afternoon (sun west of south in N hemisphere)
        // FARSITE convention: azimuth = 0 = N, increasing CW
        // Geographic azimuth from N: morning=E-side (az < 180°), afternoon=W-side
        if (HRA_rad > amrex::Real(0.0)) {
            azimuth_rad = amrex::Real(2.0) * WildfireConst::PI - az_rad;
        } else {
            azimuth_rad = az_rad;
        }
    } else {
        // Sun at or below horizon: azimuth undefined; set to 0
        azimuth_rad = amrex::Real(0.0);
    }

    SolarPosition sp;
    sp.elevation_rad = elev_rad;
    sp.azimuth_rad   = azimuth_rad;
    return sp;
}

// ============================================================================
// compute_shade_fraction_mf
//
// Fills a single-component MultiFab (shade_mf) with per-cell shade fractions
// in [0, 1].  A value of 1.0 means fully shaded; 0.0 means fully insolated.
//
// Terrain shading: a cell is shaded when the sun is below the terrain horizon
// (elevation ≤ 0) or when the cosine of the incidence angle is ≤ 0:
//
//   cos_i = cos(θ_s)·sin(β)·cos(φ_s − α) + sin(θ_s)·cos(β)
//
// Canopy shading (optional): shade_canopy = canopy_cover [%] / 100.
//
// Cloud-cover shading (optional): cloud_cover [0,1] uniformly increases the
// effective shade fraction:
//   f_total = max(f_terrain_canopy, cloud_cover)
// When cloud_cover = 1 the entire domain is fully overcast (f = 1 everywhere).
//
// Topographic horizon shading (optional): the pre-computed 8-direction horizon
// MultiFab is consulted.  The horizon elevation angle for the solar azimuth is
// obtained by linear interpolation between the two bounding compass directions
// (each covers 45°).  If the solar elevation angle is below that interpolated
// horizon angle, the cell is in the shadow of a distant ridge.
//
// Parameters:
//   shade_mf     – output MultiFab (1 component), filled in-place
//   slope_mf     – terrain slope [deg] (1 component, FARSITE convention)
//   aspect_mf    – terrain aspect [deg] (1 component, FARSITE 0=N, CW)
//   cc_mf        – canopy cover [%] (1 component); nullptr = ignore canopy
//   sun          – solar position (elevation_rad, azimuth_rad)
//   use_canopy   – 1 to include canopy_cover shading, 0 to ignore
//   cloud_cover  – domain-uniform cloud fraction [0, 1] (default 0 = clear sky)
//   horizon_mf   – 8-comp horizon-angle MultiFab from
//                  compute_topographic_horizon_angles(); nullptr = skip
// ============================================================================
inline void compute_shade_fraction_mf(
    MultiFab&           shade_mf,
    const MultiFab&     slope_mf,
    const MultiFab&     aspect_mf,
    const MultiFab*     cc_mf,
    const SolarPosition& sun,
    int                 use_canopy,
    amrex::Real         cloud_cover = amrex::Real(0.0),
    const MultiFab*     horizon_mf  = nullptr)
{
    constexpr amrex::Real DEG2RAD = WildfireConst::PI / amrex::Real(180.0);

    const amrex::Real sun_elev = sun.elevation_rad;
    const amrex::Real sun_az   = sun.azimuth_rad;

    // Clamp cloud_cover to [0, 1]
    const amrex::Real cc_frac = amrex::max(amrex::Real(0.0),
                                           amrex::min(amrex::Real(1.0), cloud_cover));

    // When the sun is at or below the horizon, every cell is fully shaded.
    if (sun_elev <= amrex::Real(0.0)) {
        shade_mf.setVal(amrex::Real(1.0));
        return;
    }

    const amrex::Real cos_elev = std::cos(sun_elev);
    const amrex::Real sin_elev = std::sin(sun_elev);
    const bool has_cc          = (cc_mf != nullptr && use_canopy == 1);
    const bool has_hz          = (horizon_mf != nullptr);

    // Pre-compute direction interpolation weights for the horizon check once,
    // outside the per-cell kernel (same solar azimuth for the whole domain).
    // The 8 horizon directions are evenly spaced at 45° starting from N=0°.
    // dir_lo / dir_hi are the bounding direction indices; hz_frac is the
    // fractional weight toward dir_hi.
    int         hz_dir_lo = 0;
    int         hz_dir_hi = 1;
    amrex::Real hz_frac   = amrex::Real(0.0);
    if (has_hz) {
        constexpr amrex::Real RAD2DEG = amrex::Real(180.0) / WildfireConst::PI;
        amrex::Real az_deg = sun_az * RAD2DEG;
        // Normalise to [0, 360)
        az_deg = std::fmod(az_deg, amrex::Real(360.0));
        if (az_deg < amrex::Real(0.0)) az_deg += amrex::Real(360.0);
        amrex::Real dir_f = az_deg / amrex::Real(45.0);
        hz_dir_lo = static_cast<int>(dir_f) % 8;
        hz_dir_hi = (hz_dir_lo + 1) % 8;
        hz_frac   = dir_f - std::floor(dir_f);
    }

    for (MFIter mfi(shade_mf); mfi.isValid(); ++mfi) {
        const Box& bx    = mfi.validbox();
        auto       shade  = shade_mf.array(mfi);
        auto const slope  = slope_mf.const_array(mfi);
        auto const aspect = aspect_mf.const_array(mfi);

        // Capture canopy cover array only when needed
        Array4<const Real> cc_arr;
        if (has_cc) cc_arr = cc_mf->const_array(mfi);

        // Capture horizon array only when needed
        Array4<const Real> hz_arr;
        if (has_hz) hz_arr = horizon_mf->const_array(mfi);

        const amrex::Real cloud_f    = cc_frac;
        const int         d_lo       = hz_dir_lo;
        const int         d_hi       = hz_dir_hi;
        const amrex::Real frac       = hz_frac;
        const bool        use_hz     = has_hz;

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            amrex::Real slope_deg  = slope (i, j, k);
            amrex::Real aspect_deg = aspect(i, j, k);

            amrex::Real beta_rad  = slope_deg  * DEG2RAD;
            amrex::Real alpha_rad = aspect_deg * DEG2RAD;

            // Incidence angle cosine (Finney 2004, Eq. for terrain shading)
            amrex::Real cos_i = cos_elev * std::sin(beta_rad)
                                         * std::cos(sun_az - alpha_rad)
                              + sin_elev * std::cos(beta_rad);

            // Terrain shadow: cos_i ≤ 0 → sun below local terrain horizon
            amrex::Real f = (cos_i <= amrex::Real(0.0))
                            ? amrex::Real(1.0) : amrex::Real(0.0);

            // Canopy shading: additional fraction from canopy cover [%]
            if (has_cc) {
                amrex::Real cc_pct = cc_arr(i, j, k);
                // Clamp canopy cover to [0, 100] %
                cc_pct = amrex::max(amrex::Real(0.0),
                         amrex::min(amrex::Real(100.0), cc_pct));
                amrex::Real f_canopy = cc_pct / amrex::Real(100.0);
                f = amrex::max(f, f_canopy);
            }

            // Topographic horizon shading (FARSITE far-field ridge check).
            // Interpolate horizon angle between the two bounding compass
            // directions and shade the cell if the sun is below that angle.
            if (use_hz) {
                amrex::Real hz = hz_arr(i, j, k, d_lo) * (amrex::Real(1.0) - frac)
                               + hz_arr(i, j, k, d_hi) * frac;
                if (sun_elev < hz) f = amrex::Real(1.0);
            }

            // Cloud-cover shading: domain-uniform overcast fraction
            f = amrex::max(f, cloud_f);

            shade(i, j, k) = f;
        });
    }
}

// ============================================================================
// compute_topographic_horizon_angles
//
// One-time precomputation of per-cell terrain horizon angles in 8 compass
// directions (N, NE, E, SE, S, SW, W, NW — FARSITE clockwise convention).
//
// For each cell (i, j) and each direction d, the function ray-marches outward
// step by step and tracks the maximum elevation angle to any terrain point:
//
//   θ_horizon[d](i,j) = max_{s=1..S}  atan2( z(i+s·di, j+s·dj) − z(i,j),
//                                             s · step_dist[d] )
//
// where step_dist accounts for the actual physical cell spacing (diagonal
// directions use √(dx²+dy²)).
//
// The result is stored in an 8-component MultiFab with the layout:
//   comp 0 = N   (azimuth   0°)    comp 4 = S   (azimuth 180°)
//   comp 1 = NE  (azimuth  45°)    comp 5 = SW  (azimuth 225°)
//   comp 2 = E   (azimuth  90°)    comp 6 = W   (azimuth 270°)
//   comp 3 = SE  (azimuth 135°)    comp 7 = NW  (azimuth 315°)
//
// Implementation notes
// --------------------
//  • The computation is CPU-only.  The full elevation field is gathered to
//    all MPI ranks via ParallelDescriptor::ReduceRealSum before the scan.
//  • A pinned-memory (host-accessible) MultiFab is used to move elevation
//    data between device and host in GPU builds.
//  • THIS IS EXPENSIVE.  On an N×N grid the scan is O(8 N³) in the worst
//    case (full-domain scan).  Control cost with max_scan_dist_m > 0.
//  • Callers should guard this call behind
//      inputs.solar_radiation.use_topographic_horizon == 1
//    so that users can skip it when the expense is not justified.
//
// Parameters:
//   horizon_mf        – output 8-component MultiFab (same layout as elevation_mf)
//   elevation_mf      – ground elevation [m], 1 component
//   geom              – AMReX geometry (provides cell size and domain bounds)
//   max_scan_dist_m   – maximum ray-march distance [m]; 0 = scan full domain
// ============================================================================
inline void compute_topographic_horizon_angles(
    MultiFab&           horizon_mf,
    const MultiFab&     elevation_mf,
    const Geometry&     geom,
    amrex::Real         max_scan_dist_m = amrex::Real(0.0))
{
    // ---- Domain bookkeeping ------------------------------------------------
    const Box& domain = geom.Domain();
    const int  ni     = domain.length(0);
    const int  nj     = domain.length(1);
    const int  lo_i   = domain.smallEnd(0);
    const int  lo_j   = domain.smallEnd(1);

    const amrex::Real dx_m = static_cast<amrex::Real>(geom.CellSize(0));
    const amrex::Real dy_m = static_cast<amrex::Real>(geom.CellSize(1));
    const amrex::Real diag = std::sqrt(dx_m * dx_m + dy_m * dy_m);

    // Maximum number of steps per ray
    const int max_steps = (max_scan_dist_m > amrex::Real(0.0))
        ? static_cast<int>(std::ceil(static_cast<double>(max_scan_dist_m)
                                     / std::min(dx_m, dy_m)))
        : std::max(ni, nj);

    // ---- Step vectors for the 8 compass directions -------------------------
    // FARSITE CW from N: 0=N, 1=NE, 2=E, 3=SE, 4=S, 5=SW, 6=W, 7=NW
    const int         DI[8]        = { 0,  1,  1,  1,  0, -1, -1, -1};
    const int         DJ[8]        = { 1,  1,  0, -1, -1, -1,  0,  1};
    const amrex::Real step_dist[8] = {dy_m, diag, dx_m, diag,
                                      dy_m, diag, dx_m, diag};

    // ---- Copy elevation to pinned (host-accessible) memory -----------------
    // MFInfo::SetArena(The_Pinned_Arena()) gives CPU-accessible memory in both
    // CPU and GPU builds (regular heap on CPU, pinned host memory on GPU).
    MultiFab elev_pinned(elevation_mf.boxArray(),
                         elevation_mf.DistributionMap(), 1, 0,
                         MFInfo().SetArena(amrex::The_Pinned_Arena()));
    MultiFab::Copy(elev_pinned, elevation_mf, 0, 0, 1, 0);
    amrex::Gpu::streamSynchronize();   // ensure DtoH copy is complete

    // ---- Build global flat elevation array (all ranks) ---------------------
    std::vector<amrex::Real> elev_flat(ni * nj, amrex::Real(0.0));
    for (MFIter mfi(elev_pinned); mfi.isValid(); ++mfi) {
        const Box& bx  = mfi.validbox();
        auto const arr = elev_pinned[mfi].const_array();  // host pointer
        for (int jj = bx.smallEnd(1); jj <= bx.bigEnd(1); ++jj)
            for (int ii = bx.smallEnd(0); ii <= bx.bigEnd(0); ++ii)
                elev_flat[(jj - lo_j) * ni + (ii - lo_i)] = arr(ii, jj, 0, 0);
    }
    // Sum across MPI ranks: each rank owns non-overlapping tiles, so the sum
    // reconstructs the complete domain elevation on every rank.
    amrex::ParallelDescriptor::ReduceRealSum(elev_flat.data(), ni * nj);

    // ---- Horizon scan on each owned tile -----------------------------------
    MultiFab hz_pinned(horizon_mf.boxArray(),
                       horizon_mf.DistributionMap(), 8, 0,
                       MFInfo().SetArena(amrex::The_Pinned_Arena()));

    for (MFIter mfi(hz_pinned); mfi.isValid(); ++mfi) {
        const Box& bx  = mfi.validbox();
        auto       arr = hz_pinned[mfi].array();  // host pointer

        for (int jj = bx.smallEnd(1); jj <= bx.bigEnd(1); ++jj) {
            for (int ii = bx.smallEnd(0); ii <= bx.bigEnd(0); ++ii) {
                const amrex::Real elev_src =
                    elev_flat[(jj - lo_j) * ni + (ii - lo_i)];

                for (int d = 0; d < 8; ++d) {
                    amrex::Real max_angle = amrex::Real(-0.5) * WildfireConst::PI; // -π/2 floor
                    const amrex::Real sd  = step_dist[d];

                    for (int s = 1; s <= max_steps; ++s) {
                        const int it = ii + s * DI[d] - lo_i;
                        const int jt = jj + s * DJ[d] - lo_j;
                        if (it < 0 || it >= ni || jt < 0 || jt >= nj) break;

                        const amrex::Real elev_tgt = elev_flat[jt * ni + it];
                        const amrex::Real angle =
                            std::atan2(elev_tgt - elev_src,
                                       static_cast<amrex::Real>(s) * sd);
                        if (angle > max_angle) max_angle = angle;
                    }
                    arr(ii, jj, 0, d) = max_angle;
                }
            }
        }
    }

    // ---- Upload results to device (no-op on CPU builds) --------------------
    MultiFab::Copy(horizon_mf, hz_pinned, 0, 0, 8, 0);
    amrex::Gpu::streamSynchronize();

    amrex::Print() << "Topographic horizon angles precomputed for all 8 directions.\n";
}


// Computes per-cell shade-adjusted equilibrium moisture content (EMC) using
// the Simard (1968) equations with a solar-heating-corrected fuel temperature,
// then stores the result in components 0-2 of spatial_moisture_mf
// (M_d1, M_d10, M_d100).  Live-fuel moistures (components 3-4) are unchanged.
//
// Physics:
//   T_fuel(i,j,k) = T_air_C + solar_heating_C · (1 − shade(i,j,k))
//   EMC_pct       = Simard_1968(T_fuel, RH_pct)        [%]
//   M_d1          = EMC_pct / 100                       [fraction]
//   M_d10         = EMC_pct · 1.10 / 100               [fraction, Nelson 2000]
//   M_d100        = EMC_pct · 1.30 / 100               [fraction, Nelson 2000]
//
// Parameters:
//   spatial_moisture_mf – in/out FMS_MOISTURE_NCOMP-component MultiFab
//   shade_mf            – per-cell shade fraction [0, 1]
//   T_air_C             – ambient air temperature [°C] at current timestep
//   RH_pct              – relative humidity [%] at current timestep
//   solar_heating_C     – maximum solar heating increment [°C] (FARSITE default ≈ 17)
// ============================================================================
inline void apply_solar_emc_to_spatial_moisture(
    MultiFab&       spatial_moisture_mf,
    const MultiFab& shade_mf,
    amrex::Real     T_air_C,
    amrex::Real     RH_pct,
    amrex::Real     solar_heating_C)
{
    // Clamp ambient temperature and RH to physically meaningful ranges
    T_air_C = amrex::max(amrex::Real(-40.0), amrex::min(amrex::Real(60.0), T_air_C));
    RH_pct  = amrex::max(amrex::Real(1.0),   amrex::min(amrex::Real(100.0), RH_pct));

    const amrex::Real T_air   = T_air_C;
    const amrex::Real RH      = RH_pct;
    const amrex::Real dT_sol  = solar_heating_C;

    for (MFIter mfi(spatial_moisture_mf); mfi.isValid(); ++mfi) {
        const Box& bx  = mfi.validbox();
        auto       sm  = spatial_moisture_mf.array(mfi);
        auto const sh  = shade_mf.const_array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            // Per-cell fuel temperature (shade reduces solar heating)
            amrex::Real T_fuel = T_air + dT_sol * (amrex::Real(1.0) - sh(i, j, k));
            T_fuel = amrex::max(amrex::Real(-40.0),
                     amrex::min(amrex::Real(80.0), T_fuel));

            amrex::Real rh = RH / amrex::Real(100.0);  // [0, 1] fractional

            // Simard (1968) EMC equations [%]
            amrex::Real emc_pct;
            if (RH < amrex::Real(10.0)) {
                emc_pct = amrex::Real(0.03229)
                        + amrex::Real(0.281073) * rh
                        - amrex::Real(0.000578) * T_fuel * rh;
            } else if (RH < amrex::Real(50.0)) {
                emc_pct = amrex::Real(2.22749)
                        + amrex::Real(0.160107) * rh
                        - amrex::Real(0.014784) * T_fuel;
            } else {
                emc_pct = amrex::Real(21.0606)
                        + amrex::Real(0.005565) * rh * rh
                        - amrex::Real(0.00035)  * T_fuel * rh
                        - amrex::Real(0.483199) * rh;
            }
            emc_pct = amrex::max(amrex::Real(0.5), emc_pct);

            // Store timelag-corrected dead-fuel moistures (Nelson 2000)
            sm(i, j, k, 0) = emc_pct / amrex::Real(100.0);           // 1-hr
            sm(i, j, k, 1) = emc_pct * amrex::Real(1.10) / amrex::Real(100.0); // 10-hr
            sm(i, j, k, 2) = emc_pct * amrex::Real(1.30) / amrex::Real(100.0); // 100-hr
            // Components 3 (M_lh) and 4 (M_lw) are unchanged.
        });
    }
}

// ============================================================================
// advance_solar_date
//
// Advances the calendar date by the given number of elapsed seconds,
// updating year, month, day, and returning the fractional hour offset.
// Handles month boundaries and leap years.
//
// Parameters (all modified in-place):
//   year, month, day  – starting calendar date; updated on return
//   sim_start_hour    – local clock hour at simulation start
//   elapsed_s         – elapsed simulation time [s]
//
// Returns: current local clock time [decimal hours]
// ============================================================================
inline amrex::Real advance_solar_date(
    int&        year,
    int&        month,
    int&        day,
    amrex::Real sim_start_hour,
    amrex::Real elapsed_s)
{
    // Total elapsed hours from simulation start
    double total_h = static_cast<double>(sim_start_hour)
                   + static_cast<double>(elapsed_s) / 3600.0;

    // Use std::floor to correctly compute whole days elapsed,
    // handling both positive and negative total_h.
    double extra_days_d = (total_h >= 0.0) ? std::floor(total_h / 24.0) : 0.0;
    int    extra_days   = static_cast<int>(extra_days_d);
    // Fractional remaining hour within the day
    double hour_of_day = std::fmod(total_h, 24.0);
    if (hour_of_day < 0.0) hour_of_day += 24.0;  // fmod can return negative

    // Advance calendar by extra_days using per-month day counts
    for (int d = 0; d < extra_days; ++d) {
        int dim[12];
        solar_days_in_months(year, dim);
        ++day;
        if (day > dim[month - 1]) {
            day = 1;
            ++month;
            if (month > 12) {
                month = 1;
                ++year;
            }
        }
    }

    return static_cast<amrex::Real>(hour_of_day);
}

// ============================================================================
// apply_solar_radiation_step
//
// Computes the solar position for the given simulation time, derives per-cell
// terrain and canopy shade fractions, and (when the diurnal moisture model is
// active) writes shade-adjusted dead-fuel moistures into spatial_moisture_mf
// (components 0-2: M_d1, M_d10, M_d100).
//
// Called once at t=0 (elapsed_s = 0) for the initial state and once per
// timestep inside the main time loop.
//
// Parameters:
//   inputs            – full input parameter set (solar_radiation, diurnal_moisture)
//   shade_fraction_mf – output per-cell shade fraction [0,1] (1 component)
//   slope_mf          – terrain slope [deg], FARSITE convention
//   aspect_mf         – terrain aspect [deg], FARSITE 0=N CW
//   has_spatial_crown – true when the canopy-cover MultiFab is populated
//   cc_mf             – canopy-cover [%] MultiFab (used only when has_spatial_crown)
//   spatial_moisture_mf – in/out: dead-fuel moistures updated when diurnal_moisture active
//   elapsed_s         – elapsed simulation time [s] (0 at start)
//   horizon_mf        – 8-comp horizon-angle MultiFab from
//                       compute_topographic_horizon_angles(); nullptr = skip
//   print_position    – when true, print the computed solar elevation/azimuth
// ============================================================================
inline void apply_solar_radiation_step(
    const InputParameters& inputs,
    MultiFab&              shade_fraction_mf,
    const MultiFab&        slope_mf,
    const MultiFab&        aspect_mf,
    bool                   has_spatial_crown,
    const MultiFab&        cc_mf,
    MultiFab&              spatial_moisture_mf,
    amrex::Real            elapsed_s,
    const MultiFab*        horizon_mf    = nullptr,
    bool                   print_position = false)
{
    // Advance calendar date by elapsed simulation time
    int sol_year  = inputs.solar_radiation.year;
    int sol_month = inputs.solar_radiation.month;
    int sol_day   = inputs.solar_radiation.day;
    amrex::Real hour_local = advance_solar_date(
        sol_year, sol_month, sol_day,
        inputs.solar_radiation.sim_start_hour,
        elapsed_s);

    // Compute solar position
    SolarPosition sun = compute_solar_position(
        inputs.solar_radiation.latitude,
        inputs.solar_radiation.longitude,
        sol_year, sol_month, sol_day,
        hour_local,
        inputs.solar_radiation.timezone_offset);

    if (print_position) {
        amrex::Print() << "Solar position at t=0: elevation="
                       << sun.elevation_rad * (180.0 / WildfireConst::PI) << " deg"
                       << "  azimuth="
                       << sun.azimuth_rad   * (180.0 / WildfireConst::PI) << " deg\n";
    }

    // Per-cell terrain + canopy shade fraction
    const MultiFab* cc_for_shade = (inputs.solar_radiation.use_canopy_shading == 1
                                    && has_spatial_crown) ? &cc_mf : nullptr;
    compute_shade_fraction_mf(shade_fraction_mf,
                              slope_mf, aspect_mf,
                              cc_for_shade,
                              sun,
                              inputs.solar_radiation.use_canopy_shading,
                              static_cast<amrex::Real>(inputs.solar_radiation.cloud_cover),
                              horizon_mf);

    // Per-cell shade-adjusted EMC (requires diurnal T_air / RH)
    if (inputs.diurnal_moisture.enable == 1) {
        const double phase = WildfireConst::TWO_PI *
            (inputs.diurnal_moisture.t_start_s +
             static_cast<double>(elapsed_s) -
             inputs.diurnal_moisture.t_T_peak_s)
            / WildfireConst::DAY_S;
        const double T_mean  = 0.5 * (inputs.diurnal_moisture.T_max +
                                      inputs.diurnal_moisture.T_min);
        const double A_T     = 0.5 * (inputs.diurnal_moisture.T_max -
                                      inputs.diurnal_moisture.T_min);
        const double RH_mean = 0.5 * (inputs.diurnal_moisture.RH_max +
                                      inputs.diurnal_moisture.RH_min);
        const double A_RH    = 0.5 * (inputs.diurnal_moisture.RH_max -
                                      inputs.diurnal_moisture.RH_min);
        double T_C  = T_mean  + A_T  * std::sin(phase);
        double RH_p = RH_mean - A_RH * std::sin(phase);
        T_C  = std::max(-40.0, std::min(60.0,  T_C));
        RH_p = std::max(1.0,   std::min(100.0, RH_p));

        apply_solar_emc_to_spatial_moisture(
            spatial_moisture_mf,
            shade_fraction_mf,
            static_cast<amrex::Real>(T_C),
            static_cast<amrex::Real>(RH_p),
            inputs.solar_radiation.solar_heating_C);
    }
}

#endif // SOLAR_RADIATION_H
