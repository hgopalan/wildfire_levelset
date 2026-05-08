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
// Parameters:
//   shade_mf     – output MultiFab (1 component), filled in-place
//   slope_mf     – terrain slope [deg] (1 component, FARSITE convention)
//   aspect_mf    – terrain aspect [deg] (1 component, FARSITE 0=N, CW)
//   cc_mf        – canopy cover [%] (1 component); nullptr = ignore canopy
//   sun          – solar position (elevation_rad, azimuth_rad)
//   use_canopy   – 1 to include canopy_cover shading, 0 to ignore
//   cloud_cover  – domain-uniform cloud fraction [0, 1] (default 0 = clear sky)
// ============================================================================
inline void compute_shade_fraction_mf(
    MultiFab&           shade_mf,
    const MultiFab&     slope_mf,
    const MultiFab&     aspect_mf,
    const MultiFab*     cc_mf,
    const SolarPosition& sun,
    int                 use_canopy,
    amrex::Real         cloud_cover = amrex::Real(0.0))
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

    for (MFIter mfi(shade_mf); mfi.isValid(); ++mfi) {
        const Box& bx    = mfi.validbox();
        auto       shade  = shade_mf.array(mfi);
        auto const slope  = slope_mf.const_array(mfi);
        auto const aspect = aspect_mf.const_array(mfi);

        // Capture canopy cover array only when needed
        Array4<const Real> cc_arr;
        if (has_cc) cc_arr = cc_mf->const_array(mfi);

        const amrex::Real cloud_f = cc_frac;

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

            // Cloud-cover shading: domain-uniform overcast fraction
            f = amrex::max(f, cloud_f);

            shade(i, j, k) = f;
        });
    }
}

// ============================================================================
// apply_solar_emc_to_spatial_moisture
//
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

#endif // SOLAR_RADIATION_H
