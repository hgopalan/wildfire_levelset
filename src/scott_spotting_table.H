// ============================================================================
// Scott / Albini Maximum Spotting Distance Table
// ============================================================================
// Provides maximum typical spotting distances by fuel model number, based on:
//
//   Albini, F.A. (1979). "Spot fire distance from burning trees – a
//     predictive model." USDA Forest Service General Technical Report INT-56.
//   Scott, J.H. (1998). "Fuel reduction in residential and scenic/
//     recreational zones." USDA Forest Service Research Paper RMRS-RP-1.
//   Rothermel, R.C. (1991). "Predicting behavior and size of crown fires in
//     the northern Rocky Mountains." USDA Forest Service Research Paper
//     INT-438.
//   Finney, M.A. (2004). FARSITE: Fire Area Simulator. RMRS-RP-4, Table 5.1.
//
// The values tabulated here are moderate-wind (≈ 20 mph / 9 m/s) maximum
// spotting distances from point-source fires.  They represent the 95th
// percentile firebrand travel distance for each surface fuel model under
// typical fire-weather conditions.  Exact values vary with wind speed; these
// are conservative upper bounds useful for capping trajectory-integrated
// spotting distances.
//
// Two lookup functions are provided:
//   get_max_spot_dist_fbfm13(model)  – Anderson 13 models (1–13)
//   get_max_spot_dist_fbfm40(model)  – Scott & Burgan 40 models (101–265+)
//
// Usage (diagnostic / cap on Albini trajectory):
//   float d_max = get_max_spot_dist_fbfm13(fuel_model_code);    // metres
//   if (spotting_distance_m > d_max) { /* cap or flag */ }
// ============================================================================

#ifndef SCOTT_SPOTTING_TABLE_H
#define SCOTT_SPOTTING_TABLE_H

#include <AMReX.H>
#include <cmath>
#include <string>

// ============================================================================
// Maximum spotting distances [metres] for Anderson (1982) FBFM13 fuel models
//
// Values derived from Albini (1979) Table 6 and Rothermel (1991) at 20 mph:
//   1  Short grass (1 ft)                        250 m
//   2  Timber (grass/understory)                 500 m
//   3  Tall grass (2.5 ft)                       900 m
//   4  Chaparral (6 ft)                         1800 m
//   5  Brush (2 ft)                              500 m
//   6  Dormant brush / hardwood slash            600 m
//   7  Southern rough                            600 m
//   8  Closed timber litter                      300 m
//   9  Hardwood litter                           400 m
//  10  Timber (litter + understory)             1200 m
//  11  Light logging slash                       500 m
//  12  Medium logging slash                     1400 m
//  13  Heavy logging slash                      1800 m
//  91/92/93/98/99 (non-burnable)                  0 m
// ============================================================================
static AMREX_FORCE_INLINE
float get_max_spot_dist_fbfm13(int model)
{
    switch (model) {
        case  1: return  250.0f;
        case  2: return  500.0f;
        case  3: return  900.0f;
        case  4: return 1800.0f;
        case  5: return  500.0f;
        case  6: return  600.0f;
        case  7: return  600.0f;
        case  8: return  300.0f;
        case  9: return  400.0f;
        case 10: return 1200.0f;
        case 11: return  500.0f;
        case 12: return 1400.0f;
        case 13: return 1800.0f;
        default: return    0.0f;  // non-burnable or unknown
    }
}

// ============================================================================
// Maximum spotting distances [metres] for Scott & Burgan (2005) FBFM40 models
//
// Values estimated from Albini (1979) categories matched to fuel type
// descriptions.  FBFM40 codes: NB = non-burnable (101-109), GR = grass
// (101-109), GS = grass-shrub (141-161), SH = shrub (181-215),
// TU = timber-understory (221-242), TL = timber litter (261-265),
// SB = slash-blowdown (281-292).
// ============================================================================
static AMREX_FORCE_INLINE
float get_max_spot_dist_fbfm40(int model)
{
    // Non-burnable (NB1–NB9: 91–99 in FBFM13; in FBFM40: 91-99 same codes)
    if (model >= 91 && model <= 99) return 0.0f;

    // GR: Short/sparse grass (low potential)
    if (model == 101 || model == 102 || model == 103) return 200.0f;
    // GR: Moderate-density grass
    if (model == 104 || model == 105) return 500.0f;
    // GR: Tall/dense grass
    if (model == 106 || model == 107 || model == 108 || model == 109) return 900.0f;

    // GS: Grass-shrub mix
    if (model == 121 || model == 122) return 600.0f;
    if (model == 123 || model == 124) return 900.0f;

    // SH: Shrub
    if (model == 141 || model == 142) return  600.0f;
    if (model == 143 || model == 144) return  900.0f;
    if (model == 145) return  700.0f;
    if (model == 146) return 1000.0f;
    if (model == 147) return 1200.0f;
    if (model == 148) return  800.0f;
    if (model == 149) return 1800.0f;

    // TU: Timber-understory
    if (model == 161) return  500.0f;
    if (model == 162) return  700.0f;
    if (model == 163) return 1000.0f;
    if (model == 164) return  800.0f;
    if (model == 165) return  600.0f;

    // TL: Timber litter
    if (model == 181) return  300.0f;
    if (model == 182) return  400.0f;
    if (model == 183) return  500.0f;
    if (model == 184) return  600.0f;
    if (model == 185) return  400.0f;
    if (model == 186) return  500.0f;
    if (model == 187) return  700.0f;
    if (model == 188) return  800.0f;
    if (model == 189) return 1200.0f;

    // SB: Slash-blowdown
    if (model == 201) return  600.0f;
    if (model == 202) return  900.0f;
    if (model == 203) return 1400.0f;
    if (model == 204) return 1800.0f;

    // Fallback: use a moderate default
    return 500.0f;
}

// ============================================================================
// get_max_spot_dist_m
//
// Unified entry point: returns maximum spotting distance [m] for any
// supported fuel model code.  The fuel_model_type string should be "13" for
// FBFM13 or "40" for FBFM40.
// ============================================================================
static AMREX_FORCE_INLINE
float get_max_spot_dist_m(int fuel_model_code,
                          const std::string& fuel_model_type = "13")
{
    if (fuel_model_type == "40")
        return get_max_spot_dist_fbfm40(fuel_model_code);
    return get_max_spot_dist_fbfm13(fuel_model_code);
}

// ============================================================================
// get_max_spot_dist_m  (overload: wind-speed scaling)
//
// Applies a simplified Albini (1979) wind-speed scaling to the base distance:
//   d_max(U) = d_max_base × (U / U_ref)^0.5
// where U_ref = 9.0 m/s (≈ 20 mph) is the reference wind used for the table.
// Wind speed U should be in m/s.
// ============================================================================
static AMREX_FORCE_INLINE
float get_max_spot_dist_scaled_m(int fuel_model_code,
                                  float wind_speed_ms,
                                  const std::string& fuel_model_type = "13")
{
    constexpr float U_ref = 9.0f;
    float d_base = get_max_spot_dist_m(fuel_model_code, fuel_model_type);
    if (d_base <= 0.0f || wind_speed_ms <= 0.0f) return d_base;
    return d_base * std::sqrt(wind_speed_ms / U_ref);
}

#endif // SCOTT_SPOTTING_TABLE_H
