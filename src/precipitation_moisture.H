// ============================================================================
// Precipitation-Driven Dead Fuel Moisture Wetting
// ============================================================================
// When precipitation falls during a simulation, dead fuel moisture increases
// towards saturation and then dries exponentially once rain stops.
//
// Physics (Nelson 2000 / Rothermel 1983 wetting / drying model):
//
//  Wetting phase (rain_rate_mm_hr > precip_threshold_mm_hr):
//    M_new = M_sat * (1 - exp(-dt / tau_wet))
//          + M_old * exp(-dt / tau_wet)
//
//  Drying phase (no rain):
//    M_new = M_emc + (M_old - M_emc) * exp(-dt / tau_dry)
//
// where:
//   M_sat   = saturation moisture content [fraction] (default: 1.20 = 120%)
//   M_emc   = equilibrium moisture content (from Nelson 2000 / DiurnalMoisture)
//   tau_wet = wetting time constant [s] (default: 3600 = 1 hr)
//   tau_dry = drying time constant [s] (default: 86400 = 24 hr for 1-hr fuel)
//
// Per size-class time constants follow Rothermel (1983) NFDRS conventions:
//   1-hr  fuel: tau_wet = 3600 s, tau_dry = 14400 s   (4 hr to 63% drying)
//   10-hr fuel: tau_wet = 7200 s, tau_dry = 86400 s   (24 hr)
//   100-hr fuel: tau_wet = 14400 s, tau_dry = 864000 s (10 days)
//
// The wetting model is called from the apply_fmd_moisture lambda in main.cpp
// when diurnal_moisture.enable = 1 AND diurnal_moisture.precip_mm_hr > 0.
//
// References:
//   Nelson, R.M. Jr. (2000). Prediction of diurnal change in 10-h fuel stick
//     moisture content. Canadian Journal of Forest Research, 30(7), 1071-1087.
//   Rothermel, R.C. (1983). How to Predict the Spread and Intensity of Forest
//     and Range Fires. Gen. Tech. Rep. INT-143. USDA Forest Service.
// ============================================================================

#ifndef PRECIPITATION_MOISTURE_H
#define PRECIPITATION_MOISTURE_H

#include <cmath>
#include "fuel_moisture_scheduler.H"   // for RothermelMoistures

// ============================================================================
// PrecipState
//
// Carries the per-class moisture values as they evolve under wetting / drying.
// Initialised from the global Rothermel moisture at t=0, then updated each
// timestep by apply_precipitation_moisture().
// ============================================================================
struct PrecipState {
    float M_d1   = 0.08f;   // current 1-hr   dead moisture [fraction]
    float M_d10  = 0.08f;   // current 10-hr  dead moisture [fraction]
    float M_d100 = 0.08f;   // current 100-hr dead moisture [fraction]
    float M_d1000 = 0.15f;  // current 1000-hr dead moisture [fraction]
    bool  initialized = false;
};

// ============================================================================
// compute_mcarthur_tau_scale
//
// McArthur (1967) / Noble et al. (1980) temperature-dependent time constant
// scaling for fine fuel moisture response. Higher temperatures increase drying
// rate (reduce tau); higher relative humidity slows drying (increases tau).
//
// Empirical formula:
//   f_scale = exp(-k_T * (T_C - T_ref) / T_ref) * exp(k_RH * RH_frac)
//
// where:
//   T_C      – ambient temperature [°C]
//   RH_frac  – relative humidity [0-1]
//   T_ref    – reference temperature [°C] (default: 20°C)
//   k_T      – temperature sensitivity coefficient (default: 0.05)
//   k_RH     – humidity sensitivity coefficient (default: 0.3)
//
// References:
//   McArthur, A.G. (1967). Fire behaviour in eucalypt forests.
//     Forestry & Timber Bureau, Canberra, Australia.
//   Noble, I.R., Bary, G.A.V., & Gill, A.M. (1980). McArthur's fire-danger
//     meters expressed as equations. Australian Journal of Ecology, 5(2):201-203.
// ============================================================================
inline float compute_mcarthur_tau_scale(
    float T_C,
    float RH_frac,
    float T_ref = 20.0f,
    float k_T   = 0.05f,
    float k_RH  = 0.3f)
{
    // Temperature effect: higher T → faster drying → smaller tau → smaller scale factor
    const float temp_factor = std::exp(-k_T * (T_C - T_ref) / std::max(T_ref, 1.0f));
    
    // Humidity effect: higher RH → slower drying → larger tau → larger scale factor
    const float humid_factor = std::exp(k_RH * RH_frac);
    
    const float scale = temp_factor * humid_factor;
    
    // Apply reasonable bounds to prevent extreme values
    return std::max(0.1f, std::min(scale, 10.0f));
}

// ============================================================================
// apply_precipitation_moisture
//
// Update dead fuel moisture for one timestep dt_s [s].
//
// Parameters:
//   state           – evolving moisture state (updated in-place)
//   m_emc           – equilibrium moisture from DiurnalMoisture / Nelson 2000
//                     (used as the drying target when not raining)
//   rain_rate_mm_hr – current rain rate [mm/hr]; 0 = dry
//   dt_s            – timestep size [s]
//   precip_threshold_mm_hr – minimum rain rate to trigger wetting (default 0.25 mm/hr)
//   M_sat           – saturation moisture content [fraction] (default: 1.20)
//   use_mcarthur_scaling – 1 to enable McArthur temperature/RH scaling (default: 0)
//   T_C             – ambient temperature [°C] for McArthur scaling (default: 20.0)
//   RH_frac         – relative humidity [0-1] for McArthur scaling (default: 0.5)
// ============================================================================
inline void apply_precipitation_moisture(
    PrecipState&               state,
    const RothermelMoistures&  m_emc,
    float                      rain_rate_mm_hr,
    float                      dt_s,
    float                      precip_threshold_mm_hr = 0.25f,
    float                      M_sat                  = 1.20f,
    int                        use_mcarthur_scaling   = 0,
    float                      T_C                    = 20.0f,
    float                      RH_frac                = 0.5f)
{
    // Per size-class wetting / drying time constants [s]
    constexpr float TAU_WET_1HR    =  3600.0f;
    constexpr float TAU_DRY_1HR    = 14400.0f;
    constexpr float TAU_WET_10HR   =  7200.0f;
    constexpr float TAU_DRY_10HR   = 86400.0f;
    constexpr float TAU_WET_100HR  = 14400.0f;
    constexpr float TAU_DRY_100HR  = 864000.0f;
    // 1000-hr fuel: very slow wetting (8 hr) and very slow drying (25 days)
    // Following Nelson (2000) extension of the Rothermel (1983) time-lag approach
    constexpr float TAU_WET_1000HR = 28800.0f;
    constexpr float TAU_DRY_1000HR = 2160000.0f;

    // Apply McArthur temperature/humidity scaling if enabled
    float tau_scale = 1.0f;
    if (use_mcarthur_scaling == 1) {
        tau_scale = compute_mcarthur_tau_scale(T_C, RH_frac);
    }

    const bool raining = (rain_rate_mm_hr > precip_threshold_mm_hr);

    auto update = [&](float M_old, float M_emc_val,
                      float tau_wet, float tau_dry) -> float {
        // Apply scaling to drying time constant only (wetting is less sensitive)
        const float tau_dry_scaled = tau_dry * tau_scale;
        
        if (raining) {
            // Exponential wetting towards saturation
            const float w = std::exp(-dt_s / tau_wet);
            return M_sat * (1.0f - w) + M_old * w;
        } else {
            // Exponential drying towards EMC (with optional McArthur scaling)
            const float d = std::exp(-dt_s / tau_dry_scaled);
            return M_emc_val + (M_old - M_emc_val) * d;
        }
    };

    state.M_d1   = update(state.M_d1,   m_emc.M_d1,
                          TAU_WET_1HR,   TAU_DRY_1HR);
    state.M_d10  = update(state.M_d10,  m_emc.M_d10,
                          TAU_WET_10HR,  TAU_DRY_10HR);
    state.M_d100 = update(state.M_d100, m_emc.M_d100,
                          TAU_WET_100HR, TAU_DRY_100HR);
    // 1000-hr fuel: use m_emc.M_d1000 (already scaled ≈1.5×M_d100) as the
    // equilibrium target.  When the EMC model has not set M_d1000 explicitly,
    // it is populated by compute_diurnal_emc / get_moisture_at_time using the
    // M_D1000_SCALE empirical factor from fuel_moisture_scheduler.H.
    // This is an approximation; actual 1000-hr EMC equations are not yet
    // standardised in FARSITE but the time-lag approach is consistent with
    // Rothermel (1983) large-log wetting guidelines.
    const float M_emc_1000 = m_emc.M_d1000;
    state.M_d1000 = update(state.M_d1000, M_emc_1000,
                           TAU_WET_1000HR, TAU_DRY_1000HR);
}

#endif // PRECIPITATION_MOISTURE_H
