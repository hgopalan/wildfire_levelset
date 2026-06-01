// ============================================================================
// Diurnal Weather Cycles
// ============================================================================
// Implements diurnal (daily) variation of temperature, relative humidity,
// and wind speed for realistic multi-day fire simulations.
//
// Fire weather exhibits strong diurnal patterns driven by solar heating:
//   - Temperature: minimum at sunrise, maximum in mid-afternoon
//   - Relative humidity: maximum at sunrise, minimum in mid-afternoon
//   - Wind speed: minimum at night, maximum in afternoon
//
// These cycles significantly affect fire behavior, causing:
//   - Nighttime fire spread reduction (cooler temps, higher RH, lighter winds)
//   - Afternoon fire intensification (hot, dry, windy conditions)
//   - Morning transition periods with changing fire activity
//
// Model Formulations:
//
// **Temperature (°C):**
//   T(t) = T_mean + T_amp × sin(2π × (t - t_min) / 24)
//
//   where T_mean = (T_max + T_min) / 2, T_amp = (T_max - T_min) / 2
//
// **Relative Humidity (%):**
//   RH(t) = RH_mean - RH_amp × sin(2π × (t - t_min) / 24)
//
//   where RH is anti-phase with temperature
//
// **Wind Speed (m/s):**
//   U(t) = U_mean × [1 + α × sin(2π × (t - 6) / 24)]
//
//   where α is diurnal amplitude factor (default: 0.4 for 40% variation)
//
// Time Convention:
//   - t = local time in hours (0-24)
//   - t_min = time of minimum temperature (default: 6:00 = sunrise)
//   - t_max = time of maximum temperature (default: 15:00 = 3:00 PM)
//
// References:
//   Lawson, B.D., Frandsen, W.H., Hawkes, B.C., & Dalrymple, G.N. (1996).
//     "Probability of sustained ignition as a function of nature of the
//     ignition source." International Journal of Wildland Fire, 6:29-36.
//   Tymstra, C., et al. (2010). "Development and structure of Prometheus:
//     the Canadian wildland fire growth simulation model."
//     Natural Resources Canada, Information Report NOR-X-417.
//   Simard, A.J. (1968). "The moisture content of forest fuels."
//     Canadian Forestry Service, Publication 1236.
// ============================================================================

#ifndef DIURNAL_WEATHER_H
#define DIURNAL_WEATHER_H

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_GpuLaunch.H>
#include <cmath>

using namespace amrex;

// ============================================================================
// Diurnal Weather Constants
// ============================================================================
namespace DiurnalConst {
    constexpr Real PI = Real(3.14159265358979323846);
    constexpr Real HOURS_PER_DAY = Real(24.0);
    
    // Default times for temperature extrema (local solar time)
    constexpr Real T_MIN_HOUR = Real(6.0);   // Sunrise (6:00 AM)
    constexpr Real T_MAX_HOUR = Real(15.0);  // Afternoon (3:00 PM)
    
    // Wind speed diurnal amplitude factor
    constexpr Real WIND_AMPLITUDE = Real(0.4);  // ±40% variation
    constexpr Real WIND_PEAK_HOUR = Real(15.0); // Peak wind at 3:00 PM
}

// ============================================================================
// compute_diurnal_temperature
//
// Computes temperature at a given time using sinusoidal diurnal variation.
//
// Formulation:
//   T(t) = T_mean + T_amp × sin(2π × (t - t_min) / 24)
//
// Parameters:
//   hour_of_day - local time [hours] (0-24)
//   T_min       - daily minimum temperature [°C]
//   T_max       - daily maximum temperature [°C]
//   t_min_hour  - time of minimum temperature [hours] (default: 6:00 AM)
//
// Returns: temperature [°C]
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real compute_diurnal_temperature(Real hour_of_day,
                                  Real T_min,
                                  Real T_max,
                                  Real t_min_hour = DiurnalConst::T_MIN_HOUR)
{
    // Mean and amplitude
    Real T_mean = (T_max + T_min) / Real(2.0);
    Real T_amp = (T_max - T_min) / Real(2.0);
    
    // Wrap hour to [0, 24)
    hour_of_day = std::fmod(hour_of_day, DiurnalConst::HOURS_PER_DAY);
    if (hour_of_day < Real(0.0)) hour_of_day += DiurnalConst::HOURS_PER_DAY;
    
    // Phase angle (radians)
    Real phase = Real(2.0) * DiurnalConst::PI * (hour_of_day - t_min_hour) / DiurnalConst::HOURS_PER_DAY;
    
    // Sinusoidal temperature
    Real T = T_mean + T_amp * std::sin(phase);
    
    return T;
}

// ============================================================================
// compute_diurnal_rh
//
// Computes relative humidity using inverse relationship with temperature.
//
// Formulation:
//   RH(t) = RH_mean - RH_amp × sin(2π × (t - t_min) / 24)
//
// RH is anti-phase with temperature (maximum at sunrise, minimum at afternoon).
//
// Parameters:
//   hour_of_day - local time [hours] (0-24)
//   RH_min      - daily minimum RH [%]
//   RH_max      - daily maximum RH [%]
//   t_min_hour  - time of minimum RH (= time of max T) [hours] (default: 15:00)
//
// Returns: relative humidity [%]
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real compute_diurnal_rh(Real hour_of_day,
                        Real RH_min,
                        Real RH_max,
                        Real t_min_rh_hour = DiurnalConst::T_MAX_HOUR)
{
    // Mean and amplitude
    Real RH_mean = (RH_max + RH_min) / Real(2.0);
    Real RH_amp = (RH_max - RH_min) / Real(2.0);
    
    // Wrap hour to [0, 24)
    hour_of_day = std::fmod(hour_of_day, DiurnalConst::HOURS_PER_DAY);
    if (hour_of_day < Real(0.0)) hour_of_day += DiurnalConst::HOURS_PER_DAY;
    
    // Phase angle (radians) - RH is anti-phase with T
    Real phase = Real(2.0) * DiurnalConst::PI * (hour_of_day - t_min_rh_hour) / DiurnalConst::HOURS_PER_DAY;
    
    // Sinusoidal RH (minus sign for anti-phase)
    Real RH = RH_mean - RH_amp * std::sin(phase);
    
    // Clamp to physical range [0, 100]
    return amrex::max(Real(0.0), amrex::min(RH, Real(100.0)));
}

// ============================================================================
// compute_diurnal_wind_speed
//
// Computes wind speed diurnal scaling factor.
//
// Formulation:
//   U(t) = U_mean × [1 + α × sin(2π × (t - 6) / 24)]
//
// where α is the amplitude factor (default: 0.4 for ±40% variation).
//
// Wind speed is typically minimum at night/early morning (6:00 AM) and
// maximum in mid-afternoon (3:00 PM).
//
// Parameters:
//   hour_of_day  - local time [hours] (0-24)
//   U_mean       - daily mean wind speed [m/s]
//   amplitude    - diurnal amplitude factor (default: 0.4)
//   peak_hour    - time of maximum wind [hours] (default: 15:00)
//
// Returns: wind speed [m/s]
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real compute_diurnal_wind_speed(Real hour_of_day,
                                 Real U_mean,
                                 Real amplitude = DiurnalConst::WIND_AMPLITUDE,
                                 Real peak_hour = DiurnalConst::WIND_PEAK_HOUR)
{
    // Wrap hour to [0, 24)
    hour_of_day = std::fmod(hour_of_day, DiurnalConst::HOURS_PER_DAY);
    if (hour_of_day < Real(0.0)) hour_of_day += DiurnalConst::HOURS_PER_DAY;
    
    // Phase angle (radians) - minimum at sunrise, maximum at peak_hour
    Real t_min = peak_hour - Real(9.0);  // Minimum 9 hours before peak
    if (t_min < Real(0.0)) t_min += DiurnalConst::HOURS_PER_DAY;
    
    Real phase = Real(2.0) * DiurnalConst::PI * (hour_of_day - t_min) / DiurnalConst::HOURS_PER_DAY;
    
    // Sinusoidal wind speed
    Real U = U_mean * (Real(1.0) + amplitude * std::sin(phase));
    
    // Ensure non-negative wind speed
    return amrex::max(Real(0.0), U);
}

// ============================================================================
// compute_hour_of_day
//
// Extracts local hour of day from simulation time.
//
// Parameters:
//   sim_time    - simulation time [seconds]
//   time_offset - local time offset from UTC [hours] (default: 0)
//
// Returns: hour of day [hours] (0-24)
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real compute_hour_of_day(Real sim_time, Real time_offset = Real(0.0))
{
    // Convert simulation time to hours
    Real hours = sim_time / Real(3600.0);
    
    // Add local time offset
    hours += time_offset;
    
    // Wrap to [0, 24)
    Real hour_of_day = std::fmod(hours, DiurnalConst::HOURS_PER_DAY);
    if (hour_of_day < Real(0.0)) hour_of_day += DiurnalConst::HOURS_PER_DAY;
    
    return hour_of_day;
}

// ============================================================================
// apply_diurnal_weather
//
// Applies diurnal cycles to temperature, RH, and wind speed fields.
//
// Parameters:
//   temperature_mf - temperature field [°C] (modified in-place)
//   rh_mf          - relative humidity field [%] (modified in-place)
//   wind_speed_mf  - wind speed field [m/s] (modified in-place)
//   sim_time       - current simulation time [seconds]
//   T_min          - daily minimum temperature [°C]
//   T_max          - daily maximum temperature [°C]
//   RH_min         - daily minimum RH [%]
//   RH_max         - daily maximum RH [%]
//   U_mean         - daily mean wind speed [m/s]
//   time_offset    - local time offset from UTC [hours] (default: 0)
// ============================================================================
inline void apply_diurnal_weather(MultiFab& temperature_mf,
                                   MultiFab& rh_mf,
                                   MultiFab& wind_speed_mf,
                                   Real sim_time,
                                   Real T_min,
                                   Real T_max,
                                   Real RH_min,
                                   Real RH_max,
                                   Real U_mean,
                                   Real time_offset = Real(0.0))
{
    // Compute current hour of day
    Real hour = compute_hour_of_day(sim_time, time_offset);
    
    // Compute diurnal values
    Real T_current = compute_diurnal_temperature(hour, T_min, T_max);
    Real RH_current = compute_diurnal_rh(hour, RH_min, RH_max);
    Real U_current = compute_diurnal_wind_speed(hour, U_mean);
    
    // Update fields (uniform across domain for now)
    temperature_mf.setVal(T_current);
    rh_mf.setVal(RH_current);
    
    // Scale wind speed field
    for (MFIter mfi(wind_speed_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        Array4<Real> const& wind = wind_speed_mf.array(mfi);
        Real scale = U_current / U_mean;
        
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            wind(i,j,k) *= scale;
        });
    }
}

#endif // DIURNAL_WEATHER_H
