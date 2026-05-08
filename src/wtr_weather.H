// ============================================================================
// FARSITE Weather File (.wtr) Parser
// ============================================================================
// Reads FARSITE's hourly weather file format and populates a WindDirSchedule,
// diurnal moisture arrays, and precipitation state from a single input file,
// replacing the need for separate wind_dir_schedule_file, fmd_file, and manual
// diurnal_moisture parameters.
//
// FARSITE .wtr file format (whitespace or comma delimited; # or ! lines are
// comments):
//
//   MONTH  DAY  HOUR  TEMP_F  RH_PCT  PRECIP_IN  WIND_SPEED_MPH  WIND_DIR_DEG
//   7      1    1000  85      20      0.00        10              270
//   7      1    1100  90      18      0.00        12              265
//   ...
//
// Fields:
//   MONTH         – calendar month (1-12)
//   DAY           – calendar day (1-31)
//   HOUR          – local clock time in HHMM (0 = midnight, 1400 = 14:00)
//   TEMP_F        – air temperature [°F]
//   RH_PCT        – relative humidity [%]
//   PRECIP_IN     – hourly precipitation [inches]
//   WIND_SPEED_MPH – 20-ft open wind speed [mph]
//   WIND_DIR_DEG  – meteorological wind direction [degrees from North]
//
// Output:
//   The parser populates:
//   • WindDirSchedule – (time_s, speed_m/s, dir_deg) ready for get_wind_at_time()
//   • WtrRecord vector – (time_s, T_C, RH_pct, precip_mm_hr) for moisture models
//
// Usage:
//   WtrWeatherData wtr;
//   load_wtr_weather("fire.wtr", wtr);
//   // Wind: delegate to get_wind_at_time(wtr.wind_sched, sim_time_s);
//   // Moisture: query wtr.get_TRH_at_time(sim_time_s) for (T_C, RH_pct)
//   //           and wtr.get_precip_at_time(sim_time_s) for rain_mm_hr.
//
// References:
//   Finney, M.A. (2004). FARSITE: Fire Area Simulator. USDA FS RMRS-RP-4.
//   FARSITE Weather and Wind file formats (Sec. 3.2 of the FARSITE manual).
// ============================================================================

#ifndef WTR_WEATHER_H
#define WTR_WEATHER_H

#include <AMReX.H>
#include "wind_dir_schedule.H"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <utility>

// One parsed row from the .wtr file
struct WtrRecord {
    double time_s;        // seconds from simulation start
    double T_C;           // air temperature [°C]
    double RH_pct;        // relative humidity [%]
    double precip_mm_hr;  // precipitation rate [mm/hr]
};

// Container for all .wtr data
struct WtrWeatherData {
    WindDirSchedule        wind_sched;  // wind direction/speed schedule
    std::vector<WtrRecord> weather;     // T, RH, precip time series
    bool empty() const { return weather.empty(); }

    // Linearly interpolated (T_C, RH_pct) at simulation time t_s
    std::pair<double,double> get_TRH_at_time(double t_s) const {
        if (weather.empty()) return {25.0, 30.0};
        if (t_s <= weather.front().time_s)
            return {weather.front().T_C, weather.front().RH_pct};
        if (t_s >= weather.back().time_s)
            return {weather.back().T_C, weather.back().RH_pct};
        auto it = std::upper_bound(weather.begin(), weather.end(), t_s,
            [](double t, const WtrRecord& r){ return t < r.time_s; });
        auto prev = std::prev(it);
        const double alpha = (t_s - prev->time_s) / (it->time_s - prev->time_s);
        return {prev->T_C   + alpha * (it->T_C   - prev->T_C),
                prev->RH_pct + alpha * (it->RH_pct - prev->RH_pct)};
    }

    // Linearly interpolated precipitation rate [mm/hr] at simulation time t_s
    double get_precip_at_time(double t_s) const {
        if (weather.empty()) return 0.0;
        if (t_s <= weather.front().time_s) return weather.front().precip_mm_hr;
        if (t_s >= weather.back().time_s)  return weather.back().precip_mm_hr;
        auto it = std::upper_bound(weather.begin(), weather.end(), t_s,
            [](double t, const WtrRecord& r){ return t < r.time_s; });
        auto prev = std::prev(it);
        const double alpha = (t_s - prev->time_s) / (it->time_s - prev->time_s);
        return prev->precip_mm_hr + alpha * (it->precip_mm_hr - prev->precip_mm_hr);
    }
};

// ============================================================================
// load_wtr_weather
//
// Parse a FARSITE .wtr file into a WtrWeatherData structure.
//
// The start date/hour are used to convert .wtr calendar timestamps to
// simulation time [s].  When start_year/month/day/hour are all 0, the
// first record defines t = 0.
// ============================================================================
inline void load_wtr_weather(const std::string& filename,
                              WtrWeatherData&    wtr,
                              int start_year  = 0,
                              int start_month = 0,
                              int start_day   = 0,
                              int start_hour  = 0)
{
    wtr.wind_sched.records.clear();
    wtr.weather.clear();

    std::ifstream f(filename);
    if (!f.is_open())
        amrex::Abort("WtrWeather: cannot open file: " + filename);

    // Helper: day number relative to a fixed epoch (modified from Gregorian day count).
    // Formula produces a monotone integer day counter; value 719591 aligns
    // the reference to 1 Jan 1900 (common epoch for FARSITE date arithmetic).
    // Algorithm: Adapted from the Gregorian proleptic calendar formula by
    //   Richards (2013), in Explanatory Supplement to the Astronomical Almanac.
    auto days_since_epoch = [](int y, int m, int d) -> long {
        if (m < 3) { m += 12; --y; }
        return 365L * y + y/4 - y/100 + y/400
               + (306L*(m+1))/10 + d - 719591L;
    };

    long   ref_s = 0;
    bool   ref_set = (start_year != 0);

    if (ref_set) {
        ref_s = days_since_epoch(start_year, start_month, start_day) * 86400L
              + static_cast<long>(start_hour / 100) * 3600L
              + static_cast<long>(start_hour % 100) * 60L;
    }

    // Unit conversion constants
    constexpr double F_TO_C  = 5.0 / 9.0;   // °F offset handled below
    constexpr double MPH_TO_MS = 0.44704;    // 1 mph = 0.44704 m/s
    constexpr double IN_PER_HR_TO_MM_PER_HR = 25.4;

    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#' || line[0] == '!') continue;
        for (char& c : line) if (c == ',') c = ' ';
        std::istringstream iss(line);
        int month, day, hour_hhmm;
        double temp_F, rh_pct, precip_in, wspd_mph, wdir_deg;
        if (!(iss >> month >> day >> hour_hhmm
                  >> temp_F >> rh_pct >> precip_in >> wspd_mph >> wdir_deg)) continue;
        if (month < 1 || month > 12 || day < 1 || day > 31) continue;

        // Convert HHMM → seconds into day
        int hh = hour_hhmm / 100;
        int mm = hour_hhmm % 100;
        long block_day_s = days_since_epoch(2000, month, day) * 86400L
                         + static_cast<long>(hh) * 3600L
                         + static_cast<long>(mm) * 60L;

        if (!ref_set) {
            // First record defines simulation t = 0.  Year 2000 is used as the
            // nominal calendar year since FARSITE .wtr files often omit the year;
            // only the relative time differences matter for schedule interpolation.
            ref_s   = block_day_s;
            ref_set = true;
        }
        double t_s = static_cast<double>(block_day_s - ref_s);

        // Temperature: °F → °C
        double T_C = (temp_F - 32.0) * F_TO_C;
        // Clamp to physical range
        T_C   = std::max(-40.0, std::min(60.0,  T_C));
        rh_pct = std::max(1.0,  std::min(100.0, rh_pct));
        // Precipitation: in/hr → mm/hr
        double precip_mm_hr = std::max(0.0, precip_in) * IN_PER_HR_TO_MM_PER_HR;
        // Wind: mph → m/s
        double spd_ms = std::max(0.0, wspd_mph) * MPH_TO_MS;
        // Wind direction: meteorological degrees (0 = N)
        wdir_deg = std::fmod(std::max(0.0, wdir_deg), 360.0);

        // Store weather record
        WtrRecord wrec;
        wrec.time_s        = t_s;
        wrec.T_C           = T_C;
        wrec.RH_pct        = rh_pct;
        wrec.precip_mm_hr  = precip_mm_hr;
        wtr.weather.push_back(wrec);

        // Store wind record
        WindDirRecord wdrec;
        wdrec.time_s   = t_s;
        wdrec.speed_ms = spd_ms;
        wdrec.dir_deg  = wdir_deg;
        wtr.wind_sched.records.push_back(wdrec);
    }

    if (wtr.weather.empty())
        amrex::Abort("WtrWeather: no valid records found in: " + filename);

    // Sort by time
    std::sort(wtr.weather.begin(), wtr.weather.end(),
              [](const WtrRecord& a, const WtrRecord& b){ return a.time_s < b.time_s; });
    std::sort(wtr.wind_sched.records.begin(), wtr.wind_sched.records.end(),
              [](const WindDirRecord& a, const WindDirRecord& b){ return a.time_s < b.time_s; });

    amrex::Print() << "WtrWeather: loaded " << wtr.weather.size()
                   << " hourly weather records from " << filename << "\n";
}

#endif // WTR_WEATHER_H
