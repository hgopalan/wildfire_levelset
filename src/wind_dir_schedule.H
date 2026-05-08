// ============================================================================
// Compact Wind Direction / Speed Schedule
// ============================================================================
// Reads a compact CSV schedule file that specifies time-varying 20-ft wind
// speed and direction without requiring a full spatial X/Y wind-field CSV per
// snapshot.
//
// File format (whitespace or comma separated; # lines are comments):
//   time_s  speed_ms  dir_deg
//   0       3.0       270
//   3600    4.5       255
//   7200    5.0       240
//   ...
//
// Fields:
//   time_s    – simulation time [s] at which this wind becomes effective
//   speed_ms  – 20-ft wind speed [m/s]
//   dir_deg   – meteorological wind direction [degrees from North, where wind
//               comes FROM; e.g. 270 = westerly wind blowing East]
//
// The schedule is linearly interpolated between rows.  Outside the schedule
// range the nearest endpoint value is used.
//
// Usage:
//   WindDirSchedule wds;
//   load_wind_dir_schedule("wind_schedule.csv", wds);
//   // at each timestep:
//   auto [ux, uy] = get_wind_at_time(wds, sim_time_s);
//   // then: init_velocity_constant(vel, geom, ux, uy, 0.0);
//
// Relation to time-dependent wind files:
//   This schedule avoids having to write N full spatial wind CSVs when the
//   wind is spatially uniform.  It is the compact equivalent of the multi-file
//   use_time_dependent_wind = 1 mode for uniform-wind domains.
//
// Reference:
//   Finney, M.A. (2004). FARSITE: Fire Area Simulator. USDA FS RMRS-RP-4.
//   FARSITE accepts per-weather-station time-stamped (speed, direction) without
//   spatial grids; this schedule mirrors that input style for the solver.
// ============================================================================

#ifndef WIND_DIR_SCHEDULE_H
#define WIND_DIR_SCHEDULE_H

#include <AMReX.H>
#include <AMReX_Math.H>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

// Portable pi constant via AMReX (avoids relying on non-standard M_PI)
static constexpr double WDSC_PI = amrex::Math::pi<double>();

// One entry in the wind direction schedule
struct WindDirRecord {
    double time_s;     // simulation time [s]
    double speed_ms;   // wind speed [m/s]
    double dir_deg;    // meteorological direction [degrees from North]
};

// The full schedule (sorted by time_s)
struct WindDirSchedule {
    std::vector<WindDirRecord> records;
    bool empty() const { return records.empty(); }
};

// ============================================================================
// load_wind_dir_schedule
//
// Parse a CSV schedule file into a WindDirSchedule.
// ============================================================================
inline void load_wind_dir_schedule(const std::string& filename, WindDirSchedule& sched)
{
    sched.records.clear();
    std::ifstream f(filename);
    if (!f.is_open()) {
        amrex::Abort("WindDirSchedule: cannot open file: " + filename);
    }
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#' || line[0] == '!') continue;
        for (char& c : line) if (c == ',') c = ' ';
        std::istringstream iss(line);
        WindDirRecord r;
        if (iss >> r.time_s >> r.speed_ms >> r.dir_deg) {
            sched.records.push_back(r);
        }
    }
    if (sched.records.empty()) {
        amrex::Abort("WindDirSchedule: no valid records in: " + filename);
    }
    std::sort(sched.records.begin(), sched.records.end(),
              [](const WindDirRecord& a, const WindDirRecord& b){
                  return a.time_s < b.time_s; });
    amrex::Print() << "WindDirSchedule: loaded " << sched.records.size()
                   << " records from " << filename << "\n";
}

// ============================================================================
// get_wind_at_time
//
// Linearly interpolate wind (speed, direction) at simulation time t_s.
// Returns (Ux, Uy) in m/s using East=+x, North=+y convention.
//
// Note: wind direction is interpolated as a circular quantity to avoid wrap-
// around artefacts (e.g. 355° and 5° correctly interpolate through 0°).
// ============================================================================
inline std::pair<double, double> get_wind_at_time(const WindDirSchedule& sched,
                                                   double t_s)
{
    if (sched.empty()) return {0.0, 0.0};

    const auto& recs = sched.records;

    // Clamp to schedule endpoints
    if (t_s <= recs.front().time_s) {
        const double spd = recs.front().speed_ms;
        const double d   = recs.front().dir_deg;
        const double ang = (270.0 - d) * WDSC_PI / 180.0;
        return {spd * std::cos(ang), spd * std::sin(ang)};
    }
    if (t_s >= recs.back().time_s) {
        const double spd = recs.back().speed_ms;
        const double d   = recs.back().dir_deg;
        const double ang = (270.0 - d) * WDSC_PI / 180.0;
        return {spd * std::cos(ang), spd * std::sin(ang)};
    }

    // Binary search for bracketing interval
    auto it = std::upper_bound(recs.begin(), recs.end(), t_s,
                               [](double t, const WindDirRecord& r){ return t < r.time_s; });
    auto prev = std::prev(it);
    const double alpha = (t_s - prev->time_s) / (it->time_s - prev->time_s);

    // Linearly interpolate speed
    const double spd = prev->speed_ms + alpha * (it->speed_ms - prev->speed_ms);

    // Circular interpolation for direction (shortest arc)
    double d0 = prev->dir_deg;
    double d1 = it->dir_deg;
    double dd = d1 - d0;
    // Wrap dd to [-180, 180]
    while (dd >  180.0) dd -= 360.0;
    while (dd < -180.0) dd += 360.0;
    const double dir = d0 + alpha * dd;

    // Convert to (Ux, Uy): meteorological convention → math angle
    const double ang = (270.0 - dir) * WDSC_PI / 180.0;
    return {spd * std::cos(ang), spd * std::sin(ang)};
}

#endif // WIND_DIR_SCHEDULE_H
