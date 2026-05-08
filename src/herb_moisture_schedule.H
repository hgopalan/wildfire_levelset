// ============================================================================
// Live Herbaceous Moisture / Curing Schedule
// ============================================================================
// Provides a date-driven live herbaceous moisture content (M_lh) schedule that
// reflects the seasonal curing of grasses and forbs.  FARSITE uses this to
// represent the transition from green (high M_lh) to cured (low M_lh) herbage
// as the season progresses, which significantly increases surface fire spread
// in grass-dominated fuel models.
//
// Two modes are supported:
//
//  1. CSV schedule file (``herb_moisture_schedule.file``):
//     Two-column CSV (day_of_year [1–365], m_lh_pct [%]).
//     Values are linearly interpolated between rows; clamped outside range.
//
//     Example:
//       # doy, m_lh_pct
//       1,  30
//       90, 120
//       150, 180
//       240, 100
//       300, 30
//       365, 30
//
//  2. Built-in parametric curing curve (``use_curing_curve = 1``):
//     Simple seasonal model:
//       - Winter dormant (doy < spring_start):          M_lh = M_lh_min
//       - Green-up (spring_start ≤ doy < summer_peak): linear rise to M_lh_max
//       - Summer (summer_peak ≤ doy < fall_start):     M_lh = M_lh_max
//       - Curing (fall_start ≤ doy ≤ fall_end):        linear fall to M_lh_min
//       - Post-curing (doy > fall_end):                 M_lh = M_lh_min
//
//     Default breakpoints match typical temperate-grassland phenology:
//       spring_start = 90  (early April)
//       summer_peak  = 150 (late May)
//       fall_start   = 200 (late July – onset of curing)
//       fall_end     = 270 (late September – fully cured)
//       M_lh_min     = 30  [%]  (fully cured)
//       M_lh_max     = 180 [%]  (peak green)
//
// Usage:
//   HerbMoistureSchedule sched;
//   load_herb_moisture_schedule("herb.csv", sched, 1 /*start_doy*/);
//   // at each timestep:
//   double m_lh_pct = get_herb_moisture_at_time(sched, sim_time_s);
//   inputs.rothermel.M_lh = m_lh_pct / 100.0;   // convert % → fraction
//
// References:
//   Cheney, N.P. & Gould, J.S. (1995). "Fire Growth in Grassland Fuels."
//     International Journal of Wildland Fire, 5(4), 237-247.
//   Finney, M.A. (1998). FARSITE: Fire Area Simulator — Model Development and
//     Evaluation. USDA Forest Service Research Paper RMRS-RP-4.
// ============================================================================

#ifndef HERB_MOISTURE_SCHEDULE_H
#define HERB_MOISTURE_SCHEDULE_H

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <AMReX.H>

// One entry in the live herbaceous moisture schedule
struct HerbMoistureRecord {
    double doy;      // day of year [1–365]
    double m_lh_pct; // live herbaceous moisture content [%]
};

// The full schedule (sorted by doy)
struct HerbMoistureSchedule {
    std::vector<HerbMoistureRecord> records;
    int    start_doy;    // day-of-year at simulation t = 0 (default: 1)
    bool   empty() const { return records.empty(); }
};

// ============================================================================
// load_herb_moisture_schedule
//
// Read a two-column CSV file (doy, m_lh_pct) into sched.
// Lines starting with # or ! are treated as comments.
// ============================================================================
inline void load_herb_moisture_schedule(const std::string& filename,
                                         HerbMoistureSchedule& sched,
                                         int start_doy = 1)
{
    sched.records.clear();
    sched.start_doy = start_doy;

    std::ifstream f(filename);
    if (!f.is_open()) {
        amrex::Abort("HerbMoistureSchedule: cannot open file: " + filename);
    }

    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#' || line[0] == '!') continue;
        for (char& c : line) if (c == ',') c = ' ';
        std::istringstream iss(line);
        HerbMoistureRecord r;
        if (iss >> r.doy >> r.m_lh_pct) {
            sched.records.push_back(r);
        }
    }

    if (sched.records.empty()) {
        amrex::Abort("HerbMoistureSchedule: no valid records in file: " + filename);
    }

    std::sort(sched.records.begin(), sched.records.end(),
              [](const HerbMoistureRecord& a, const HerbMoistureRecord& b){
                  return a.doy < b.doy;
              });

    amrex::Print() << "Herb moisture schedule: loaded " << sched.records.size()
                   << " records from " << filename
                   << "  start_doy=" << start_doy << "\n";
}

// ============================================================================
// build_herb_curing_curve
//
// Populate sched from the built-in parametric seasonal curing curve.
// ============================================================================
inline void build_herb_curing_curve(HerbMoistureSchedule& sched,
                                     int    start_doy    = 1,
                                     int    spring_start = 90,
                                     int    summer_peak  = 150,
                                     int    fall_start   = 200,
                                     int    fall_end     = 270,
                                     double m_lh_min     = 30.0,
                                     double m_lh_max     = 180.0)
{
    sched.records.clear();
    sched.start_doy = start_doy;

    for (int doy = 1; doy <= 365; doy += 10) {
        double mlh;
        if (doy < spring_start) {
            mlh = m_lh_min;
        } else if (doy < summer_peak) {
            double t = double(doy - spring_start) / double(summer_peak - spring_start);
            mlh = m_lh_min + t * (m_lh_max - m_lh_min);
        } else if (doy < fall_start) {
            mlh = m_lh_max;
        } else if (doy <= fall_end) {
            double t = double(doy - fall_start) / double(fall_end - fall_start);
            mlh = m_lh_max - t * (m_lh_max - m_lh_min);
        } else {
            mlh = m_lh_min;
        }
        sched.records.push_back({double(doy), mlh});
    }
    sched.records.push_back({365.0, m_lh_min});

    amrex::Print() << "Herb moisture schedule: built parametric curing curve"
                   << "  spring_start=" << spring_start
                   << "  summer_peak=" << summer_peak
                   << "  fall_start=" << fall_start
                   << "  fall_end=" << fall_end
                   << "  M_lh=[" << m_lh_min << "," << m_lh_max << "]%\n";
}

// ============================================================================
// get_herb_moisture_at_time
//
// Return the live herbaceous moisture [%] at simulation time t_s [s].
// Linear interpolation; clamped outside schedule range.
// ============================================================================
inline double get_herb_moisture_at_time(const HerbMoistureSchedule& sched,
                                         double t_s)
{
    if (sched.empty()) return 100.0;  // sensible default if not loaded

    const double seconds_per_day = 86400.0;
    double doy = sched.start_doy + t_s / seconds_per_day;
    doy = std::fmod(doy - 1.0, 365.0);
    if (doy < 0.0) doy += 365.0;
    doy += 1.0;  // restore to [1, 365]

    const auto& recs = sched.records;
    if (doy <= recs.front().doy) return recs.front().m_lh_pct;
    if (doy >= recs.back().doy)  return recs.back().m_lh_pct;

    auto it = std::upper_bound(recs.begin(), recs.end(), doy,
                               [](double d, const HerbMoistureRecord& r){
                                   return d < r.doy;
                               });
    if (it == recs.begin()) return recs.front().m_lh_pct;
    auto prev = std::prev(it);

    const double t = (doy - prev->doy) / (it->doy - prev->doy);
    return prev->m_lh_pct + t * (it->m_lh_pct - prev->m_lh_pct);
}

#endif // HERB_MOISTURE_SCHEDULE_H
