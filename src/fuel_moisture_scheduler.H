// ============================================================================
// Time-Varying Fuel Moisture Scheduler
// ============================================================================
// Reads a FARSITE-compatible fuel moisture schedule file (.fmd format) and
// provides per-simulation-time moisture values that are used to update the
// global RothermelParams at each timestep.
//
// FMD file format (whitespace-delimited; # or ! lines are comments):
//
//   MONTH  DAY  HOUR  PRECIP  NUM_MODELS
//   MODEL  1HR%  10HR%  100HR%  LHERB%  LWOOD%
//   ...                                 (one row per fuel model)
//   [next timestamp block]
//   ...
//
// HOUR is 0-2300 (HHMM format, e.g. 1400 = 14:00).
// Moisture values are integer percentages (8 → 0.08 fraction).
// PRECIP is precipitation in hundredths of an inch (informational, not used
// here; FARSITE uses it to adjust moistures internally).
//
// The scheduler:
//   • Converts each timestamp to seconds-since-start-of-simulation using
//     an optional reference start date (year/month/day/hour from inputs).
//     When no start date is given, the first record's timestamp is t = 0.
//   • Provides get_moisture_at_time() which returns linearly interpolated
//     moisture values for any simulation time [s].
//   • A fuel_model of 0 means "apply to all models"; specific model numbers
//     override the global entry.
//
// Usage:
//   FuelMoistureSchedule sched;
//   load_fuel_moisture_schedule("fire.fmd", sched);
//   // at each timestep:
//   RothermelMoistures m = get_moisture_at_time(sched, sim_time_s, model_id);
//   inputs.rothermel.M_d1 = m.M_d1;  // etc.
//
// References:
//   Finney, M.A. (2004). FARSITE: Fire Area Simulator — Model Development
//     and Evaluation. USDA Forest Service Research Paper RMRS-RP-4.
// ============================================================================

#ifndef FUEL_MOISTURE_SCHEDULER_H
#define FUEL_MOISTURE_SCHEDULER_H

#include "constants.H"
#include <AMReX.H>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "parse_inputs.H"

// Per-fuel-model moisture snapshot at one timestamp
struct MoistureSnapshot {
    double time_s;   // seconds since simulation start
    int    fuel_model;
    float  M_d1;     // 1-hr   dead [fraction]
    float  M_d10;    // 10-hr  dead [fraction]
    float  M_d100;   // 100-hr dead [fraction]
    float  M_lh;     // live herbaceous [fraction]
    float  M_lw;     // live woody [fraction]
};

// The full schedule: sorted by time_s, one entry per (timestamp, fuel_model)
struct FuelMoistureSchedule {
    std::vector<MoistureSnapshot> records;  // sorted by time_s
    bool empty() const { return records.empty(); }
};

// Returned by get_moisture_at_time()
struct RothermelMoistures {
    float M_d1   = 0.08f;
    float M_d10  = 0.08f;
    float M_d100 = 0.08f;
    float M_d1000 = 0.15f;  // 1000-hr dead moisture (Nelson 2000 / Rothermel 1983 ext.)
    float M_lh   = 0.90f;
    float M_lw   = 1.20f;
};

// Empirical scaling factor used when approximating 1000-hr dead fuel moisture
// as a multiple of the 100-hr value.  Based on the Rothermel (1983) time-lag
// progressions where 1000-hr fuels equilibrate ~1.5× more slowly than 100-hr fuels.
// Reference: Rothermel, R.C. (1983), How to Predict the Spread and Intensity of
// Forest and Range Fires, Gen. Tech. Rep. INT-143, USDA Forest Service, Table 1
// (timelag / stick moisture relationships for dead fuels).
static constexpr float M_D1000_SCALE = 1.5f;

// ============================================================================
// load_fuel_moisture_schedule
//
// Parses the FMD file into a FuelMoistureSchedule.
// start_year/month/day/hour: if all are 0, the first record's clock time is
// treated as t = 0 [s].  Otherwise timestamps are converted to seconds from
// the given reference.
// ============================================================================
inline void load_fuel_moisture_schedule(
    const std::string&   filename,
    FuelMoistureSchedule& sched,
    int start_year  = 0,
    int start_month = 0,
    int start_day   = 0,
    int start_hour  = 0)
{
    sched.records.clear();

    std::ifstream f(filename);
    if (!f.is_open())
        amrex::Abort("Cannot open fuel moisture schedule file: " + filename);

    // Helper: compute days since epoch (1 Jan 1900) – simple Gregorian
    auto days_since_epoch = [](int y, int m, int d) -> long {
        // Zeller-like computation (sufficient for ordering; not ISO epoch)
        if (m < 3) { m += 12; --y; }
        long J = 365L * y + y / 4 - y / 100 + y / 400
                + (306 * (m + 1)) / 10 + d - 719591;
        return J;
    };

    long   ref_day_s = 0;
    bool   ref_set   = (start_year != 0);

    if (ref_set) {
        ref_day_s = days_since_epoch(start_year, start_month, start_day) * 86400L
                  + start_hour * 3600L;
    }

    std::string line;
    int line_num = 0;

    while (std::getline(f, line)) {
        ++line_num;
        size_t pos = line.find_first_not_of(" \t\r\n");
        if (pos == std::string::npos) continue;
        line = line.substr(pos);
        if (line[0] == '#' || line[0] == '!') continue;

        std::istringstream iss(line);
        // Try to read header line: MONTH DAY HOUR PRECIP NUM_MODELS
        int month, day, hour_hhmm, precip, nmodels;
        if (!(iss >> month >> day >> hour_hhmm >> precip >> nmodels)) continue;

        // Silently skip clearly invalid values
        if (month < 1 || month > 12 || day < 1 || day > 31 || nmodels < 1) continue;

        // Convert HHMM → seconds into day
        int hh = hour_hhmm / 100;
        int mm = hour_hhmm % 100;
        long block_s = static_cast<long>(hh) * 3600L + static_cast<long>(mm) * 60L;

        // For year we use a running year estimate (FMD files often omit the year)
        // Use a nominal year 2000; ordering/interpolation only needs relative times.
        long block_day_s =
            days_since_epoch(2000, month, day) * 86400L + block_s;

        if (!ref_set) {
            // First block defines t = 0
            ref_day_s = block_day_s;
            ref_set   = true;
        }

        double t_s = static_cast<double>(block_day_s - ref_day_s);

        // Read nmodels moisture rows
        for (int m_idx = 0; m_idx < nmodels; ++m_idx) {
            std::string mline;
            while (std::getline(f, mline)) {
                ++line_num;
                size_t mpos = mline.find_first_not_of(" \t\r\n");
                if (mpos == std::string::npos) continue;
                mline = mline.substr(mpos);
                if (mline[0] == '#' || mline[0] == '!') continue;
                break;
            }
            std::istringstream miss(mline);
            int fm;
            float m1, m10, m100, mlh, mlw;
            if (!(miss >> fm >> m1 >> m10 >> m100 >> mlh >> mlw)) continue;

            MoistureSnapshot snap;
            snap.time_s    = t_s;
            snap.fuel_model = fm;
            snap.M_d1   = m1   / 100.0f;
            snap.M_d10  = m10  / 100.0f;
            snap.M_d100 = m100 / 100.0f;
            snap.M_lh   = mlh  / 100.0f;
            snap.M_lw   = mlw  / 100.0f;
            sched.records.push_back(snap);
        }
    }

    // Sort by time, then fuel_model
    std::stable_sort(sched.records.begin(), sched.records.end(),
        [](const MoistureSnapshot& a, const MoistureSnapshot& b) {
            if (a.time_s != b.time_s) return a.time_s < b.time_s;
            return a.fuel_model < b.fuel_model;
        });

    amrex::Print() << "Loaded " << sched.records.size()
                   << " fuel moisture records from " << filename << "\n";
}

// ============================================================================
// get_moisture_at_time
//
// Returns linearly interpolated moisture values for a given simulation time
// [s] and a specific fuel model code.
//
// Search order: exact fuel_model match, then fallback to fuel_model == 0
// (global / "applies to all" entry).  Linear interpolation between the two
// nearest timestamps; extrapolation clamps to the first / last record.
//
// If the schedule is empty the function returns the supplied defaults
// (default_m) unchanged.
// ============================================================================
inline RothermelMoistures get_moisture_at_time(
    const FuelMoistureSchedule& sched,
    double sim_time_s,
    int    fuel_model_code,
    const RothermelMoistures& defaults = {})
{
    if (sched.empty()) return defaults;

    // Collect all distinct timestamps for this fuel_model (or 0)
    std::vector<const MoistureSnapshot*> matches;
    for (const auto& r : sched.records) {
        if (r.fuel_model == fuel_model_code || r.fuel_model == 0)
            matches.push_back(&r);
    }
    if (matches.empty()) return defaults;

    // Find bracketing snapshots by time
    // Use last record whose time_s <= sim_time_s as lo; first > sim_time_s as hi
    const MoistureSnapshot* lo = matches.front();
    const MoistureSnapshot* hi = matches.back();

    for (const auto* r : matches) {
        if (r->time_s <= sim_time_s) lo = r;
        if (r->time_s > sim_time_s) { hi = r; break; }
    }

    if (lo->time_s >= hi->time_s || lo == hi) {
        // Clamp / exact match
        RothermelMoistures r;
        r.M_d1    = lo->M_d1;
        r.M_d10   = lo->M_d10;
        r.M_d100  = lo->M_d100;
        r.M_d1000 = lo->M_d100 * M_D1000_SCALE;  // proxy: 1000-hr ≈ 1.5× 100-hr
        r.M_lh    = lo->M_lh;
        r.M_lw    = lo->M_lw;
        return r;
    }

    // Linear interpolation
    double alpha = (sim_time_s - lo->time_s) / (hi->time_s - lo->time_s);
    alpha = std::max(0.0, std::min(1.0, alpha));
    float a = static_cast<float>(alpha);

    RothermelMoistures m;
    m.M_d1    = lo->M_d1   + a * (hi->M_d1   - lo->M_d1);
    m.M_d10   = lo->M_d10  + a * (hi->M_d10  - lo->M_d10);
    m.M_d100  = lo->M_d100 + a * (hi->M_d100 - lo->M_d100);
    m.M_d1000 = m.M_d100 * M_D1000_SCALE;  // proxy: 1000-hr ≈ 1.5× 100-hr
    m.M_lh    = lo->M_lh   + a * (hi->M_lh   - lo->M_lh);
    m.M_lw    = lo->M_lw   + a * (hi->M_lw   - lo->M_lw);
    return m;
}

// ============================================================================
// compute_diurnal_emc
//
// Returns dead fuel moisture content values computed from a sinusoidal diurnal
// cycle of temperature and relative humidity using the Nelson (2000) /
// Simard (1968) EMC equations.
//
// The diurnal cycle is:
//   T(t)  = T_mean + A_T  * sin(2π * (t + t_peak - t_ref) / 86400)
//   RH(t) = RH_mean - A_RH * sin(2π * (t + t_peak - t_ref) / 86400)
//
// where t_peak is the clock time [s from midnight] when T is maximum,
// t_ref is the clock time [s from midnight] at simulation start.
//
// EMC equations (Simard 1968):
//   Region 1 (RH < 10%): EMC = 0.03229 + 0.281073*rh - 0.000578*T*rh
//   Region 2 (10% ≤ RH < 50%): EMC = 2.22749 + 0.160107*rh - 0.014784*T
//   Region 3 (RH ≥ 50%): EMC = 21.0606 + 0.005565*rh² - 0.00035*T*rh - 0.483199*rh
// where rh = RH/100 and T is in °C, result in %.
//
// Timelag corrections (Nelson 2000):
//   M_1hr   = EMC [fraction]
//   M_10hr  = EMC * 1.10 [fraction]
//   M_100hr = EMC * 1.30 [fraction]
//
// Parameters:
//   dp         – DiurnalMoistureParams (T_min/max, RH_min/max, t_start_s, t_T_peak_s)
//   sim_time_s – current simulation time [s] (relative to the simulation start)
//   defaults   – fallback moisture values (returned unchanged if not enabled)
// ============================================================================
inline RothermelMoistures compute_diurnal_emc(
    const InputParameters::DiurnalMoistureParams& dp,
    double sim_time_s,
    const RothermelMoistures& defaults = {})
{
    if (dp.enable != 1) return defaults;

    constexpr double TWO_PI = WildfireConst::TWO_PI;
    constexpr double DAY_S  = WildfireConst::DAY_S;

    // Phase: how far into the day we are when the simulation starts,
    // shifted so that the peak of the sine occurs at t_T_peak_s.
    // phase = (t_start_s + sim_time_s - t_T_peak_s) / 86400 * 2π
    double phase = TWO_PI * (dp.t_start_s + sim_time_s - dp.t_T_peak_s) / DAY_S;

    double T_mean  = 0.5 * (dp.T_max + dp.T_min);
    double A_T     = 0.5 * (dp.T_max - dp.T_min);
    double RH_mean = 0.5 * (dp.RH_max + dp.RH_min);
    double A_RH    = 0.5 * (dp.RH_max - dp.RH_min);

    double T_C  = T_mean  + A_T  * std::sin(phase);
    double RH_p = RH_mean - A_RH * std::sin(phase);  // RH in % (anti-phase with T)

    // Clamp to physically meaningful ranges
    T_C  = std::max(-40.0, std::min(60.0,  T_C));
    RH_p = std::max(1.0,   std::min(100.0, RH_p));

    double rh = RH_p / 100.0;

    // Nelson (2000) / Simard (1968) EMC [%]
    double emc_pct;
    if (RH_p < 10.0) {
        emc_pct = 0.03229 + 0.281073 * rh - 0.000578 * T_C * rh;
    } else if (RH_p < 50.0) {
        emc_pct = 2.22749 + 0.160107 * rh - 0.014784 * T_C;
    } else {
        emc_pct = 21.0606 + 0.005565 * rh * rh - 0.00035 * T_C * rh - 0.483199 * rh;
    }
    emc_pct = std::max(0.5, emc_pct);  // physical lower bound

    RothermelMoistures m;
    m.M_d1    = static_cast<float>(emc_pct / 100.0);
    m.M_d10   = static_cast<float>(emc_pct * 1.10 / 100.0);
    m.M_d100  = static_cast<float>(emc_pct * 1.30 / 100.0);
    m.M_d1000 = static_cast<float>(emc_pct * M_D1000_SCALE / 100.0);  // 1000-hr approx (Nelson 2000 ext.)
    // Live fuel moisture is not well-modelled by Nelson EMC; keep defaults
    m.M_lh   = defaults.M_lh;
    m.M_lw   = defaults.M_lw;
    return m;
}

#endif // FUEL_MOISTURE_SCHEDULER_H
