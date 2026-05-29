// ============================================================================
// Foliar Moisture Content (FMC) Seasonal / Phenological Schedule
// ============================================================================
// Provides a date-driven FMC value that reflects the seasonal green-up and
// brown-off cycle of forest canopy foliage.  FARSITE uses this same
// phenological curve to adjust the Van Wagner (1977) crown fire initiation
// threshold: higher FMC → higher critical surface fire intensity needed to
// initiate crown fire.
//
// Two modes are supported:
//
//  1. CSV schedule file (``fmc_schedule_file``):
//     Two-column CSV (day_of_year [1–365], fmc_pct [%]).
//     Simulation time is mapped to a day-of-year using start_doy (default: 1).
//     Values are linearly interpolated between rows; clamped outside range.
//
//     Example:
//       # doy, fmc_pct
//       1,  85
//       80, 85
//       120, 110
//       180, 140
//       240, 120
//       300, 85
//       365, 85
//
//  2. Built-in FARSITE-style parametric curve (``fmc_use_farsite_curve = 1``):
//     Based on Rothermel (1983) / Finney (1998) seasonal phenology:
//       - Winter dormant phase (doy < spring_start):        FMC = FMC_min
//       - Green-up   (spring_start ≤ doy < summer_peak):   linear rise to FMC_max
//       - Summer maximum (summer_peak ≤ doy < fall_start): FMC = FMC_max
//       - Curing     (fall_start ≤ doy ≤ fall_end):        linear fall to FMC_min
//       - Post-curing (doy > fall_end):                     FMC = FMC_min
//
//     Default breakpoints match Southern California chaparral phenology
//     (Rothermel 1983):
//       spring_start = 90  (early April)
//       summer_peak  = 150 (late May)
//       fall_start   = 240 (late August)
//       fall_end     = 300 (late October)
//       FMC_min      = 85  [%]  (dormant / cured foliage)
//       FMC_max      = 140 [%]  (peak green foliage)
//
// Usage:
//   FMCSchedule fmc_sched;
//   load_fmc_schedule("fmc.csv", fmc_sched);
//   // or: build_farsite_fmc_curve(fmc_sched, ...);
//   // at each timestep:
//   Real fmc_pct = get_fmc_at_time(fmc_sched, sim_time_s);
//   inputs.crown.FMC = fmc_pct;
//
// References:
//   Rothermel, R.C. (1983). How to Predict the Spread and Intensity of Forest
//     and Range Fires. Gen. Tech. Rep. INT-143. USDA Forest Service.
//   Finney, M.A. (1998). FARSITE: Fire Area Simulator — Model Development and
//     Evaluation. USDA Forest Service Research Paper RMRS-RP-4.
//   Van Wagner, C.E. (1977). Conditions for the start and spread of crown fire.
//     Canadian Journal of Forest Research, 7(1), 23-34.
// ============================================================================

#ifndef FMC_SCHEDULE_H
#define FMC_SCHEDULE_H

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <AMReX.H>

// One entry in the FMC schedule
struct FMCRecord {
    double doy;      // day of year [1–365]
    double fmc_pct;  // foliar moisture content [%]
};

// The full FMC schedule (sorted by doy)
struct FMCSchedule {
    std::vector<FMCRecord> records;
    int    start_doy;    // day-of-year corresponding to simulation t = 0 (default: 1)
    bool   empty() const { return records.empty(); }
};

// ============================================================================
// load_fmc_schedule
//
// Read a two-column CSV file (doy, fmc_pct) into fmc_sched.
// Lines starting with # or ! are treated as comments.
// ============================================================================
inline void load_fmc_schedule(const std::string& filename, FMCSchedule& sched,
                               int start_doy = 1)
{
    sched.records.clear();
    sched.start_doy = start_doy;

    std::ifstream f(filename);
    if (!f.is_open()) {
        amrex::Abort("FMCSchedule: cannot open file: " + filename);
    }

    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#' || line[0] == '!') continue;
        // Replace commas with spaces for uniform parsing
        for (char& c : line) if (c == ',') c = ' ';
        std::istringstream iss(line);
        FMCRecord r;
        if (iss >> r.doy >> r.fmc_pct) {
            sched.records.push_back(r);
        }
    }

    if (sched.records.empty()) {
        amrex::Abort("FMCSchedule: no valid records in file: " + filename);
    }

    // Sort by doy
    std::sort(sched.records.begin(), sched.records.end(),
              [](const FMCRecord& a, const FMCRecord& b){ return a.doy < b.doy; });

    amrex::Print() << "FMC schedule: loaded " << sched.records.size()
                   << " records from " << filename
                   << "  start_doy=" << start_doy << "\n";
}

// ============================================================================
// build_farsite_fmc_curve
//
// Populate fmc_sched from the built-in parametric phenological curve.
// ============================================================================
inline void build_farsite_fmc_curve(FMCSchedule& sched,
                                    int    start_doy    = 1,
                                    int    spring_start = 90,
                                    int    summer_peak  = 150,
                                    int    fall_start   = 240,
                                    int    fall_end     = 300,
                                    double fmc_min      = 85.0,
                                    double fmc_max      = 140.0)
{
    sched.records.clear();
    sched.start_doy = start_doy;

    // Build curve at 10-day intervals (enough resolution for interpolation)
    for (int doy = 1; doy <= 365; doy += 10) {
        double fmc;
        if (doy < spring_start) {
            fmc = fmc_min;
        } else if (doy < summer_peak) {
            double t = double(doy - spring_start) / double(summer_peak - spring_start);
            fmc = fmc_min + t * (fmc_max - fmc_min);
        } else if (doy < fall_start) {
            fmc = fmc_max;
        } else if (doy <= fall_end) {
            double t = double(doy - fall_start) / double(fall_end - fall_start);
            fmc = fmc_max - t * (fmc_max - fmc_min);
        } else {
            fmc = fmc_min;
        }
        sched.records.push_back({double(doy), fmc});
    }
    // Ensure we cover doy=365
    sched.records.push_back({365.0, fmc_min});

    amrex::Print() << "FMC schedule: built parametric FARSITE phenology curve"
                   << "  spring_start=" << spring_start
                   << "  summer_peak=" << summer_peak
                   << "  fall_start=" << fall_start
                   << "  fall_end=" << fall_end
                   << "  FMC=[" << fmc_min << "," << fmc_max << "]%\n";
}

// ============================================================================
// get_fmc_at_time
//
// Return the FMC [%] at simulation time t_s [s], using the start_doy anchor
// to convert elapsed time to a day of year.
// Linear interpolation; clamped outside [records.front, records.back].
// ============================================================================
inline double get_fmc_at_time(const FMCSchedule& sched, double t_s)
{
    if (sched.empty()) return 100.0;  // sensible default if schedule not loaded

    // Convert simulation seconds to a day of year (wrapping annually).
    // Use modulo arithmetic for O(1) wrapping regardless of elapsed time.
    const double seconds_per_day = 86400.0;
    double doy = sched.start_doy + t_s / seconds_per_day;
    // Wrap into [1, 365] with a single modulo operation
    doy = std::fmod(doy - 1.0, 365.0);
    if (doy < 0.0) doy += 365.0;
    doy += 1.0;  // restore to [1, 365] range

    const auto& recs = sched.records;
    if (doy <= recs.front().doy) return recs.front().fmc_pct;
    if (doy >= recs.back().doy)  return recs.back().fmc_pct;

    // Binary search for bracketing interval
    auto it = std::upper_bound(recs.begin(), recs.end(), doy,
                               [](double d, const FMCRecord& r){ return d < r.doy; });
    if (it == recs.begin()) return recs.front().fmc_pct;
    auto prev = std::prev(it);

    const double t = (doy - prev->doy) / (it->doy - prev->doy);
    return prev->fmc_pct + t * (it->fmc_pct - prev->fmc_pct);
}

// ============================================================================
// build_sinusoidal_fmc_phenology
//
// Build FMC schedule using a sinusoidal phenology curve based on day-of-year.
// This provides a smooth seasonal variation suitable for temperate and
// subtropical climates.
//
// Formula:
//   FMC(doy) = FMC_min + (FMC_max - FMC_min) *
//              sin²(π * (doy - doy_offset) / 365)
//
// where doy_offset controls the phase (default: 0 = peak at mid-year ~day 182).
//
// References:
//   Jolly, W.M., et al. (2005). A generalized, bioclimatic index to predict
//     foliar phenology in response to climate. Global Change Biology, 11(4):619-632.
//   Dennison, P.E., et al. (2008). Use of normalized difference water index for
//     monitoring live fuel moisture. International J. of Remote Sensing, 29(8):2369-2375.
// ============================================================================
inline void build_sinusoidal_fmc_phenology(FMCSchedule& sched,
                                           int    start_doy    = 1,
                                           double fmc_min      = 85.0,
                                           double fmc_max      = 140.0,
                                           double doy_offset   = 0.0)
{
    sched.records.clear();
    sched.start_doy = start_doy;

    const double pi = 3.141592653589793;
    
    // Build curve at 5-day intervals for smooth interpolation
    for (int doy = 1; doy <= 365; doy += 5) {
        const double phase = pi * double(doy - doy_offset) / 365.0;
        const double sin_val = std::sin(phase);
        const double fmc = fmc_min + (fmc_max - fmc_min) * sin_val * sin_val;
        sched.records.push_back({double(doy), fmc});
    }
    // Ensure we cover doy=365
    const double phase_365 = pi * double(365 - doy_offset) / 365.0;
    const double sin_365 = std::sin(phase_365);
    const double fmc_365 = fmc_min + (fmc_max - fmc_min) * sin_365 * sin_365;
    sched.records.push_back({365.0, fmc_365});

    amrex::Print() << "FMC schedule: built sinusoidal phenology curve"
                   << "  FMC=[" << fmc_min << "," << fmc_max << "]%"
                   << "  doy_offset=" << doy_offset << "\n";
}

// ============================================================================
// GDDState
//
// Growing Degree Day (GDD) accumulator for phenology modeling.
// Tracks cumulative heat units to predict foliar development.
// ============================================================================
struct GDDState {
    double gdd_cumulative = 0.0;  // Accumulated growing degree days [°C·day]
    bool initialized = false;
};

// ============================================================================
// compute_gdd_fmc
//
// Compute FMC based on Growing Degree Days (GDD) accumulation.
// Simple empirical model where FMC increases sigmoidally with accumulated GDD.
//
// GDD accumulation (daily average temp T_avg):
//   dGDD = max(0, T_avg - T_base)
//
// FMC(GDD):
//   FMC = FMC_min + (FMC_max - FMC_min) / (1 + exp(-k * (GDD - GDD_mid)))
//
// where:
//   T_base   – base temperature for growth [°C] (default: 5.0°C)
//   GDD_mid  – midpoint GDD for green-up [°C·day] (default: 500.0)
//   k        – logistic growth rate (default: 0.01)
//
// References:
//   Jolly, W.M., et al. (2005). A generalized, bioclimatic index to predict
//     foliar phenology in response to climate. Global Change Biology, 11(4):619-632.
// ============================================================================
inline double compute_gdd_fmc(GDDState& gdd_state,
                              double T_avg_C,
                              double dt_days,
                              double fmc_min   = 85.0,
                              double fmc_max   = 140.0,
                              double T_base    = 5.0,
                              double GDD_mid   = 500.0,
                              double k         = 0.01)
{
    // Accumulate GDD
    const double dGDD = std::max(0.0, T_avg_C - T_base) * dt_days;
    gdd_state.gdd_cumulative += dGDD;
    
    // Reset GDD at end of year (simple annual cycle)
    // In practice this would be triggered by doy==1 or similar
    if (gdd_state.gdd_cumulative > 2000.0) {
        gdd_state.gdd_cumulative = std::fmod(gdd_state.gdd_cumulative, 2000.0);
    }
    
    // Logistic FMC response
    const double exponent = -k * (gdd_state.gdd_cumulative - GDD_mid);
    const double sigmoid = 1.0 / (1.0 + std::exp(exponent));
    const double fmc = fmc_min + (fmc_max - fmc_min) * sigmoid;
    
    return fmc;
}

#endif // FMC_SCHEDULE_H
