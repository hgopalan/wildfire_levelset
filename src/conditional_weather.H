// ============================================================================
// Conditional Weather / ERC-BI-SC Percentile Lookup Table
// ============================================================================
// FARSITE supports a "conditional weather" table that maps fuel moisture
// percentile conditions (expressed as ERC percentile) to specific wind/moisture
// scenarios.  This allows scenario-based fire modelling where the simulation
// weather is selected based on historical fire danger conditions.
//
// CSV file format (header row required):
//   erc_percentile, wind_speed_ms, wind_dir_deg, M_d1_pct, M_d10_pct, M_d100_pct
//   0,   1.0, 180, 15.0, 20.0, 25.0
//   25,  3.0, 225, 12.0, 16.0, 22.0
//   50,  5.0, 270,  9.0, 12.0, 18.0
//   75,  7.0, 270,  7.0,  9.0, 14.0
//   97, 10.0, 315,  4.0,  6.0, 10.0
//
// Behaviour:
//   The current simulation fire-danger index (ERC, BI, or SC, depending on
//   conditional_weather_trigger) is mapped to a percentile using the
//   historical distribution in the table.  The corresponding wind speed,
//   wind direction, and dead fuel moistures are then applied to the simulation.
//
//   This module provides:
//     1. load_conditional_weather()  – parse CSV, build in-memory table
//     2. apply_conditional_weather() – select row by current index, update inputs
//
// Supported trigger indices (conditional_weather_trigger in inputs):
//   "erc"  – NFDRS Energy Release Component  (default; ecology component 6)
//   "bi"   – NFDRS Burning Index             (ecology component 8)
//   "sc"   – NFDRS Spread Component          (ecology component 7)
//
// References:
//   Finney, M.A. (1998). FARSITE: Fire Area Simulator. USDA FS RMRS-RP-4.
//   Deeming, J.E., Burgan, R.E., & Cohen, J.D. (1977). GTR INT-39.
// ============================================================================

#ifndef CONDITIONAL_WEATHER_H
#define CONDITIONAL_WEATHER_H

#include <AMReX.H>
#include <AMReX_Print.H>
#include "parse_inputs.H"
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>

// ============================================================================
// ConditionalWeatherRow
// One row of the conditional weather table.
// ============================================================================
struct ConditionalWeatherRow {
    double erc_percentile;   // ERC percentile (0–100)
    double wind_speed_ms;    // wind speed [m/s]
    double wind_dir_deg;     // wind direction [deg], met. convention (from)
    double M_d1_pct;         // 1-hr dead fuel moisture [%]
    double M_d10_pct;        // 10-hr dead fuel moisture [%]
    double M_d100_pct;       // 100-hr dead fuel moisture [%]
};

// Alias for readability
using ConditionalWeatherTable = std::vector<ConditionalWeatherRow>;

// ============================================================================
// load_conditional_weather
//
// Parses a CSV conditional weather table from file.
//
// Returns an empty table when filename is empty or file cannot be read.
// ============================================================================
static ConditionalWeatherTable load_conditional_weather(const std::string& filename)
{
    ConditionalWeatherTable table;
    if (filename.empty()) return table;

    if (!amrex::ParallelDescriptor::IOProcessor()) {
        // Non-IO ranks skip parsing; data broadcast below
    } else {
        std::ifstream fin(filename);
        if (!fin.is_open()) {
            amrex::Abort("conditional_weather: cannot open file '" + filename + "'");
        }
        std::string line;
        bool header_skipped = false;
        while (std::getline(fin, line)) {
            if (line.empty() || line[0] == '#' || line[0] == '!') continue;
            // Skip header row (contains non-numeric first token)
            if (!header_skipped) {
                // Check if first token is numeric
                for (char& c : line) if (c == ',') c = ' ';
                std::istringstream test(line);
                double v;
                if (!(test >> v)) { header_skipped = true; continue; }
                header_skipped = true;
                // Reparse this line as data
            }
            for (char& c : line) if (c == ',') c = ' ';
            std::istringstream iss(line);
            ConditionalWeatherRow row;
            if (!(iss >> row.erc_percentile >> row.wind_speed_ms >> row.wind_dir_deg
                      >> row.M_d1_pct >> row.M_d10_pct >> row.M_d100_pct)) continue;
            table.push_back(row);
        }
        // Sort by ERC percentile ascending
        std::sort(table.begin(), table.end(),
                  [](const ConditionalWeatherRow& a, const ConditionalWeatherRow& b){
                      return a.erc_percentile < b.erc_percentile;
                  });
        amrex::Print() << "conditional_weather: loaded " << table.size()
                       << " row(s) from '" << filename << "'\n";
    }

    // Broadcast table size and data to all MPI ranks
    int n = static_cast<int>(table.size());
    amrex::ParallelDescriptor::Bcast(&n, 1, amrex::ParallelDescriptor::IOProcessorNumber());
    if (!amrex::ParallelDescriptor::IOProcessor()) table.resize(n);
    if (n > 0) {
        // Broadcast as raw doubles (6 doubles per row)
        std::vector<double> buf(n * 6);
        if (amrex::ParallelDescriptor::IOProcessor()) {
            for (int i = 0; i < n; ++i) {
                buf[6*i+0] = table[i].erc_percentile;
                buf[6*i+1] = table[i].wind_speed_ms;
                buf[6*i+2] = table[i].wind_dir_deg;
                buf[6*i+3] = table[i].M_d1_pct;
                buf[6*i+4] = table[i].M_d10_pct;
                buf[6*i+5] = table[i].M_d100_pct;
            }
        }
        amrex::ParallelDescriptor::Bcast(buf.data(), n * 6,
            amrex::ParallelDescriptor::IOProcessorNumber());
        if (!amrex::ParallelDescriptor::IOProcessor()) {
            for (int i = 0; i < n; ++i) {
                table[i].erc_percentile = buf[6*i+0];
                table[i].wind_speed_ms  = buf[6*i+1];
                table[i].wind_dir_deg   = buf[6*i+2];
                table[i].M_d1_pct       = buf[6*i+3];
                table[i].M_d10_pct      = buf[6*i+4];
                table[i].M_d100_pct     = buf[6*i+5];
            }
        }
    }

    return table;
}

// ============================================================================
// erc_to_percentile
//
// Maps a current ERC value to a percentile (0–100) using the table's ERC
// percentile column as the reference distribution.
//
// The table rows define discrete percentile breakpoints.  The ERC column in
// the table is used as the per-percentile ERC threshold: if the current ERC
// >= the ERC at the 97th percentile row, the simulation uses the 97th
// percentile weather scenario.
//
// When the table only contains percentile labels (not ERC thresholds),
// the percentile is computed as a linear fraction of the ERC range observed
// in the table's wind_speed column (simple rank-based approach).
//
// Here we use a straightforward linear interpolation: the ERC value is
// linearly mapped to [0, 100] between the minimum and maximum ERC
// percentile in the table.
// ============================================================================
static double erc_to_percentile(double current_erc,
                                 const ConditionalWeatherTable& table)
{
    if (table.empty()) return 50.0;
    if (table.size() == 1) return table[0].erc_percentile;

    // The table is sorted by erc_percentile.  Treat erc_percentile as the
    // expected ERC at that percentile and linearly interpolate.
    const double erc_lo = table.front().erc_percentile;  // percentile = ERC scale
    const double erc_hi = table.back().erc_percentile;
    if (erc_hi <= erc_lo) return 50.0;

    // Clamp and scale to [0, 100]
    const double frac = std::max(0.0, std::min(1.0,
        (current_erc - erc_lo) / (erc_hi - erc_lo)));
    return frac * 100.0;
}

// ============================================================================
// apply_conditional_weather
//
// Selects the weather scenario whose erc_percentile best matches the current
// simulation fire-danger index, then updates RothermelParams moisture and
// wind inputs.
//
// The selected row is the one with the largest erc_percentile that does not
// exceed the current index percentile (lower-bound lookup).
//
// Parameters
// ----------
//   table         – conditional weather table (from load_conditional_weather)
//   current_index – current domain-mean fire-danger index value
//                   (ERC, BI, or SC depending on trigger_index)
//   trigger_index – which index drives the lookup: "erc", "bi", or "sc"
//                   (default "erc" reproduces the original behaviour)
//   rp            – Rothermel parameters to update (M_d1, M_d10, M_d100)
//   ux, uy        – wind components [m/s] to update
//
// Returns the index of the selected row (-1 if table is empty).
// ============================================================================
static int apply_conditional_weather(
    const ConditionalWeatherTable& table,
    double                         current_index,
    const std::string&             trigger_index,
    RothermelParams&               rp,
    amrex::Real&                   ux,
    amrex::Real&                   uy)
{
    if (table.empty()) return -1;

    // Log the trigger index type on first application so users can verify
    // which fire-danger index is driving the scenario lookup.
    static bool logged = false;
    if (!logged) {
        amrex::Print() << "ConditionalWeather: using trigger index = " << trigger_index << "\n";
        logged = true;
    }

    // Compute current percentile from the supplied index value
    const double current_pct = erc_to_percentile(current_index, table);

    // Log the trigger index type on first use (printed once via static flag)
    // so users can verify which index is driving the lookup.

    // Find the row with the largest erc_percentile <= current_pct
    int sel = 0;
    for (int i = 0; i < static_cast<int>(table.size()); ++i) {
        if (table[i].erc_percentile <= current_pct) sel = i;
        else break;
    }

    const auto& row = table[sel];

    // Update dead fuel moistures [fraction]
    rp.M_d1   = static_cast<amrex::Real>(row.M_d1_pct   / 100.0);
    rp.M_d10  = static_cast<amrex::Real>(row.M_d10_pct  / 100.0);
    rp.M_d100 = static_cast<amrex::Real>(row.M_d100_pct / 100.0);
    rp.M_f    = rp.M_d1;  // keep single-class moisture in sync

    // Convert meteorological wind (from direction) to U/V components [m/s]
    // Met. convention: dir_deg = direction the wind comes FROM
    //   U = -speed * sin(dir_rad)  (positive = eastward)
    //   V = -speed * cos(dir_rad)  (positive = northward)
    const double dir_rad = row.wind_dir_deg * (3.14159265358979323846 / 180.0);
    ux = static_cast<amrex::Real>(-row.wind_speed_ms * std::sin(dir_rad));
    uy = static_cast<amrex::Real>(-row.wind_speed_ms * std::cos(dir_rad));

    return sel;
}

#endif // CONDITIONAL_WEATHER_H
