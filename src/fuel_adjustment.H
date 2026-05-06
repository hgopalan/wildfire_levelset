// ============================================================================
// FARSITE Fuel Adjustment File (.adj) Reader
// ============================================================================
// A FARSITE fuel adjustment file specifies a per-fuel-model rate-of-spread
// multiplier that scales R₀ (the no-wind, no-slope ROS) for each fuel model
// in the landscape.  This allows regional calibration without modifying the
// standard fuel database.
//
// File format (whitespace-delimited, # or ! lines are comments):
//
//   [optional: <num_fuel_models>]
//   <fuel_model_number>  <adj_factor>
//   ...
//
// Examples:
//   1  1.0       (no change for model 1)
//   4  1.35      (35% increase for chaparral)
//   10 0.80      (20% reduction for timber/understory)
//
// The adj_factor is applied as:
//   R_adjusted = R_original × adj_factor
//
// Usage:
//   auto adjs = parse_fuel_adjustment_file("fire.adj");
//   apply_fuel_adjustment_to_table(h_fuel_table, adjs);
//
// References:
//   Finney, M.A. (2004). FARSITE: Fire Area Simulator — Model Development
//     and Evaluation. USDA Forest Service Research Paper RMRS-RP-4.
// ============================================================================

#ifndef FUEL_ADJUSTMENT_H
#define FUEL_ADJUSTMENT_H

#include "rothermel_model.H"   // RothermelComputed

#include <AMReX.H>

#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

// Per-fuel-model ROS adjustment entry
struct FuelAdjustment {
    int   fuel_model;    // FBFM13 or FBFM40 fuel model number
    float adj_factor;    // multiplicative ROS scale (default 1.0)
};

// ============================================================================
// parse_fuel_adjustment_file
//
// Reads a FARSITE .adj file and returns one FuelAdjustment per valid entry.
// Non-positive adjustment factors are clamped to 0 and a warning is printed.
// ============================================================================
inline std::vector<FuelAdjustment>
parse_fuel_adjustment_file(const std::string& filename)
{
    std::ifstream f(filename);
    if (!f.is_open())
        amrex::Abort("Cannot open fuel adjustment file: " + filename);

    std::vector<FuelAdjustment> result;
    std::string line;
    int line_num = 0;
    int declared_count = -1;   // from optional first numeric line
    bool count_line_consumed = false;

    while (std::getline(f, line)) {
        ++line_num;
        // Strip leading whitespace
        size_t first = line.find_first_not_of(" \t\r\n");
        if (first == std::string::npos) continue;
        line = line.substr(first);
        // Skip comments
        if (line[0] == '#' || line[0] == '!') continue;

        std::istringstream iss(line);
        int   a;
        float b;

        if (!(iss >> a)) continue;  // cannot parse

        if (!(iss >> b)) {
            // Only one token: treat as the optional count line
            if (!count_line_consumed) {
                declared_count = a;
                count_line_consumed = true;
            }
            continue;
        }

        // Validate
        if (b < 0.0f) {
            amrex::Print() << "WARNING fuel_adjustment: negative adj_factor "
                           << b << " for fuel model " << a
                           << " (line " << line_num << ") — clamped to 0\n";
            b = 0.0f;
        }
        result.push_back({a, b});
    }

    amrex::Print() << "Parsed " << result.size()
                   << " fuel adjustment entries from " << filename << "\n";
    if (declared_count > 0 && static_cast<int>(result.size()) != declared_count) {
        amrex::Print() << "  (file declared " << declared_count
                       << " entries; found " << result.size() << ")\n";
    }
    return result;
}

// ============================================================================
// apply_fuel_adjustment_to_table
//
// Scales R0 in every element of a RothermelComputed lookup table whose
// 1-based index matches a fuel model number in the adjustment list.
//
// The table is indexed as: table[fuel_model_code] where fuel_model_code is the
// integer returned by build_fuel_rothermel_table (codes 0..N-1; code 0 is
// reserved for non-burnable; actual FBFM numbers are stored as the index).
//
// If the table size is smaller than the highest adjusted fuel model number,
// the extra entries are silently ignored.
// ============================================================================
inline void apply_fuel_adjustment_to_table(
    std::vector<RothermelComputed>& table,
    const std::vector<FuelAdjustment>& adjs)
{
    if (table.empty() || adjs.empty()) return;

    // Build a fast lookup map: fuel_model_code → adj_factor
    std::map<int, float> adj_map;
    for (const auto& a : adjs)
        adj_map[a.fuel_model] = a.adj_factor;

    int applied = 0;
    for (const auto& kv : adj_map) {
        int code = kv.first;
        if (code < 0 || code >= static_cast<int>(table.size())) continue;
        float f = kv.second;
        // Scale R0 (no-wind, no-slope ROS)
        table[code].R0 *= static_cast<amrex::Real>(f);
        ++applied;
    }
    amrex::Print() << "Applied " << applied
                   << " fuel adjustments to Rothermel lookup table\n";
}

// ============================================================================
// apply_fuel_adjustment_to_params
//
// Convenience overload for when there is no per-cell fuel table (single global
// fuel model).  Looks up the adjustment for the current global fuel model code
// and scales rp.w0 (and hence R0 at run time) by the factor.
//
// fuel_model_code should be the integer code used in the landscape, or 0 to
// apply the first adjustment that has a factor != 1.
// ============================================================================
inline void apply_fuel_adjustment_to_params(
    RothermelParams& rp,
    const std::vector<FuelAdjustment>& adjs,
    int fuel_model_code)
{
    if (adjs.empty()) return;
    for (const auto& a : adjs) {
        if (a.fuel_model == fuel_model_code) {
            // Scale the oven-dry total fuel load – this propagates through
            // the Rothermel reaction intensity and R0 consistently.
            rp.w0    *= static_cast<amrex::Real>(a.adj_factor);
            rp.w_d1  *= static_cast<amrex::Real>(a.adj_factor);
            rp.w_d10 *= static_cast<amrex::Real>(a.adj_factor);
            rp.w_d100*= static_cast<amrex::Real>(a.adj_factor);
            rp.w_lh  *= static_cast<amrex::Real>(a.adj_factor);
            rp.w_lw  *= static_cast<amrex::Real>(a.adj_factor);
            amrex::Print() << "Applied fuel adjustment " << a.adj_factor
                           << " to global RothermelParams (fuel model "
                           << fuel_model_code << ")\n";
            return;
        }
    }
}

#endif // FUEL_ADJUSTMENT_H
