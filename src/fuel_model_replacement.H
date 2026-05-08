// ============================================================================
// Post-Fire Fuel Model Replacement
// ============================================================================
// After the fire front passes through a cell (phi flips from positive to
// negative), the per-cell fuel model code in fuel_model_mf is replaced with
// a post-fire code from a user-supplied lookup table.
//
// Physical motivation:
//   Burned vegetation has fundamentally different fuel properties than unburned
//   fuel.  FARSITE supports a post-fire fuel model substitution table so that
//   re-entry burns (e.g., spot fires re-burning partially consumed fuel) and
//   subsequent fire spread through burned areas use appropriate fuel parameters.
//   A classic example: Anderson Model 4 (chaparral) → Model 2 (moderate fuel)
//   after the initial burn.
//
// CSV file format (whitespace or comma delimited; # or ! comment lines):
//   original_code   replacement_code
//   4               2
//   5               1
//   10              3
//
// Usage:
//   At each time step, newly burned cells (phi < 0 AND arrival_time set for the
//   first time at this step) have their fuel model code updated.
//
// References:
//   Finney, M.A. (1998). FARSITE: Fire Area Simulator. USDA FS RMRS-RP-4.
//   Scott, J.H. (2006). Comparison of Crown Fire Modelling Systems Used in
//     Three Fire Management Applications. USDA FS RMRS-RP-58.
// ============================================================================

#ifndef FUEL_MODEL_REPLACEMENT_H
#define FUEL_MODEL_REPLACEMENT_H

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_Gpu.H>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>

using namespace amrex;

// ============================================================================
// FuelModelReplacementTable
// Stores original → replacement code mappings for use on the GPU.
// ============================================================================
struct FuelModelReplacementTable {
    // Parallel arrays: original[i] → replacement[i]
    amrex::Gpu::DeviceVector<int> d_original;
    amrex::Gpu::DeviceVector<int> d_replacement;
    int n_entries = 0;

    bool empty() const { return n_entries == 0; }
};

// ============================================================================
// load_post_fire_fuel_map
//
// Parses the post-fire fuel model CSV and builds a GPU-accessible lookup table.
//
// Parameters
// ----------
//   filename  – path to the post-fire fuel map CSV
//
// Returns a FuelModelReplacementTable; empty table when file is empty/missing.
// ============================================================================
static FuelModelReplacementTable load_post_fire_fuel_map(const std::string& filename)
{
    FuelModelReplacementTable table;
    if (filename.empty()) return table;

    std::vector<int> h_original, h_replacement;

    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ifstream fin(filename);
        if (!fin.is_open()) {
            amrex::Abort("post_fire_fuel_map: cannot open file '" + filename + "'");
        }
        std::string line;
        int n_loaded = 0;
        while (std::getline(fin, line)) {
            if (line.empty() || line[0] == '#' || line[0] == '!') continue;
            for (char& c : line) if (c == ',') c = ' ';
            std::istringstream iss(line);
            int orig, repl;
            if (!(iss >> orig >> repl)) continue;
            if (orig == repl) {
                amrex::Print() << "  post_fire_fuel_map: skipping identity mapping "
                               << orig << " → " << repl << "\n";
                continue;
            }
            h_original.push_back(orig);
            h_replacement.push_back(repl);
            ++n_loaded;
        }
        amrex::Print() << "post_fire_fuel_map: loaded "
                       << n_loaded << " code replacement(s) from '"
                       << filename << "'\n";
    }

    // Broadcast to all MPI ranks
    int n = static_cast<int>(h_original.size());
    amrex::ParallelDescriptor::Bcast(&n, 1, amrex::ParallelDescriptor::IOProcessorNumber());
    h_original.resize(n);
    h_replacement.resize(n);
    if (n > 0) {
        amrex::ParallelDescriptor::Bcast(h_original.data(), n,
            amrex::ParallelDescriptor::IOProcessorNumber());
        amrex::ParallelDescriptor::Bcast(h_replacement.data(), n,
            amrex::ParallelDescriptor::IOProcessorNumber());
    }

    table.n_entries = n;
    if (n > 0) {
        table.d_original.resize(n);
        table.d_replacement.resize(n);
        amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                         h_original.begin(), h_original.end(),
                         table.d_original.begin());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                         h_replacement.begin(), h_replacement.end(),
                         table.d_replacement.begin());
    }

    return table;
}

// ============================================================================
// apply_post_fire_fuel_replacement
//
// Replaces fuel model codes in fuel_model_mf for newly burned cells.
// A cell is "newly burned" when phi < 0 AND arrival_time was just set at
// current_time (arrival_time is within one time-step of current_time).
//
// This function is called once per time step after phi has been updated.
//
// Parameters
// ----------
//   fuel_model_mf – per-cell fuel model code MultiFab (modified in-place)
//   phi           – level-set field (phi < 0 = burned)
//   arrival_time  – simulation time of first ignition per cell (-1 = unburned)
//   table         – replacement code lookup table
//   current_time  – current simulation time [s]
//   dt            – current time step size [s]
// ============================================================================
inline void apply_post_fire_fuel_replacement(
    MultiFab&                              fuel_model_mf,
    const MultiFab&                        phi,
    const MultiFab&                        arrival_time,
    const FuelModelReplacementTable&       table,
    amrex::Real                            current_time,
    amrex::Real                            dt)
{
    if (table.empty()) return;

    const int n_entries = table.n_entries;
    const int* d_orig   = table.d_original.data();
    const int* d_repl   = table.d_replacement.data();

    // A cell is newly burned when its arrival_time falls within the current step
    const Real t_lo = current_time - Real(1.01) * dt;
    const Real t_hi = current_time + Real(0.01) * dt;

    for (MFIter mfi(fuel_model_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto       fm  = fuel_model_mf.array(mfi);
        auto const p   = phi.const_array(mfi);
        auto const at  = arrival_time.const_array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            if (p(i, j, k) >= Real(0.0)) return;       // unburned
            const Real ar = at(i, j, k);
            if (ar < t_lo || ar > t_hi) return;        // not newly burned

            const int code = static_cast<int>(fm(i, j, k));
            for (int e = 0; e < n_entries; ++e) {
                if (d_orig[e] == code) {
                    fm(i, j, k) = static_cast<Real>(d_repl[e]);
                    return;
                }
            }
        });
    }
}

#endif // FUEL_MODEL_REPLACEMENT_H
