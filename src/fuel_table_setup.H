// fuel_table_setup.H
// Builds all GPU-side lookup tables (Rothermel, burnout-time, Balbi) and
// pre-computes the Cruz (2005) crown fire coefficients.  These tables are
// computed once at startup from the landscape fuel map and fire-spread model
// settings.
#ifndef FUEL_TABLE_SETUP_H
#define FUEL_TABLE_SETUP_H

// All required headers provided by wildfire_includes.H (the umbrella header).
// fuel_table_setup.H additionally requires AMReX_GpuContainers.H for
// Gpu::DeviceVector and Gpu::copy.
#include <AMReX_GpuContainers.H>

using namespace amrex;

// All GPU lookup tables and pre-computed model coefficients produced at startup.
struct FuelTableData {
    // Rothermel per-cell lookup table
    Gpu::DeviceVector<RothermelComputed> d_fuel_table;
    const RothermelComputed* d_fuel_table_ptr = nullptr;
    int fuel_table_size = 0;

    // Per-fuel burnout (residence) time table
    Gpu::DeviceVector<FuelResidenceTime> d_tau_table;
    const FuelResidenceTime* d_tau_table_ptr = nullptr;
    int tau_table_size = 0;

    // Balbi per-cell lookup table
    Gpu::DeviceVector<BalbiComputed> d_balbi_table;
    const BalbiComputed* d_balbi_table_ptr = nullptr;
    int balbi_table_size = 0;
    BalbiComputed bc_global_default;

    // Cruz, Alexander & Wakimoto (2005) crown fire coefficients
    CruzCrownComputed ccc_global;
    const CruzCrownComputed* ccc_ptr = nullptr;

    FuelTableData() { bc_global_default.A_coeff = 0.0; bc_global_default.v_b = 1.0; }
};

// Build all lookup tables and pre-computed coefficients.
// inputs may be modified (fuel adjustment applied to global params when no LCP).
inline FuelTableData setup_fuel_tables(InputParameters& inputs,
                                       const MultiFab& fuel_model_mf)
{
    FuelTableData ftd;

    // ---- Rothermel per-cell table ----
    if (!inputs.rothermel.landscape_file.empty()) {
        std::vector<RothermelComputed> h_table =
            build_fuel_rothermel_table(inputs.rothermel,
                                       inputs.rothermel.landscape_fuel_type);
        if (!inputs.fuel_adj_file.empty()) {
            auto adjs = parse_fuel_adjustment_file(inputs.fuel_adj_file);
            apply_fuel_adjustment_to_table(h_table, adjs);
        }
        ftd.fuel_table_size = static_cast<int>(h_table.size());
        ftd.d_fuel_table.resize(ftd.fuel_table_size);
        Gpu::copy(Gpu::hostToDevice,
                  h_table.begin(), h_table.end(), ftd.d_fuel_table.begin());
        ftd.d_fuel_table_ptr = ftd.d_fuel_table.data();
        amrex::Print() << "Built per-cell Rothermel lookup table: "
                       << ftd.fuel_table_size << " entries ("
                       << (inputs.rothermel.landscape_fuel_type == "40"
                           ? "FBFM40" : "FBFM13") << ")\n";
    } else if (!inputs.fuel_adj_file.empty()) {
        auto adjs = parse_fuel_adjustment_file(inputs.fuel_adj_file);
        apply_fuel_adjustment_to_params(inputs.rothermel, adjs,
                                        inputs.fuel_adj_model);
    }

    // ---- Per-fuel burnout (residence) time table ----
    if (ftd.d_fuel_table_ptr != nullptr && ftd.fuel_table_size > 0) {
        std::vector<FuelResidenceTime> h_tau =
            build_fuel_tau_table(ftd.fuel_table_size);
        ftd.tau_table_size = static_cast<int>(h_tau.size());
        ftd.d_tau_table.resize(ftd.tau_table_size);
        Gpu::copy(Gpu::hostToDevice,
                  h_tau.begin(), h_tau.end(), ftd.d_tau_table.begin());
        ftd.d_tau_table_ptr = ftd.d_tau_table.data();
        amrex::Print() << "Built per-fuel burnout time table: "
                       << ftd.tau_table_size << " entries\n";
    }

    // ---- Cruz crown coefficients ----
    ftd.ccc_global = compute_cruz_crown_params(inputs.cruz_crown);
    ftd.ccc_ptr    = (inputs.crown.use_cruz_crown == 1) ? &ftd.ccc_global : nullptr;

    // ---- Balbi per-cell table ----
    if (inputs.fire_spread_model == "balbi") {
        print_balbi_fuel_table(inputs.rothermel, inputs.balbi,
                               inputs.rothermel.landscape_fuel_type);
        ftd.bc_global_default =
            compute_balbi_params(inputs.rothermel, inputs.balbi);
        if (!inputs.rothermel.landscape_file.empty()) {
            std::vector<BalbiComputed> h_balbi =
                build_fuel_balbi_table(inputs.rothermel, inputs.balbi,
                                       inputs.rothermel.landscape_fuel_type);
            ftd.balbi_table_size = static_cast<int>(h_balbi.size());
            ftd.d_balbi_table.resize(ftd.balbi_table_size);
            Gpu::copy(Gpu::hostToDevice,
                      h_balbi.begin(), h_balbi.end(),
                      ftd.d_balbi_table.begin());
            ftd.d_balbi_table_ptr = ftd.d_balbi_table.data();
            amrex::Print() << "Built per-cell Balbi lookup table: "
                           << ftd.balbi_table_size << " entries ("
                           << (inputs.rothermel.landscape_fuel_type == "40"
                               ? "FBFM40" : "FBFM13") << ")\n";
        }
    }

    return ftd;
}

#endif // FUEL_TABLE_SETUP_H
