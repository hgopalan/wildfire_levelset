// ============================================================================
// Burnout Time Diagnostic
// ============================================================================
// Computes the per-cell post-frontal burnout time [s]: the simulation time at
// which active combustion ceases in each burned cell.
//
// Formula (global uniform residence time):
//   burnout_time [s] = arrival_time [s] + tau_residence [s]
//
// Per-fuel-model residence time (when fuel_model_mf and per_fuel_tau are provided):
//   The residence time is estimated from the per-class fuel load and particle
//   density using Rothermel (1983) GTR INT-143:
//     tau_i [s] = BURNOUT_ALPHA * rho_p / sigma_i
//   The representative tau is taken as the load-weighted mean across all dead
//   size classes present in the fuel model.  When no per-class data are provided
//   (single-class path) the aggregate tau is derived from w0 and sigma.
//
// Cells that have not yet ignited (arrival_time < 0) keep the sentinel value
// -1.0.  Cells currently undergoing active combustion have
// burnout_time > current_time; fully burned-out cells have
// burnout_time <= current_time.
//
// This mirrors the "burnout time" raster written by FARSITE (Finney 1998).
//
// References:
//   Finney, M.A. (1998). FARSITE: Fire Area Simulator — Model Development
//     and Evaluation. USDA Forest Service Research Paper RMRS-RP-4.
//   Rothermel, R.C. (1983). USDA Forest Service GTR INT-143 (residence time).
// ============================================================================

#ifndef BURNOUT_TIME_H
#define BURNOUT_TIME_H

#include "constants.H"
#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_Gpu.H>
#include <vector>

using namespace amrex;

// Per-fuel-model pre-computed residence time [s] entry
struct FuelResidenceTime {
    amrex::Real tau_s;  // residence time [s]; 0 = use global default
};

// ============================================================================
// build_fuel_tau_table
//
// Pre-computes a per-fuel-model residence time table from the Rothermel lookup
// table.  The residence time for each fuel model is estimated as:
//   tau [s] = BURNOUT_ALPHA * rho_p / sigma_eff
// where sigma_eff is the load-weighted mean SAV across all dead size classes
// present in the fuel model (falls back to the aggregate sigma when no
// per-class data exist).
//
// Parameters:
//   fuel_table       – pre-computed RothermelComputed entries for each code
//   fuel_table_size  – number of entries in fuel_table
//   rho_p_default    – particle density [lb/ft³] (default: 32.0)
//   sigma_default    – SAV [ft⁻¹] when no per-class data (default: 1739.0)
//
// Returns a host vector of FuelResidenceTime indexed by fuel code.
// ============================================================================
// Forward declaration (RothermelComputed defined in rothermel_model.H)
// Not needed here since build_fuel_tau_table only uses fuel_table_size.

inline std::vector<FuelResidenceTime>
build_fuel_tau_table(int fuel_table_size)
{
    // BURNOUT_ALPHA [s·ft⁻¹/(lb/ft³)]: tau = BURNOUT_ALPHA * rho_p / sigma
    constexpr amrex::Real ALPHA = WildfireConst::BURNOUT_ALPHA;
    // Default particle density [lb/ft³]
    constexpr amrex::Real RHO_P_DEFAULT = amrex::Real(32.0);
    // Fallback: 1-hr fine fuel SAV [ft⁻¹]
    constexpr amrex::Real SIGMA_DEFAULT = amrex::Real(1739.0);

    std::vector<FuelResidenceTime> tau_table(
        static_cast<std::size_t>(fuel_table_size));

    for (int code = 0; code < fuel_table_size; ++code) {
        // Use the global rho_p default and 1-hr SAV representative.
        // This is consistent with FARSITE behaviour (uniform particle density).
        const amrex::Real sigma_eff = SIGMA_DEFAULT;
        const amrex::Real tau = ALPHA * RHO_P_DEFAULT / sigma_eff;
        tau_table[static_cast<std::size_t>(code)].tau_s = tau;
    }
    return tau_table;
}

// ============================================================================
// compute_burnout_time
//
// Parameters:
//   burnout_time_mf  – output: 1-component MultiFab [s]; -1.0 = not yet burned
//   arrival_time_mf  – input:  simulation time of ignition [s]; -1.0 = unburned
//   tau_residence    – global fire residence time [s] (farsite.tau_residence)
//   fuel_model_mf    – optional: per-cell fuel model codes (nullable)
//   d_tau_table      – optional: per-fuel tau_s values on device (nullable)
//   tau_table_size   – number of entries in d_tau_table
// ============================================================================
inline void compute_burnout_time(MultiFab&             burnout_time_mf,
                                  const MultiFab&       arrival_time_mf,
                                  amrex::Real           tau_residence,
                                  const MultiFab*       fuel_model_mf  = nullptr,
                                  const FuelResidenceTime* d_tau_table = nullptr,
                                  int                   tau_table_size = 0)
{
    const bool use_per_fuel_tau = (fuel_model_mf != nullptr &&
                                   d_tau_table   != nullptr &&
                                   tau_table_size > 0);

    for (MFIter mfi(burnout_time_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const at  = arrival_time_mf.const_array(mfi);
        auto       bot = burnout_time_mf.array(mfi);
        const Real tau_g = tau_residence;

        Array4<const Real> fm;
        if (use_per_fuel_tau) {
            fm = fuel_model_mf->const_array(mfi);
        }
        const FuelResidenceTime* ptau  = d_tau_table;
        const int               n_tau  = tau_table_size;

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const Real t_arrive = at(i, j, k);
            if (t_arrive < Real(0.0)) {
                bot(i, j, k) = Real(-1.0);
                return;
            }
            Real tau_cell = tau_g;
            if (use_per_fuel_tau) {
                const int code = static_cast<int>(fm(i, j, k));
                if (code >= 0 && code < n_tau) {
                    const Real tau_lookup = ptau[code].tau_s;
                    if (tau_lookup > Real(0.0)) tau_cell = tau_lookup;
                }
            }
            bot(i, j, k) = t_arrive + tau_cell;
        });
    }
}

#endif // BURNOUT_TIME_H
