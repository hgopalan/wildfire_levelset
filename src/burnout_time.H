// ============================================================================
// Burnout Time Diagnostic
// ============================================================================
// Computes the per-cell post-frontal burnout time [s]: the simulation time at
// which active combustion ceases in each burned cell.
//
// Formula:
//   burnout_time [s] = arrival_time [s] + tau_residence [s]
//
// where
//   arrival_time  = simulation time when the cell first ignited (phi < 0)
//                   (-1.0 sentinel for cells not yet burned)
//   tau_residence = fire residence time [s] (farsite.tau_residence)
//
// Cells that have not yet ignited (arrival_time < 0) keep the sentinel value
// -1.0.  Cells currently undergoing active combustion have
// burnout_time > current_time; fully burned-out cells have
// burnout_time <= current_time.
//
// This mirrors the "burnout time" raster written by FARSITE (Finney 1998)
// and is useful for computing the duration of active combustion per cell and
// for visualising the post-frontal burn pattern.
//
// References:
//   Finney, M.A. (1998). FARSITE: Fire Area Simulator — Model Development
//     and Evaluation. USDA Forest Service Research Paper RMRS-RP-4.
//   Rothermel, R.C. (1983). USDA Forest Service GTR INT-143 (residence time).
// ============================================================================

#ifndef BURNOUT_TIME_H
#define BURNOUT_TIME_H

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_GpuLaunch.H>

using namespace amrex;

// ============================================================================
// compute_burnout_time
//
// Parameters:
//   burnout_time_mf  – output: 1-component MultiFab [s]; -1.0 = not yet burned
//   arrival_time_mf  – input:  simulation time of ignition [s]; -1.0 = unburned
//   tau_residence    – fire residence time [s] (farsite.tau_residence)
// ============================================================================
inline void compute_burnout_time(MultiFab&       burnout_time_mf,
                                  const MultiFab& arrival_time_mf,
                                  amrex::Real     tau_residence)
{
    for (MFIter mfi(burnout_time_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const at  = arrival_time_mf.const_array(mfi);
        auto       bot = burnout_time_mf.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const Real t_arrive = at(i, j, k);
            // Keep sentinel -1.0 for cells that have not ignited yet
            bot(i, j, k) = (t_arrive >= Real(0.0))
                           ? t_arrive + tau_residence
                           : Real(-1.0);
        });
    }
}

#endif // BURNOUT_TIME_H
