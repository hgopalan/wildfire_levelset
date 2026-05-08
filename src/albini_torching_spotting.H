// ============================================================================
// Albini (1979) Torching-Tree Spotting Model
// ============================================================================
// Models long-range firebrand spotting originating from torching trees — trees
// whose crowns are ignited by the passing fire front.  This is distinct from
// the existing Albini (1983) surface-fire spotting model: the 1979 model uses
// flame height above the canopy top (driven by fireline intensity and stand
// height) to estimate the firebrand lofting height and maximum spotting
// distance.
//
// Physical model (simplified from Albini 1979 / Rothermel 1991 / FARSITE):
//
//   1.  Torching condition: crown_activity >= min_crown_activity AND
//       fireline_intensity >= I_B_min [kW/m].
//
//   2.  Effective lofting height:
//         H_f  [m] = flame_length  (Byram 1959 from fireline_intensity_mf)
//         H_eff [m] = max(H_f, H_canopy)          [m]
//       When flame penetrates above the canopy top, the firebrand reaches a
//       greater height and thus a larger spotting distance.
//
//   3.  Maximum spotting distance (simplified Albini 1979):
//         d_spot [m] = k_torch × U [m/s] × sqrt(H_eff [m])
//       The default k_torch = 4.24 ≈ sqrt(2/g) × (empirical scaling) is
//       consistent with a ballistic trajectory: a particle rising to H_eff
//       and falling at terminal velocity while carried by wind speed U.
//
//   4.  Spot ignition: with probability P_base, a new fire is placed at
//       distance d_spot downwind (plus ±45° random lateral scatter), within
//       a Gaussian spot of radius spot_radius [m].
//
// References:
//   Albini, F.A. (1979). "Spot fire distance from burning trees — a predictive
//       model." USDA Forest Service General Technical Report INT-56.
//   Rothermel, R.C. (1991). "Predicting behavior and size of crown fires in
//       the northern Rocky Mountains." USDA Forest Service Research Paper INT-438.
//   Finney, M.A. (1998). FARSITE: Fire Area Simulator. RMRS-RP-4.
// ============================================================================

#ifndef ALBINI_TORCHING_SPOTTING_H
#define ALBINI_TORCHING_SPOTTING_H

#include "constants.H"
#include "parse_inputs.H"
#include <AMReX.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <cmath>
#include <random>
#include <vector>

using namespace amrex;

// ============================================================================
// compute_albini_torching_spots
//
// Identifies torching-tree cells, computes their spotting distances, and
// (stochastically) places new spot-fire ignitions in phi.
//
// Parameters:
//   phi                    – level-set field (modified in place: spotted cells
//                            are set to -1)
//   ecology_mf             – fire ecology MultiFab; component 3 = crown_activity
//                            (0=surface, 1=passive, 2=active)
//   fireline_intensity_mf  – Byram fireline intensity [kW/m]
//   flame_length_mf        – Byram flame length [m]; used as H_f
//   canopy_height_mf       – stand (canopy) height [m]; used for H_eff clamp
//   vel                    – wind velocity [m/s]; used for wind direction and
//                            magnitude
//   geom                   – AMReX geometry (cell sizes, domain lo)
//   tp                     – torching-spotting parameters
//   step                   – current time step (used to seed RNG)
// ============================================================================
inline void compute_albini_torching_spots(
    MultiFab&       phi,
    const MultiFab& ecology_mf,
    const MultiFab& fireline_intensity_mf,
    const MultiFab& flame_length_mf,
    const MultiFab& canopy_height_mf,
    const MultiFab& vel,
    const Geometry& geom,
    const InputParameters::TorchingSpottingParams& tp,
    int step)
{
    if (tp.enable != 1) return;

    const auto dx     = geom.CellSize();
    const Real *prob_lo = geom.ProbLo();
#if (AMREX_SPACEDIM == 3)
    const GpuArray<Real,3> gdx{dx[0], dx[1], dx[2]};
#else
    const GpuArray<Real,2> gdx{dx[0], dx[1]};
#endif

    const Real k_torch      = tp.k_torch;
    const Real I_B_min      = tp.I_B_min;
    const Real spot_radius  = tp.spot_radius;
    const Real P_base       = tp.P_base;
    const int  min_ca       = tp.min_crown_activity;

    // RNG: seed from step + torching seed.
    // Multipliers (1000, 31, 7) spread seeds across different Lehmer multiplier
    // sub-sequences so this model does not produce correlated draws with the
    // Albini 1983 and stochastic spotting models that use step-seeded RNGs.
    std::mt19937 rng(static_cast<unsigned>(tp.random_seed * 1000 + step * 31 + 7));
    std::uniform_real_distribution<Real> uniform_01(Real(0.0), Real(1.0));
    std::uniform_real_distribution<Real> lateral_angle(-Real(0.785398), Real(0.785398)); // ±45°

    // Collect torching cells in serial (CPU-only, matches existing spotting pattern)
    struct TorchingCell { Real x, y, d_spot, wind_dir; };
    std::vector<TorchingCell> cells;

    for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
        if (!ParallelDescriptor::IOProcessor() &&
            ParallelDescriptor::NProcs() > 1) continue; // serial on IO rank only

        const Box& bx = mfi.validbox();
        auto const p   = phi.const_array(mfi);
        auto const eco = ecology_mf.const_array(mfi);
        auto const ib  = fireline_intensity_mf.const_array(mfi);
        auto const fl  = flame_length_mf.const_array(mfi);
        auto const cht = canopy_height_mf.const_array(mfi);
        auto const v   = vel.const_array(mfi);

        amrex::LoopOnCpu(bx, [&](int i, int j, int k) {
            // Only consider cells near or inside fire (phi <= 0)
            if (p(i, j, k) > Real(0.1)) return;

            // Crown activity threshold
            const Real ca = eco(i, j, k, 3);
            if (ca < static_cast<Real>(min_ca) - Real(0.5)) return;

            // Fireline intensity threshold
            const Real IB = ib(i, j, k);
            if (IB < I_B_min) return;

            // Stochastic gate
            if (uniform_01(rng) > P_base) return;

            // Effective lofting height [m]: max(flame_length, canopy_height)
            const Real H_f      = fl(i, j, k);
            const Real H_canopy = cht(i, j, k);
            const Real H_eff    = amrex::max(H_f, amrex::max(H_canopy, Real(1.0)));

            // Wind speed and direction
            const Real ux = v(i, j, k, 0);
            const Real uy = v(i, j, k, 1);
            const Real U  = std::sqrt(ux*ux + uy*uy);
            if (U < Real(1.0e-6)) return;  // no wind → no directional spotting

            // d_spot [m] = k_torch × U [m/s] × sqrt(H_eff [m])
            const Real d_spot = k_torch * U * std::sqrt(H_eff);
            if (d_spot < Real(1.0e-3)) return;

            // Wind direction angle [rad]
            const Real wind_dir = std::atan2(uy, ux);

            // Cell centre position
            const Real cx = prob_lo[0] + (i + Real(0.5)) * gdx[0];
            const Real cy = prob_lo[1] + (j + Real(0.5)) * gdx[1];

            cells.push_back({cx, cy, d_spot, wind_dir});
        });
    }

    if (cells.empty()) return;

    // Place spot fires
    for (const auto& c : cells) {
        // Random lateral scatter in ±45°
        const Real scatter = lateral_angle(rng);
        const Real angle   = c.wind_dir + scatter;
        const Real x_spot  = c.x + c.d_spot * std::cos(angle);
        const Real y_spot  = c.y + c.d_spot * std::sin(angle);

        // Stamp a Gaussian spot fire in phi
        const int radius_cells = static_cast<int>(std::ceil(spot_radius / amrex::min(dx[0], dx[1])));

        for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            auto p_arr = phi.array(mfi);
            const int i0 = static_cast<int>((x_spot - prob_lo[0]) / dx[0] - Real(0.5));
            const int j0 = static_cast<int>((y_spot - prob_lo[1]) / dx[1] - Real(0.5));
            const int ilo = amrex::max(i0 - radius_cells, bx.smallEnd(0));
            const int ihi = amrex::min(i0 + radius_cells, bx.bigEnd(0));
            const int jlo = amrex::max(j0 - radius_cells, bx.smallEnd(1));
            const int jhi = amrex::min(j0 + radius_cells, bx.bigEnd(1));
#if (AMREX_SPACEDIM == 3)
            const int kfix = (bx.smallEnd(2) + bx.bigEnd(2)) / 2;
#else
            const int kfix = 0;
#endif
            const Real sig2 = spot_radius * spot_radius;
            for (int ii = ilo; ii <= ihi; ++ii) {
                for (int jj = jlo; jj <= jhi; ++jj) {
                    if (!bx.contains(IntVect(AMREX_D_DECL(ii, jj, kfix)))) continue;
                    if (p_arr(ii, jj, kfix) < Real(0.0)) continue; // already burned
                    const Real px = prob_lo[0] + (ii + Real(0.5)) * dx[0];
                    const Real py = prob_lo[1] + (jj + Real(0.5)) * dx[1];
                    const Real r2 = (px - x_spot)*(px - x_spot) + (py - y_spot)*(py - y_spot);
                    if (r2 <= sig2) {
                        p_arr(ii, jj, kfix) = Real(-1.0);
                    }
                }
            }
        }
    }

    if (!cells.empty()) {
        amrex::Print() << "  Torching-tree spotting: " << cells.size()
                       << " spot(s) placed this step\n";
    }
}

#endif // ALBINI_TORCHING_SPOTTING_H
