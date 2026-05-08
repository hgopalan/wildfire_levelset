// ============================================================================
// Aerial Retardant Suppression Model
// ============================================================================
// Models aerial retardant drops as time-scheduled polygon zones where the
// rate of spread is reduced by an effectiveness factor that decays
// exponentially back to ambient over hours.
//
// Physical model:
//   When a retardant drop covers a cell, the effective ROS multiplier is:
//     f_ret(t) = 1 - effectiveness * exp(-(t - t_drop) / tau_decay)
//   where:
//     effectiveness ∈ [0,1]  initial fractional reduction (0 = no effect,
//                            1 = full suppression)
//     tau_decay [s]          half-life time constant for decay (typical: 4-8 hrs)
//     t - t_drop [s]         age of the drop
//
//   ROS_effective = ROS_original * f_ret(t)
//
//   Simultaneously, the retardant raises dead fuel moisture, simulating the
//   wetting effect.  A simplified approach scales M_d1 and M_d10 by a moisture
//   boost factor that also decays over time:
//     M_boost(t) = moisture_boost * exp(-(t - t_drop) / tau_decay)
//   (applied externally through the rothermel moisture update; not done here)
//
// File format (whitespace or comma delimited; # or ! lines are comments):
//   DROP_TIME_S  EFFECTIVENESS  HALF_LIFE_S  N_VERTICES
//   X1  Y1
//   X2  Y2
//   ...                   (polygon vertices, one per line)
//   [next drop]
//
// Usage:
//   RetardantDropList drops;
//   load_retardant_drops("retardant.csv", drops);
//   // at each timestep:
//   apply_retardant_to_ros(R_mf, drops, geom, current_time);
//
// References:
//   Finney, M.A. (2004). FARSITE: Fire Area Simulator. USDA FS RMRS-RP-4.
//   National Interagency Fire Center (NIFC) aerial retardant effectiveness
//     guidelines (NFES 2106, 2011).
// ============================================================================

#ifndef RETARDANT_DROP_H
#define RETARDANT_DROP_H

#include "parse_inputs.H"
#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_Gpu.H>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace amrex;

// One retardant drop zone
struct RetardantDrop {
    double              drop_time_s;     // time when drop occurs [s]
    double              effectiveness;   // initial fractional ROS reduction [0-1]
    double              tau_decay_s;     // exponential decay time constant [s]
    std::vector<double> vx;              // polygon vertex X [m]
    std::vector<double> vy;              // polygon vertex Y [m]
};

// Container for all drops
struct RetardantDropList {
    std::vector<RetardantDrop> drops;
    bool empty() const { return drops.empty(); }
};

// ============================================================================
// load_retardant_drops
//
// Parse a retardant schedule file into a RetardantDropList.
// ============================================================================
inline void load_retardant_drops(const std::string& filename, RetardantDropList& list)
{
    list.drops.clear();
    std::ifstream f(filename);
    if (!f.is_open())
        amrex::Abort("RetardantDrop: cannot open file: " + filename);

    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#' || line[0] == '!') continue;
        for (char& c : line) if (c == ',') c = ' ';
        std::istringstream iss(line);
        double t_drop, eff, tau;
        int nv;
        if (!(iss >> t_drop >> eff >> tau >> nv)) continue;
        if (nv < 3) {
            amrex::Print() << "RetardantDrop: skipping drop with < 3 vertices\n";
            continue;
        }

        RetardantDrop drop;
        drop.drop_time_s   = t_drop;
        drop.effectiveness = std::max(0.0, std::min(1.0, eff));
        drop.tau_decay_s   = std::max(1.0, tau);
        drop.vx.resize(nv);
        drop.vy.resize(nv);
        bool ok = true;
        for (int vi = 0; vi < nv && ok; ++vi) {
            std::string vline;
            while (std::getline(f, vline)) {
                if (vline.empty() || vline[0] == '#' || vline[0] == '!') continue;
                for (char& c : vline) if (c == ',') c = ' ';
                break;
            }
            std::istringstream viss(vline);
            if (!(viss >> drop.vx[vi] >> drop.vy[vi])) {
                amrex::Print() << "RetardantDrop: failed to parse vertex " << vi
                               << " of drop at t=" << t_drop
                               << " s; skipping this drop.\n";
                ok = false;
                break;
            }
        }
        if (ok) list.drops.push_back(std::move(drop));
    }
    amrex::Print() << "RetardantDrop: loaded " << list.drops.size()
                   << " drop zones from " << filename << "\n";
}

// ============================================================================
// point_in_polygon_2d
//
// Simple winding-number test (GPU-safe).
// Returns true when (px, py) lies inside the polygon defined by (vx, vy).
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
static bool retardant_point_in_polygon(amrex::Real px, amrex::Real py,
                                       const amrex::Real* vx, const amrex::Real* vy,
                                       int n)
{
    int winding = 0;
    for (int i = 0; i < n; ++i) {
        int j = (i + 1) % n;
        if (vy[i] <= py) {
            if (vy[j] > py) {
                amrex::Real cross = (vx[j]-vx[i])*(py-vy[i]) - (px-vx[i])*(vy[j]-vy[i]);
                if (cross > amrex::Real(0.0)) ++winding;
            }
        } else {
            if (vy[j] <= py) {
                amrex::Real cross = (vx[j]-vx[i])*(py-vy[i]) - (px-vx[i])*(vy[j]-vy[i]);
                if (cross < amrex::Real(0.0)) --winding;
            }
        }
    }
    return (winding != 0);
}

// ============================================================================
// apply_retardant_to_ros
//
// For each active drop zone (drop has occurred, not yet fully decayed), scale
// R_mf by (1 - effectiveness * exp(-(current_time - drop_time) / tau_decay))
// inside the drop polygon.
//
// Parameters:
//   R_mf         – rate of spread field (modified in place)
//   drops        – retardant drop list (unchanged)
//   geom         – AMReX geometry
//   current_time – current simulation time [s]
// ============================================================================
inline void apply_retardant_to_ros(MultiFab&               R_mf,
                                   const RetardantDropList& drops,
                                   const Geometry&          geom,
                                   amrex::Real              current_time)
{
    if (drops.empty()) return;

    const auto dx  = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();

    for (const auto& drop : drops.drops) {
        const amrex::Real t_drop = static_cast<amrex::Real>(drop.drop_time_s);
        if (current_time < t_drop) continue;  // drop hasn't happened yet

        amrex::Real age = current_time - t_drop;
        amrex::Real tau = static_cast<amrex::Real>(drop.tau_decay_s);
        amrex::Real eff = static_cast<amrex::Real>(drop.effectiveness);
        amrex::Real f_active = eff * std::exp(-age / tau);
        if (f_active < amrex::Real(1.0e-6)) continue;  // negligible effect

        const int nv = static_cast<int>(drop.vx.size());

        // Copy polygon vertices to device
        amrex::Gpu::DeviceVector<amrex::Real> d_vx(nv), d_vy(nv);
        for (int vi = 0; vi < nv; ++vi) {
            d_vx[vi] = static_cast<amrex::Real>(drop.vx[vi]);
            d_vy[vi] = static_cast<amrex::Real>(drop.vy[vi]);
        }
        const amrex::Real* pvx = d_vx.data();
        const amrex::Real* pvy = d_vy.data();

        const amrex::Real f_ret_val = amrex::Real(1.0) - f_active;
        const int nv_g = nv;

        for (MFIter mfi(R_mf); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            auto       R  = R_mf.array(mfi);
            const auto gdx  = dx;
            const auto gplo = plo;

            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                amrex::Real x = gplo[0] + (i + amrex::Real(0.5)) * gdx[0];
                amrex::Real y = gplo[1] + (j + amrex::Real(0.5)) * gdx[1];
                if (retardant_point_in_polygon(x, y, pvx, pvy, nv_g)) {
                    R(i, j, k) = amrex::max(amrex::Real(0.0), R(i, j, k) * f_ret_val);
                }
            });
        }
    }
}

// ============================================================================
// apply_retardant_to_spotting_probability
//
// Inside active retardant drop zones, scale the spotting probability field
// (component 0 of spotting_prob_mf) by the same (1 - f_active) factor used
// for ROS suppression.  This prevents firebrand spot fires from starting in
// retardant-covered cells where sustained ignition is less likely.
//
// Parameters:
//   spotting_prob_mf – spotting data MultiFab (comp 0 = probability, modified in place)
//   drops            – retardant drop list (unchanged)
//   geom             – AMReX geometry
//   current_time     – current simulation time [s]
// ============================================================================
inline void apply_retardant_to_spotting_probability(
    MultiFab&               spotting_prob_mf,
    const RetardantDropList& drops,
    const Geometry&          geom,
    amrex::Real              current_time)
{
    if (drops.empty()) return;

    const auto dx  = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();

    for (const auto& drop : drops.drops) {
        const amrex::Real t_drop = static_cast<amrex::Real>(drop.drop_time_s);
        if (current_time < t_drop) continue;

        amrex::Real age = current_time - t_drop;
        amrex::Real tau = static_cast<amrex::Real>(drop.tau_decay_s);
        amrex::Real eff = static_cast<amrex::Real>(drop.effectiveness);
        amrex::Real f_active = eff * std::exp(-age / tau);
        if (f_active < amrex::Real(1.0e-6)) continue;

        const int nv = static_cast<int>(drop.vx.size());

        amrex::Gpu::DeviceVector<amrex::Real> d_vx(nv), d_vy(nv);
        for (int vi = 0; vi < nv; ++vi) {
            d_vx[vi] = static_cast<amrex::Real>(drop.vx[vi]);
            d_vy[vi] = static_cast<amrex::Real>(drop.vy[vi]);
        }
        const amrex::Real* pvx = d_vx.data();
        const amrex::Real* pvy = d_vy.data();

        const amrex::Real f_ret_val = amrex::Real(1.0) - f_active;
        const int nv_g = nv;

        for (MFIter mfi(spotting_prob_mf); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            auto       sp = spotting_prob_mf.array(mfi);
            const auto gdx  = dx;
            const auto gplo = plo;

            // Scale only component 0 (spotting probability)
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                amrex::Real x = gplo[0] + (i + amrex::Real(0.5)) * gdx[0];
                amrex::Real y = gplo[1] + (j + amrex::Real(0.5)) * gdx[1];
                if (retardant_point_in_polygon(x, y, pvx, pvy, nv_g)) {
                    sp(i, j, k, 0) = amrex::max(amrex::Real(0.0), sp(i, j, k, 0) * f_ret_val);
                }
            });
        }
    }
}

#endif // RETARDANT_DROP_H
