// ============================================================================
// Minimum Travel Time (MTT) Fire Propagation
// ============================================================================
//
// MTT computes, once at startup, the minimum time for fire to reach every
// grid cell from the initial burn region.  The level-set field phi is then
// simply updated at each simulation time step T as:
//
//     phi(i,j,k) = arrival_time(i,j,k) - T
//
// giving:
//   phi < 0  →  burned  (fire arrived before current time)
//   phi = 0  →  fire front  (arriving exactly now)
//   phi > 0  →  unburned  (fire not yet arrived)
//
// The arrival times are computed using a Dijkstra-like fast-marching scheme
// on the Cartesian grid.  For each cell the travel time from a burning
// neighbour is:
//
//     t_new = t_neighbour + ds / R_interface
//
// where ds is the cell-centre-to-cell-centre distance (dx or dy) and
// R_interface is the harmonic mean of the ROS values at the two cell centres.
// (Harmonic mean avoids division by small ROS values blowing up travel times.)
//
// The algorithm is CPU-only (priority queue cannot efficiently run on GPU);
// it is called once before the time loop and the result stored in the
// provided arrival_time_mf MultiFab (same field used by the rest of the solver
// for diagnostics and checkpoint I/O).
//
// Reference:
//   Finney, M.A. (2002). "Fire growth using minimum travel time methods."
//   Canadian Journal of Forest Research, 32(8), 1420–1424.
// ============================================================================

#ifndef MTT_PROPAGATION_H
#define MTT_PROPAGATION_H

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <queue>
#include <vector>
#include <limits>
#include <cmath>

using namespace amrex;

// ============================================================================
// compute_mtt_arrival_times
//
// Parameters
// ----------
//   arrival_time_mf  – output: per-cell minimum arrival time [s].
//                      Initialised to -1.0 on input for unburned cells and
//                      to 0.0 for already-burning cells.
//                      On return, every reachable cell has arrival_time >= 0.
//                      Unreachable cells (isolated by zero-ROS walls) retain
//                      a large sentinel value (1e30).
//   phi              – level-set field at t=0; cells with phi < 0 are ignition
//                      sources (arrival_time = 0).
//   R_mf             – rate-of-spread field [m/s], computed by the selected
//                      fire spread model; must have at least 0 ghost cells.
//   geom             – AMReX Geometry describing the Cartesian grid.
// ============================================================================
inline void compute_mtt_arrival_times(MultiFab& arrival_time_mf,
                                      const MultiFab& phi,
                                      const MultiFab& R_mf,
                                      const Geometry& geom)
{
    const auto dx_arr = geom.CellSize();
    const Real dx = static_cast<Real>(dx_arr[0]);
    const Real dy = static_cast<Real>(dx_arr[1]);
#if (AMREX_SPACEDIM == 3)
    const Real dz = static_cast<Real>(dx_arr[2]);
#endif

    // -----------------------------------------------------------------
    // 1. Copy phi, R, and initial arrival_time to host vectors for the
    //    priority-queue algorithm.
    // -----------------------------------------------------------------
    // We work with the domain box (no ghosts) on the I/O processor's portion
    // of the grid.  Since MTT is a global-connectivity algorithm we need
    // the full domain.  Copy all data to rank-0 (host) arrays via AllGather.

    const Box& domain = geom.Domain();
    const IntVect lo   = domain.smallEnd();
    const IntVect hi   = domain.bigEnd();
    const int nx = hi[0] - lo[0] + 1;
    const int ny = hi[1] - lo[1] + 1;
#if (AMREX_SPACEDIM == 3)
    const int nz = hi[2] - lo[2] + 1;
#else
    const int nz = 1;
#endif
    const long ncells = static_cast<long>(nx) * ny * nz;

    // Helper: flat index from (i,j,k) relative to lo
    // In 2-D builds AMREX_SPACEDIM == 2 so lo/hi are 2-element IntVects;
    // accessing lo[2] would trigger an assertion.  Use #if to guard.
    auto idx = [&](int i, int j, int k) -> long {
#if (AMREX_SPACEDIM == 3)
        return static_cast<long>(i - lo[0]) +
               static_cast<long>(nx) * (j - lo[1]) +
               static_cast<long>(nx) * static_cast<long>(ny) * (k - lo[2]);
#else
        amrex::ignore_unused(k);
        return static_cast<long>(i - lo[0]) +
               static_cast<long>(nx) * (j - lo[1]);
#endif
    };

    std::vector<Real> h_phi(ncells, Real(1.0));
    std::vector<Real> h_R(ncells, Real(0.0));
    std::vector<Real> h_t(ncells, Real(1e30));

    // Fill from MultiFabs (each rank only owns part of the data)
    for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const p_arr = phi.const_array(mfi);
        auto const r_arr = R_mf.const_array(mfi);
        auto const t_arr = arrival_time_mf.const_array(mfi);
        amrex::LoopOnCpu(bx, [&](int i, int j, int k) {
            long n = idx(i, j, k);
            h_phi[n] = p_arr(i, j, k);
            h_R[n]   = r_arr(i, j, k);
            // Use existing arrival_time initialisation (0 for ignition, -1 for rest)
            Real at = t_arr(i, j, k);
            h_t[n]  = (at >= Real(0.0)) ? at : Real(1e30);
        });
    }

    // Reduce across MPI ranks so every rank has the full arrays
    amrex::ParallelDescriptor::ReduceRealSum(h_phi.data(), ncells);
    amrex::ParallelDescriptor::ReduceRealSum(h_R.data(), ncells);
    // For arrival_time: we want the minimum nonzero value → use custom reduce
    {
        // Replace 1e30 with 0 for summation-based "min" (each rank owns 
        // exactly the cells in its validbox; sum gives the correct values).
        std::vector<Real> h_t_local(ncells, Real(0.0));
        for (MFIter mfi(arrival_time_mf); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            auto const t_arr = arrival_time_mf.const_array(mfi);
            amrex::LoopOnCpu(bx, [&](int i, int j, int k) {
                long n = idx(i, j, k);
                Real at = t_arr(i, j, k);
                h_t_local[n] = (at >= Real(0.0)) ? at : Real(0.0);
            });
        }
        amrex::ParallelDescriptor::ReduceRealSum(h_t_local.data(), ncells);
        // Any cell with 0 in h_t_local and phi >= 0 → unvisited → sentinel 1e30
        // Any cell with h_t_local > 0 is an ignition cell
        for (long n = 0; n < ncells; ++n) {
            if (h_t_local[n] > Real(0.0))
                h_t[n] = h_t_local[n];   // ignition cell at t = h_t_local[n]
        }
    }

    // Seed: mark ignition cells (initial phi < 0) as arrival_time = 0
    for (long n = 0; n < ncells; ++n) {
        if (h_phi[n] < Real(0.0))
            h_t[n] = Real(0.0);
    }

    // -----------------------------------------------------------------
    // 2. Dijkstra priority-queue propagation
    // -----------------------------------------------------------------
    // Min-heap entry: (arrival_time, flat_index)
    using Entry = std::pair<Real, long>;
    std::priority_queue<Entry, std::vector<Entry>, std::greater<Entry>> pq;

    for (long n = 0; n < ncells; ++n) {
        if (h_t[n] < Real(1e29))
            pq.push({h_t[n], n});
    }

    // Neighbour offsets: ±x, ±y, ±z directions
    struct Neighbour { int di, dj, dk; Real ds; };
    std::vector<Neighbour> neighbours;
    neighbours.push_back({+1,  0,  0, dx});
    neighbours.push_back({-1,  0,  0, dx});
    neighbours.push_back({ 0, +1,  0, dy});
    neighbours.push_back({ 0, -1,  0, dy});
#if (AMREX_SPACEDIM == 3)
    neighbours.push_back({ 0,  0, +1, dz});
    neighbours.push_back({ 0,  0, -1, dz});
#endif

    while (!pq.empty()) {
        auto [t_curr, n_curr] = pq.top();
        pq.pop();

        // Skip if this is a stale entry
        if (t_curr > h_t[n_curr] + Real(1e-12)) continue;

        // Decode flat index to (i,j,k)
        long tmp = n_curr;
        int i = lo[0] + static_cast<int>(tmp % nx);
        tmp /= nx;
        int j = lo[1] + static_cast<int>(tmp % ny);
        int k = lo[2] + static_cast<int>(tmp / ny);

        const Real R_curr = h_R[n_curr];

        for (const auto& nb : neighbours) {
            int ni = i + nb.di;
            int nj = j + nb.dj;
            int nk = k + nb.dk;

            // Domain bounds check
            if (ni < lo[0] || ni > hi[0]) continue;
            if (nj < lo[1] || nj > hi[1]) continue;
#if (AMREX_SPACEDIM == 3)
            if (nk < lo[2] || nk > hi[2]) continue;
#endif

            long n_nb = idx(ni, nj, nk);
            const Real R_nb = h_R[n_nb];

            // Skip zero-ROS barriers (no fire can pass through)
            if (R_curr <= Real(0.0) || R_nb <= Real(0.0)) continue;

            // Harmonic mean ROS at the interface
            const Real R_interface = Real(2.0) * R_curr * R_nb / (R_curr + R_nb);

            const Real t_new = h_t[n_curr] + nb.ds / R_interface;

            if (t_new < h_t[n_nb] - Real(1e-12)) {
                h_t[n_nb] = t_new;
                pq.push({t_new, n_nb});
            }
        }
    }

    // -----------------------------------------------------------------
    // 3. Write results back to the arrival_time MultiFab
    // -----------------------------------------------------------------
    for (MFIter mfi(arrival_time_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto t_arr = arrival_time_mf.array(mfi);
        amrex::LoopOnCpu(bx, [&](int i, int j, int k) {
            long n = idx(i, j, k);
            Real at = h_t[n];
            // Keep sentinel as -1 for truly unreachable cells (consistent with
            // the rest of the solver which uses -1 for "not yet burned")
            t_arr(i, j, k) = (at < Real(1e29)) ? at : Real(-1.0);
        });
    }

    amrex::Print() << "MTT: arrival times computed for " << ncells << " cells.\n";
}

// ============================================================================
// apply_mtt_phi_update
//
// At simulation time T, set phi = arrival_time - T so the zero-contour
// tracks the fire front.
//
// Parameters
// ----------
//   phi             – level-set field to update in-place.
//   arrival_time_mf – MTT arrival times [s] (computed by compute_mtt_arrival_times).
//   time            – current simulation time [s].
// ============================================================================
inline void apply_mtt_phi_update(MultiFab& phi,
                                 const MultiFab& arrival_time_mf,
                                 Real time)
{
    for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto p_arr = phi.array(mfi);
        auto const t_arr = arrival_time_mf.const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const Real at = t_arr(i, j, k);
            if (at < Real(0.0)) {
                // Unreachable cell: keep a large positive phi (unburned)
                p_arr(i, j, k) = Real(1e30);
            } else {
                p_arr(i, j, k) = at - time;
            }
        });
    }
}

#endif // MTT_PROPAGATION_H
