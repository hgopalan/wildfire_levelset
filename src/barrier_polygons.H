// ============================================================================
// FARSITE Barrier Polygon / Firebreak Support
// ============================================================================
//
// Reads one or more CSV files, each defining a polyline or closed polygon of
// barrier vertices (X Y coordinates, one per line).  For each vertex the
// nearest grid cell is located.  At every call to apply_barrier_polygons()
// any of those cells that are currently burning (phi < 0) are switched to
// not-burning by setting phi = +0.5 * min(dx,dy).
//
// This models firebreaks, fuel breaks, roads, and other barriers that prevent
// fire spread without using the Embedded Boundary (EB) machinery.
//
// Input format (each file)
// -------------------------
// Lines beginning with '#' or empty lines are ignored.
// Each data line must contain at least two whitespace- or comma-separated
// floating-point values:   X  Y  [Z]
// The coordinates must be in the same units as the simulation domain (metres
// for UTM; model units otherwise).
//
// Usage
// -----
// From the input file:
//
//   barrier_files = barrier1.csv  barrier2.csv
//
// or with a single file:
//
//   barrier_files = road_break.csv
//
// The barrier is applied every time step.  The applied cells are determined
// once at initialisation (load_barrier_cells) and reused every time step.
// ============================================================================

#ifndef BARRIER_POLYGONS_H
#define BARRIER_POLYGONS_H

#include <AMReX.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_Gpu.H>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>

using namespace amrex;

// ============================================================================
// load_barrier_cells
//
// Reads a list of CSV barrier files and returns the set of cell indices
// (i, j, k) that lie nearest to any barrier vertex.
//
// Parameters
// ----------
//   barrier_files – list of CSV file paths (may be empty → returns empty set)
//   geom          – AMReX Geometry
//
// Returns a vector of IntVect, one entry per unique barrier cell (domain-clipped).
// ============================================================================
static std::vector<IntVect> load_barrier_cells(
    const std::vector<std::string>& barrier_files,
    const Geometry& geom)
{
    if (barrier_files.empty()) return {};

    const auto plo  = geom.ProbLo();
    const auto dx   = geom.CellSize();
    const Box& domain = geom.Domain();
    const IntVect lo  = domain.smallEnd();
    const IntVect hi  = domain.bigEnd();

    // Use a simple sorted-unique approach to de-duplicate cells
    std::vector<IntVect> cells;

    // Read on IOProcessor then broadcast
    std::vector<Real> all_x, all_y, all_z;

    if (amrex::ParallelDescriptor::IOProcessor()) {
        for (const auto& filename : barrier_files) {
            std::ifstream fin(filename);
            if (!fin.is_open()) {
                amrex::Warning("barrier_polygons: cannot open file '" + filename + "' – skipping");
                continue;
            }
            int n_pts = 0;
            std::string line;
            while (std::getline(fin, line)) {
                if (line.empty() || line[0] == '#') continue;
                for (char& c : line) if (c == ',') c = ' ';
                std::istringstream iss(line);
                Real x = 0.0, y = 0.0, z = 0.0;
                if (!(iss >> x >> y)) continue;
                iss >> z; // optional Z
                all_x.push_back(x);
                all_y.push_back(y);
                all_z.push_back(z);
                ++n_pts;
            }
            amrex::Print() << "barrier_polygons: loaded " << n_pts
                           << " vertex/vertices from '" << filename << "'\n";
        }
    }

    // Broadcast to all MPI ranks
    {
        int n = static_cast<int>(all_x.size());
        amrex::ParallelDescriptor::Bcast(&n, 1, amrex::ParallelDescriptor::IOProcessorNumber());
        all_x.resize(n); all_y.resize(n); all_z.resize(n);
        amrex::ParallelDescriptor::Bcast(all_x.data(), n, amrex::ParallelDescriptor::IOProcessorNumber());
        amrex::ParallelDescriptor::Bcast(all_y.data(), n, amrex::ParallelDescriptor::IOProcessorNumber());
        amrex::ParallelDescriptor::Bcast(all_z.data(), n, amrex::ParallelDescriptor::IOProcessorNumber());
    }

    if (all_x.empty()) return {};

    // Convert each vertex to the nearest cell index
    for (int n = 0; n < static_cast<int>(all_x.size()); ++n) {
        int i = lo[0] + static_cast<int>(std::floor((all_x[n] - plo[0]) / dx[0]));
        int j = lo[1] + static_cast<int>(std::floor((all_y[n] - plo[1]) / dx[1]));
#if (AMREX_SPACEDIM == 3)
        int k = lo[2] + static_cast<int>(std::floor((all_z[n] - plo[2]) / dx[2]));
#else
        int k = 0; // always 0 in 2-D (single z-layer; lo[2] would be out-of-bounds)
#endif
        // Clamp to domain
        i = std::max(lo[0], std::min(hi[0], i));
        j = std::max(lo[1], std::min(hi[1], j));
#if (AMREX_SPACEDIM == 3)
        k = std::max(lo[2], std::min(hi[2], k));
#else
        amrex::ignore_unused(k);
#endif
        cells.push_back(IntVect(AMREX_D_DECL(i, j, k)));
    }

    // De-duplicate
    std::sort(cells.begin(), cells.end(),
              [](const IntVect& a, const IntVect& b) {
                  for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                      if (a[d] != b[d]) return a[d] < b[d];
                  }
                  return false;
              });
    cells.erase(std::unique(cells.begin(), cells.end()), cells.end());

    amrex::Print() << "barrier_polygons: " << cells.size()
                   << " unique barrier cell(s) registered.\n";
    return cells;
}

// ============================================================================
// apply_barrier_polygons
//
// For each barrier cell that is currently burning (phi < 0), set
// phi = +half_cell so it becomes unburned.
//
// Parameters
// ----------
//   phi           – level-set field, modified in-place.
//   barrier_cells – list of barrier cell indices returned by load_barrier_cells.
//   geom          – AMReX Geometry (used to obtain cell spacing).
// ============================================================================
inline void apply_barrier_polygons(MultiFab& phi,
                                   const std::vector<IntVect>& barrier_cells,
                                   const Geometry& geom)
{
    if (barrier_cells.empty()) return;

    const auto dx  = geom.CellSize();
#if (AMREX_SPACEDIM == 3)
    const Real half_cell = Real(0.5) * std::min({dx[0], dx[1], dx[2]});
#else
    const Real half_cell = Real(0.5) * std::min(dx[0], dx[1]);
#endif

    // Copy barrier cells to a GPU-accessible device vector
    amrex::Gpu::DeviceVector<IntVect> d_barrier(barrier_cells.size());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                     barrier_cells.begin(), barrier_cells.end(),
                     d_barrier.begin());

    const int n_barrier = static_cast<int>(barrier_cells.size());
    const IntVect* barrier_ptr = d_barrier.data();

    for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto p_arr = phi.array(mfi);

        ParallelFor(n_barrier, [=] AMREX_GPU_DEVICE (int b) noexcept {
            const IntVect iv = barrier_ptr[b];
            // Only process if this cell belongs to the current MFIter box
            if (!bx.contains(iv)) return;
            int i = iv[0], j = iv[1];
#if (AMREX_SPACEDIM == 3)
            int k = iv[2];
#else
            int k = 0;
#endif
            // Switch to not-burning if currently burning
            if (p_arr(i, j, k) < Real(0.0)) {
                p_arr(i, j, k) = half_cell;
            }
        });
    }
}

#endif // BARRIER_POLYGONS_H
