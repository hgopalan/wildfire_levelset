// ============================================================================
// Ground Suppression Resources: Dozer Lines and Hand Crew Lines
// ============================================================================
// Models time-scheduled firebreak construction by ground suppression resources
// (dozer lines, hand crew lines).  A suppression line becomes an active
// firebreak once construction is complete (current_time >= line_end_s).
//
// This mirrors FARSITE's Resource File (*.rsc) which defines lines with:
//   – start time (when construction begins),
//   – completion time (when the line is active as a firebreak),
//   – line type (dozer = 0, hand crew = 1; both treated as total firebreak),
//   – polyline vertices (X Y coordinates).
//
// CSV file format (whitespace or comma delimited; # or ! comment lines):
//   LINE_START_S   LINE_END_S   LINE_TYPE   N_VERTICES
//   X1  Y1
//   X2  Y2
//   ...   (N_VERTICES lines of vertex coordinates)
//   [next LINE_START_S LINE_END_S LINE_TYPE N_VERTICES block]
//
// Example:
//   # Ground crew hand line (completed in 2 hours)
//   0.0   7200.0   1   3
//   330100.0  3775100.0
//   330200.0  3775200.0
//   330300.0  3775200.0
//
// Behaviour:
//   At each time step, cells along any line whose LINE_END_S <= current_time
//   are treated as firebreaks: phi is set to +half_cell if currently burning
//   (phi < 0).  Cells on lines not yet complete (LINE_END_S > current_time)
//   are left unchanged.
//
// References:
//   Finney, M.A. (1998). FARSITE: Fire Area Simulator. USDA FS RMRS-RP-4.
// ============================================================================

#ifndef SUPPRESSION_RESOURCES_H
#define SUPPRESSION_RESOURCES_H

#include <AMReX.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_Gpu.H>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>

using namespace amrex;

// ============================================================================
// SuppressionLine
// Represents one time-scheduled suppression line with its activation window
// and the set of grid cells it covers.
// ============================================================================
struct SuppressionLine {
    double start_s;              // construction start time [s]
    double end_s;                // completion time [s] – line active when current_time >= end_s
    int    line_type;            // 0 = dozer, 1 = hand crew (informational only)
    std::vector<IntVect> cells;  // grid cells along the line
};

// ============================================================================
// load_suppression_lines
//
// Parses the suppression CSV file and converts polyline vertices to grid cells.
//
// Parameters
// ----------
//   filename  – path to suppression schedule CSV
//   geom      – AMReX Geometry (for coordinate → cell mapping)
//
// Returns a vector of SuppressionLine, one per line block in the CSV.
// ============================================================================
static std::vector<SuppressionLine> load_suppression_lines(
    const std::string& filename,
    const Geometry& geom)
{
    std::vector<SuppressionLine> lines;
    if (filename.empty()) return lines;

    const auto plo = geom.ProbLo();
    const auto dx  = geom.CellSize();
    const Box& dom = geom.Domain();
    const IntVect lo = dom.smallEnd();
    const IntVect hi = dom.bigEnd();

    // Read entire file on IOProcessor, then broadcast
    std::vector<double> h_start, h_end;
    std::vector<int>    h_type;
    // Flat list: groups of (start, end, type, n) followed by n (x,y) pairs
    // We use a simple serial parser.
    struct RawLine {
        double start_s, end_s;
        int    type;
        std::vector<double> xs, ys;
    };
    std::vector<RawLine> raw_lines;

    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ifstream fin(filename);
        if (!fin.is_open()) {
            amrex::Abort("suppression_resources: cannot open file '" + filename + "'");
        }
        std::string line;
        while (std::getline(fin, line)) {
            if (line.empty() || line[0] == '#' || line[0] == '!') continue;
            for (char& c : line) if (c == ',') c = ' ';
            std::istringstream iss(line);
            double start_s, end_s;
            int type, n_verts;
            if (!(iss >> start_s >> end_s >> type >> n_verts)) continue;
            if (n_verts < 2) {
                amrex::Warning("suppression_resources: line with < 2 vertices skipped");
                // Skip n_verts lines
                for (int v = 0; v < n_verts; ++v) std::getline(fin, line);
                continue;
            }
            RawLine rl;
            rl.start_s = start_s;
            rl.end_s   = end_s;
            rl.type    = type;
            for (int v = 0; v < n_verts; ++v) {
                std::string vline;
                while (std::getline(fin, vline)) {
                    if (!vline.empty() && vline[0] != '#' && vline[0] != '!') break;
                }
                for (char& c : vline) if (c == ',') c = ' ';
                std::istringstream vss(vline);
                double vx, vy;
                if (!(vss >> vx >> vy)) {
                    amrex::Warning("suppression_resources: malformed vertex line; skipping");
                    break;
                }
                rl.xs.push_back(vx);
                rl.ys.push_back(vy);
            }
            if (static_cast<int>(rl.xs.size()) >= 2)
                raw_lines.push_back(rl);
        }
        amrex::Print() << "suppression_resources: loaded "
                       << raw_lines.size() << " suppression line(s) from '"
                       << filename << "'\n";
    }

    // Broadcast number of lines
    int n_lines = static_cast<int>(raw_lines.size());
    amrex::ParallelDescriptor::Bcast(&n_lines, 1,
        amrex::ParallelDescriptor::IOProcessorNumber());
    if (n_lines == 0) return lines;

    // For simplicity, broadcast each line's data sequentially
    // (suppression lines are parsed once at startup, so cost is negligible)
    if (!amrex::ParallelDescriptor::IOProcessor()) {
        raw_lines.resize(n_lines);
    }
    for (int li = 0; li < n_lines; ++li) {
        auto& rl = raw_lines[li];
        amrex::ParallelDescriptor::Bcast(&rl.start_s, 1,
            amrex::ParallelDescriptor::IOProcessorNumber());
        amrex::ParallelDescriptor::Bcast(&rl.end_s, 1,
            amrex::ParallelDescriptor::IOProcessorNumber());
        amrex::ParallelDescriptor::Bcast(&rl.type, 1,
            amrex::ParallelDescriptor::IOProcessorNumber());
        int nv = static_cast<int>(rl.xs.size());
        amrex::ParallelDescriptor::Bcast(&nv, 1,
            amrex::ParallelDescriptor::IOProcessorNumber());
        rl.xs.resize(nv);
        rl.ys.resize(nv);
        amrex::ParallelDescriptor::Bcast(rl.xs.data(), nv,
            amrex::ParallelDescriptor::IOProcessorNumber());
        amrex::ParallelDescriptor::Bcast(rl.ys.data(), nv,
            amrex::ParallelDescriptor::IOProcessorNumber());
    }

    // Convert vertices to grid cells by rasterising each polyline segment
    lines.resize(n_lines);
    for (int li = 0; li < n_lines; ++li) {
        const auto& rl = raw_lines[li];
        lines[li].start_s  = rl.start_s;
        lines[li].end_s    = rl.end_s;
        lines[li].line_type = rl.type;

        std::vector<IntVect>& cells = lines[li].cells;
        const int nv = static_cast<int>(rl.xs.size());

        // For each segment, trace cells using Bresenham-style rasterisation
        for (int v = 0; v < nv; ++v) {
            // Add vertex cell
            int ci = lo[0] + static_cast<int>(std::floor((rl.xs[v] - plo[0]) / dx[0]));
            int cj = lo[1] + static_cast<int>(std::floor((rl.ys[v] - plo[1]) / dx[1]));
            ci = std::max(lo[0], std::min(hi[0], ci));
            cj = std::max(lo[1], std::min(hi[1], cj));
            cells.push_back(IntVect(AMREX_D_DECL(ci, cj, lo[2])));

            if (v + 1 >= nv) break;

            // Rasterise segment v → v+1 by stepping in cell-width increments
            double x0 = rl.xs[v],   y0 = rl.ys[v];
            double x1 = rl.xs[v+1], y1 = rl.ys[v+1];
            double seg_len = std::sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0));
            int n_steps = std::max(1, static_cast<int>(std::ceil(seg_len / std::min(dx[0], dx[1]))));
            for (int s = 1; s <= n_steps; ++s) {
                double frac = static_cast<double>(s) / n_steps;
                double sx = x0 + frac * (x1 - x0);
                double sy = y0 + frac * (y1 - y0);
                int si = lo[0] + static_cast<int>(std::floor((sx - plo[0]) / dx[0]));
                int sj = lo[1] + static_cast<int>(std::floor((sy - plo[1]) / dx[1]));
                si = std::max(lo[0], std::min(hi[0], si));
                sj = std::max(lo[1], std::min(hi[1], sj));
                cells.push_back(IntVect(AMREX_D_DECL(si, sj, lo[2])));
            }
        }

        // De-duplicate cells
        std::sort(cells.begin(), cells.end(),
                  [](const IntVect& a, const IntVect& b) {
                      for (int d = 0; d < AMREX_SPACEDIM; ++d)
                          if (a[d] != b[d]) return a[d] < b[d];
                      return false;
                  });
        cells.erase(std::unique(cells.begin(), cells.end()), cells.end());

        amrex::Print() << "  Suppression line " << li
                       << " type=" << rl.type
                       << " active_at=" << rl.end_s << " s"
                       << " (" << cells.size() << " cells)\n";
    }

    return lines;
}

// ============================================================================
// apply_suppression_lines
//
// For each suppression line whose construction is complete
// (current_time >= line.end_s), sets phi = +half_cell for any burning cell
// (phi < 0) that lies on the line.
//
// Parameters
// ----------
//   phi           – level-set field, modified in-place
//   lines         – suppression lines (from load_suppression_lines)
//   current_time  – current simulation time [s]
//   geom          – AMReX Geometry
// ============================================================================
inline void apply_suppression_lines(MultiFab& phi,
                                    const std::vector<SuppressionLine>& lines,
                                    amrex::Real current_time,
                                    const Geometry& geom)
{
    if (lines.empty()) return;

    const auto dx  = geom.CellSize();
#if (AMREX_SPACEDIM == 3)
    const Real half_cell = Real(0.5) * std::min({dx[0], dx[1], dx[2]});
#else
    const Real half_cell = Real(0.5) * std::min(dx[0], dx[1]);
#endif

    for (const auto& sl : lines) {
        // Only apply completed lines
        if (static_cast<amrex::Real>(sl.end_s) > current_time) continue;
        if (sl.cells.empty()) continue;

        // Copy this line's cells to a GPU-accessible device vector
        amrex::Gpu::DeviceVector<IntVect> d_cells(sl.cells.size());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                         sl.cells.begin(), sl.cells.end(),
                         d_cells.begin());
        const IntVect* d_ptr = d_cells.data();
        const int n_cells    = static_cast<int>(sl.cells.size());

        for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            auto p = phi.array(mfi);
            const Real hc = half_cell;

            ParallelFor(n_cells, [=] AMREX_GPU_DEVICE (int idx) noexcept {
                const IntVect& iv = d_ptr[idx];
                if (!bx.contains(iv)) return;
                if (p(iv) < Real(0.0)) {
                    p(iv) = hc;
                }
            });
        }
    }
}

#endif // SUPPRESSION_RESOURCES_H
