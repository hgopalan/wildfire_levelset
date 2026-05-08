// ============================================================================
// Polygon / Polyline Fire Ignition Rasterizer
// ============================================================================
// Initialises the level-set field phi from a fire ignition polygon or polyline
// described as a sequence of (X, Y [, Z]) vertices in a CSV file.
//
// Two modes:
//
//  1. **Closed polygon** (source_type = "polygon"):
//     The interior of the polygon (inside the closed boundary) is set to
//     phi = -half_cell (burned), the exterior to the signed distance from the
//     polygon boundary (approximate Euclidean distance, capped at a multiple
//     of the cell size).  Polygon closure is handled automatically — the last
//     vertex is connected back to the first.
//
//  2. **Polyline / line fire** (source_type = "polyline"):
//     Cells within `polyline_width` [physical units] of any segment of the
//     polyline are set to phi = -half_cell; other cells receive the distance
//     to the nearest segment.
//
// CSV file format (whitespace or comma separated, comment lines start with #):
//   X1, Y1[, Z1]
//   X2, Y2[, Z2]
//   ...
//
// Input parameters:
//   fire_polygon_file   – path to CSV vertex file
//   source_type         – "polygon" or "polyline"
//   polyline_width      – half-width [m] of line-ignition zone (default: 10.0)
//   fire_polygon_level  – z-coordinate for the SDF burn zone in 3-D runs;
//                         cells with |z - fire_polygon_level| > dz are not set
//                         (default: 0.5, i.e. the first vertical layer)
//
// Algorithm (polygon SDF approximation):
//   For each grid cell, compute the minimum distance to any polygon edge.
//   The sign is determined by the winding number: cells inside the polygon
//   receive a negative sign (burned), cells outside receive a positive sign.
//   This gives a signed distance field consistent with the level-set convention
//   used throughout the solver.
//
// Reference:
//   Finney, M.A. (1998). FARSITE: Fire Area Simulator. USDA FS RMRS-RP-4.
//   "An ignition polygon defines the initial fire perimeter at the start of the
//    simulation period."
// ============================================================================

#ifndef POLYGON_IGNITION_H
#define POLYGON_IGNITION_H

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuContainers.H>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace amrex;

// ---------------------------------------------------------------------------
// read_polygon_vertices
//
// Parse a CSV/whitespace-delimited vertex file into vectors of X and Y coords.
// Lines starting with # or ! are comments.
// ---------------------------------------------------------------------------
inline void read_polygon_vertices(const std::string& filename,
                                  std::vector<Real>& xs,
                                  std::vector<Real>& ys)
{
    xs.clear();
    ys.clear();
    std::ifstream f(filename);
    if (!f.is_open()) {
        amrex::Abort("PolygonIgnition: cannot open file: " + filename);
    }
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#' || line[0] == '!') continue;
        for (char& c : line) if (c == ',') c = ' ';
        std::istringstream iss(line);
        Real x, y;
        if (iss >> x >> y) {
            xs.push_back(x);
            ys.push_back(y);
        }
    }
    if (xs.size() < 2) {
        amrex::Abort("PolygonIgnition: need at least 2 vertices, got "
                     + std::to_string(xs.size()) + " in: " + filename);
    }
    amrex::Print() << "PolygonIgnition: loaded " << xs.size()
                   << " vertices from " << filename << "\n";
}

// ---------------------------------------------------------------------------
// point_in_polygon (winding number test, 2-D)
// Returns 1 if (px, py) is strictly inside the polygon, 0 otherwise.
// ---------------------------------------------------------------------------
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int polygon_winding_number(Real px, Real py,
                           const Real* vx, const Real* vy, int n)
{
    int wn = 0;
    for (int i = 0; i < n; ++i) {
        int j = (i + 1) % n;
        Real x0 = vx[i], y0 = vy[i];
        Real x1 = vx[j], y1 = vy[j];
        if (y0 <= py) {
            if (y1 > py) {
                // Upward crossing
                Real cross = (x1 - x0) * (py - y0) - (px - x0) * (y1 - y0);
                if (cross > Real(0.0)) ++wn;
            }
        } else {
            if (y1 <= py) {
                // Downward crossing
                Real cross = (x1 - x0) * (py - y0) - (px - x0) * (y1 - y0);
                if (cross < Real(0.0)) --wn;
            }
        }
    }
    return wn;
}

// ---------------------------------------------------------------------------
// segment_distance_sq (2-D squared distance from point to line segment)
// ---------------------------------------------------------------------------
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real segment_dist_sq(Real px, Real py,
                     Real ax, Real ay,
                     Real bx, Real by)
{
    const Real dx = bx - ax, dy = by - ay;
    const Real len2 = dx*dx + dy*dy;
    // Use epsilon appropriate for metre-scale UTM coordinates
    if (len2 < Real(1.0e-10)) {
        const Real ex = px-ax, ey = py-ay;
        return ex*ex + ey*ey;
    }
    Real t = ((px-ax)*dx + (py-ay)*dy) / len2;
    t = amrex::max(Real(0.0), amrex::min(Real(1.0), t));
    const Real cx = ax + t*dx - px;
    const Real cy = ay + t*dy - py;
    return cx*cx + cy*cy;
}

// ============================================================================
// init_phi_from_polygon
//
// Rasterise a closed polygon ignition into phi as a signed-distance function.
//
// Parameters:
//   phi       – level-set MultiFab (1 component, modified in-place)
//   geom      – AMReX geometry (provides cell centres)
//   xs, ys    – polygon vertex coordinates [physical units]
//   z_level   – Z centre of the ignition layer in 3-D runs (ignored in 2-D)
// ============================================================================
inline void init_phi_from_polygon(MultiFab&              phi,
                                  const Geometry&        geom,
                                  const std::vector<Real>& xs,
                                  const std::vector<Real>& ys,
                                  Real                   z_level = 0.5)
{
    const int n = static_cast<int>(xs.size());

    // Copy vertices to device
    amrex::Gpu::DeviceVector<Real> d_vx(xs.size());
    amrex::Gpu::DeviceVector<Real> d_vy(ys.size());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, xs.begin(), xs.end(), d_vx.begin());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, ys.begin(), ys.end(), d_vy.begin());
    const Real* vx = d_vx.data();
    const Real* vy = d_vy.data();

    const Real* prob_lo = geom.ProbLo();
    const Real* dx_arr  = geom.CellSize();

    for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto p = phi.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real cx = prob_lo[0] + (i + Real(0.5)) * dx_arr[0];
            Real cy = prob_lo[1] + (j + Real(0.5)) * dx_arr[1];
#if (AMREX_SPACEDIM == 3)
            Real cz = prob_lo[2] + (k + Real(0.5)) * dx_arr[2];
            // Only modify cells in the ignition z-layer
            if (std::abs(cz - z_level) > dx_arr[2]) return;
#else
            (void)z_level;
#endif
            // Minimum distance to any polygon edge
            Real min_dist2 = Real(1.0e30);
            for (int seg = 0; seg < n; ++seg) {
                int nxt = (seg + 1) % n;
                Real d2 = segment_dist_sq(cx, cy, vx[seg], vy[seg], vx[nxt], vy[nxt]);
                if (d2 < min_dist2) min_dist2 = d2;
            }
            Real dist = std::sqrt(min_dist2);

            // Sign from winding number: inside = negative (burned)
            const int wn = polygon_winding_number(cx, cy, vx, vy, n);
            p(i, j, k) = (wn != 0) ? -dist : dist;
        });
    }
    amrex::Print() << "PolygonIgnition: rasterised closed polygon ("
                   << n << " vertices)\n";
}

// ============================================================================
// init_phi_from_polyline
//
// Rasterise a polyline (line fire) into phi.  Cells within half_width of any
// segment receive phi = −half_width (inside); others receive the distance to
// the nearest segment (always positive).
// ============================================================================
inline void init_phi_from_polyline(MultiFab&               phi,
                                   const Geometry&         geom,
                                   const std::vector<Real>& xs,
                                   const std::vector<Real>& ys,
                                   Real half_width = 10.0,
                                   Real z_level    = 0.5)
{
    const int n = static_cast<int>(xs.size());
    const int nseg = n - 1;

    amrex::Gpu::DeviceVector<Real> d_vx(xs.size());
    amrex::Gpu::DeviceVector<Real> d_vy(ys.size());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, xs.begin(), xs.end(), d_vx.begin());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, ys.begin(), ys.end(), d_vy.begin());
    const Real* vx  = d_vx.data();
    const Real* vy  = d_vy.data();

    const Real* prob_lo = geom.ProbLo();
    const Real* dx_arr  = geom.CellSize();
    const Real hw = half_width;

    for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto p = phi.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real cx = prob_lo[0] + (i + Real(0.5)) * dx_arr[0];
            Real cy = prob_lo[1] + (j + Real(0.5)) * dx_arr[1];
#if (AMREX_SPACEDIM == 3)
            Real cz = prob_lo[2] + (k + Real(0.5)) * dx_arr[2];
            if (std::abs(cz - z_level) > dx_arr[2]) return;
#else
            (void)z_level;
#endif
            Real min_dist2 = Real(1.0e30);
            for (int seg = 0; seg < nseg; ++seg) {
                Real d2 = segment_dist_sq(cx, cy, vx[seg], vy[seg], vx[seg+1], vy[seg+1]);
                if (d2 < min_dist2) min_dist2 = d2;
            }
            Real dist = std::sqrt(min_dist2);
            p(i, j, k) = (dist <= hw) ? -hw : dist;
        });
    }
    amrex::Print() << "PolygonIgnition: rasterised polyline ("
                   << nseg << " segments, half-width=" << half_width << " m)\n";
}

#endif // POLYGON_IGNITION_H
