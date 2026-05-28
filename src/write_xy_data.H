// write_xy_data.H - Functions for writing x-y data files for phi analysis
#ifndef WRITE_XY_DATA_H
#define WRITE_XY_DATA_H

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <limits>
#include <cmath>

using namespace amrex;

// ======================= Write all negative phi cells as x-y file ================
// Writes coordinates (x, y) of all cells where phi < 0 to a file
// ==================================================================================
inline void write_negative_phi_xy(const MultiFab& phi, 
                                   const Geometry& geom, 
                                   const std::string& filename)
{
    // Collect all negative phi cell centers
    std::vector<std::pair<Real, Real>> negative_points;
    
    const auto dx = geom.CellSizeArray();
    const auto problo = geom.ProbLoArray();
    
    // Ensure any pending GPU kernels that wrote to phi have completed
    // before reading its values on the host.
    Gpu::streamSynchronize();

    for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const phi_arr = phi.array(mfi);
        
        // Loop over valid cells
#if (AMREX_SPACEDIM == 3)
        for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k) {
#else
        int k = 0;
        {
#endif
            for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
                for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
                    if (phi_arr(i,j,k,0) < 0.0) {
                        // Compute cell center coordinates
                        Real x = problo[0] + (i + Real(0.5)) * dx[0];
                        Real y = problo[1] + (j + Real(0.5)) * dx[1];
                        negative_points.push_back({x, y});
                    }
                }
            }
        }
    }
    
    // Gather all points from all processors (if running in parallel)
#ifdef AMREX_USE_MPI
    std::vector<int> counts(ParallelDescriptor::NProcs());
    int local_count = negative_points.size();
    MPI_Allgather(&local_count, 1, MPI_INT, counts.data(), 1, MPI_INT, 
                  ParallelDescriptor::Communicator());
    
    std::vector<int> displacements(ParallelDescriptor::NProcs());
    int total_count = 0;
    for (int i = 0; i < ParallelDescriptor::NProcs(); ++i) {
        displacements[i] = total_count;
        total_count += counts[i];
        counts[i] *= 2; // Each point has 2 values (x, y)
        displacements[i] *= 2;
    }
    
    std::vector<Real> local_data;
    for (const auto& pt : negative_points) {
        local_data.push_back(pt.first);
        local_data.push_back(pt.second);
    }
    
    std::vector<Real> all_data;
    if (ParallelDescriptor::IOProcessor()) {
        all_data.resize(total_count * 2);
    }
    
    MPI_Gatherv(local_data.data(), local_data.size(), 
                ParallelDescriptor::Mpi_typemap<Real>::type(),
                all_data.data(), counts.data(), displacements.data(), 
                ParallelDescriptor::Mpi_typemap<Real>::type(),
                ParallelDescriptor::IOProcessorNumber(),
                ParallelDescriptor::Communicator());
    
    if (ParallelDescriptor::IOProcessor()) {
        negative_points.clear();
        for (size_t i = 0; i < all_data.size(); i += 2) {
            negative_points.push_back({all_data[i], all_data[i+1]});
        }
    }
#endif
    
    // Write to file (only on IO processor)
    if (ParallelDescriptor::IOProcessor()) {
        std::ofstream outfile(filename);
        if (!outfile.is_open()) {
            amrex::Abort("Cannot open output file: " + filename);
        }
        
        outfile << "# X Y coordinates of cells where phi < 0\n";
        outfile << "# Total points: " << negative_points.size() << "\n";
        
        for (const auto& pt : negative_points) {
            outfile << pt.first << " " << pt.second << "\n";
        }
        
        outfile.close();
        amrex::Print() << "Wrote " << negative_points.size() << " negative phi points to " 
                       << filename << "\n";
    }
}

// ======================= Convex Hull Helper Functions =============================
// Cross product of vectors OA and OB where O is the origin point
static Real cross_product(const std::pair<Real, Real>& O, 
                         const std::pair<Real, Real>& A, 
                         const std::pair<Real, Real>& B)
{
    return (A.first - O.first) * (B.second - O.second) - 
           (A.second - O.second) * (B.first - O.first);
}

// Compute convex hull using Graham scan algorithm
static std::vector<std::pair<Real, Real>> compute_convex_hull(
    std::vector<std::pair<Real, Real>>& points)
{
    if (points.size() < 3) {
        return points;
    }
    
    // Find the point with lowest y-coordinate (and leftmost if tie)
    size_t min_idx = 0;
    for (size_t i = 1; i < points.size(); ++i) {
        if (points[i].second < points[min_idx].second || 
            (points[i].second == points[min_idx].second && 
             points[i].first < points[min_idx].first)) {
            min_idx = i;
        }
    }
    std::swap(points[0], points[min_idx]);
    
    const auto& pivot = points[0];
    
    // Sort points by polar angle with respect to pivot
    std::sort(points.begin() + 1, points.end(), 
        [&pivot](const std::pair<Real, Real>& a, const std::pair<Real, Real>& b) {
            Real cp = cross_product(pivot, a, b);
            if (std::abs(cp) < 1e-10) {
                // Collinear points - sort by distance
                Real dist_a = (a.first - pivot.first) * (a.first - pivot.first) + 
                             (a.second - pivot.second) * (a.second - pivot.second);
                Real dist_b = (b.first - pivot.first) * (b.first - pivot.first) + 
                             (b.second - pivot.second) * (b.second - pivot.second);
                return dist_a < dist_b;
            }
            return cp > 0;
        });
    
    // Build convex hull
    std::vector<std::pair<Real, Real>> hull;
    hull.push_back(points[0]);
    hull.push_back(points[1]);
    
    for (size_t i = 2; i < points.size(); ++i) {
        // Remove points that make clockwise turn
        while (hull.size() > 1 && 
               cross_product(hull[hull.size()-2], hull[hull.size()-1], points[i]) <= 0) {
            hull.pop_back();
        }
        hull.push_back(points[i]);
    }
    
    return hull;
}

// ======================= Write convex hull of negative phi values =================
// Writes the convex hull of all cells where phi < 0 as x-y file.
// For speed, only the leftmost and rightmost negative-phi cell per grid row is
// collected; interior cells cannot lie on the convex hull so discarding them has
// no effect on the result while drastically reducing the point count.
// ==================================================================================
inline void write_negative_phi_convex_hull(const MultiFab& phi, 
                                           const Geometry& geom, 
                                           const std::string& filename)
{
    // For each row j, track the minimum and maximum i index of negative phi cells.
    // Using only these extremal cells reduces the candidate count from O(N_cells)
    // to O(N_rows) before the convex hull computation.
    std::map<int, std::pair<int,int>> row_extremes; // j -> {min_i, max_i}

    const auto dx = geom.CellSizeArray();
    const auto problo = geom.ProbLoArray();

    // Ensure any pending GPU kernels that wrote to phi have completed
    // before reading its values on the host.
    Gpu::streamSynchronize();

    for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const phi_arr = phi.array(mfi);

        // Loop over valid cells
#if (AMREX_SPACEDIM == 3)
        for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k) {
#else
        int k = 0;
        {
#endif
            for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
                int min_i = std::numeric_limits<int>::max();
                int max_i = std::numeric_limits<int>::lowest();
                for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
                    if (phi_arr(i,j,k,0) < 0.0) {
                        if (i < min_i) min_i = i;
                        if (i > max_i) max_i = i;
                    }
                }
                if (min_i != std::numeric_limits<int>::max()) {
                    auto [it, inserted] = row_extremes.try_emplace(j, min_i, max_i);
                    if (!inserted) {
                        if (min_i < it->second.first)  it->second.first  = min_i;
                        if (max_i > it->second.second) it->second.second = max_i;
                    }
                }
            }
        }
    }

    // Convert extremal indices to (x, y) coordinates
    std::vector<std::pair<Real, Real>> negative_points;
    negative_points.reserve(row_extremes.size() * 2);
    for (const auto& kv : row_extremes) {
        int j = kv.first;
        Real y = problo[1] + (j + Real(0.5)) * dx[1];
        Real x_left  = problo[0] + (kv.second.first  + Real(0.5)) * dx[0];
        Real x_right = problo[0] + (kv.second.second + Real(0.5)) * dx[0];
        negative_points.push_back({x_left, y});
        if (kv.second.first != kv.second.second) {
            negative_points.push_back({x_right, y});
        }
    }
    
    // Gather all points from all processors (if running in parallel)
#ifdef AMREX_USE_MPI
    std::vector<int> counts(ParallelDescriptor::NProcs());
    int local_count = negative_points.size();
    MPI_Allgather(&local_count, 1, MPI_INT, counts.data(), 1, MPI_INT, 
                  ParallelDescriptor::Communicator());
    
    std::vector<int> displacements(ParallelDescriptor::NProcs());
    int total_count = 0;
    for (int i = 0; i < ParallelDescriptor::NProcs(); ++i) {
        displacements[i] = total_count;
        total_count += counts[i];
        counts[i] *= 2; // Each point has 2 values (x, y)
        displacements[i] *= 2;
    }
    
    std::vector<Real> local_data;
    for (const auto& pt : negative_points) {
        local_data.push_back(pt.first);
        local_data.push_back(pt.second);
    }
    
    std::vector<Real> all_data;
    if (ParallelDescriptor::IOProcessor()) {
        all_data.resize(total_count * 2);
    }
    
    MPI_Gatherv(local_data.data(), local_data.size(), 
                ParallelDescriptor::Mpi_typemap<Real>::type(),
                all_data.data(), counts.data(), displacements.data(), 
                ParallelDescriptor::Mpi_typemap<Real>::type(),
                ParallelDescriptor::IOProcessorNumber(),
                ParallelDescriptor::Communicator());
    
    if (ParallelDescriptor::IOProcessor()) {
        negative_points.clear();
        for (size_t i = 0; i < all_data.size(); i += 2) {
            negative_points.push_back({all_data[i], all_data[i+1]});
        }
    }
#endif
    
    // Compute and write convex hull (only on IO processor)
    if (ParallelDescriptor::IOProcessor()) {
        if (negative_points.empty()) {
            amrex::Print() << "No negative phi points found. Skipping convex hull output.\n";
            return;
        }
        
        // Compute convex hull
        std::vector<std::pair<Real, Real>> hull = compute_convex_hull(negative_points);
        
        std::ofstream outfile(filename);
        if (!outfile.is_open()) {
            amrex::Abort("Cannot open output file: " + filename);
        }
        
        outfile << "# X Y coordinates of convex hull of phi < 0 region\n";
        outfile << "# Total hull points: " << hull.size() << "\n";
        outfile << "# Total negative phi points: " << negative_points.size() << "\n";
        
        for (const auto& pt : hull) {
            outfile << pt.first << " " << pt.second << "\n";
        }
        
        // Optionally close the hull by adding the first point at the end
        if (!hull.empty()) {
            outfile << hull[0].first << " " << hull[0].second << "\n";
        }
        
        outfile.close();
        amrex::Print() << "Wrote " << hull.size() << " convex hull points (from " 
                       << negative_points.size() << " negative phi points) to " 
                       << filename << "\n";
    }
}

// ============================================================================
// write_fire_perimeter_csv
// ============================================================================
// Writes the fire perimeter cells (cells where phi < 0 with at least one
// phi >= 0 neighbour) as a CSV file with columns:  X, Y
//
// These are the cell-centre coordinates in simulation physical units [m].
// When sorting = true (default) the boundary points are sorted by angle
// around the centroid to produce a traversal-ordered perimeter.
// ============================================================================
inline void write_fire_perimeter_csv(const MultiFab& phi,
                                     const Geometry& geom,
                                     const std::string& filename)
{
    const auto dx     = geom.CellSizeArray();
    const auto problo = geom.ProbLoArray();

    Gpu::streamSynchronize();

    std::vector<std::pair<Real, Real>> perim_pts;

    for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
        const Box& bx  = mfi.validbox();
        auto const phi_arr = phi.array(mfi);
        const IntVect lo   = bx.smallEnd();
        const IntVect hi   = bx.bigEnd();

        amrex::LoopOnCpu(bx, [&](int i, int j, int k) {
            if (phi_arr(i, j, k, 0) >= Real(0.0)) return;
            bool on_perim = false;
            if (i > lo[0]) on_perim |= (phi_arr(i-1, j, k, 0) >= Real(0.0));
            if (i < hi[0]) on_perim |= (phi_arr(i+1, j, k, 0) >= Real(0.0));
            if (j > lo[1]) on_perim |= (phi_arr(i, j-1, k, 0) >= Real(0.0));
            if (j < hi[1]) on_perim |= (phi_arr(i, j+1, k, 0) >= Real(0.0));
            if (on_perim) {
                Real x = problo[0] + (i + Real(0.5)) * dx[0];
                Real y = problo[1] + (j + Real(0.5)) * dx[1];
                perim_pts.push_back({x, y});
            }
        });
    }

#ifdef AMREX_USE_MPI
    {
        std::vector<int> counts(ParallelDescriptor::NProcs());
        int local_n = static_cast<int>(perim_pts.size());
        MPI_Allgather(&local_n, 1, MPI_INT, counts.data(), 1, MPI_INT,
                      ParallelDescriptor::Communicator());
        std::vector<int> displ(ParallelDescriptor::NProcs());
        int total = 0;
        for (int r = 0; r < ParallelDescriptor::NProcs(); ++r) {
            displ[r] = total; total += counts[r];
            counts[r] *= 2;  displ[r] *= 2;
        }
        std::vector<Real> ldata;
        ldata.reserve(2 * local_n);
        for (auto& pt : perim_pts) { ldata.push_back(pt.first); ldata.push_back(pt.second); }
        std::vector<Real> all;
        if (ParallelDescriptor::IOProcessor()) all.resize(total * 2);
        MPI_Gatherv(ldata.data(), ldata.size(),
                    ParallelDescriptor::Mpi_typemap<Real>::type(),
                    all.data(), counts.data(), displ.data(),
                    ParallelDescriptor::Mpi_typemap<Real>::type(),
                    ParallelDescriptor::IOProcessorNumber(),
                    ParallelDescriptor::Communicator());
        if (ParallelDescriptor::IOProcessor()) {
            perim_pts.clear();
            for (size_t k = 0; k + 1 < all.size(); k += 2)
                perim_pts.push_back({all[k], all[k+1]});
        }
    }
#endif

    if (!ParallelDescriptor::IOProcessor()) return;
    if (perim_pts.empty()) {
        amrex::Print() << "write_fire_perimeter_csv: no perimeter cells found\n";
        return;
    }

    // Sort by angle around centroid for a cleaner perimeter polygon
    Real cx = 0, cy = 0;
    for (auto& p : perim_pts) { cx += p.first; cy += p.second; }
    cx /= static_cast<Real>(perim_pts.size());
    cy /= static_cast<Real>(perim_pts.size());
    std::sort(perim_pts.begin(), perim_pts.end(),
        [cx, cy](const std::pair<Real,Real>& a, const std::pair<Real,Real>& b) {
            return std::atan2(a.second - cy, a.first - cx) <
                   std::atan2(b.second - cy, b.first - cx);
        });

    std::ofstream out(filename);
    if (!out.is_open()) amrex::Abort("Cannot open perimeter CSV: " + filename);
    out << "# Fire perimeter cell centres (phi < 0, boundary cells)\n";
    out << "X,Y\n";
    for (auto& p : perim_pts)
        out << p.first << "," << p.second << "\n";
    out.close();
    amrex::Print() << "Wrote " << perim_pts.size()
                   << " perimeter points to " << filename << "\n";
}

// ============================================================================
// write_fire_perimeter_geojson
// ============================================================================
// Writes the fire perimeter as a GeoJSON FeatureCollection containing:
//   - A Polygon feature whose coordinates are the angle-sorted boundary cells
//     (same set as write_fire_perimeter_csv, first = last to close the ring).
//   - Properties: step, time_s, n_perimeter_cells.
//
// Coordinates are in simulation physical units [m] (easting/northing when the
// domain is in a UTM projection, or arbitrary for idealised runs).
// To use with mapping tools, specify the CRS separately or pre-project.
// ============================================================================
inline void write_fire_perimeter_geojson(const MultiFab& phi,
                                         const Geometry& geom,
                                         const std::string& filename,
                                         int  step,
                                         Real time_s)
{
    const auto dx     = geom.CellSizeArray();
    const auto problo = geom.ProbLoArray();

    Gpu::streamSynchronize();

    std::vector<std::pair<Real, Real>> perim_pts;

    for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
        const Box& bx  = mfi.validbox();
        auto const phi_arr = phi.array(mfi);
        const IntVect lo   = bx.smallEnd();
        const IntVect hi   = bx.bigEnd();

        amrex::LoopOnCpu(bx, [&](int i, int j, int k) {
            if (phi_arr(i, j, k, 0) >= Real(0.0)) return;
            bool on_perim = false;
            if (i > lo[0]) on_perim |= (phi_arr(i-1, j, k, 0) >= Real(0.0));
            if (i < hi[0]) on_perim |= (phi_arr(i+1, j, k, 0) >= Real(0.0));
            if (j > lo[1]) on_perim |= (phi_arr(i, j-1, k, 0) >= Real(0.0));
            if (j < hi[1]) on_perim |= (phi_arr(i, j+1, k, 0) >= Real(0.0));
            if (on_perim) {
                Real x = problo[0] + (i + Real(0.5)) * dx[0];
                Real y = problo[1] + (j + Real(0.5)) * dx[1];
                perim_pts.push_back({x, y});
            }
        });
    }

#ifdef AMREX_USE_MPI
    {
        std::vector<int> counts(ParallelDescriptor::NProcs());
        int local_n = static_cast<int>(perim_pts.size());
        MPI_Allgather(&local_n, 1, MPI_INT, counts.data(), 1, MPI_INT,
                      ParallelDescriptor::Communicator());
        std::vector<int> displ(ParallelDescriptor::NProcs());
        int total = 0;
        for (int r = 0; r < ParallelDescriptor::NProcs(); ++r) {
            displ[r] = total; total += counts[r];
            counts[r] *= 2;  displ[r] *= 2;
        }
        std::vector<Real> ldata;
        ldata.reserve(2 * local_n);
        for (auto& pt : perim_pts) { ldata.push_back(pt.first); ldata.push_back(pt.second); }
        std::vector<Real> all;
        if (ParallelDescriptor::IOProcessor()) all.resize(total * 2);
        MPI_Gatherv(ldata.data(), ldata.size(),
                    ParallelDescriptor::Mpi_typemap<Real>::type(),
                    all.data(), counts.data(), displ.data(),
                    ParallelDescriptor::Mpi_typemap<Real>::type(),
                    ParallelDescriptor::IOProcessorNumber(),
                    ParallelDescriptor::Communicator());
        if (ParallelDescriptor::IOProcessor()) {
            perim_pts.clear();
            for (size_t k = 0; k + 1 < all.size(); k += 2)
                perim_pts.push_back({all[k], all[k+1]});
        }
    }
#endif

    if (!ParallelDescriptor::IOProcessor()) return;
    if (perim_pts.empty()) {
        amrex::Print() << "write_fire_perimeter_geojson: no perimeter cells\n";
        return;
    }

    // Sort by angle around centroid
    Real cx = 0, cy = 0;
    for (auto& p : perim_pts) { cx += p.first; cy += p.second; }
    cx /= static_cast<Real>(perim_pts.size());
    cy /= static_cast<Real>(perim_pts.size());
    std::sort(perim_pts.begin(), perim_pts.end(),
        [cx, cy](const std::pair<Real,Real>& a, const std::pair<Real,Real>& b) {
            return std::atan2(a.second - cy, a.first - cx) <
                   std::atan2(b.second - cy, b.first - cx);
        });

    std::ofstream out(filename);
    if (!out.is_open()) amrex::Abort("Cannot open perimeter GeoJSON: " + filename);

    out << "{\n";
    out << "  \"type\": \"FeatureCollection\",\n";
    out << "  \"features\": [\n";
    out << "    {\n";
    out << "      \"type\": \"Feature\",\n";
    out << "      \"properties\": {\n";
    out << "        \"step\": " << step << ",\n";
    out << "        \"time_s\": " << time_s << ",\n";
    out << "        \"n_perimeter_cells\": " << perim_pts.size() << "\n";
    out << "      },\n";
    out << "      \"geometry\": {\n";
    out << "        \"type\": \"Polygon\",\n";
    out << "        \"coordinates\": [[\n";

    // Write ring: close by repeating first point
    for (size_t k = 0; k < perim_pts.size(); ++k) {
        out << "          [" << perim_pts[k].first
            << ", "          << perim_pts[k].second << "]";
        out << ",\n";
    }
    // Close the ring
    out << "          [" << perim_pts[0].first
        << ", "          << perim_pts[0].second << "]\n";
    out << "        ]]\n";
    out << "      }\n";
    out << "    }\n";
    out << "  ]\n";
    out << "}\n";
    out.close();
    amrex::Print() << "Wrote fire perimeter GeoJSON (" << perim_pts.size()
                   << " points) to " << filename << "\n";
}

// ============================================================================
// UTM → WGS-84 conversion (inverse Transverse Mercator)
//
// Converts UTM easting/northing [m] to WGS-84 geographic coordinates.
// Returns {lat_deg, lon_deg}.
//
// zone_number : UTM zone number 1-60.
// northern    : true for northern hemisphere, false for southern.
//
// Based on the standard WGS-84 inverse TM algorithm (Redfearn 1948 series).
// Accurate to < 1 mm within a UTM zone.
// ============================================================================
namespace {

inline std::pair<double,double> utm_to_latlon_wgs84(
    double easting, double northing, int zone_number, bool northern)
{
    // WGS-84 ellipsoid constants — all double precision for accuracy
    constexpr double WGS84_PI = 3.14159265358979323846;
    constexpr double a    = 6378137.0;
    constexpr double f    = 1.0 / 298.257223563;
    constexpr double b    = a * (1.0 - f);
    constexpr double e2   = 1.0 - (b / a) * (b / a);
    constexpr double ep2  = e2 / (1.0 - e2);
    constexpr double k0   = 0.9996;
    constexpr double E0   = 500000.0;
    constexpr double N0S  = 10000000.0;  // southern-hemisphere false northing

    const double x = easting - E0;
    const double y = northern ? northing : northing - N0S;

    const double e4 = e2 * e2;
    const double e6 = e4 * e2;

    // Compute footpoint latitude
    const double m  = y / k0;
    const double mu = m / (a * (1.0 - e2 / 4.0 - 3.0 * e4 / 64.0 - 5.0 * e6 / 256.0));

    const double e1 = (1.0 - std::sqrt(1.0 - e2)) / (1.0 + std::sqrt(1.0 - e2));
    const double lat1 = mu
        + (3.0 * e1 / 2.0 - 27.0 * e1 * e1 * e1 / 32.0) * std::sin(2.0 * mu)
        + (21.0 * e1 * e1 / 16.0 - 55.0 * e1 * e1 * e1 * e1 / 32.0) * std::sin(4.0 * mu)
        + (151.0 * e1 * e1 * e1 / 96.0) * std::sin(6.0 * mu)
        + (1097.0 * e1 * e1 * e1 * e1 / 512.0) * std::sin(8.0 * mu);

    const double sinL1 = std::sin(lat1);
    const double cosL1 = std::cos(lat1);
    const double tanL1 = sinL1 / cosL1;

    const double N1 = a / std::sqrt(1.0 - e2 * sinL1 * sinL1);
    const double T1 = tanL1 * tanL1;
    const double C1 = ep2 * cosL1 * cosL1;
    const double R1 = a * (1.0 - e2) / std::pow(1.0 - e2 * sinL1 * sinL1, 1.5);
    const double D  = x / (N1 * k0);

    const double lat = lat1 - (N1 * tanL1 / R1) * (
          D * D / 2.0
        - (5.0 + 3.0 * T1 + 10.0 * C1 - 4.0 * C1 * C1 - 9.0 * ep2)
          * D * D * D * D / 24.0
        + (61.0 + 90.0 * T1 + 298.0 * C1 + 45.0 * T1 * T1 - 252.0 * ep2
           - 3.0 * C1 * C1)
          * D * D * D * D * D * D / 720.0);

    const double lon0 = ((zone_number - 1) * 6.0 - 180.0 + 3.0) * (WGS84_PI / 180.0);
    const double lon  = lon0 + (
          D
        - (1.0 + 2.0 * T1 + C1) * D * D * D / 6.0
        + (5.0 - 2.0 * C1 + 28.0 * T1 - 3.0 * C1 * C1 + 8.0 * ep2 + 24.0 * T1 * T1)
          * D * D * D * D * D / 120.0) / cosL1;

    return {lat * 180.0 / WGS84_PI, lon * 180.0 / WGS84_PI};
}

} // anonymous namespace

// ============================================================================
// write_fire_perimeter_kml
// ============================================================================
// Writes the fire perimeter as a KML (Keyhole Markup Language) document
// suitable for opening in Google Earth or other GIS tools that support KML.
// Coordinates are optionally converted from UTM [m] to WGS-84 lon/lat using
// the utm_zone and utm_northern parameters.
//
// Parameters:
//   phi           – level-set field (perimeter cells are phi < 0 on boundary)
//   geom          – AMReX Geometry
//   filename      – output file path (e.g. "perimeter_0010.kml")
//   step          – timestep index (written as KML description)
//   time_s        – simulation time [s] (written as KML description)
//   utm_zone      – UTM zone number 1-60; 0 = write raw UTM coords (not WGS-84)
//   utm_northern  – 1 = northern hemisphere, 0 = southern
// ============================================================================
inline void write_fire_perimeter_kml(const MultiFab& phi,
                                     const Geometry& geom,
                                     const std::string& filename,
                                     int  step,
                                     Real time_s,
                                     int  utm_zone,
                                     int  utm_northern)
{
    const auto dx     = geom.CellSizeArray();
    const auto problo = geom.ProbLoArray();

    Gpu::streamSynchronize();

    std::vector<std::pair<Real, Real>> perim_pts;

    for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
        const Box& bx  = mfi.validbox();
        auto const phi_arr = phi.array(mfi);
        const IntVect lo   = bx.smallEnd();
        const IntVect hi   = bx.bigEnd();

        amrex::LoopOnCpu(bx, [&](int i, int j, int k) {
            if (phi_arr(i, j, k, 0) >= Real(0.0)) return;
            bool on_perim = false;
            if (i > lo[0]) on_perim |= (phi_arr(i-1, j, k, 0) >= Real(0.0));
            if (i < hi[0]) on_perim |= (phi_arr(i+1, j, k, 0) >= Real(0.0));
            if (j > lo[1]) on_perim |= (phi_arr(i, j-1, k, 0) >= Real(0.0));
            if (j < hi[1]) on_perim |= (phi_arr(i, j+1, k, 0) >= Real(0.0));
            if (on_perim) {
                Real x = problo[0] + (i + Real(0.5)) * dx[0];
                Real y = problo[1] + (j + Real(0.5)) * dx[1];
                perim_pts.push_back({x, y});
            }
        });
    }

#ifdef AMREX_USE_MPI
    {
        std::vector<int> counts(ParallelDescriptor::NProcs());
        int local_n = static_cast<int>(perim_pts.size());
        MPI_Allgather(&local_n, 1, MPI_INT, counts.data(), 1, MPI_INT,
                      ParallelDescriptor::Communicator());
        std::vector<int> displ(ParallelDescriptor::NProcs());
        int total = 0;
        for (int r = 0; r < ParallelDescriptor::NProcs(); ++r) {
            displ[r] = total; total += counts[r];
            counts[r] *= 2;  displ[r] *= 2;
        }
        std::vector<Real> ldata;
        ldata.reserve(2 * local_n);
        for (auto& pt : perim_pts) { ldata.push_back(pt.first); ldata.push_back(pt.second); }
        std::vector<Real> all;
        if (ParallelDescriptor::IOProcessor()) all.resize(total * 2);
        MPI_Gatherv(ldata.data(), ldata.size(),
                    ParallelDescriptor::Mpi_typemap<Real>::type(),
                    all.data(), counts.data(), displ.data(),
                    ParallelDescriptor::Mpi_typemap<Real>::type(),
                    ParallelDescriptor::IOProcessorNumber(),
                    ParallelDescriptor::Communicator());
        if (ParallelDescriptor::IOProcessor()) {
            perim_pts.clear();
            for (size_t k = 0; k + 1 < all.size(); k += 2)
                perim_pts.push_back({all[k], all[k+1]});
        }
    }
#endif

    if (!ParallelDescriptor::IOProcessor()) return;
    if (perim_pts.empty()) {
        amrex::Print() << "write_fire_perimeter_kml: no perimeter cells\n";
        return;
    }

    // Sort by angle around centroid
    Real cx = 0, cy = 0;
    for (auto& p : perim_pts) { cx += p.first; cy += p.second; }
    cx /= static_cast<Real>(perim_pts.size());
    cy /= static_cast<Real>(perim_pts.size());
    std::sort(perim_pts.begin(), perim_pts.end(),
        [cx, cy](const std::pair<Real,Real>& a, const std::pair<Real,Real>& b) {
            return std::atan2(a.second - cy, a.first - cx) <
                   std::atan2(b.second - cy, b.first - cx);
        });

    // Convert UTM → WGS-84 when zone is specified
    const bool do_convert = (utm_zone >= 1 && utm_zone <= 60);
    const bool northern   = (utm_northern != 0);

    std::ofstream out(filename);
    if (!out.is_open()) amrex::Abort("Cannot open perimeter KML: " + filename);

    out << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    out << "<kml xmlns=\"http://www.opengis.net/kml/2.2\">\n";
    out << "<Document>\n";
    out << "  <name>Fire Perimeter</name>\n";
    out << "  <description>Step " << step
        << " | Time " << std::fixed << std::setprecision(1)
        << static_cast<double>(time_s) / 3600.0 << " h</description>\n";
    out << "  <Style id=\"firePerimeter\">\n";
    out << "    <LineStyle><color>ff0000ff</color><width>3</width></LineStyle>\n";
    out << "    <PolyStyle><color>660000ff</color></PolyStyle>\n";
    out << "  </Style>\n";
    out << "  <Placemark>\n";
    out << "    <styleUrl>#firePerimeter</styleUrl>\n";
    out << "    <name>Fire Perimeter t=" << std::fixed << std::setprecision(2)
        << static_cast<double>(time_s) / 3600.0 << " h</name>\n";
    out << "    <Polygon><outerBoundaryIs><LinearRing><coordinates>\n";

    for (size_t k = 0; k < perim_pts.size(); ++k) {
        double lon, lat;
        if (do_convert) {
            auto [la, lo] = utm_to_latlon_wgs84(
                static_cast<double>(perim_pts[k].first),
                static_cast<double>(perim_pts[k].second),
                utm_zone, northern);
            lat = la; lon = lo;
        } else {
            lon = static_cast<double>(perim_pts[k].first);
            lat = static_cast<double>(perim_pts[k].second);
        }
        out << "      " << std::fixed << std::setprecision(7)
            << lon << "," << lat << ",0\n";
    }
    // Close the ring with the first point
    {
        double lon, lat;
        if (do_convert) {
            auto [la, lo] = utm_to_latlon_wgs84(
                static_cast<double>(perim_pts[0].first),
                static_cast<double>(perim_pts[0].second),
                utm_zone, northern);
            lat = la; lon = lo;
        } else {
            lon = static_cast<double>(perim_pts[0].first);
            lat = static_cast<double>(perim_pts[0].second);
        }
        out << "      " << std::fixed << std::setprecision(7)
            << lon << "," << lat << ",0\n";
    }

    out << "    </coordinates></LinearRing></outerBoundaryIs></Polygon>\n";
    out << "  </Placemark>\n";
    out << "</Document>\n";
    out << "</kml>\n";
    out.close();
    amrex::Print() << "Wrote fire perimeter KML (" << perim_pts.size()
                   << " points) to " << filename << "\n";
}

// ============================================================================
// write_fire_stats_header
//
// Writes the CSV header row for the burned area / perimeter time series.
// Should be called once at simulation startup before the first append.
// ============================================================================
inline void write_fire_stats_header(const std::string& filename)
{
    if (!ParallelDescriptor::IOProcessor()) return;
    std::ofstream out(filename, std::ios::out | std::ios::trunc);
    if (!out.is_open()) {
        amrex::Print() << "WARNING: Cannot open fire_stats file: " << filename << "\n";
        return;
    }
    out << "step,time_s,burned_area_ha,perimeter_km,"
           "active_front_cells,head_ros_ms,max_flame_length_m,"
           "total_co2_kg_m2,total_co_kg_m2,total_pm25_kg_m2\n";
    out.close();
}

// ============================================================================
// append_fire_stats
//
// Computes burned area [ha], perimeter length [km], active front cell count,
// head ROS [m/s] (maximum ROS in perimeter cells), maximum flame length [m],
// and domain-mean emission totals, then appends one row to the CSV file.
//
// Burned area and perimeter are based on arrival_time_mf >= 0 (every cell
// that has ever been ignited, including the burned interior).  In the FARSITE
// propagation mode phi is reset to 0 after each step so phi < 0 only marks
// the thin current fire-front ring; arrival_time correctly accumulates the
// full burned region for all propagation modes.
//
// emissions_mf is optional (pass nullptr to skip emission stats).
// R_mf is optional (pass nullptr to skip head ROS tracking).
// flame_length_mf is optional (pass nullptr to skip max flame length tracking).
// ============================================================================

// ============================================================================
// Helper: Calculate fire ellipse major and minor axes from perimeter and area
// ============================================================================
// Approximates the fire as an ellipse given its perimeter and burned area.
// Uses the relationships:
//   - Area = π * a * b  (where a = semi-major axis, b = semi-minor axis)
//   - Perimeter ≈ π * (a + b)  (Ramanujan's first approximation)
//
// Returns {major_axis_m, minor_axis_m} (full axes, not semi-axes)
// ============================================================================
inline std::pair<Real, Real> compute_fire_ellipse_axes(Real area_m2, Real perimeter_m)
{
    // Avoid division by zero
    if (perimeter_m <= 0.0 || area_m2 <= 0.0) {
        return {0.0, 0.0};
    }

    constexpr Real PI = M_PI;
    
    // From Area = π * a * b  and  Perimeter ≈ π * (a + b)
    // we get: a + b ≈ Perimeter / π
    // and:    a * b = Area / π
    // These are the sum and product of the semi-axes.
    // Solving the quadratic: x² - (a+b)x + ab = 0
    
    Real sum_ab    = perimeter_m / PI;
    Real prod_ab   = area_m2 / PI;
    
    // Clamp to avoid numerical issues
    if (sum_ab <= 0.0 || prod_ab <= 0.0) {
        return {0.0, 0.0};
    }
    
    // Solve for semi-major and semi-minor axes using quadratic formula
    Real discriminant = sum_ab * sum_ab - 4.0 * prod_ab;
    if (discriminant < 0.0) {
        // No real solution; return a circle approximation
        Real radius = std::sqrt(area_m2 / PI);
        return {2.0 * radius, 2.0 * radius};
    }
    
    Real sqrt_disc = std::sqrt(discriminant);
    Real semi_a = (sum_ab + sqrt_disc) / 2.0;  // semi-major axis
    Real semi_b = (sum_ab - sqrt_disc) / 2.0;  // semi-minor axis
    
    // Return full axes (2 * semi-axes)
    return {2.0 * semi_a, 2.0 * semi_b};
}

// ============================================================================
// Helper: Read previous area and time from fire_stats CSV for growth rate
// ============================================================================
// Reads the last row of the CSV file to get the previous burned area and time.
// Returns {prev_area_ha, prev_time_s} or {0, 0} if file doesn't exist or is empty.
// ============================================================================
inline std::pair<Real, Real> get_previous_fire_stats(const std::string& filename)
{
    std::ifstream in(filename);
    if (!in.is_open()) {
        return {0.0, 0.0};
    }
    
    std::string line;
    Real last_time = 0.0;
    Real last_area = 0.0;
    bool found = false;
    
    // Skip header line
    std::getline(in, line);
    
    // Read through all lines to get the last one
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        
        // Parse: step, time_s, burned_area_ha, ...
        std::istringstream iss(line);
        int step_tmp;
        Real time_tmp, area_tmp;
        char comma;
        
        if (iss >> step_tmp >> comma >> time_tmp >> comma >> area_tmp) {
            last_time = time_tmp;
            last_area = area_tmp;
            found = true;
        }
    }
    in.close();
    
    return found ? std::make_pair(last_area, last_time) : std::make_pair(0.0, 0.0);
}

// ============================================================================
// Write initial header for fire_stats.csv if file is empty
// ============================================================================
inline void write_fire_stats_header_if_needed(const std::string& filename)
{
    std::ifstream check(filename);
    bool file_exists = check.good();
    check.close();
    
    if (!file_exists) {
        std::ofstream out(filename, std::ios::out);
        if (out.is_open()) {
            out << "step,time_s,burned_area_ha,perimeter_km,active_front_cells,"
                << "head_ros_ms,max_flame_length_m,total_co2_kg_m2,total_co_kg_m2,"
                << "total_pm25_kg_m2,growth_rate_ha_min,major_axis_m,minor_axis_m\n";
            out.close();
        }
    }
}

inline void append_fire_stats(const MultiFab& phi,
                               const Geometry& geom,
                               const MultiFab* emissions_mf,
                               int step, Real time_s,
                               const std::string& filename,
                               const MultiFab& arrival_time_mf,
                               const MultiFab* R_mf = nullptr,
                               const MultiFab* flame_length_mf = nullptr)
{
    // Ensure CSV header exists
    if (step == 0 || !ParallelDescriptor::IOProcessor()) {
        if (step == 0 && ParallelDescriptor::IOProcessor()) {
            write_fire_stats_header_if_needed(filename);
        }
        ParallelDescriptor::Barrier();
    }
    const auto dxx = geom.CellSize();
    Real cell_area = dxx[0] * dxx[1];

    // Burned area and perimeter use arrival_time_mf >= 0 (cumulative burned
    // region) via the edge-crossing method: for every cell face separating a
    // burned cell (at >= 0) from an unburned cell (at < 0), the face length
    // is added to the total perimeter.  phi is still used for n_active (cells
    // near the live fire front).
    amrex::Long n_burned = 0, n_active = 0;
    Real perim_m = Real(0.0);
    Real sum_co2 = 0, sum_co = 0, sum_pm25 = 0;
    Real head_ros = Real(0.0);  // maximum ROS in perimeter cells [m/s]
    Real max_fl = Real(0.0);     // maximum flame length [m]

    for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const p  = phi.const_array(mfi);
        auto const at = arrival_time_mf.const_array(mfi);
        const IntVect lo = bx.smallEnd();
        const IntVect hi = bx.bigEnd();

        // Head ROS: optional per-cell R_mf array for this MFIter
        Array4<const Real> R_arr;
        if (R_mf) R_arr = R_mf->const_array(mfi);
        
        // Flame length: optional per-cell flame_length_mf array
        Array4<const Real> FL_arr;
        if (flame_length_mf) FL_arr = flame_length_mf->const_array(mfi);

        amrex::LoopOnCpu(bx, [&](int i, int j, int k) {
            bool burned = (at(i, j, k, 0) >= Real(0.0));
            Real pv     = p(i, j, k, 0);
            if (burned) {
                ++n_burned;
                // Track head ROS and flame length for perimeter cells
                // (any burned cell adjacent to an unburned neighbor qualifies as a perimeter cell).
                bool on_perim = false;
                if (i > lo[0]) on_perim |= (at(i-1,j,k,0) < Real(0.0));
                if (i < hi[0]) on_perim |= (at(i+1,j,k,0) < Real(0.0));
                if (j > lo[1]) on_perim |= (at(i,j-1,k,0) < Real(0.0));
                if (j < hi[1]) on_perim |= (at(i,j+1,k,0) < Real(0.0));
                
                if (on_perim) {
                    if (R_mf) {
                        Real r = R_arr(i, j, k, 0);
                        if (r > head_ros) head_ros = r;
                    }
                    if (flame_length_mf) {
                        Real fl = FL_arr(i, j, k, 0);
                        if (fl > max_fl) max_fl = fl;
                    }
                }
            } else if (pv < Real(0.1) && pv >= Real(0.0)) {
                ++n_active;  // unburned cells within phi_thresh of the front
            }
            // Edge-crossing perimeter (arrival_time boundary):
            // x-face between (i,j) and (i+1,j): face length = dxx[1]
            if (i < hi[0]) {
                bool burned_n = (at(i+1, j, k, 0) >= Real(0.0));
                if (burned != burned_n) perim_m += dxx[1];
            }
            // y-face between (i,j) and (i,j+1): face length = dxx[0]
            if (j < hi[1]) {
                bool burned_n = (at(i, j+1, k, 0) >= Real(0.0));
                if (burned != burned_n) perim_m += dxx[0];
            }
        });

        if (emissions_mf) {
            auto const em = emissions_mf->const_array(mfi);
            amrex::LoopOnCpu(bx, [&](int i, int j, int k) {
                sum_co2  += em(i, j, k, 0);
                sum_co   += em(i, j, k, 1);
                sum_pm25 += em(i, j, k, 2);
            });
        }
    }

    ParallelDescriptor::ReduceLongSum(n_burned);
    ParallelDescriptor::ReduceLongSum(n_active);
    ParallelDescriptor::ReduceRealSum(perim_m);
    if (R_mf) {
        ParallelDescriptor::ReduceRealMax(head_ros);
    }
    if (flame_length_mf) {
        ParallelDescriptor::ReduceRealMax(max_fl);
    }
    if (emissions_mf) {
        ParallelDescriptor::ReduceRealSum(sum_co2);
        ParallelDescriptor::ReduceRealSum(sum_co);
        ParallelDescriptor::ReduceRealSum(sum_pm25);
    }

    Real burned_ha    = static_cast<Real>(n_burned) * cell_area / Real(1.0e4);
    Real perim_km     = perim_m / Real(1000.0);
    // Convert domain-sum emission [kg/m²·cell] to mean [kg/m²]
    amrex::Long total_cells  = geom.Domain().numPts();
    Real mean_co2  = (total_cells > 0) ? sum_co2  / static_cast<Real>(total_cells) : Real(0.0);
    Real mean_co   = (total_cells > 0) ? sum_co   / static_cast<Real>(total_cells) : Real(0.0);
    Real mean_pm25 = (total_cells > 0) ? sum_pm25 / static_cast<Real>(total_cells) : Real(0.0);
    
    // Calculate fire ellipse major/minor axes from area and perimeter
    Real area_m2 = burned_ha * Real(1.0e4);  // convert ha to m²
    auto [major_axis_m, minor_axis_m] = compute_fire_ellipse_axes(area_m2, perim_m);
    
    // Calculate growth rate (ha/min)
    Real growth_rate_ha_min = 0.0;
    if (step > 0 && time_s > 0.0) {
        auto [prev_area_ha, prev_time_s] = get_previous_fire_stats(filename);
        Real dt_min = (time_s - prev_time_s) / 60.0;
        if (dt_min > 1e-6) {  // avoid division by very small numbers
            growth_rate_ha_min = (burned_ha - prev_area_ha) / dt_min;
        }
    }

    if (!ParallelDescriptor::IOProcessor()) return;

    std::ofstream out(filename, std::ios::out | std::ios::app);
    if (!out.is_open()) {
        amrex::Print() << "WARNING: Cannot append to fire_stats file: " << filename << "\n";
        return;
    }
    out << step       << ","
        << time_s     << ","
        << burned_ha  << ","
        << perim_km   << ","
        << n_active   << ","
        << head_ros   << ","
        << max_fl     << ","
        << mean_co2   << ","
        << mean_co    << ","
        << mean_pm25  << ","
        << growth_rate_ha_min << ","
        << major_axis_m << ","
        << minor_axis_m << "\n";
    out.close();
}

// ============================================================================
// FARSITE Fire Spread Atlas (.fsa) – Cumulative Perimeter Archive
// ============================================================================
// Appends a time-stamped fire perimeter snapshot to a single cumulative file.
// This mirrors FARSITE's Fire Spread Atlas which collects all perimeter
// polygons across the simulation into one file for post-processing.
//
// File format (text, human-readable):
//   # Fire Spread Atlas – wildfire_levelset
//   # Each snapshot block:
//   SNAPSHOT step  time_s  n_points
//   x1 y1
//   x2 y2
//   ...
//   END_SNAPSHOT
//
// Call write_fsa_header() once at startup, then append_fire_spread_atlas()
// at each plotfile step.
//
// References:
//   Finney, M.A. (1998). FARSITE: Fire Area Simulator. USDA FS RMRS-RP-4.
// ============================================================================

inline void write_fsa_header(const std::string& filename)
{
    if (!ParallelDescriptor::IOProcessor()) return;
    if (filename.empty()) return;
    std::ofstream out(filename, std::ios::out | std::ios::trunc);
    if (!out.is_open()) {
        amrex::Print() << "WARNING: Cannot open .fsa file: " << filename << "\n";
        return;
    }
    out << "# Fire Spread Atlas - wildfire_levelset\n"
        << "# Each snapshot: SNAPSHOT <step> <time_s> <n_points>\n"
        << "#                <x> <y>  (one pair per line)\n"
        << "#                END_SNAPSHOT\n";
    out.close();
}

inline void append_fire_spread_atlas(const MultiFab& phi,
                                     const Geometry& geom,
                                     int step, Real time_s,
                                     const std::string& filename)
{
    if (filename.empty()) return;

    const auto dx     = geom.CellSizeArray();
    const auto problo = geom.ProbLoArray();

    Gpu::streamSynchronize();

    std::vector<std::pair<Real, Real>> perim_pts;

    for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
        const Box& bx  = mfi.validbox();
        auto const phi_arr = phi.array(mfi);
        const IntVect lo   = bx.smallEnd();
        const IntVect hi   = bx.bigEnd();

        amrex::LoopOnCpu(bx, [&](int i, int j, int k) {
            if (phi_arr(i, j, k, 0) >= Real(0.0)) return;
            bool on_perim = false;
            if (i > lo[0]) on_perim |= (phi_arr(i-1, j, k, 0) >= Real(0.0));
            if (i < hi[0]) on_perim |= (phi_arr(i+1, j, k, 0) >= Real(0.0));
            if (j > lo[1]) on_perim |= (phi_arr(i, j-1, k, 0) >= Real(0.0));
            if (j < hi[1]) on_perim |= (phi_arr(i, j+1, k, 0) >= Real(0.0));
            if (on_perim) {
                Real x = problo[0] + (i + Real(0.5)) * dx[0];
                Real y = problo[1] + (j + Real(0.5)) * dx[1];
                perim_pts.push_back({x, y});
            }
        });
    }

#ifdef AMREX_USE_MPI
    {
        int local_n = static_cast<int>(perim_pts.size());
        std::vector<int> counts(ParallelDescriptor::NProcs());
        MPI_Allgather(&local_n, 1, MPI_INT, counts.data(), 1, MPI_INT,
                      ParallelDescriptor::Communicator());
        std::vector<int> displ(ParallelDescriptor::NProcs());
        int total = 0;
        for (int r = 0; r < ParallelDescriptor::NProcs(); ++r) {
            displ[r] = total; total += counts[r];
            counts[r] *= 2;  displ[r] *= 2;
        }
        std::vector<Real> ldata;
        ldata.reserve(2 * local_n);
        for (auto& pt : perim_pts) { ldata.push_back(pt.first); ldata.push_back(pt.second); }
        std::vector<Real> all;
        if (ParallelDescriptor::IOProcessor()) all.resize(total * 2);
        MPI_Gatherv(ldata.data(), (int)ldata.size(),
                    ParallelDescriptor::Mpi_typemap<Real>::type(),
                    all.data(), counts.data(), displ.data(),
                    ParallelDescriptor::Mpi_typemap<Real>::type(),
                    ParallelDescriptor::IOProcessorNumber(),
                    ParallelDescriptor::Communicator());
        if (ParallelDescriptor::IOProcessor()) {
            perim_pts.clear();
            for (size_t kk = 0; kk + 1 < all.size(); kk += 2)
                perim_pts.push_back({all[kk], all[kk+1]});
        }
    }
#endif

    if (!ParallelDescriptor::IOProcessor()) return;

    // Sort by angle around centroid for a consistent ordering
    if (!perim_pts.empty()) {
        Real cx = 0, cy = 0;
        for (auto& p : perim_pts) { cx += p.first; cy += p.second; }
        cx /= static_cast<Real>(perim_pts.size());
        cy /= static_cast<Real>(perim_pts.size());
        std::sort(perim_pts.begin(), perim_pts.end(),
            [cx, cy](const std::pair<Real,Real>& a, const std::pair<Real,Real>& b) {
                return std::atan2(a.second - cy, a.first - cx) <
                       std::atan2(b.second - cy, b.first - cx);
            });
    }

    std::ofstream out(filename, std::ios::out | std::ios::app);
    if (!out.is_open()) {
        amrex::Print() << "WARNING: Cannot append to .fsa file: " << filename << "\n";
        return;
    }
    out << "SNAPSHOT " << step << " " << time_s << " " << perim_pts.size() << "\n";
    for (auto& p : perim_pts)
        out << std::fixed << std::setprecision(2) << p.first << " " << p.second << "\n";
    out << "END_SNAPSHOT\n";
    out.close();
}

// ============================================================================
// FARSITE Post-processing Statistics (.pst)
// ============================================================================
// Writes per-step fire statistics to a FARSITE-compatible .pst text file.
// Columns extend the fire_stats CSV with maximum fireline intensity and
// maximum flame length observed at the fire front, and maximum spotting
// distance when spotting is active.
//
// File format:
//   # FARSITE Post-processing Statistics – wildfire_levelset
//   # step time_h burned_ha perimeter_km active_cells max_intensity_kW_m
//   #       max_flame_len_m max_spot_dist_m
//
// References:
//   Finney, M.A. (1998). FARSITE: Fire Area Simulator. USDA FS RMRS-RP-4.
// ============================================================================

inline void write_pst_header(const std::string& filename)
{
    if (!ParallelDescriptor::IOProcessor()) return;
    if (filename.empty()) return;
    std::ofstream out(filename, std::ios::out | std::ios::trunc);
    if (!out.is_open()) {
        amrex::Print() << "WARNING: Cannot open .pst file: " << filename << "\n";
        return;
    }
    out << "# FARSITE Post-processing Statistics - wildfire_levelset\n";
    out << "step,time_h,burned_ha,perimeter_km,active_cells,"
           "max_fireline_intensity_kW_m,max_flame_length_m,max_spotting_dist_m\n";
    out.close();
}

inline void append_pst_stats(const MultiFab& phi,
                              const Geometry& geom,
                              const MultiFab& arrival_time_mf,
                              const MultiFab* fireline_intensity_mf,
                              const MultiFab* flame_length_mf,
                              const MultiFab* spotting_data_mf,
                              int step, Real time_s,
                              const std::string& filename)
{
    if (filename.empty()) return;

    const auto dxx = geom.CellSize();
    Real cell_area = dxx[0] * dxx[1];

    amrex::Long n_burned = 0, n_active = 0;
    Real perim_m = Real(0.0);
    Real max_fi  = Real(0.0);  // max fireline intensity [kW/m]
    Real max_fl  = Real(0.0);  // max flame length [m]
    Real max_sd  = Real(0.0);  // max spotting distance [m]

    for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const at  = arrival_time_mf.const_array(mfi);
        auto const pv  = phi.const_array(mfi);
        const IntVect hi = bx.bigEnd();

        Array4<const Real> fi_arr, fl_arr, sp_arr;
        if (fireline_intensity_mf) fi_arr = fireline_intensity_mf->const_array(mfi);
        if (flame_length_mf)       fl_arr = flame_length_mf->const_array(mfi);
        if (spotting_data_mf)      sp_arr = spotting_data_mf->const_array(mfi);

        amrex::LoopOnCpu(bx, [&](int i, int j, int k) {
            bool burned = (at(i, j, k, 0) >= Real(0.0));
            if (burned) {
                ++n_burned;
                if (fireline_intensity_mf) {
                    Real v = fi_arr(i, j, k, 0);
                    if (v > max_fi) max_fi = v;
                }
                if (flame_length_mf) {
                    Real v = fl_arr(i, j, k, 0);
                    if (v > max_fl) max_fl = v;
                }
            } else if (pv(i, j, k, 0) < Real(0.1) && pv(i, j, k, 0) >= Real(0.0)) {
                ++n_active;
            }
            // Spotting: comp 2 = landing distance [m]
            if (spotting_data_mf) {
                Real v = sp_arr(i, j, k, 2);
                if (v > max_sd) max_sd = v;
            }
            // Edge-crossing perimeter
            if (i < hi[0]) {
                bool burned_n = (at(i+1, j, k, 0) >= Real(0.0));
                if (burned != burned_n) perim_m += dxx[1];
            }
            if (j < hi[1]) {
                bool burned_n = (at(i, j+1, k, 0) >= Real(0.0));
                if (burned != burned_n) perim_m += dxx[0];
            }
        });
    }

    ParallelDescriptor::ReduceLongSum(n_burned);
    ParallelDescriptor::ReduceLongSum(n_active);
    ParallelDescriptor::ReduceRealSum(perim_m);
    if (fireline_intensity_mf) ParallelDescriptor::ReduceRealMax(max_fi);
    if (flame_length_mf)       ParallelDescriptor::ReduceRealMax(max_fl);
    if (spotting_data_mf)      ParallelDescriptor::ReduceRealMax(max_sd);

    if (!ParallelDescriptor::IOProcessor()) return;

    Real burned_ha = static_cast<Real>(n_burned) * cell_area / Real(1.0e4);
    Real perim_km  = perim_m / Real(1000.0);
    Real time_h    = time_s / Real(3600.0);

    std::ofstream out(filename, std::ios::out | std::ios::app);
    if (!out.is_open()) {
        amrex::Print() << "WARNING: Cannot append to .pst file: " << filename << "\n";
        return;
    }
    out << step        << ","
        << time_h      << ","
        << burned_ha   << ","
        << perim_km    << ","
        << n_active    << ","
        << max_fi      << ","
        << max_fl      << ","
        << max_sd      << "\n";
    out.close();
}

#endif // WRITE_XY_DATA_H