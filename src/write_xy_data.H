// write_xy_data.H - Functions for writing x-y data files for phi analysis
#ifndef WRITE_XY_DATA_H
#define WRITE_XY_DATA_H

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <fstream>
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
           "active_front_cells,head_ros_ms,"
           "total_co2_kg_m2,total_co_kg_m2,total_pm25_kg_m2\n";
    out.close();
}

// ============================================================================
// append_fire_stats
//
// Computes burned area [ha], perimeter length [km], active front cell count,
// head ROS [m/s] (maximum ROS in perimeter cells), and domain-mean emission
// totals, then appends one row to the CSV file.
//
// emissions_mf is optional (pass nullptr to skip emission stats).
// R_mf is optional (pass nullptr to skip head ROS tracking).
// ============================================================================
inline void append_fire_stats(const MultiFab& phi,
                               const Geometry& geom,
                               const MultiFab* emissions_mf,
                               int step, Real time_s,
                               const std::string& filename,
                               const MultiFab* R_mf = nullptr)
{
    const auto dxx = geom.CellSize();
    Real cell_area = dxx[0] * dxx[1];

    amrex::Long n_burned = 0, n_perim = 0, n_active = 0;
    Real sum_co2 = 0, sum_co = 0, sum_pm25 = 0;
    Real head_ros = Real(0.0);  // maximum ROS in perimeter cells [m/s]

    for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const p  = phi.const_array(mfi);
        const IntVect lo = bx.smallEnd();
        const IntVect hi = bx.bigEnd();

        // Head ROS: optional per-cell R_mf array for this MFIter
        Array4<const Real> R_arr;
        if (R_mf) R_arr = R_mf->const_array(mfi);

        amrex::LoopOnCpu(bx, [&](int i, int j, int k) {
            Real pv = p(i, j, k, 0);
            if (pv < Real(0.0)) {
                ++n_burned;
                bool on_perim = false;
                if (i > lo[0]) on_perim |= (p(i-1,j,k,0) >= Real(0.0));
                if (i < hi[0]) on_perim |= (p(i+1,j,k,0) >= Real(0.0));
                if (j > lo[1]) on_perim |= (p(i,j-1,k,0) >= Real(0.0));
                if (j < hi[1]) on_perim |= (p(i,j+1,k,0) >= Real(0.0));
                if (on_perim) {
                    ++n_perim;
                    if (R_mf) {
                        Real r = R_arr(i, j, k, 0);
                        if (r > head_ros) head_ros = r;
                    }
                }
            } else if (pv < Real(0.1)) {
                ++n_active;  // cells within phi_thresh = 0.1 of the front
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
    ParallelDescriptor::ReduceLongSum(n_perim);
    ParallelDescriptor::ReduceLongSum(n_active);
    if (R_mf) {
        ParallelDescriptor::ReduceRealMax(head_ros);
    }
    if (emissions_mf) {
        ParallelDescriptor::ReduceRealSum(sum_co2);
        ParallelDescriptor::ReduceRealSum(sum_co);
        ParallelDescriptor::ReduceRealSum(sum_pm25);
    }

    Real burned_ha    = static_cast<Real>(n_burned) * cell_area / Real(1.0e4);
    Real perim_km     = static_cast<Real>(n_perim)  * std::sqrt(cell_area) / Real(1000.0);
    // Convert domain-sum emission [kg/m²·cell] to mean [kg/m²]
    amrex::Long total_cells  = geom.Domain().numPts();
    Real mean_co2  = (total_cells > 0) ? sum_co2  / static_cast<Real>(total_cells) : Real(0.0);
    Real mean_co   = (total_cells > 0) ? sum_co   / static_cast<Real>(total_cells) : Real(0.0);
    Real mean_pm25 = (total_cells > 0) ? sum_pm25 / static_cast<Real>(total_cells) : Real(0.0);

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
        << mean_co2   << ","
        << mean_co    << ","
        << mean_pm25  << "\n";
    out.close();
}

#endif // WRITE_XY_DATA_H