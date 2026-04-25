// write_xy_data.H - Functions for writing x-y data files for phi analysis
#ifndef WRITE_XY_DATA_H
#define WRITE_XY_DATA_H

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <algorithm>
#include <cmath>

using namespace amrex;

// ======================= Write all negative phi cells as x-y file ================
// Writes coordinates (x, y) of all cells where phi < 0 to a file
// ==================================================================================
static void write_negative_phi_xy(const MultiFab& phi, 
                                   const Geometry& geom, 
                                   const std::string& filename)
{
    // Collect all negative phi cell centers
    std::vector<std::pair<Real, Real>> negative_points;
    
    const auto dx = geom.CellSizeArray();
    const auto problo = geom.ProbLoArray();
    
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
// Writes the convex hull of all cells where phi < 0 as x-y file
// ==================================================================================
static void write_negative_phi_convex_hull(const MultiFab& phi, 
                                           const Geometry& geom, 
                                           const std::string& filename)
{
    // Collect all negative phi cell centers
    std::vector<std::pair<Real, Real>> negative_points;
    
    const auto dx = geom.CellSizeArray();
    const auto problo = geom.ProbLoArray();
    
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

#endif // WRITE_XY_DATA_H