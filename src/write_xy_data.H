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
        for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k) {
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

// ======================= Write envelope of negative phi values ====================
// Writes the boundary/envelope of the region where phi < 0 as x-y file
// The envelope is defined as cells where phi < 0 and at least one neighbor has phi >= 0
// ==================================================================================
static void write_negative_phi_envelope(const MultiFab& phi, 
                                        const Geometry& geom, 
                                        const std::string& filename)
{
    // Collect envelope points (cells on the boundary of phi < 0 region)
    std::vector<std::pair<Real, Real>> envelope_points;
    
    const auto dx = geom.CellSizeArray();
    const auto problo = geom.ProbLoArray();
    
    for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const phi_arr = phi.array(mfi);
        
        // Loop over valid cells
        for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k) {
            for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
                for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
                    if (phi_arr(i,j,k,0) < 0.0) {
                        // Check if this cell is on the envelope by checking neighbors
                        // A cell is on the envelope if:
                        // 1. It has phi < 0 AND
                        // 2. At least one neighbor has phi >= 0 OR it's at domain boundary
                        bool on_envelope = false;
                        
                        // Check 4 neighbors in 2D (or 6 in 3D, but we only need x-y for output)
                        // Left neighbor
                        if (i > bx.smallEnd(0)) {
                            if (phi_arr(i-1,j,k,0) >= 0.0) on_envelope = true;
                        } else {
                            on_envelope = true; // At domain boundary
                        }
                        // Right neighbor
                        if (i < bx.bigEnd(0)) {
                            if (phi_arr(i+1,j,k,0) >= 0.0) on_envelope = true;
                        } else {
                            on_envelope = true; // At domain boundary
                        }
                        // Bottom neighbor
                        if (j > bx.smallEnd(1)) {
                            if (phi_arr(i,j-1,k,0) >= 0.0) on_envelope = true;
                        } else {
                            on_envelope = true; // At domain boundary
                        }
                        // Top neighbor
                        if (j < bx.bigEnd(1)) {
                            if (phi_arr(i,j+1,k,0) >= 0.0) on_envelope = true;
                        } else {
                            on_envelope = true; // At domain boundary
                        }
                        
                        if (on_envelope) {
                            // Compute cell center coordinates
                            Real x = problo[0] + (i + Real(0.5)) * dx[0];
                            Real y = problo[1] + (j + Real(0.5)) * dx[1];
                            envelope_points.push_back({x, y});
                        }
                    }
                }
            }
        }
    }
    
    // Gather all points from all processors (if running in parallel)
#ifdef AMREX_USE_MPI
    std::vector<int> counts(ParallelDescriptor::NProcs());
    int local_count = envelope_points.size();
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
    for (const auto& pt : envelope_points) {
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
        envelope_points.clear();
        for (size_t i = 0; i < all_data.size(); i += 2) {
            envelope_points.push_back({all_data[i], all_data[i+1]});
        }
    }
#endif
    
    // Write to file (only on IO processor)
    if (ParallelDescriptor::IOProcessor()) {
        std::ofstream outfile(filename);
        if (!outfile.is_open()) {
            amrex::Abort("Cannot open output file: " + filename);
        }
        
        outfile << "# X Y coordinates of envelope (boundary) of phi < 0 region\n";
        outfile << "# Total points: " << envelope_points.size() << "\n";
        
        for (const auto& pt : envelope_points) {
            outfile << pt.first << " " << pt.second << "\n";
        }
        
        outfile.close();
        amrex::Print() << "Wrote " << envelope_points.size() << " envelope points to " 
                       << filename << "\n";
    }
}

#endif // WRITE_XY_DATA_H
