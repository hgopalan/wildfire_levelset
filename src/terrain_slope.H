// terrain_slope.H - Terrain slope calculation from X Y Z data
#ifndef TERRAIN_SLOPE_H
#define TERRAIN_SLOPE_H

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_Gpu.H>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "spatial_grid.H"

using namespace amrex;

// ======================= Read terrain data and compute slopes ====================
// Reads terrain file in X Y Z format and computes slopes (tan(angle)) for each cell
// Slopes are stored in a MultiFab with 2 components: slope_x and slope_y
// ==================================================================================
inline void compute_slopes_from_terrain(MultiFab& slopes, 
                                        const Geometry& geom, 
                                        const std::string& filename)
{
    // Read terrain data from file (X Y Z format)
    std::vector<Real> x_data, y_data, z_data;
    
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        amrex::Abort("Cannot open terrain file: " + filename);
    }
    
    std::string line;
    while (std::getline(infile, line)) {
        // Skip empty lines and comments
        if (line.empty() || line[0] == '#') continue;
        
        // Replace commas with spaces to handle CSV format
        std::replace(line.begin(), line.end(), ',', ' ');
        
        std::istringstream iss(line);
        Real x, y, z;
        if (iss >> x >> y >> z) {
            x_data.push_back(x);
            y_data.push_back(y);
            z_data.push_back(z);
        }
    }
    infile.close();
    
    if (x_data.empty()) {
        amrex::Abort("No data read from terrain file: " + filename);
    }

    amrex::Print() << "Read " << x_data.size() << " terrain data points from " << filename << "\n";

    int n_points = static_cast<int>(x_data.size());

    // Characteristic data spacing: h ≈ sqrt(area / N)
    Real x_min_d = *std::min_element(x_data.begin(), x_data.end());
    Real x_max_d = *std::max_element(x_data.begin(), x_data.end());
    Real y_min_d = *std::min_element(y_data.begin(), y_data.end());
    Real y_max_d = *std::max_element(y_data.begin(), y_data.end());
    Real area    = (x_max_d - x_min_d) * (y_max_d - y_min_d);
    Real h = (area > Real(0.0) && n_points > 1)
                 ? std::sqrt(area / static_cast<Real>(n_points))
                 : Real(1.0);

    // Build cell-list spatial index (cutoff R = 3h, bin size = R).
    // Points within R of a query are found in O(K) rather than O(N).
    Gpu::DeviceVector<Real> d_sx, d_sy, d_sz;
    Gpu::DeviceVector<int>  d_bstart, d_bcount;
    CellListParams cl;

    std::vector<int> sidx = build_cell_list(x_data, y_data, h,
                                            d_sx, d_sy,
                                            d_bstart, d_bcount, cl);

    // Reorder z to match sorted layout and copy to device
    std::vector<Real> z_sorted = reorder_values(z_data, sidx);
    d_sz.resize(n_points);
    Gpu::copy(Gpu::hostToDevice, z_sorted.begin(), z_sorted.end(), d_sz.begin());

    // Set device pointers in params
    cl.sx     = d_sx.data();
    cl.sy     = d_sy.data();
    cl.bstart = d_bstart.data();
    cl.bcount = d_bcount.data();

    Real const* z_ptr = d_sz.data();

    // Get domain bounds and cell size
    const auto dx     = geom.CellSizeArray();
    const auto problo = geom.ProbLoArray();

    // Compute slopes for each cell using central finite differences.
    // Elevation at stencil points is interpolated via Wendland C² RBF with
    // the cell-list index – O(K) per cell instead of O(N).
    for (MFIter mfi(slopes); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const s  = slopes.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real x = problo[0] + (i + Real(0.5)) * dx[0];
            Real y = problo[1] + (j + Real(0.5)) * dx[1];

            Real z_xplus  = cell_list_interp(cl, z_ptr, x + dx[0], y);
            Real z_xminus = cell_list_interp(cl, z_ptr, x - dx[0], y);
            Real z_yplus  = cell_list_interp(cl, z_ptr, x,          y + dx[1]);
            Real z_yminus = cell_list_interp(cl, z_ptr, x,          y - dx[1]);

            // slope_x = ∂z/∂x,  slope_y = ∂z/∂y  (central differences)
            s(i,j,k,0) = (z_xplus - z_xminus) / (Real(2.0) * dx[0]);
            s(i,j,k,1) = (z_yplus - z_yminus) / (Real(2.0) * dx[1]);
        });
    }

    amrex::Print() << "Computed terrain slopes from unstructured data\n";
}

// ======================= Read terrain data and interpolate elevation ====================
// Reads terrain file in X Y Z format and interpolates elevation (Z) to each cell
// using a Gaussian RBF kernel.  The result is stored in a single-component
// MultiFab (component 0 = elevation in the same units as the terrain file).
// ==================================================================================
inline void compute_elevation_from_terrain(MultiFab& elevation,
                                           const Geometry& geom,
                                           const std::string& filename)
{
    // Read terrain data from file (X Y Z format)
    std::vector<Real> x_data, y_data, z_data;

    std::ifstream infile(filename);
    if (!infile.is_open()) {
        amrex::Abort("Cannot open terrain file: " + filename);
    }

    std::string line;
    while (std::getline(infile, line)) {
        if (line.empty() || line[0] == '#') continue;
        
        // Replace commas with spaces to handle CSV format
        std::replace(line.begin(), line.end(), ',', ' ');
        
        std::istringstream iss(line);
        Real x, y, z;
        if (iss >> x >> y >> z) {
            x_data.push_back(x);
            y_data.push_back(y);
            z_data.push_back(z);
        }
    }
    infile.close();

    if (x_data.empty()) {
        amrex::Abort("No data read from terrain file: " + filename);
    }

    int n_points = static_cast<int>(x_data.size());

    // Characteristic data spacing
    Real x_min_d = *std::min_element(x_data.begin(), x_data.end());
    Real x_max_d = *std::max_element(x_data.begin(), x_data.end());
    Real y_min_d = *std::min_element(y_data.begin(), y_data.end());
    Real y_max_d = *std::max_element(y_data.begin(), y_data.end());
    Real area    = (x_max_d - x_min_d) * (y_max_d - y_min_d);
    Real h = (area > Real(0.0) && n_points > 1)
                 ? std::sqrt(area / static_cast<Real>(n_points))
                 : Real(1.0);

    // Build cell-list spatial index
    Gpu::DeviceVector<Real> d_sx, d_sy, d_sz;
    Gpu::DeviceVector<int>  d_bstart, d_bcount;
    CellListParams cl;

    std::vector<int> sidx = build_cell_list(x_data, y_data, h,
                                            d_sx, d_sy,
                                            d_bstart, d_bcount, cl);

    std::vector<Real> z_sorted = reorder_values(z_data, sidx);
    d_sz.resize(n_points);
    Gpu::copy(Gpu::hostToDevice, z_sorted.begin(), z_sorted.end(), d_sz.begin());

    cl.sx     = d_sx.data();
    cl.sy     = d_sy.data();
    cl.bstart = d_bstart.data();
    cl.bcount = d_bcount.data();

    Real const* z_ptr = d_sz.data();

    const auto dx     = geom.CellSizeArray();
    const auto problo = geom.ProbLoArray();

    for (MFIter mfi(elevation); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const elev = elevation.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real x = problo[0] + (i + Real(0.5)) * dx[0];
            Real y = problo[1] + (j + Real(0.5)) * dx[1];
            elev(i, j, k) = cell_list_interp(cl, z_ptr, x, y);
        });
    }

    amrex::Print() << "Interpolated terrain elevation from " << n_points << " data points\n";
}

// ======================= Derive slope (degrees) and aspect (degrees) ====================
// Converts a 2-component slopes MultiFab (slope_x = tan*cos, slope_y = tan*sin) into
// slope magnitude in degrees and aspect in degrees (FARSITE convention: 0=North, CW).
// ========================================================================================
inline void compute_slope_aspect_from_slopes(const MultiFab& slopes,
                                             MultiFab& slope_out,
                                             MultiFab& aspect_out)
{
    for (MFIter mfi(slopes); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const s  = slopes.const_array(mfi);
        auto const sd = slope_out.array(mfi);
        auto const ad = aspect_out.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real sx = s(i, j, k, 0); // tan(slope) * cos(aspect_math)
            Real sy = s(i, j, k, 1); // tan(slope) * sin(aspect_math)

            Real tan_slope = std::sqrt(sx * sx + sy * sy);
            // slope in degrees
            sd(i, j, k) = std::atan(tan_slope) * Real(180.0) / Real(M_PI);

            // aspect: recover mathematical angle, then convert to FARSITE compass bearing
            // aspect_math = atan2(sy, sx),  aspect_farsite = 90 - aspect_math
            Real aspect_math_deg = std::atan2(sy, sx) * Real(180.0) / Real(M_PI);
            Real aspect_farsite  = Real(90.0) - aspect_math_deg;
            // Normalise to [0, 360)
            aspect_farsite = std::fmod(aspect_farsite, Real(360.0));
            if (aspect_farsite < Real(0.0)) aspect_farsite += Real(360.0);
            ad(i, j, k) = aspect_farsite;
        });
    }
}

#endif // TERRAIN_SLOPE_H
