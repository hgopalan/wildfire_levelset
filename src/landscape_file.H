// landscape_file.H - FARSITE landscape file (.lcp) reader for ASCII format
// Supports reading elevation, slope, aspect, and fuel model data
#ifndef LANDSCAPE_FILE_H
#define LANDSCAPE_FILE_H

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

using namespace amrex;

// ============================================================================
// Read FARSITE landscape file (ASCII format)
// Format: X Y ELEVATION SLOPE ASPECT FUEL_MODEL
// - X, Y: coordinates in meters
// - ELEVATION: elevation in meters
// - SLOPE: slope in degrees (0-90)
// - ASPECT: aspect in degrees (0-360, 0=North, 90=East, 180=South, 270=West)
// - FUEL_MODEL: fuel model number (optional)
// ============================================================================
static void read_landscape_file(const std::string& filename,
                                std::vector<Real>& x_data,
                                std::vector<Real>& y_data,
                                std::vector<Real>& elevation_data,
                                std::vector<Real>& slope_data,
                                std::vector<Real>& aspect_data,
                                std::vector<int>& fuel_model_data)
{
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        amrex::Abort("Cannot open landscape file: " + filename);
    }
    
    std::string line;
    int line_num = 0;
    
    while (std::getline(infile, line)) {
        line_num++;
        
        // Skip empty lines and comments
        if (line.empty() || line[0] == '#') continue;
        
        std::istringstream iss(line);
        Real x, y, elevation, slope, aspect;
        int fuel_model = 0; // default fuel model
        
        // Read required fields
        if (!(iss >> x >> y >> elevation >> slope >> aspect)) {
            amrex::Print() << "WARNING: Skipping malformed line " << line_num << " in landscape file\n";
            continue;
        }
        
        // Read optional fuel model field
        iss >> fuel_model;
        
        x_data.push_back(x);
        y_data.push_back(y);
        elevation_data.push_back(elevation);
        slope_data.push_back(slope);
        aspect_data.push_back(aspect);
        fuel_model_data.push_back(fuel_model);
    }
    infile.close();
    
    if (x_data.empty()) {
        amrex::Abort("No data read from landscape file: " + filename);
    }
    
    amrex::Print() << "Read " << x_data.size() << " landscape data points from " << filename << "\n";
}

// ============================================================================
// Compute slopes from landscape file data and store in MultiFab
// The landscape file contains slope and aspect, which we convert to
// slope_x and slope_y (tangent of slope angle in x and y directions)
//
// The landscape data is assumed to lie on a regular rectangular grid.  The
// code detects the grid spacing automatically and uses bilinear interpolation
// (O(1) per cell) instead of brute-force IDW (O(n_points) per cell), which
// would hang for large datasets.
// ============================================================================
static void compute_slopes_from_landscape(MultiFab& slopes, 
                                         const Geometry& geom, 
                                         const std::string& filename)
{
    // Read landscape data from file
    std::vector<Real> x_data, y_data, elevation_data, slope_data, aspect_data;
    std::vector<int> fuel_model_data;
    
    read_landscape_file(filename, x_data, y_data, elevation_data, 
                       slope_data, aspect_data, fuel_model_data);

    // ------------------------------------------------------------------
    // Detect regular grid structure
    // ------------------------------------------------------------------
    std::vector<Real> xs = x_data;
    std::vector<Real> ys = y_data;
    std::sort(xs.begin(), xs.end());
    xs.erase(std::unique(xs.begin(), xs.end()), xs.end());
    std::sort(ys.begin(), ys.end());
    ys.erase(std::unique(ys.begin(), ys.end()), ys.end());

    int nx_land = static_cast<int>(xs.size());
    int ny_land = static_cast<int>(ys.size());

    if (nx_land * ny_land != static_cast<int>(x_data.size())) {
        amrex::Abort("Landscape file does not appear to be a regular rectangular grid.");
    }

    if (nx_land < 2 || ny_land < 2) {
        amrex::Abort("Landscape grid must have at least 2 points in each direction.");
    }

    Real x_min_land = xs.front();
    Real y_min_land = ys.front();
    Real dx_land    = (xs.back() - xs.front()) / (nx_land - 1);
    Real dy_land    = (ys.back() - ys.front()) / (ny_land - 1);

    // Build 2-D arrays for slope and aspect indexed as [ix * ny_land + iy].
    std::vector<Real> slope_grid(nx_land * ny_land, Real(0.0));
    std::vector<Real> aspect_grid(nx_land * ny_land, Real(0.0));
    for (int n = 0; n < static_cast<int>(x_data.size()); ++n) {
        int ix = static_cast<int>(std::round((x_data[n] - x_min_land) / dx_land));
        int iy = static_cast<int>(std::round((y_data[n] - y_min_land) / dy_land));
        // Clamp to guard against floating-point precision at the boundary
        ix = std::max(0, std::min(ix, nx_land - 1));
        iy = std::max(0, std::min(iy, ny_land - 1));
        slope_grid [ix * ny_land + iy] = slope_data [n];
        aspect_grid[ix * ny_land + iy] = aspect_data[n];
    }

    // Copy grids to device memory.
    Gpu::DeviceVector<Real> d_slope_grid(slope_grid.size());
    Gpu::DeviceVector<Real> d_aspect_grid(aspect_grid.size());
    Gpu::copy(Gpu::hostToDevice, slope_grid.begin(),  slope_grid.end(),  d_slope_grid.begin());
    Gpu::copy(Gpu::hostToDevice, aspect_grid.begin(), aspect_grid.end(), d_aspect_grid.begin());
    Real const* slope_ptr  = d_slope_grid.data();
    Real const* aspect_ptr = d_aspect_grid.data();

    // Get domain bounds and cell size
    const auto dx = geom.CellSizeArray();
    const auto problo = geom.ProbLoArray();

    // ------------------------------------------------------------------
    // Bilinear interpolation helper (O(1) per query point)
    // ------------------------------------------------------------------
    // Returns the value at landscape coordinates (px, py) using bilinear
    // interpolation on the regular grid.  Queries outside the grid are
    // clamped to the nearest grid edge.
    auto interpolate_grid = [=] AMREX_GPU_DEVICE (Real const* grid_ptr, Real px, Real py) -> Real
    {
        Real fx = (px - x_min_land) / dx_land;
        Real fy = (py - y_min_land) / dy_land;

        fx = amrex::max(Real(0.0), amrex::min(fx, Real(nx_land - 1)));
        fy = amrex::max(Real(0.0), amrex::min(fy, Real(ny_land - 1)));

        int ix0 = static_cast<int>(fx);
        int iy0 = static_cast<int>(fy);
        int ix1 = amrex::min(ix0 + 1, nx_land - 1);
        int iy1 = amrex::min(iy0 + 1, ny_land - 1);

        Real tx = fx - ix0;
        Real ty = fy - iy0;

        Real v00 = grid_ptr[ix0 * ny_land + iy0];
        Real v10 = grid_ptr[ix1 * ny_land + iy0];
        Real v01 = grid_ptr[ix0 * ny_land + iy1];
        Real v11 = grid_ptr[ix1 * ny_land + iy1];

        return (Real(1.0) - tx) * (Real(1.0) - ty) * v00
             +              tx  * (Real(1.0) - ty) * v10
             + (Real(1.0) - tx) *              ty  * v01
             +              tx  *              ty  * v11;
    };

    // Compute slopes for each cell
    // slopes has 2 components: slope_x (∂z/∂x) and slope_y (∂z/∂y)
    for (MFIter mfi(slopes); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const s = slopes.array(mfi);
        
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            // Compute cell center coordinates
            Real x = problo[0] + (i + Real(0.5)) * dx[0];
            Real y = problo[1] + (j + Real(0.5)) * dx[1];

            // Bilinear interpolation of slope and aspect at this cell center
            Real result_slope  = interpolate_grid(slope_ptr,  x, y);
            Real result_aspect = interpolate_grid(aspect_ptr, x, y);

            // Convert slope (degrees) and aspect (degrees) to slope_x and slope_y
            // slope_x = tan(slope) * cos(aspect_radians_from_east)
            // slope_y = tan(slope) * sin(aspect_radians_from_east)
            // 
            // FARSITE aspect convention: 0=North, 90=East, 180=South, 270=West
            // Mathematical convention: 0=East, 90=North, 180=West, 270=South
            // Convert: aspect_math = 90 - aspect_farsite
            
            Real slope_radians  = result_slope  * Real(M_PI) / Real(180.0);
            Real aspect_math    = Real(90.0) - result_aspect;
            Real aspect_radians = aspect_math   * Real(M_PI) / Real(180.0);
            
            Real tan_slope = std::tan(slope_radians);
            s(i, j, k, 0) = tan_slope * std::cos(aspect_radians); // slope_x
            s(i, j, k, 1) = tan_slope * std::sin(aspect_radians); // slope_y
        });
    }
    
    amrex::Print() << "Computed slopes from landscape file on " << nx_land << "x" << ny_land
                   << " regular grid (spacing " << dx_land << " x " << dy_land << ")\n";
}

#endif // LANDSCAPE_FILE_H
