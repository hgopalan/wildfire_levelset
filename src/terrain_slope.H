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

using namespace amrex;

// ======================= Read terrain data and compute slopes ====================
// Reads terrain file in X Y Z format and computes slopes (tan(angle)) for each cell
// Slopes are stored in a MultiFab with 2 components: slope_x and slope_y
//
// The terrain data is assumed to lie on a regular rectangular grid.  The code
// detects the grid spacing automatically and uses bilinear interpolation (O(1)
// per cell) instead of brute-force IDW (O(n_points) per cell), which would hang
// for large datasets such as the 251×251 Gaussian-hill terrain file.
// ==================================================================================
static void compute_slopes_from_terrain(MultiFab& slopes, 
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
        // Skip empty lines and comment/header lines
        if (line.empty() || line[0] == '#') continue;
        
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

    // ------------------------------------------------------------------
    // Detect regular grid structure
    // ------------------------------------------------------------------
    // Collect sorted unique x and y coordinates.
    std::vector<Real> xs = x_data;
    std::vector<Real> ys = y_data;
    std::sort(xs.begin(), xs.end());
    xs.erase(std::unique(xs.begin(), xs.end()), xs.end());
    std::sort(ys.begin(), ys.end());
    ys.erase(std::unique(ys.begin(), ys.end()), ys.end());

    int nx_terrain = static_cast<int>(xs.size());
    int ny_terrain = static_cast<int>(ys.size());

    if (nx_terrain * ny_terrain != static_cast<int>(x_data.size())) {
        amrex::Abort("Terrain file does not appear to be a regular rectangular grid.");
    }

    Real x_min_terrain = xs.front();
    Real y_min_terrain = ys.front();
    Real dx_terrain    = (nx_terrain > 1) ? (xs.back() - xs.front()) / (nx_terrain - 1) : Real(1.0);
    Real dy_terrain    = (ny_terrain > 1) ? (ys.back() - ys.front()) / (ny_terrain - 1) : Real(1.0);

    // Build a 2-D elevation array indexed as z_grid[ix * ny_terrain + iy].
    std::vector<Real> z_grid(nx_terrain * ny_terrain, Real(0.0));
    for (int n = 0; n < static_cast<int>(x_data.size()); ++n) {
        int ix = static_cast<int>(std::round((x_data[n] - x_min_terrain) / dx_terrain));
        int iy = static_cast<int>(std::round((y_data[n] - y_min_terrain) / dy_terrain));
        z_grid[ix * ny_terrain + iy] = z_data[n];
    }

    // Copy the elevation grid to device memory.
    Gpu::DeviceVector<Real> d_z_grid(z_grid.size());
    Gpu::copy(Gpu::hostToDevice, z_grid.begin(), z_grid.end(), d_z_grid.begin());
    Real const* z_ptr = d_z_grid.data();

    // Get domain bounds and cell size
    const auto dx = geom.CellSizeArray();
    const auto problo = geom.ProbLoArray();

    // ------------------------------------------------------------------
    // Bilinear interpolation helper (O(1) per query point)
    // ------------------------------------------------------------------
    // Returns the elevation at terrain coordinates (px, py) using bilinear
    // interpolation on the regular grid.  Queries outside the grid are
    // clamped to the nearest grid edge.
    auto bilinear_z = [=] AMREX_GPU_DEVICE (Real px, Real py) -> Real
    {
        Real fx = (px - x_min_terrain) / dx_terrain;
        Real fy = (py - y_min_terrain) / dy_terrain;

        // Clamp to valid index range
        fx = amrex::max(Real(0.0), amrex::min(fx, Real(nx_terrain - 1)));
        fy = amrex::max(Real(0.0), amrex::min(fy, Real(ny_terrain - 1)));

        int ix0 = static_cast<int>(fx);
        int iy0 = static_cast<int>(fy);
        int ix1 = amrex::min(ix0 + 1, nx_terrain - 1);
        int iy1 = amrex::min(iy0 + 1, ny_terrain - 1);

        Real tx = fx - ix0;
        Real ty = fy - iy0;

        Real z00 = z_ptr[ix0 * ny_terrain + iy0];
        Real z10 = z_ptr[ix1 * ny_terrain + iy0];
        Real z01 = z_ptr[ix0 * ny_terrain + iy1];
        Real z11 = z_ptr[ix1 * ny_terrain + iy1];

        return (Real(1.0) - tx) * (Real(1.0) - ty) * z00
             +              tx  * (Real(1.0) - ty) * z10
             + (Real(1.0) - tx) *              ty  * z01
             +              tx  *              ty  * z11;
    };

    // Compute slopes for each cell using finite differences
    // slopes has 2 components: slope_x (∂z/∂x) and slope_y (∂z/∂y)
    for (MFIter mfi(slopes); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const s = slopes.array(mfi);
        
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            // Compute cell center coordinates
            Real x = problo[0] + (i + Real(0.5)) * dx[0];
            Real y = problo[1] + (j + Real(0.5)) * dx[1];
            
            // Interpolate elevation at neighboring points using bilinear interp
            Real z_xplus  = bilinear_z(x + dx[0], y);
            Real z_xminus = bilinear_z(x - dx[0], y);
            Real z_yplus  = bilinear_z(x, y + dx[1]);
            Real z_yminus = bilinear_z(x, y - dx[1]);
            
            // Compute slopes using central differences
            // slope_x = ∂z/∂x = tan(angle_x)
            // slope_y = ∂z/∂y = tan(angle_y)
            s(i,j,k,0) = (z_xplus - z_xminus) / (Real(2.0) * dx[0]);
            s(i,j,k,1) = (z_yplus - z_yminus) / (Real(2.0) * dx[1]);
        });
    }
    
    amrex::Print() << "Computed terrain slopes from " << nx_terrain << "x" << ny_terrain
                   << " regular grid (spacing " << dx_terrain << " x " << dy_terrain << ")\n";
}

#endif // TERRAIN_SLOPE_H
