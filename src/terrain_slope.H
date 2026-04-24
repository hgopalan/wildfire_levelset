// terrain_slope.H - Terrain slope calculation from X Y Z data
#ifndef TERRAIN_SLOPE_H
#define TERRAIN_SLOPE_H

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_Gpu.H>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace amrex;

// ======================= Read terrain data and compute slopes ====================
// Reads terrain file in X Y Z format and computes slopes (tan(angle)) for each cell
// Slopes are stored in a MultiFab with 2 components: slope_x and slope_y
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
        // Skip empty lines and comments
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
    
    // Copy data to device-accessible memory
    Gpu::DeviceVector<Real> d_x(x_data.size());
    Gpu::DeviceVector<Real> d_y(y_data.size());
    Gpu::DeviceVector<Real> d_z(z_data.size());
    
    Gpu::copy(Gpu::hostToDevice, x_data.begin(), x_data.end(), d_x.begin());
    Gpu::copy(Gpu::hostToDevice, y_data.begin(), y_data.end(), d_y.begin());
    Gpu::copy(Gpu::hostToDevice, z_data.begin(), z_data.end(), d_z.begin());
    
    Real const* x_ptr = d_x.data();
    Real const* y_ptr = d_y.data();
    Real const* z_ptr = d_z.data();
    int n_points = static_cast<int>(x_data.size());
    
    // Get domain bounds and cell size
    const auto dx = geom.CellSizeArray();
    const auto problo = geom.ProbLoArray();
    
    // Compute slopes for each cell using finite differences
    // slopes has 2 components: slope_x (∂z/∂x) and slope_y (∂z/∂y)
    for (MFIter mfi(slopes); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const s = slopes.array(mfi);
        
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            // Compute cell center coordinates
            Real x = problo[0] + (i + Real(0.5)) * dx[0];
            Real y = problo[1] + (j + Real(0.5)) * dx[1];
            
            // Interpolate elevation at this cell center and neighboring points
            // using inverse distance weighting
            Real epsilon = Real(1.0e-10);
            
            // Helper function for IDW interpolation
            auto interpolate_z = [&](Real px, Real py) -> Real {
                Real sum_weights = Real(0.0);
                Real sum_z = Real(0.0);
                bool exact_match = false;
                Real result = Real(0.0);
                
                for (int n = 0; n < n_points; ++n) {
                    Real dx_val = px - x_ptr[n];
                    Real dy_val = py - y_ptr[n];
                    Real dist_sq = dx_val * dx_val + dy_val * dy_val;
                    
                    // If very close to a data point, use its value directly
                    if (dist_sq < epsilon * epsilon) {
                        result = z_ptr[n];
                        exact_match = true;
                        break;
                    }
                    
                    // For power=2, weight = 1/dist^2
                    Real weight = Real(1.0) / dist_sq;
                    sum_weights += weight;
                    sum_z += weight * z_ptr[n];
                }
                
                if (!exact_match) {
                    result = sum_z / sum_weights;
                }
                return result;
            };
            
            // Interpolate elevation at cell center and neighbors
            // Note: z_center is computed but not used in the slope calculation
            // Real z_center = interpolate_z(x, y);
            Real z_xplus  = interpolate_z(x + dx[0], y);
            Real z_xminus = interpolate_z(x - dx[0], y);
            Real z_yplus  = interpolate_z(x, y + dx[1]);
            Real z_yminus = interpolate_z(x, y - dx[1]);
            
            // Compute slopes using central differences
            // slope_x = ∂z/∂x = tan(angle_x)
            // slope_y = ∂z/∂y = tan(angle_y)
            s(i,j,k,0) = (z_xplus - z_xminus) / (Real(2.0) * dx[0]);
            s(i,j,k,1) = (z_yplus - z_yminus) / (Real(2.0) * dx[1]);
        });
    }
    
    amrex::Print() << "Computed terrain slopes from unstructured data\n";
}

#endif // TERRAIN_SLOPE_H
