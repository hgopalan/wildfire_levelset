// landscape_file.H - FARSITE landscape file (.lcp) reader for ASCII format
// Supports reading elevation, slope, aspect, and fuel model data
#ifndef LANDSCAPE_FILE_H
#define LANDSCAPE_FILE_H

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_Gpu.H>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>

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
    
    // Copy data to device-accessible memory
    Gpu::DeviceVector<Real> d_x(x_data.size());
    Gpu::DeviceVector<Real> d_y(y_data.size());
    Gpu::DeviceVector<Real> d_slope(slope_data.size());
    Gpu::DeviceVector<Real> d_aspect(aspect_data.size());
    
    Gpu::copy(Gpu::hostToDevice, x_data.begin(), x_data.end(), d_x.begin());
    Gpu::copy(Gpu::hostToDevice, y_data.begin(), y_data.end(), d_y.begin());
    Gpu::copy(Gpu::hostToDevice, slope_data.begin(), slope_data.end(), d_slope.begin());
    Gpu::copy(Gpu::hostToDevice, aspect_data.begin(), aspect_data.end(), d_aspect.begin());
    
    Real const* x_ptr = d_x.data();
    Real const* y_ptr = d_y.data();
    Real const* slope_ptr = d_slope.data();
    Real const* aspect_ptr = d_aspect.data();
    int n_points = static_cast<int>(x_data.size());
    
    // Get domain bounds and cell size
    const auto dx = geom.CellSizeArray();
    const auto problo = geom.ProbLoArray();
    
    // Compute slopes for each cell using inverse distance weighting interpolation
    // slopes has 2 components: slope_x (∂z/∂x) and slope_y (∂z/∂y)
    for (MFIter mfi(slopes); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const s = slopes.array(mfi);
        
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            // Compute cell center coordinates
            Real x = problo[0] + (i + Real(0.5)) * dx[0];
            Real y = problo[1] + (j + Real(0.5)) * dx[1];
            
            // Interpolate slope and aspect at this cell center using inverse distance weighting
            Real epsilon = Real(1.0e-10);
            Real sum_weights = Real(0.0);
            Real sum_slope = Real(0.0);
            Real sum_aspect = Real(0.0);
            bool exact_match = false;
            Real result_slope = Real(0.0);
            Real result_aspect = Real(0.0);
            
            for (int n = 0; n < n_points; ++n) {
                Real dx_val = x - x_ptr[n];
                Real dy_val = y - y_ptr[n];
                Real dist_sq = dx_val * dx_val + dy_val * dy_val;
                
                // If very close to a data point, use its value directly
                if (dist_sq < epsilon * epsilon) {
                    result_slope = slope_ptr[n];
                    result_aspect = aspect_ptr[n];
                    exact_match = true;
                    break;
                }
                
                // Inverse distance weighting using sqrt(dist_sq) for better numerical stability
                // This prevents extreme weights when points are very close
                Real dist = std::sqrt(dist_sq);
                Real weight = Real(1.0) / (dist + epsilon);
                sum_weights += weight;
                sum_slope += weight * slope_ptr[n];
                sum_aspect += weight * aspect_ptr[n];
            }
            
            if (!exact_match && sum_weights > epsilon) {
                result_slope = sum_slope / sum_weights;
                result_aspect = sum_aspect / sum_weights;
            }
            
            // Convert slope (degrees) and aspect (degrees) to slope_x and slope_y
            // slope_x = tan(slope) * cos(aspect_radians_from_east)
            // slope_y = tan(slope) * sin(aspect_radians_from_east)
            // 
            // FARSITE aspect convention: 0=North, 90=East, 180=South, 270=West
            // Mathematical convention: 0=East, 90=North, 180=West, 270=South
            // Convert: aspect_math = 90 - aspect_farsite
            
            Real slope_radians = result_slope * Real(M_PI) / Real(180.0);
            Real aspect_farsite = result_aspect;
            Real aspect_math = Real(90.0) - aspect_farsite; // Convert to math convention
            Real aspect_radians = aspect_math * Real(M_PI) / Real(180.0);
            
            Real tan_slope = std::tan(slope_radians);
            s(i, j, k, 0) = tan_slope * std::cos(aspect_radians); // slope_x
            s(i, j, k, 1) = tan_slope * std::sin(aspect_radians); // slope_y
        });
    }
    
    amrex::Print() << "Computed slopes from landscape file with " << n_points << " data points\n";
}

#endif // LANDSCAPE_FILE_H
