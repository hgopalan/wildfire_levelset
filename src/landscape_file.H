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

#include "spatial_grid.H"

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
    Gpu::DeviceVector<Real> d_sx, d_sy, d_sslope, d_saspect;
    Gpu::DeviceVector<int>  d_bstart, d_bcount;
    CellListParams cl;

    // Characteristic data spacing
    int n_points = static_cast<int>(x_data.size());
    Real x_min_d = *std::min_element(x_data.begin(), x_data.end());
    Real x_max_d = *std::max_element(x_data.begin(), x_data.end());
    Real y_min_d = *std::min_element(y_data.begin(), y_data.end());
    Real y_max_d = *std::max_element(y_data.begin(), y_data.end());
    Real area    = (x_max_d - x_min_d) * (y_max_d - y_min_d);
    Real h = (area > Real(0.0) && n_points > 1)
                 ? std::sqrt(area / static_cast<Real>(n_points))
                 : Real(1.0);

    std::vector<int> sidx = build_cell_list(x_data, y_data, h,
                                            d_sx, d_sy,
                                            d_bstart, d_bcount, cl);

    std::vector<Real> slope_sorted  = reorder_values(slope_data,  sidx);
    std::vector<Real> aspect_sorted = reorder_values(aspect_data, sidx);
    d_sslope.resize(n_points);
    d_saspect.resize(n_points);
    Gpu::copy(Gpu::hostToDevice, slope_sorted.begin(),  slope_sorted.end(),  d_sslope.begin());
    Gpu::copy(Gpu::hostToDevice, aspect_sorted.begin(), aspect_sorted.end(), d_saspect.begin());

    cl.sx     = d_sx.data();
    cl.sy     = d_sy.data();
    cl.bstart = d_bstart.data();
    cl.bcount = d_bcount.data();

    Real const* slope_ptr  = d_sslope.data();
    Real const* aspect_ptr = d_saspect.data();
    
    // Get domain bounds and cell size
    const auto dx     = geom.CellSizeArray();
    const auto problo = geom.ProbLoArray();
    
    // Compute slopes for each cell using Wendland C² RBF with cell-list index
    for (MFIter mfi(slopes); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const s  = slopes.array(mfi);
        
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real x = problo[0] + (i + Real(0.5)) * dx[0];
            Real y = problo[1] + (j + Real(0.5)) * dx[1];

            Real result_slope  = cell_list_interp(cl, slope_ptr,  x, y);
            Real result_aspect = cell_list_interp(cl, aspect_ptr, x, y);

            // Convert slope (degrees) and aspect (degrees) to slope_x and slope_y
            // slope_x = tan(slope) * cos(aspect_math),  slope_y = tan(slope) * sin(aspect_math)
            //
            // FARSITE aspect convention: 0°=North, increases clockwise (compass bearing).
            // Mathematical (trig) convention: 0°=East, increases counter-clockwise.
            // Relationship: aspect_math = 90° - aspect_farsite.
            // This is exact because the cos/sin decomposition of a compass bearing α onto
            // Cartesian (East, North) axes gives (sin α, cos α), which equals
            // (cos(90-α), sin(90-α)) in standard trig convention.
            Real slope_radians  = result_slope * Real(M_PI) / Real(180.0);
            Real aspect_farsite = result_aspect;
            Real aspect_math    = Real(90.0) - aspect_farsite;  // compass → trig
            Real aspect_radians = aspect_math * Real(M_PI) / Real(180.0);

            Real tan_slope = std::tan(slope_radians);
            s(i, j, k, 0) = tan_slope * std::cos(aspect_radians); // slope_x
            s(i, j, k, 1) = tan_slope * std::sin(aspect_radians); // slope_y
        });
    }
    
    amrex::Print() << "Computed slopes from landscape file with " << n_points << " data points\n";
}

// ============================================================================
// Interpolate elevation from landscape file data to a MultiFab.
// Uses inverse distance weighting on the elevation column of the file.
// ============================================================================
static void compute_elevation_from_landscape(MultiFab& elevation,
                                             const Geometry& geom,
                                             const std::string& filename)
{
    std::vector<Real> x_data, y_data, elevation_data, slope_data, aspect_data;
    std::vector<int>  fuel_model_data;

    read_landscape_file(filename, x_data, y_data, elevation_data,
                        slope_data, aspect_data, fuel_model_data);

    // Copy coordinates and elevation to device memory using cell-list index
    Gpu::DeviceVector<Real> d_sx, d_sy, d_selev;
    Gpu::DeviceVector<int>  d_bstart, d_bcount;
    CellListParams cl;

    int n_points = static_cast<int>(x_data.size());
    Real x_min_d = *std::min_element(x_data.begin(), x_data.end());
    Real x_max_d = *std::max_element(x_data.begin(), x_data.end());
    Real y_min_d = *std::min_element(y_data.begin(), y_data.end());
    Real y_max_d = *std::max_element(y_data.begin(), y_data.end());
    Real area    = (x_max_d - x_min_d) * (y_max_d - y_min_d);
    Real h = (area > Real(0.0) && n_points > 1)
                 ? std::sqrt(area / static_cast<Real>(n_points))
                 : Real(1.0);

    std::vector<int> sidx = build_cell_list(x_data, y_data, h,
                                            d_sx, d_sy,
                                            d_bstart, d_bcount, cl);

    std::vector<Real> elev_sorted = reorder_values(elevation_data, sidx);
    d_selev.resize(n_points);
    Gpu::copy(Gpu::hostToDevice, elev_sorted.begin(), elev_sorted.end(), d_selev.begin());

    cl.sx     = d_sx.data();
    cl.sy     = d_sy.data();
    cl.bstart = d_bstart.data();
    cl.bcount = d_bcount.data();

    Real const* elev_ptr = d_selev.data();

    const auto dx     = geom.CellSizeArray();
    const auto problo = geom.ProbLoArray();

    for (MFIter mfi(elevation); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const elev = elevation.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real x = problo[0] + (i + Real(0.5)) * dx[0];
            Real y = problo[1] + (j + Real(0.5)) * dx[1];
            elev(i, j, k) = cell_list_interp(cl, elev_ptr, x, y);
        });
    }

    amrex::Print() << "Interpolated landscape elevation from " << n_points << " data points\n";
}

// ============================================================================
// Interpolate fuel model from landscape file data to a MultiFab.
// Used when a terrain file provides elevation/slope/aspect but fuel
// information is still needed from the landscape file.
// ============================================================================
static void compute_fuel_model_from_landscape(MultiFab& fuel_model_out,
                                              const Geometry& geom,
                                              const std::string& filename)
{
    std::vector<Real> x_data, y_data, elevation_data, slope_data, aspect_data;
    std::vector<int>  fuel_model_data;

    read_landscape_file(filename, x_data, y_data, elevation_data,
                        slope_data, aspect_data, fuel_model_data);

    int n_points = static_cast<int>(x_data.size());

    Real x_min_d = *std::min_element(x_data.begin(), x_data.end());
    Real x_max_d = *std::max_element(x_data.begin(), x_data.end());
    Real y_min_d = *std::min_element(y_data.begin(), y_data.end());
    Real y_max_d = *std::max_element(y_data.begin(), y_data.end());
    Real area    = (x_max_d - x_min_d) * (y_max_d - y_min_d);
    Real h = (area > Real(0.0) && n_points > 1)
                 ? std::sqrt(area / static_cast<Real>(n_points))
                 : Real(1.0);

    Gpu::DeviceVector<Real> d_sx, d_sy;
    Gpu::DeviceVector<int>  d_bstart, d_bcount;
    CellListParams cl;

    std::vector<int> sidx = build_cell_list(x_data, y_data, h,
                                            d_sx, d_sy,
                                            d_bstart, d_bcount, cl);

    std::vector<Real> fuel_real(fuel_model_data.begin(), fuel_model_data.end());
    std::vector<Real> fuel_sorted = reorder_values(fuel_real, sidx);

    Gpu::DeviceVector<Real> d_sfuel(n_points);
    Gpu::copy(Gpu::hostToDevice, fuel_sorted.begin(), fuel_sorted.end(), d_sfuel.begin());

    cl.sx     = d_sx.data();
    cl.sy     = d_sy.data();
    cl.bstart = d_bstart.data();
    cl.bcount = d_bcount.data();

    Real const* fuel_ptr = d_sfuel.data();

    const auto dx     = geom.CellSizeArray();
    const auto problo = geom.ProbLoArray();

    for (MFIter mfi(fuel_model_out); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const f  = fuel_model_out.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real x = problo[0] + (i + Real(0.5)) * dx[0];
            Real y = problo[1] + (j + Real(0.5)) * dx[1];
            // Round to nearest integer: fuel model is a discrete category
            f(i, j, k) = std::round(cell_list_interp(cl, fuel_ptr, x, y));
        });
    }

    amrex::Print() << "Interpolated fuel model from landscape file with " << n_points << " data points\n";
}

// ============================================================================
// Interpolate slope (degrees), aspect (degrees), and fuel model from
// landscape file data to MultiFabs.  All three MultiFabs must have 1
// component and no ghost cells.
// ============================================================================
static void compute_landscape_slope_aspect_fuel(MultiFab& slope_out,
                                                MultiFab& aspect_out,
                                                MultiFab& fuel_model_out,
                                                const Geometry& geom,
                                                const std::string& filename)
{
    std::vector<Real> x_data, y_data, elevation_data, slope_data, aspect_data;
    std::vector<int>  fuel_model_data;

    read_landscape_file(filename, x_data, y_data, elevation_data,
                        slope_data, aspect_data, fuel_model_data);

    int n_points = static_cast<int>(x_data.size());

    Real x_min_d = *std::min_element(x_data.begin(), x_data.end());
    Real x_max_d = *std::max_element(x_data.begin(), x_data.end());
    Real y_min_d = *std::min_element(y_data.begin(), y_data.end());
    Real y_max_d = *std::max_element(y_data.begin(), y_data.end());
    Real area    = (x_max_d - x_min_d) * (y_max_d - y_min_d);
    Real h = (area > Real(0.0) && n_points > 1)
                 ? std::sqrt(area / static_cast<Real>(n_points))
                 : Real(1.0);

    Gpu::DeviceVector<Real> d_sx, d_sy;
    Gpu::DeviceVector<int>  d_bstart, d_bcount;
    CellListParams cl;

    std::vector<int> sidx = build_cell_list(x_data, y_data, h,
                                            d_sx, d_sy,
                                            d_bstart, d_bcount, cl);

    std::vector<Real> slope_sorted  = reorder_values(slope_data,  sidx);
    std::vector<Real> aspect_sorted = reorder_values(aspect_data, sidx);
    // Convert int fuel model to Real for interpolation
    std::vector<Real> fuel_real(fuel_model_data.begin(), fuel_model_data.end());
    std::vector<Real> fuel_sorted   = reorder_values(fuel_real, sidx);

    Gpu::DeviceVector<Real> d_sslope(n_points), d_saspect(n_points), d_sfuel(n_points);
    Gpu::copy(Gpu::hostToDevice, slope_sorted.begin(),  slope_sorted.end(),  d_sslope.begin());
    Gpu::copy(Gpu::hostToDevice, aspect_sorted.begin(), aspect_sorted.end(), d_saspect.begin());
    Gpu::copy(Gpu::hostToDevice, fuel_sorted.begin(),   fuel_sorted.end(),   d_sfuel.begin());

    cl.sx     = d_sx.data();
    cl.sy     = d_sy.data();
    cl.bstart = d_bstart.data();
    cl.bcount = d_bcount.data();

    Real const* slope_ptr  = d_sslope.data();
    Real const* aspect_ptr = d_saspect.data();
    Real const* fuel_ptr   = d_sfuel.data();

    const auto dx     = geom.CellSizeArray();
    const auto problo = geom.ProbLoArray();

    for (MFIter mfi(slope_out); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const s  = slope_out.array(mfi);
        auto const a  = aspect_out.array(mfi);
        auto const f  = fuel_model_out.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real x = problo[0] + (i + Real(0.5)) * dx[0];
            Real y = problo[1] + (j + Real(0.5)) * dx[1];
            s(i, j, k) = cell_list_interp(cl, slope_ptr,  x, y);
            a(i, j, k) = cell_list_interp(cl, aspect_ptr, x, y);
            // Round to nearest integer: fuel model is a discrete category
            f(i, j, k) = std::round(cell_list_interp(cl, fuel_ptr, x, y));
        });
    }

    amrex::Print() << "Interpolated landscape slope/aspect/fuel from " << n_points << " data points\n";
}

#endif // LANDSCAPE_FILE_H
