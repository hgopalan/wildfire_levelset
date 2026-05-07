// landscape_file.H - FARSITE landscape file (.lcp) reader
// Supports both binary LCP and ASCII (X Y ELEV SLOPE ASPECT FUEL) formats.
// Binary format also provides canopy cover, stand height, crown base height,
// and crown bulk density when the file was created with crown fuel layers.
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
#include "lcp_binary_reader.H"

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
inline void read_landscape_file(const std::string& filename,
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
inline void compute_slopes_from_landscape(MultiFab& slopes, 
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
inline void compute_elevation_from_landscape(MultiFab& elevation,
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
inline void compute_fuel_model_from_landscape(MultiFab& fuel_model_out,
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
            // Fuel model is a discrete category: use nearest-neighbour lookup
            // to avoid non-physical fractional blending across fuel boundaries.
            f(i, j, k) = std::round(cell_list_nearest(cl, fuel_ptr, x, y));
        });
    }

    amrex::Print() << "Interpolated fuel model from landscape file with " << n_points << " data points\n";
}

// ============================================================================
// Interpolate slope (degrees), aspect (degrees), and fuel model from
// landscape file data to MultiFabs.  All three MultiFabs must have 1
// component and no ghost cells.
// ============================================================================
inline void compute_landscape_slope_aspect_fuel(MultiFab& slope_out,
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
            // Fuel model is a discrete category: use nearest-neighbour lookup
            // to avoid non-physical fractional blending across fuel boundaries.
            f(i, j, k) = std::round(cell_list_nearest(cl, fuel_ptr, x, y));
        });
    }

    amrex::Print() << "Interpolated landscape slope/aspect/fuel from " << n_points << " data points\n";
}

// ============================================================================
// Helper: build spatial index and interpolate data from a flat array of
// (x,y,value) points on to a MultiFab.  Identical kernel to the ASCII path.
// ============================================================================
inline void interp_flat_to_multifab(
    MultiFab& mf,
    const Geometry& geom,
    const std::vector<double>& x_data,
    const std::vector<double>& y_data,
    const std::vector<float>&  values,
    bool nearest_neighbor = false)
{
    int n = static_cast<int>(x_data.size());
    if (n == 0) return;

    double x_min = *std::min_element(x_data.begin(), x_data.end());
    double x_max = *std::max_element(x_data.begin(), x_data.end());
    double y_min = *std::min_element(y_data.begin(), y_data.end());
    double y_max = *std::max_element(y_data.begin(), y_data.end());
    double area  = (x_max - x_min) * (y_max - y_min);
    Real h = (area > 0.0 && n > 1)
                 ? static_cast<Real>(std::sqrt(area / static_cast<double>(n)))
                 : Real(1.0);

    // Convert to Real vectors for cell-list builder
    std::vector<Real> rx(x_data.begin(), x_data.end());
    std::vector<Real> ry(y_data.begin(), y_data.end());
    std::vector<Real> rv(values.begin(), values.end());

    Gpu::DeviceVector<Real> d_sx, d_sy, d_sv;
    Gpu::DeviceVector<int>  d_bstart, d_bcount;
    CellListParams cl;

    std::vector<int> sidx = build_cell_list(rx, ry, h,
                                             d_sx, d_sy,
                                             d_bstart, d_bcount, cl);
    std::vector<Real> v_sorted = reorder_values(rv, sidx);
    d_sv.resize(n);
    Gpu::copy(Gpu::hostToDevice, v_sorted.begin(), v_sorted.end(), d_sv.begin());

    cl.sx     = d_sx.data();
    cl.sy     = d_sy.data();
    cl.bstart = d_bstart.data();
    cl.bcount = d_bcount.data();

    Real const* vptr = d_sv.data();
    const auto dx    = geom.CellSizeArray();
    const auto plo   = geom.ProbLoArray();
    const bool nn    = nearest_neighbor;

    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto out = mf.array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real x = plo[0] + (i + Real(0.5)) * dx[0];
            Real y = plo[1] + (j + Real(0.5)) * dx[1];
            out(i, j, k) = nn ? std::round(cell_list_nearest(cl, vptr, x, y))
                               : cell_list_interp(cl, vptr, x, y);
        });
    }
}

// ============================================================================
// read_landscape_binary
//
// Reads a binary FARSITE LCP file and produces cell-centred coordinates that
// can be passed to the interpolation helpers above.  Crown layer data vectors
// (cbh_out, cbd_out, cc_out, ht_out) are populated only when
// grid_info.has_crown == true; otherwise they remain empty.
//
// Coordinate convention: the LCP stores data row 0 = northernmost row.
// We convert to (x, y) physical coordinates where x = easting, y = northing
// (using the LCP's XLL/YLL corner and CellSize).
// ============================================================================
static bool read_landscape_binary(
    const std::string&    filename,
    LcpGridInfo&          grid_info,
    std::vector<double>&  x_out,
    std::vector<double>&  y_out,
    std::vector<float>&   elev_out,
    std::vector<float>&   slope_out,
    std::vector<float>&   aspect_out,
    std::vector<float>&   fuel_float_out,    // fuel model as float for interp
    std::vector<int>&     fuel_int_out,      // fuel model as int
    std::vector<float>&   cc_out,            // canopy cover [%]
    std::vector<float>&   ht_out,            // stand height [m]
    std::vector<float>&   cbh_out,           // crown base height [m]
    std::vector<float>&   cbd_out)           // crown bulk density [kg/m³]
{
    std::vector<float> cc_raw, ht_raw, cbh_raw, cbd_raw;
    if (!read_lcp_binary(filename, grid_info, elev_out, slope_out, aspect_out,
                         fuel_int_out, cc_raw, ht_raw, cbh_raw, cbd_raw)) {
        return false;
    }

    const int R = grid_info.numrows;
    const int C = grid_info.numcols;
    const double cs   = grid_info.cellsize;
    const double xll  = grid_info.xllcorner;
    const double yll  = grid_info.yllcorner;

    // Convert foot-based grid to metres if needed
    const double cs_m = (grid_info.gridunits == 1) ? cs * 0.30480 : cs;
    const double xll_m= (grid_info.gridunits == 1) ? xll * 0.30480 : xll;
    const double yll_m= (grid_info.gridunits == 1) ? yll * 0.30480 : yll;

    int ncells = R * C;
    x_out.resize(ncells);
    y_out.resize(ncells);
    fuel_float_out.resize(ncells);

    for (int r = 0; r < R; ++r) {
        for (int c = 0; c < C; ++c) {
            int idx = r * C + c;
            // Cell centre: row 0 is northernmost → northing decreases with r
            x_out[idx] = xll_m + (c + 0.5) * cs_m;
            y_out[idx] = yll_m + (R - r - 0.5) * cs_m;
            fuel_float_out[idx] = static_cast<float>(fuel_int_out[idx]);
        }
    }

    cc_out  = cc_raw;
    ht_out  = ht_raw;
    cbh_out = cbh_raw;
    cbd_out = cbd_raw;
    return true;
}

// ============================================================================
// compute_slopes_from_lcp
//
// Reads a landscape file (binary or ASCII) and populates a 2-component
// MultiFab with (slope_x, slope_y) = tan(slope) × (cos(aspect), sin(aspect)).
// Dispatches to the binary reader when the file is recognised as binary LCP;
// falls back to the existing ASCII path otherwise.
//
// This is a drop-in replacement for compute_slopes_from_landscape() that
// also handles binary LCP files.
// ============================================================================
[[maybe_unused]] inline void compute_slopes_from_lcp(MultiFab& slopes,
                                    const Geometry& geom,
                                    const std::string& filename)
{
    if (is_binary_lcp(filename)) {
        LcpGridInfo info;
        std::vector<double> xd, yd;
        std::vector<float>  elev, slope_f, aspect_f, fuel_f, cc, ht, cbh, cbd;
        std::vector<int>    fuel_i;

        if (!read_landscape_binary(filename, info, xd, yd, elev, slope_f,
                                   aspect_f, fuel_f, fuel_i, cc, ht, cbh, cbd)) {
            amrex::Abort("Failed to read binary LCP file: " + filename);
        }

        int n = static_cast<int>(slope_f.size());
        amrex::Print() << "Binary LCP: " << info.numrows << "×" << info.numcols
                       << " cells,  crown layers: "
                       << (info.has_crown ? "yes" : "no") << "\n";

        // Build slope_x, slope_y from slope (degrees) and aspect (degrees, FARSITE)
        // FARSITE aspect: 0°=North, CW.  Math aspect: aspect_math = 90 - aspect_farsite.
        Gpu::DeviceVector<Real> d_sx, d_sy, d_slope, d_aspect;
        Gpu::DeviceVector<int>  d_bstart, d_bcount;
        CellListParams cl;

        std::vector<Real> rx(xd.begin(), xd.end());
        std::vector<Real> ry(yd.begin(), yd.end());

        double xmin = *std::min_element(xd.begin(), xd.end());
        double xmax = *std::max_element(xd.begin(), xd.end());
        double ymin = *std::min_element(yd.begin(), yd.end());
        double ymax = *std::max_element(yd.begin(), yd.end());
        double area = (xmax - xmin) * (ymax - ymin);
        Real h = (area > 0.0 && n > 1)
                     ? static_cast<Real>(std::sqrt(area / static_cast<double>(n)))
                     : Real(1.0);

        std::vector<int> sidx = build_cell_list(rx, ry, h, d_sx, d_sy,
                                                 d_bstart, d_bcount, cl);

        std::vector<Real> slope_r(n), aspect_r(n);
        for (int i = 0; i < n; ++i) {
            slope_r[i]  = static_cast<Real>(slope_f[sidx[i]]);
            aspect_r[i] = static_cast<Real>(aspect_f[sidx[i]]);
        }
        d_slope.resize(n);  d_aspect.resize(n);
        Gpu::copy(Gpu::hostToDevice, slope_r.begin(),  slope_r.end(),  d_slope.begin());
        Gpu::copy(Gpu::hostToDevice, aspect_r.begin(), aspect_r.end(), d_aspect.begin());

        cl.sx     = d_sx.data();
        cl.sy     = d_sy.data();
        cl.bstart = d_bstart.data();
        cl.bcount = d_bcount.data();
        Real const* sp = d_slope.data();
        Real const* ap = d_aspect.data();

        const auto dx  = geom.CellSizeArray();
        const auto plo = geom.ProbLoArray();

        for (MFIter mfi(slopes); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            auto s = slopes.array(mfi);
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                Real x = plo[0] + (i + Real(0.5)) * dx[0];
                Real y = plo[1] + (j + Real(0.5)) * dx[1];
                Real sl  = cell_list_interp(cl, sp, x, y);
                Real asp = cell_list_interp(cl, ap, x, y);

                Real sl_rad  = sl  * Real(M_PI) / Real(180.0);
                Real asp_rad = (Real(90.0) - asp) * Real(M_PI) / Real(180.0);
                Real tan_sl  = std::tan(sl_rad);
                s(i, j, k, 0) = tan_sl * std::cos(asp_rad);
                s(i, j, k, 1) = tan_sl * std::sin(asp_rad);
            });
        }
        amrex::Print() << "Computed slopes from binary LCP (" << n << " cells)\n";

    } else {
        // ASCII fallback
        compute_slopes_from_landscape(slopes, geom, filename);
    }
}

// ============================================================================
// compute_elevation_from_lcp
// ============================================================================
[[maybe_unused]] inline void compute_elevation_from_lcp(MultiFab& elev_mf,
                                       const Geometry& geom,
                                       const std::string& filename)
{
    if (is_binary_lcp(filename)) {
        LcpGridInfo info;
        std::vector<double> xd, yd;
        std::vector<float>  elev_v, slope_f, aspect_f, fuel_f, cc, ht, cbh, cbd;
        std::vector<int>    fuel_i;
        if (!read_landscape_binary(filename, info, xd, yd, elev_v, slope_f,
                                   aspect_f, fuel_f, fuel_i, cc, ht, cbh, cbd))
            amrex::Abort("Failed to read binary LCP: " + filename);
        interp_flat_to_multifab(elev_mf, geom, xd, yd, elev_v, false);
        amrex::Print() << "Interpolated elevation from binary LCP\n";
    } else {
        compute_elevation_from_landscape(elev_mf, geom, filename);
    }
}

// ============================================================================
// compute_fuel_model_from_lcp
// ============================================================================
[[maybe_unused]] inline void compute_fuel_model_from_lcp(MultiFab& fuel_mf,
                                        const Geometry& geom,
                                        const std::string& filename)
{
    if (is_binary_lcp(filename)) {
        LcpGridInfo info;
        std::vector<double> xd, yd;
        std::vector<float>  elev_v, slope_f, aspect_f, fuel_f, cc, ht, cbh, cbd;
        std::vector<int>    fuel_i;
        if (!read_landscape_binary(filename, info, xd, yd, elev_v, slope_f,
                                   aspect_f, fuel_f, fuel_i, cc, ht, cbh, cbd))
            amrex::Abort("Failed to read binary LCP: " + filename);
        interp_flat_to_multifab(fuel_mf, geom, xd, yd, fuel_f, /*nearest=*/true);
        amrex::Print() << "Interpolated fuel model from binary LCP\n";
    } else {
        compute_fuel_model_from_landscape(fuel_mf, geom, filename);
    }
}

// ============================================================================
// compute_landscape_slope_aspect_fuel_lcp
//
// Populates slope [deg], aspect [deg], and fuel model MultiFabs from either
// a binary or ASCII landscape file.  Drop-in replacement for
// compute_landscape_slope_aspect_fuel() that also handles binary LCP.
// ============================================================================
[[maybe_unused]] inline void compute_landscape_slope_aspect_fuel_lcp(
    MultiFab& slope_out,
    MultiFab& aspect_out,
    MultiFab& fuel_out,
    const Geometry& geom,
    const std::string& filename)
{
    if (is_binary_lcp(filename)) {
        LcpGridInfo info;
        std::vector<double> xd, yd;
        std::vector<float>  elev_v, slope_f, aspect_f, fuel_f, cc, ht, cbh, cbd;
        std::vector<int>    fuel_i;
        if (!read_landscape_binary(filename, info, xd, yd, elev_v, slope_f,
                                   aspect_f, fuel_f, fuel_i, cc, ht, cbh, cbd))
            amrex::Abort("Failed to read binary LCP: " + filename);
        interp_flat_to_multifab(slope_out,  geom, xd, yd, slope_f,  false);
        interp_flat_to_multifab(aspect_out, geom, xd, yd, aspect_f, false);
        interp_flat_to_multifab(fuel_out,   geom, xd, yd, fuel_f,   /*nearest=*/true);
        amrex::Print() << "Interpolated slope/aspect/fuel from binary LCP\n";
    } else {
        compute_landscape_slope_aspect_fuel(slope_out, aspect_out, fuel_out,
                                            geom, filename);
    }
}

// ============================================================================
// compute_crown_layers_from_lcp
//
// Populates optional MultiFabs (all 1-component) with spatially-varying
// crown fuel properties from a binary LCP file:
//   cbh_mf          – crown base height [m]
//   cbd_mf          – crown bulk density [kg/m³]
//   cc_mf           – canopy cover [%]
//   canopy_height_mf – stand (canopy) height [m]   (layer 5 of binary LCP)
//
// Returns true if the file has crown layers (has_crown == true) and the
// MultiFabs were populated; false otherwise (MultiFabs are untouched).
//
// All MultiFab pointers may be nullptr if the caller does not need that layer.
// ============================================================================
static bool compute_crown_layers_from_lcp(
    MultiFab*          cbh_mf,
    MultiFab*          cbd_mf,
    MultiFab*          cc_mf,
    const Geometry&    geom,
    const std::string& filename,
    MultiFab*          canopy_height_mf = nullptr)
{
    if (!is_binary_lcp(filename)) {
        return false; // ASCII format has no crown layers
    }

    LcpGridInfo info;
    std::vector<double> xd, yd;
    std::vector<float>  elev_v, slope_f, aspect_f, fuel_f, cc_v, ht_v, cbh_v, cbd_v;
    std::vector<int>    fuel_i;

    if (!read_landscape_binary(filename, info, xd, yd, elev_v, slope_f,
                               aspect_f, fuel_f, fuel_i, cc_v, ht_v, cbh_v, cbd_v))
        amrex::Abort("Failed to read binary LCP: " + filename);

    // Canopy cover is always available from binary LCP (layer 4)
    if (cc_mf && !cc_v.empty())
        interp_flat_to_multifab(*cc_mf, geom, xd, yd, cc_v, false);

    if (!info.has_crown) {
        amrex::Print() << "Binary LCP: no crown fuel layers (CrownFuels = 0)\n";
        return false;
    }

    if (cbh_mf && !cbh_v.empty())
        interp_flat_to_multifab(*cbh_mf, geom, xd, yd, cbh_v, false);
    if (cbd_mf && !cbd_v.empty())
        interp_flat_to_multifab(*cbd_mf, geom, xd, yd, cbd_v, false);
    if (canopy_height_mf && !ht_v.empty())
        interp_flat_to_multifab(*canopy_height_mf, geom, xd, yd, ht_v, false);

    amrex::Print() << "Populated spatially-varying crown layers (CBH, CBD, CC, "
                      "canopy height) from binary LCP\n";
    return true;
}

#endif // LANDSCAPE_FILE_H
