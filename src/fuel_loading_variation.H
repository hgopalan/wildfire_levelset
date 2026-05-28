// ============================================================================
// Fuel Bed Bulk Density Spatial Variation Model
// ============================================================================
// Allows spatially varying fuel loading from landscape files. Current models
// use uniform fuel properties per fuel model, but real landscapes have patchy
// fuel loading that affects local ROS and intensity.
//
// Physics:
//   - Fuel loading (w0) varies spatially even within the same fuel type
//   - Patchy distribution due to terrain, moisture, past fires, etc.
//   - Higher local loading → higher intensity and potentially higher ROS
//   - Lower local loading → lower intensity and potentially lower ROS
//
// Implementation:
//   - Spatial multiplier field applied to base fuel loading
//   - w0_eff = w0_base × multiplier
//   - Can be loaded from LCP auxiliary data or separate file
//   - Simple extension to existing landscape reader
//
// References:
//   Finney, M.A. (1998). FARSITE: Fire Area Simulator.
//     USDA Forest Service RMRS-RP-4.
//   Andrews, P.L. (2018). The Rothermel surface fire spread model and associated
//     developments: A comprehensive explanation. USDA Forest Service RMRS-GTR-371.
// ============================================================================

#ifndef FUEL_LOADING_VARIATION_H
#define FUEL_LOADING_VARIATION_H

#include "constants.H"
#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_Geometry.H>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace amrex;

// Default parameters for fuel loading variation
namespace FuelLoadingConst {
    constexpr Real MULTIPLIER_MIN = 0.1;   // Minimum loading multiplier
    constexpr Real MULTIPLIER_MAX = 3.0;   // Maximum loading multiplier
    constexpr Real MULTIPLIER_DEFAULT = 1.0;  // Default (no variation)
}

// ============================================================================
// Apply spatial fuel loading multiplier
//
// w0_eff = w0_base × multiplier
//
// Parameters:
//   w0_eff_mf      – output: effective fuel loading [lb/ft²]
//   w0_base_mf     – input: base fuel loading [lb/ft²]
//   multiplier_mf  – input: spatial multiplier [dimensionless]
// ============================================================================
inline void apply_fuel_loading_multiplier(MultiFab& w0_eff_mf,
                                          const MultiFab& w0_base_mf,
                                          const MultiFab& multiplier_mf)
{
    for (MFIter mfi(w0_eff_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto w0_eff    = w0_eff_mf.array(mfi);
        auto const w0_base = w0_base_mf.const_array(mfi);
        auto const mult    = multiplier_mf.const_array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            // Apply multiplier with bounds checking
            Real m = mult(i, j, k);
            m = amrex::max(m, FuelLoadingConst::MULTIPLIER_MIN);
            m = amrex::min(m, FuelLoadingConst::MULTIPLIER_MAX);
            
            w0_eff(i, j, k) = w0_base(i, j, k) * m;
        });
    }
}

// ============================================================================
// Apply fuel loading multiplier in-place
//
// Modifies the fuel loading field directly.
//
// Parameters:
//   w0_mf         – input/output: fuel loading [lb/ft²]
//   multiplier_mf – input: spatial multiplier [dimensionless]
// ============================================================================
inline void apply_fuel_loading_multiplier_inplace(MultiFab& w0_mf,
                                                   const MultiFab& multiplier_mf)
{
    for (MFIter mfi(w0_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto w0        = w0_mf.array(mfi);
        auto const mult = multiplier_mf.const_array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real m = mult(i, j, k);
            m = amrex::max(m, FuelLoadingConst::MULTIPLIER_MIN);
            m = amrex::min(m, FuelLoadingConst::MULTIPLIER_MAX);
            
            w0(i, j, k) *= m;
        });
    }
}

// ============================================================================
// Initialize uniform fuel loading multiplier
//
// Sets all cells to a constant multiplier value (typically 1.0 for no variation).
//
// Parameters:
//   multiplier_mf  – output: spatial multiplier [dimensionless]
//   value          – uniform multiplier value
// ============================================================================
inline void init_uniform_multiplier(MultiFab& multiplier_mf,
                                    Real value = FuelLoadingConst::MULTIPLIER_DEFAULT)
{
    multiplier_mf.setVal(value);
}

// ============================================================================
// Load fuel loading multiplier from ASCII file
//
// File format: X Y multiplier (whitespace or comma-separated)
// Lines starting with '#' are ignored.
//
// Uses Wendland C² IDW interpolation (same as terrain/wind files).
//
// Parameters:
//   multiplier_mf  – output: spatial multiplier field [dimensionless]
//   geom           – AMReX Geometry object
//   filename       – path to multiplier file
//   radius         – interpolation radius [m] (default: 100.0)
//   default_value  – default multiplier for points with no data
// ============================================================================
#if (AMREX_SPACEDIM == 2)
inline void load_fuel_multiplier_from_file(MultiFab& multiplier_mf,
                                           const Geometry& geom,
                                           const std::string& filename,
                                           Real radius = Real(100.0),
                                           Real default_value = FuelLoadingConst::MULTIPLIER_DEFAULT)
{
    // Read multiplier data points from file
    std::vector<Real> x_pts, y_pts, mult_pts;
    
    std::ifstream infile(filename);
    if (!infile) {
        amrex::Print() << "Warning: Cannot open fuel multiplier file: " 
                       << filename << "\n";
        amrex::Print() << "Using default multiplier = " << default_value << "\n";
        multiplier_mf.setVal(default_value);
        return;
    }
    
    std::string line;
    while (std::getline(infile, line)) {
        // Skip comments and empty lines
        if (line.empty() || line[0] == '#') continue;
        
        // Replace commas with spaces
        std::replace(line.begin(), line.end(), ',', ' ');
        
        std::istringstream iss(line);
        Real x, y, mult;
        if (iss >> x >> y >> mult) {
            x_pts.push_back(x);
            y_pts.push_back(y);
            mult_pts.push_back(mult);
        }
    }
    infile.close();
    
    if (x_pts.empty()) {
        amrex::Print() << "Warning: No data in fuel multiplier file\n";
        multiplier_mf.setVal(default_value);
        return;
    }
    
    amrex::Print() << "Loaded " << x_pts.size() 
                   << " fuel multiplier points from " << filename << "\n";
    
    // Wendland C² kernel parameters
    const Real h = radius;
    const Real h2 = h * h;
    
    // Domain bounds
    const auto& prob_lo = geom.ProbLoArray();
    const auto& prob_hi = geom.ProbHiArray();
    const auto dx = geom.CellSizeArray();
    
    // Interpolate to grid
    for (MFIter mfi(multiplier_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto mult = multiplier_mf.array(mfi);
        
        const int n_pts = static_cast<int>(x_pts.size());
        
        // Copy data to device (for GPU)
        Gpu::DeviceVector<Real> d_x(x_pts.begin(), x_pts.end());
        Gpu::DeviceVector<Real> d_y(y_pts.begin(), y_pts.end());
        Gpu::DeviceVector<Real> d_mult(mult_pts.begin(), mult_pts.end());
        
        Real const* px = d_x.data();
        Real const* py = d_y.data();
        Real const* pm = d_mult.data();
        
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            // Cell center coordinates
            const Real x = prob_lo[0] + (Real(i) + Real(0.5)) * dx[0];
            const Real y = prob_lo[1] + (Real(j) + Real(0.5)) * dx[1];
            
            // Wendland C² IDW interpolation
            Real sum_weights = Real(0.0);
            Real sum_weighted_mult = Real(0.0);
            
            for (int p = 0; p < n_pts; ++p) {
                const Real dx_pt = x - px[p];
                const Real dy_pt = y - py[p];
                const Real r2 = dx_pt * dx_pt + dy_pt * dy_pt;
                
                if (r2 < h2) {
                    const Real r = std::sqrt(r2);
                    const Real q = r / h;
                    // Wendland C²: w = (1-q)^4 × (4q + 1)
                    const Real onemq = Real(1.0) - q;
                    const Real onemq4 = onemq * onemq * onemq * onemq;
                    const Real w = onemq4 * (Real(4.0) * q + Real(1.0));
                    
                    sum_weights += w;
                    sum_weighted_mult += w * pm[p];
                }
            }
            
            // Set multiplier: IDW result or default
            if (sum_weights > Real(1e-12)) {
                Real m = sum_weighted_mult / sum_weights;
                // Apply bounds
                m = amrex::max(m, FuelLoadingConst::MULTIPLIER_MIN);
                m = amrex::min(m, FuelLoadingConst::MULTIPLIER_MAX);
                mult(i, j, k) = m;
            } else {
                mult(i, j, k) = default_value;
            }
        });
    }
}
#endif // AMREX_SPACEDIM == 2

// ============================================================================
// Generate random fuel loading variation field
//
// Creates spatially correlated random variation for testing or ensemble runs.
//
// Parameters:
//   multiplier_mf  – output: spatial multiplier field
//   mean           – mean multiplier value (default: 1.0)
//   std_dev        – standard deviation (default: 0.2)
//   correlation_length – spatial correlation length [m]
//   geom           – geometry for spatial coordinates
//   seed           – random seed
// ============================================================================
inline void generate_random_fuel_variation(MultiFab& multiplier_mf,
                                           Real mean = Real(1.0),
                                           Real std_dev = Real(0.2),
                                           Real correlation_length = Real(50.0),
                                           const Geometry& geom = Geometry(),
                                           unsigned int seed = 12345)
{
    // Simple implementation: uniform random for now
    // In production, would use spatially correlated noise (e.g., Gaussian random field)
    
    std::srand(seed);
    
    for (MFIter mfi(multiplier_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto mult = multiplier_mf.array(mfi);
        
        const Real m = mean;
        const Real s = std_dev;
        
        // Note: This is host-only for simplicity
        // For GPU, would need cuRAND or similar
        for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k) {
        for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
        for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
            // Box-Muller transform for Gaussian
            Real u1 = Real(std::rand()) / RAND_MAX;
            Real u2 = Real(std::rand()) / RAND_MAX;
            Real z = std::sqrt(Real(-2.0) * std::log(u1)) * std::cos(Real(2.0) * M_PI * u2);
            
            Real val = m + s * z;
            val = amrex::max(val, FuelLoadingConst::MULTIPLIER_MIN);
            val = amrex::min(val, FuelLoadingConst::MULTIPLIER_MAX);
            
            mult(i, j, k) = val;
        }}}
    }
}

#endif // FUEL_LOADING_VARIATION_H
