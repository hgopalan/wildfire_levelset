// ============================================================================
// Per-Cell Live Canopy Moisture from FARSITE .fms Scenario File
// ============================================================================
// Reads a FARSITE fuel moisture scenario (.fms) file and stores per-fuel-model
// dead and live moisture values in a spatially-varying MultiFab, so that each
// cell uses the moisture appropriate for its fuel model rather than the single
// global Rothermel moisture values.
//
// FMS file format (whitespace-delimited):
//   <num_fuel_models>
//   <fuel_model>  <1hr%>  <10hr%>  <100hr%>  <live_herb%>  <live_wood%>
//   ...
//
// Moisture values are integer percentages (8 → 8 % → 0.08 fraction).
// Non-burnable fuel model codes (91–99) are silently skipped.
//
// Output MultiFab (FMS_MOISTURE_NCOMP components per cell):
//   Component 0 – M_d1   [fraction]  1-hr  dead moisture
//   Component 1 – M_d10  [fraction]  10-hr dead moisture
//   Component 2 – M_d100 [fraction]  100-hr dead moisture
//   Component 3 – M_lh   [fraction]  live herbaceous moisture
//   Component 4 – M_lw   [fraction]  live woody moisture
//
// Usage:
//   MultiFab spatial_moisture_mf(ba, dm, FMS_MOISTURE_NCOMP, 0);
//   if (!inputs.fms_file.empty())
//       load_fms_spatial_moisture(spatial_moisture_mf, fuel_model_mf,
//                                  geom, inputs.fms_file, inputs.rothermel);
//   // pass spatial_moisture_mf to compute_rothermel_R or farsite_ellipse
//
// When the spatial moisture MultiFab is populated, the per-cell moisture
// values override the global RothermelParams for any cell whose fuel model
// has a record in the FMS file.  Cells with unrecognised fuel models retain
// the global Rothermel moisture values.
//
// Reference:
//   Finney, M.A. (2004). FARSITE: Fire Area Simulator — Model Development and
//     Evaluation. USDA Forest Service Research Paper RMRS-RP-4.
// ============================================================================

#ifndef FMS_MOISTURE_H
#define FMS_MOISTURE_H

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuContainers.H>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include "parse_inputs.H"

using namespace amrex;

// Number of moisture components stored per cell
constexpr int FMS_MOISTURE_NCOMP = 5;

// Non-burnable FBFM codes
static inline bool fms_is_nonburnable(int code) {
    return (code >= 91 && code <= 99) || code == 0;
}

// ============================================================================
// load_fms_spatial_moisture
//
// Parse the .fms file and fill spatial_moisture_mf with per-cell moisture
// values based on the per-cell fuel model code in fuel_model_mf.
//
// Parameters:
//   spatial_moisture_mf – output FMS_MOISTURE_NCOMP-component MultiFab
//   fuel_model_mf       – input per-cell fuel model code (from landscape)
//   geom                – AMReX geometry (not used; for future 3-D extension)
//   fms_file            – path to FARSITE .fms file
//   rp                  – global RothermelParams (provides fallback values)
// ============================================================================
inline void load_fms_spatial_moisture(
    MultiFab&              spatial_moisture_mf,
    const MultiFab&        fuel_model_mf,
    const Geometry&        /*geom*/,
    const std::string&     fms_file,
    const RothermelParams& rp)
{
    // -----------------------------------------------------------------------
    // Parse .fms file into a host-side map
    // -----------------------------------------------------------------------
    std::map<int, std::array<float, 5>> fms_map;  // code → {M_d1, M_d10, M_d100, M_lh, M_lw}

    std::ifstream f(fms_file);
    if (!f.is_open()) {
        amrex::Abort("FMS moisture: cannot open file: " + fms_file);
    }

    std::string line;
    // First non-comment line: number of fuel models
    int nmodels = 0;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#' || line[0] == '!') continue;
        std::istringstream iss(line);
        iss >> nmodels;
        break;
    }

    int parsed = 0;
    while (std::getline(f, line) && parsed < nmodels) {
        if (line.empty() || line[0] == '#' || line[0] == '!') continue;
        for (char& c : line) if (c == ',') c = ' ';
        std::istringstream iss(line);
        int code;
        float m1, m10, m100, mlh, mlw;
        if (iss >> code >> m1 >> m10 >> m100 >> mlh >> mlw) {
            if (!fms_is_nonburnable(code)) {
                fms_map[code] = {m1 / 100.0f, m10 / 100.0f, m100 / 100.0f,
                                 mlh / 100.0f, mlw / 100.0f};
            }
            ++parsed;
        }
    }
    amrex::Print() << "FMS moisture: loaded " << fms_map.size()
                   << " fuel models from " << fms_file << "\n";

    // -----------------------------------------------------------------------
    // Build device-accessible lookup arrays
    // -----------------------------------------------------------------------
    // MAX_CODE = 261 covers the full FBFM40 range (codes 1-260; Scott & Burgan 2005)
    // plus code 0 as an unused sentinel.  Non-burnable codes (91-99) fall within the
    // range and retain the global RothermelParams fallback values set at init.
    constexpr int MAX_CODE = 261;
    std::vector<float> h_m1(MAX_CODE, static_cast<float>(rp.M_d1));
    std::vector<float> h_m10(MAX_CODE, static_cast<float>(rp.M_d10));
    std::vector<float> h_m100(MAX_CODE, static_cast<float>(rp.M_d100));
    std::vector<float> h_mlh(MAX_CODE, static_cast<float>(rp.M_lh));
    std::vector<float> h_mlw(MAX_CODE, static_cast<float>(rp.M_lw));

    for (auto& [code, vals] : fms_map) {
        if (code >= 1 && code < MAX_CODE) {
            h_m1  [code] = vals[0];
            h_m10 [code] = vals[1];
            h_m100[code] = vals[2];
            h_mlh [code] = vals[3];
            h_mlw [code] = vals[4];
        }
    }

    amrex::Gpu::DeviceVector<float> d_m1(h_m1.size());
    amrex::Gpu::DeviceVector<float> d_m10(h_m10.size());
    amrex::Gpu::DeviceVector<float> d_m100(h_m100.size());
    amrex::Gpu::DeviceVector<float> d_mlh(h_mlh.size());
    amrex::Gpu::DeviceVector<float> d_mlw(h_mlw.size());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, h_m1.begin(),   h_m1.end(),   d_m1.begin());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, h_m10.begin(),  h_m10.end(),  d_m10.begin());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, h_m100.begin(), h_m100.end(), d_m100.begin());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, h_mlh.begin(),  h_mlh.end(),  d_mlh.begin());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, h_mlw.begin(),  h_mlw.end(),  d_mlw.begin());

    const float* p_m1   = d_m1.data();
    const float* p_m10  = d_m10.data();
    const float* p_m100 = d_m100.data();
    const float* p_mlh  = d_mlh.data();
    const float* p_mlw  = d_mlw.data();

    // -----------------------------------------------------------------------
    // Fill spatial_moisture_mf on GPU
    // -----------------------------------------------------------------------
    for (MFIter mfi(spatial_moisture_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto       sm  = spatial_moisture_mf.array(mfi);
        auto const fm  = fuel_model_mf.const_array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            int code = static_cast<int>(fm(i, j, k));
            if (code < 1 || code >= MAX_CODE) {
                // Fallback to global defaults already pre-filled
                return;
            }
            sm(i, j, k, 0) = static_cast<Real>(p_m1  [code]);
            sm(i, j, k, 1) = static_cast<Real>(p_m10 [code]);
            sm(i, j, k, 2) = static_cast<Real>(p_m100[code]);
            sm(i, j, k, 3) = static_cast<Real>(p_mlh [code]);
            sm(i, j, k, 4) = static_cast<Real>(p_mlw [code]);
        });
    }
    amrex::Print() << "FMS moisture: spatial moisture MultiFab populated\n";
}

#endif // FMS_MOISTURE_H
