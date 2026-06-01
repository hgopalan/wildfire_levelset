// ============================================================================
// Grass Curing Model
// ============================================================================
// Implements grass curing percentage models for Australian and North American
// grassland fuels. Grass curing describes the fractional dead (cured) grass
// content, which significantly affects fuel moisture and fire behavior.
//
// Curing percentage (C) ranges from 0% (green, live grass) to 100% (fully
// cured/dry grass). Higher curing increases fire danger through:
//   1. Reduced effective fuel moisture (less live fuel dampening)
//   2. Increased available dead fuel load
//   3. Enhanced rate of spread in grass-dominated fuels
//
// Model Options:
//   1. **User-specified constant** (default: 100% fully cured)
//   2. **Luke & McArthur (1978) seasonal model** - sinusoidal annual cycle
//   3. **Moisture-dependent model** - curing from soil moisture deficit
//   4. **Growing degree day (GDD) model** - phenology-based curing
//
// Seasonal Model (Luke & McArthur 1978):
//   C(t) = C_mean + C_amp × cos(2π × (day_of_year - peak_day) / 365)
//
// where:
//   C_mean    - mean annual curing (default: 70%)
//   C_amp     - seasonal amplitude (default: 30%)
//   peak_day  - day of year of minimum curing (summer, default: Jan 15 = day 15)
//
// Moisture-Dependent Model (Cheney & Sullivan 2008):
//   C = 100 × (1 - exp(-k × (KBDI / 200)))
//
// where:
//   k    - calibration constant (default: 1.5)
//   KBDI - Keetch-Byram Drought Index (0-800)
//
// Effect on Rothermel Model:
//   - Adjusts live herbaceous fuel load: w_live_eff = w_live × (1 - C/100)
//   - Adjusts dead fuel load: w_dead_eff = w_dead + w_live × (C/100)
//   - Modifies live fuel moisture extinction calculation
//
// References:
//   Luke, R.H. & McArthur, A.G. (1978). "Bushfires in Australia."
//     Australian Government Publishing Service, Canberra.
//   Cheney, N.P. & Sullivan, A. (2008). "Grassfires: Fuel, Weather and Fire
//     Behaviour." CSIRO Publishing.
//   Noble, I.R., Bary, G.A.V., & Gill, A.M. (1980). "McArthur's fire-danger
//     meters expressed as equations." Australian Journal of Ecology, 5(2):201–203.
//   Anderson, H.E. (1982). "Aids to determining fuel models for estimating fire
//     behavior." USDA Forest Service General Technical Report INT-122.
// ============================================================================

#ifndef GRASS_CURING_MODEL_H
#define GRASS_CURING_MODEL_H

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_GpuLaunch.H>
#include <cmath>

using namespace amrex;

// ============================================================================
// Grass Curing Constants
// ============================================================================
namespace GrassCuringConst {
    constexpr Real C_DEFAULT = Real(100.0);     // Default fully cured [%]
    constexpr Real C_MEAN = Real(70.0);         // Mean annual curing [%]
    constexpr Real C_AMP = Real(30.0);          // Seasonal amplitude [%]
    constexpr int PEAK_GREEN_DAY = 15;          // Day of year of peak green (Jan 15)
    constexpr Real MOISTURE_K = Real(1.5);      // Moisture response calibration
    constexpr Real PI = Real(3.14159265358979323846);
}

// ============================================================================
// compute_seasonal_curing
//
// Computes grass curing percentage using a sinusoidal seasonal model.
//
// Luke & McArthur (1978) formulation:
//   C(t) = C_mean + C_amp × cos(2π × (day - peak_green_day) / 365)
//
// Parameters:
//   day_of_year - Julian day (1-365 or 1-366 for leap years)
//   C_mean      - mean annual curing percentage [%] (default: 70%)
//   C_amp       - seasonal amplitude [%] (default: 30%)
//   peak_green_day - day of minimum curing (default: 15 = mid-Jan Southern Hemisphere)
//
// Returns: curing percentage [%] (0-100)
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real compute_seasonal_curing(int day_of_year,
                              Real C_mean = GrassCuringConst::C_MEAN,
                              Real C_amp = GrassCuringConst::C_AMP,
                              int peak_green_day = GrassCuringConst::PEAK_GREEN_DAY)
{
    // Clamp day to valid range
    day_of_year = amrex::max(1, amrex::min(day_of_year, 366));
    
    // Compute phase offset in radians
    Real phase = Real(2.0) * GrassCuringConst::PI * Real(day_of_year - peak_green_day) / Real(365.0);
    
    // Sinusoidal curing model
    Real curing = C_mean + C_amp * std::cos(phase);
    
    // Clamp to physical range [0, 100]
    return amrex::max(Real(0.0), amrex::min(curing, Real(100.0)));
}

// ============================================================================
// compute_moisture_dependent_curing
//
// Computes grass curing from soil moisture deficit (KBDI).
//
// Cheney & Sullivan (2008) exponential model:
//   C = 100 × (1 - exp(-k × KBDI / 200))
//
// Parameters:
//   KBDI - Keetch-Byram Drought Index (0-800)
//   k    - calibration constant (default: 1.5)
//
// Returns: curing percentage [%] (0-100)
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real compute_moisture_dependent_curing(Real KBDI, Real k = GrassCuringConst::MOISTURE_K)
{
    // Clamp KBDI to valid range
    KBDI = amrex::max(Real(0.0), amrex::min(KBDI, Real(800.0)));
    
    // Exponential curing response to drought
    Real curing = Real(100.0) * (Real(1.0) - std::exp(-k * KBDI / Real(200.0)));
    
    // Clamp to physical range [0, 100]
    return amrex::max(Real(0.0), amrex::min(curing, Real(100.0)));
}

// ============================================================================
// compute_gdd_curing
//
// Computes grass curing from accumulated growing degree days (GDD).
//
// Phenology-based model suitable for temperate grasslands:
//   C = 100 × min(1, GDD / GDD_threshold)
//
// Parameters:
//   GDD_accum     - accumulated growing degree days since green-up [°C·days]
//   GDD_threshold - GDD for full curing (default: 1200 °C·days)
//   T_base        - base temperature for GDD calculation (default: 5°C)
//
// Returns: curing percentage [%] (0-100)
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real compute_gdd_curing(Real GDD_accum, Real GDD_threshold = Real(1200.0))
{
    // Curing increases linearly with GDD until fully cured
    Real curing = Real(100.0) * amrex::min(Real(1.0), GDD_accum / GDD_threshold);
    
    // Clamp to physical range [0, 100]
    return amrex::max(Real(0.0), amrex::min(curing, Real(100.0)));
}

// ============================================================================
// apply_curing_to_fuel_load
//
// Adjusts live and dead herbaceous fuel loads based on curing percentage.
//
// Partitioning model:
//   w_dead_eff = w_dead_base + w_herb_total × (C / 100)
//   w_live_eff = w_herb_total × (1 - C / 100)
//
// Parameters:
//   w_dead_base  - base dead fuel load (non-herbaceous) [kg/m²]
//   w_herb_total - total herbaceous fuel load [kg/m²]
//   curing_pct   - grass curing percentage [%] (0-100)
//   w_dead_out   - output: adjusted dead fuel load [kg/m²]
//   w_live_out   - output: adjusted live fuel load [kg/m²]
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void apply_curing_to_fuel_load(Real w_dead_base,
                                Real w_herb_total,
                                Real curing_pct,
                                Real& w_dead_out,
                                Real& w_live_out)
{
    // Clamp curing to valid range
    curing_pct = amrex::max(Real(0.0), amrex::min(curing_pct, Real(100.0)));
    
    // Fraction of herbaceous fuel that is cured (dead)
    Real cured_fraction = curing_pct / Real(100.0);
    
    // Partition herbaceous fuel between dead and live
    w_dead_out = w_dead_base + w_herb_total * cured_fraction;
    w_live_out = w_herb_total * (Real(1.0) - cured_fraction);
}

// ============================================================================
// update_curing_field
//
// Updates a MultiFab of grass curing values using the selected model.
//
// Parameters:
//   curing_mf    - output MultiFab to store curing [%]
//   model        - curing model: "constant", "seasonal", "moisture", "gdd"
//   day_of_year  - current Julian day (for seasonal model)
//   kbdi_mf      - KBDI field (for moisture model; can be nullptr)
//   gdd_mf       - GDD field (for gdd model; can be nullptr)
//   C_constant   - constant curing value [%] (for constant model)
// ============================================================================
inline void update_curing_field(MultiFab& curing_mf,
                                 const std::string& model,
                                 int day_of_year,
                                 const MultiFab* kbdi_mf = nullptr,
                                 const MultiFab* gdd_mf = nullptr,
                                 Real C_constant = GrassCuringConst::C_DEFAULT)
{
    for (MFIter mfi(curing_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        Array4<Real> const& curing = curing_mf.array(mfi);
        
        if (model == "seasonal") {
            // Seasonal model - uniform across domain
            Real C_seasonal = compute_seasonal_curing(day_of_year);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                curing(i,j,k) = C_seasonal;
            });
            
        } else if (model == "moisture" && kbdi_mf != nullptr) {
            // Moisture-dependent model - spatially varying from KBDI
            Array4<Real const> const& kbdi = kbdi_mf->const_array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                curing(i,j,k) = compute_moisture_dependent_curing(kbdi(i,j,k));
            });
            
        } else if (model == "gdd" && gdd_mf != nullptr) {
            // GDD-based model - spatially varying from accumulated GDD
            Array4<Real const> const& gdd = gdd_mf->const_array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                curing(i,j,k) = compute_gdd_curing(gdd(i,j,k));
            });
            
        } else {
            // Constant model (default)
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                curing(i,j,k) = C_constant;
            });
        }
    }
}

#endif // GRASS_CURING_MODEL_H
