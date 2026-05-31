// ============================================================================
// Keetch-Byram Drought Index (KBDI)
// ============================================================================
// Implements the Keetch-Byram Drought Index, a widely-used measure of soil
// moisture deficit and fire potential in the southeastern United States.
//
// The KBDI ranges from 0 (no moisture deficit) to 800 (extreme drought).
// It represents the net effect of evapotranspiration and precipitation on
// soil moisture depletion in the upper duff and soil layers.
//
// KBDI formula (Keetch & Byram 1968):
//   dQ/dt = (Q_max - Q) × DF × ET / Q_max
//
// where:
//   Q     - current KBDI value (0-800)
//   Q_max - maximum drought index = 800 (mm or 8 inches of water)
//   DF    - drought factor (function of current Q and recent rainfall)
//   ET    - daily evapotranspiration (function of temperature)
//
// Daily precipitation reduces Q:
//   Q_new = Q_old - (rainfall - 5mm)  [only if rainfall > 5mm]
//
// The KBDI is commonly used to:
//   - Assess fire danger and potential fire behavior
//   - Determine live fuel moisture content
//   - Calculate the McArthur Drought Factor for FFDI
//   - Initialize the Canadian FWI System DC code
//
// KBDI interpretation:
//     0-200:  Soil and duff layers are moist; fire potential low
//   200-400:  Soil and duff layers drying; fire potential increasing
//   400-600:  Soil moisture significantly depleted; high fire potential
//   600-800:  Severe drought; extreme fire potential
//
// References:
//   Keetch, J.J. & Byram, G.M. (1968). "A drought index for forest fire control."
//     USDA Forest Service Research Paper SE-38.
//   Janis, M.J., et al. (2002). "Monitoring the Keetch-Byram Drought Index (KBDI)
//     using the National Weather Service forecast model."
//     International Journal of Wildland Fire, 11:217-222.
//   Noble, I.R., Bary, G.A.V., & Gill, A.M. (1980). "McArthur's fire-danger
//     meters expressed as equations." Australian Journal of Ecology, 5(2):201-203.
// ============================================================================

#ifndef KBDI_INDEX_H
#define KBDI_INDEX_H

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_GpuLaunch.H>
#include <cmath>

using namespace amrex;

// ============================================================================
// KBDI State Structure
// Stores the current KBDI value and accumulated rainfall for tracking
// ============================================================================
struct KBDIState {
    Real Q;              // Current KBDI value (0-800)
    Real days_since_rain; // Days since last significant rainfall (for tracking)
};

// ============================================================================
// Default initial conditions for KBDI
// ============================================================================
namespace KBDIConst {
    constexpr Real Q_DEFAULT = 100.0;  // Typical moderate startup (0-800 range)
    constexpr Real Q_MAX = 800.0;      // Maximum KBDI value (mm of water deficit)
    constexpr Real RAIN_THRESHOLD = 5.0; // Minimum rainfall [mm] to reduce KBDI
}

// ============================================================================
// compute_kbdi_daily_increment
//
// Computes the daily increase in KBDI due to evapotranspiration.
//
// Keetch & Byram (1968) formulation:
//   dQ = (Q_max - Q) × DF × ET / Q_max
//
// where DF is a drought factor and ET is evapotranspiration.
//
// Simplified formulation (Janis et al. 2002):
//   dQ = [800 - Q] × [0.968 × exp(0.0875×T + 1.5552) - 8.3] × 10^(-3) / [1 + 10.88×exp(-0.001736×RA)]
//
// Parameters:
//   Q          - current KBDI value (0-800)
//   T_max      - daily maximum temperature [°C]
//   precip_mm  - daily precipitation [mm]
//   annual_precip_mm - mean annual precipitation [mm] (default: 1000 mm)
//
// Returns: updated KBDI value
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real compute_kbdi_daily_increment(Real Q, Real T_max, Real precip_mm, Real annual_precip_mm = Real(1000.0))
{
    // Clamp Q to valid range [0, 800]
    Q = amrex::max(Real(0.0), amrex::min(Q, Real(800.0)));
    
    // Apply rainfall reduction first (if rainfall > threshold)
    if (precip_mm > Real(5.0)) {
        // Net rainfall reduces KBDI (Keetch & Byram 1968)
        // Q_new = Q_old - (rain - 5mm), but don't go below 0
        Real reduction = precip_mm - Real(5.0);
        Q = amrex::max(Real(0.0), Q - reduction);
        // After significant rain, return current Q (no ET increment this day)
        return Q;
    }
    
    // Compute daily evapotranspiration factor
    // Based on Keetch & Byram (1968) with Janis et al. (2002) refinements
    
    // Evapotranspiration component (temperature-dependent)
    // ET factor increases exponentially with temperature
    Real ET_factor = Real(0.968) * std::exp(Real(0.0875) * T_max + Real(1.5552)) - Real(8.3);
    
    // Clamp ET_factor to positive values
    if (ET_factor < Real(0.0)) {
        ET_factor = Real(0.0);
    }
    
    // Drought factor (precipitation-dependent)
    // DF decreases with higher annual precipitation (wetter climates)
    // Using annual precipitation RA in mm
    Real DF = Real(1.0) / (Real(1.0) + Real(10.88) * std::exp(Real(-0.001736) * annual_precip_mm));
    
    // Daily KBDI increment
    // dQ = (800 - Q) × DF × ET × 10^(-3)
    Real dQ = (Real(800.0) - Q) * DF * ET_factor * Real(0.001);
    
    // Update KBDI
    Real Q_new = Q + dQ;
    
    // Clamp to valid range [0, 800]
    Q_new = amrex::max(Real(0.0), amrex::min(Q_new, Real(800.0)));
    
    return Q_new;
}

// ============================================================================
// compute_drought_factor_from_kbdi
//
// Converts KBDI to the McArthur Drought Factor (DF) used in the Australian
// Forest Fire Danger Index (FFDI).
//
// Empirical relationship (Noble et al. 1980):
//   DF = 10.5 × [1 - exp(-KBDI/800)] × [1 - 0.25 × exp(-KBDI/200)]
//
// This provides a mapping from KBDI to the 0-10 drought factor scale
// used in FFDI calculations.
//
// Parameters:
//   KBDI - Keetch-Byram Drought Index (0-800)
//
// Returns: McArthur Drought Factor (0-10)
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real compute_drought_factor_from_kbdi(Real KBDI)
{
    // Clamp KBDI to valid range
    Real Q = amrex::max(Real(0.0), amrex::min(KBDI, Real(800.0)));
    
    // Empirical McArthur Drought Factor formula
    // DF ranges from ~0 (wet) to ~10 (extreme drought)
    Real term1 = Real(1.0) - std::exp(-Q / Real(800.0));
    Real term2 = Real(1.0) - Real(0.25) * std::exp(-Q / Real(200.0));
    Real DF = Real(10.5) * term1 * term2;
    
    // Clamp to valid range [0, 10]
    DF = amrex::max(Real(0.0), amrex::min(DF, Real(10.0)));
    
    return DF;
}

// ============================================================================
// Spatially-varying KBDI computation over a MultiFab
//
// Updates KBDI field based on daily maximum temperature and precipitation.
// This function should be called once per day (not every timestep).
//
// Parameters:
//   kbdi_mf         – input/output: KBDI field (0-800)
//   temperature_mf  – input: daily maximum temperature [°C]
//   precip_mf       – input: daily precipitation [mm]
//   annual_precip   – mean annual precipitation [mm] (default: 1000)
// ============================================================================
inline void update_kbdi_field(
    MultiFab& kbdi_mf,
    const MultiFab& temperature_mf,
    const MultiFab& precip_mf,
    Real annual_precip = Real(1000.0))
{
    for (MFIter mfi(kbdi_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto kbdi = kbdi_mf.array(mfi);
        auto const T_max = temperature_mf.const_array(mfi);
        auto const precip = precip_mf.const_array(mfi);
        const Real RA = annual_precip;
        
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real Q_old = kbdi(i, j, k);
            Real T = T_max(i, j, k);
            Real P = precip(i, j, k);
            
            // Compute updated KBDI
            Real Q_new = compute_kbdi_daily_increment(Q_old, T, P, RA);
            
            kbdi(i, j, k) = Q_new;
        });
    }
}

// ============================================================================
// Uniform (scalar) KBDI computation
//
// Updates a single KBDI value based on uniform temperature and precipitation.
// Useful for small domains or when spatial variation is not needed.
//
// Parameters:
//   kbdi_state      – input/output: KBDI state structure
//   T_max           – daily maximum temperature [°C]
//   precip_mm       – daily precipitation [mm]
//   annual_precip   – mean annual precipitation [mm] (default: 1000)
// ============================================================================
inline void update_kbdi_scalar(
    KBDIState& kbdi_state,
    Real T_max,
    Real precip_mm,
    Real annual_precip = Real(1000.0))
{
    kbdi_state.Q = compute_kbdi_daily_increment(kbdi_state.Q, T_max, precip_mm, annual_precip);
    
    // Update days since rain tracker
    if (precip_mm > Real(5.0)) {
        kbdi_state.days_since_rain = Real(0.0);
    } else {
        kbdi_state.days_since_rain += Real(1.0);
    }
}

#endif // KBDI_INDEX_H
