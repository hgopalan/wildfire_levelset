// ============================================================================
// Chandler Burning Index (CBI)
// ============================================================================
// Implements the Chandler Burning Index, a fire danger rating developed for
// the U.S. National Weather Service and used in fire weather forecasts across
// the United States.
//
// The CBI is an empirical index that combines temperature, relative humidity,
// and wind speed into a single fire danger rating. It is widely used for
// public fire danger communication and fire weather watch/warning criteria.
//
// CBI Formula (Chandler et al. 1983):
//   CBI = [(110 - 1.373×RH - 0.54×(10.20 - T)] × [124 × (1 + U)^0.0001] / 60
//
// Simplified form:
//   CBI = [(110 - 1.373×RH) - 0.54×(10.20 - T)] × (1 + 0.0345×U)
//
// where:
//   T  - air temperature [°F] (converted from °C)
//   RH - relative humidity [%]
//   U  - wind speed [mph] (converted from m/s)
//
// CBI Interpretation (Chandler et al. 1983):
//     0-50:   Low fire danger
//    50-75:   Moderate fire danger
//    75-90:   High fire danger  (Fire Weather Watch criteria)
//    90-97.5: Very High fire danger (Red Flag Warning criteria)
//    97.5+:   Extreme fire danger (Extreme Red Flag Warning)
//
// Usage in Operational Fire Weather:
//   - Fire weather forecasts issue watches/warnings based on CBI thresholds
//   - Combined with other indices (Haines, KBDI) for comprehensive assessment
//   - Used for public communication of fire danger
//   - Guides fire restrictions and burn ban decisions
//
// Relationship to Other Indices:
//   - NFDRS: CBI correlates with BI (Burning Index) but uses simpler inputs
//   - KBDI: Often paired with CBI for drought + daily weather fire danger
//   - Haines: CBI (surface) + Haines (upper air) = complete fire danger picture
//
// References:
//   Chandler, C., Cheney, P., Thomas, P., Trabaud, L., & Williams, D. (1983).
//     "Fire in Forestry, Volume I: Forest Fire Behavior and Effects."
//     John Wiley & Sons, New York.
//   Schroeder, M.J. & Buck, C.C. (1970). "Fire Weather: A Guide for Application
//     of Meteorological Information to Forest Fire Control Operations."
//     USDA Forest Service Agriculture Handbook 360.
//   Fosberg, M.A. (1978). "Weather in wildland fire management: the fire weather
//     index." Proceedings Conference Sierra Nevada Meteorology, pp. 1-4.
// ============================================================================

#ifndef CHANDLER_BURNING_INDEX_H
#define CHANDLER_BURNING_INDEX_H

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_GpuLaunch.H>
#include <cmath>

using namespace amrex;

// ============================================================================
// Chandler Burning Index Constants
// ============================================================================
namespace ChandlerConst {
    // CBI formula coefficients
    constexpr Real RH_COEFF = Real(1.373);
    constexpr Real T_OFFSET = Real(10.20);  // °F
    constexpr Real T_COEFF = Real(0.54);
    constexpr Real BASE = Real(110.0);
    constexpr Real WIND_COEFF = Real(0.0345);
    
    // Unit conversions
    constexpr Real C_TO_F_SLOPE = Real(1.8);
    constexpr Real C_TO_F_OFFSET = Real(32.0);
    constexpr Real MPS_TO_MPH = Real(2.23694);
    
    // CBI thresholds for fire danger classes
    constexpr Real CBI_LOW = Real(50.0);
    constexpr Real CBI_MODERATE = Real(75.0);
    constexpr Real CBI_HIGH = Real(90.0);
    constexpr Real CBI_VERY_HIGH = Real(97.5);
}

// ============================================================================
// compute_chandler_burning_index
//
// Computes the Chandler Burning Index from temperature, RH, and wind speed.
//
// Formulation (Chandler et al. 1983):
//   CBI = [(110 - 1.373×RH) - 0.54×(10.20 - T_F)] × (1 + 0.0345×U_mph)
//
// Parameters:
//   temperature_C - air temperature [°C]
//   RH_pct        - relative humidity [%]
//   wind_speed_ms - wind speed [m/s]
//
// Returns: Chandler Burning Index (0-100+)
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real compute_chandler_burning_index(Real temperature_C,
                                     Real RH_pct,
                                     Real wind_speed_ms)
{
    // Convert temperature to Fahrenheit
    Real T_F = temperature_C * ChandlerConst::C_TO_F_SLOPE + ChandlerConst::C_TO_F_OFFSET;
    
    // Convert wind speed to mph
    Real U_mph = wind_speed_ms * ChandlerConst::MPS_TO_MPH;
    
    // Clamp RH to valid range [0, 100]
    RH_pct = amrex::max(Real(0.0), amrex::min(RH_pct, Real(100.0)));
    
    // Clamp wind to reasonable range [0, 50 mph]
    U_mph = amrex::max(Real(0.0), amrex::min(U_mph, Real(50.0)));
    
    // Temperature term
    Real temp_term = ChandlerConst::BASE - ChandlerConst::RH_COEFF * RH_pct 
                     - ChandlerConst::T_COEFF * (ChandlerConst::T_OFFSET - T_F);
    
    // Wind factor
    Real wind_factor = Real(1.0) + ChandlerConst::WIND_COEFF * U_mph;
    
    // Chandler Burning Index
    Real CBI = temp_term * wind_factor;
    
    // Clamp to reasonable range [0, 150]
    return amrex::max(Real(0.0), amrex::min(CBI, Real(150.0)));
}

// ============================================================================
// get_cbi_danger_class
//
// Returns fire danger class string based on CBI value.
//
// Parameters:
//   CBI - Chandler Burning Index value
//
// Returns: fire danger class (0-4)
//   0 - Low
//   1 - Moderate
//   2 - High (Fire Weather Watch)
//   3 - Very High (Red Flag Warning)
//   4 - Extreme (Extreme Red Flag Warning)
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int get_cbi_danger_class(Real CBI)
{
    if (CBI < ChandlerConst::CBI_LOW) {
        return 0;  // Low
    } else if (CBI < ChandlerConst::CBI_MODERATE) {
        return 1;  // Moderate
    } else if (CBI < ChandlerConst::CBI_HIGH) {
        return 2;  // High (Fire Weather Watch)
    } else if (CBI < ChandlerConst::CBI_VERY_HIGH) {
        return 3;  // Very High (Red Flag Warning)
    } else {
        return 4;  // Extreme (Extreme Red Flag Warning)
    }
}

// ============================================================================
// compute_fosberg_fire_weather_index
//
// Computes the Fosberg Fire Weather Index (FFWI), a related index.
//
// The Fosberg index is similar to CBI but uses a different formulation
// that emphasizes fuel moisture equilibrium.
//
// FFWI = η × sqrt(1 + U^2) / 0.3002
//
// where η is fuel moisture equilibrium damping.
//
// Parameters:
//   temperature_C - air temperature [°C]
//   RH_pct        - relative humidity [%]
//   wind_speed_ms - wind speed [m/s]
//
// Returns: Fosberg Fire Weather Index (0-100+)
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real compute_fosberg_fire_weather_index(Real temperature_C,
                                         Real RH_pct,
                                         Real wind_speed_ms)
{
    // Convert to mph
    Real U_mph = wind_speed_ms * ChandlerConst::MPS_TO_MPH;
    
    // Clamp RH to valid range
    RH_pct = amrex::max(Real(1.0), amrex::min(RH_pct, Real(100.0)));
    
    // Equilibrium moisture content (simplified Simard 1968)
    Real EMC = Real(0.03229) + Real(0.281073) * RH_pct 
               - Real(0.000578) * RH_pct * temperature_C;
    EMC = amrex::max(Real(1.0), amrex::min(EMC * Real(100.0), Real(50.0)));
    
    // Moisture damping
    Real eta = Real(1.0) - Real(2.0) * (EMC / Real(30.0)) 
               + Real(1.5) * std::pow(EMC / Real(30.0), 2.0) 
               - Real(0.5) * std::pow(EMC / Real(30.0), 3.0);
    eta = amrex::max(Real(0.0), eta);
    
    // Fosberg FFWI
    Real FFWI = eta * std::sqrt(Real(1.0) + U_mph * U_mph) / Real(0.3002);
    
    return amrex::max(Real(0.0), amrex::min(FFWI, Real(100.0)));
}

// ============================================================================
// compute_cbi_field
//
// Computes Chandler Burning Index for all cells in a MultiFab.
//
// Parameters:
//   cbi_mf          - output CBI field
//   temperature_mf  - temperature field [°C]
//   rh_mf           - relative humidity field [%]
//   wind_speed_mf   - wind speed field [m/s]
// ============================================================================
inline void compute_cbi_field(MultiFab& cbi_mf,
                               const MultiFab& temperature_mf,
                               const MultiFab& rh_mf,
                               const MultiFab& wind_speed_mf)
{
    for (MFIter mfi(cbi_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        Array4<Real> const& cbi = cbi_mf.array(mfi);
        Array4<Real const> const& temp = temperature_mf.const_array(mfi);
        Array4<Real const> const& rh = rh_mf.const_array(mfi);
        Array4<Real const> const& wind = wind_speed_mf.const_array(mfi);
        
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            cbi(i,j,k) = compute_chandler_burning_index(temp(i,j,k),
                                                         rh(i,j,k),
                                                         wind(i,j,k));
        });
    }
}

// ============================================================================
// compute_cbi_scalar
//
// Computes Chandler Burning Index for scalar inputs (single-point calculation).
//
// Useful for fire weather forecasts and point-based fire danger assessments.
//
// Parameters:
//   temperature_C - air temperature [°C]
//   RH_pct        - relative humidity [%]
//   wind_speed_ms - wind speed [m/s]
//   print_class   - if true, print fire danger class (default: false)
//
// Returns: Chandler Burning Index
// ============================================================================
inline Real compute_cbi_scalar(Real temperature_C,
                                Real RH_pct,
                                Real wind_speed_ms,
                                bool print_class = false)
{
    Real CBI = compute_chandler_burning_index(temperature_C, RH_pct, wind_speed_ms);
    
    if (print_class) {
        int danger_class = get_cbi_danger_class(CBI);
        const char* class_names[] = {"Low", "Moderate", "High", "Very High", "Extreme"};
        amrex::Print() << "Chandler Burning Index: " << CBI 
                       << " (" << class_names[danger_class] << " fire danger)\n";
    }
    
    return CBI;
}

#endif // CHANDLER_BURNING_INDEX_H
