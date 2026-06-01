// ============================================================================
// Elevation Temperature Lapse Rate
// ============================================================================
// Implements elevation-based temperature and pressure correction for terrain
// with significant topographic relief. Essential for accurate fire behavior
// prediction in mountainous regions.
//
// Temperature decreases with elevation due to adiabatic expansion:
//   T(z) = T_ref - Γ × (z - z_ref)
//
// where:
//   T_ref - reference temperature [°C] at reference elevation
//   z_ref - reference elevation [m]
//   z     - terrain elevation [m]
//   Γ     - environmental lapse rate [°C/m]
//
// Standard Atmosphere Lapse Rates:
//   - Dry adiabatic: Γ_d = 0.0098 °C/m (9.8 °C/km)
//   - Moist adiabatic: Γ_m = 0.0065 °C/m (6.5 °C/km) - standard atmosphere
//   - Observed mean: Γ_obs ≈ 0.006-0.007 °C/m (varies with season/humidity)
//
// Effects on Fire Behavior:
//   1. **Temperature**: Direct effect on Rothermel ROS
//   2. **Relative Humidity**: Increases with elevation (cooler air holds less moisture)
//   3. **Air Density**: Decreases with elevation (affects oxygen availability)
//   4. **Fuel Moisture**: Cooler temps at high elevation → higher fuel moisture
//
// Relative Humidity Adjustment:
//   RH(z) ≈ RH_ref × exp(0.0001 × Δz) for small Δz
//   or full Magnus formula with lapsed temperature
//
// Barometric Pressure:
//   P(z) = P_ref × exp(-g × M × Δz / (R × T_mean))
//   Simplified: P(z) = P_ref × (1 - 0.0065 × Δz / T_ref)^5.255
//
// References:
//   Wallace, J.M. & Hobbs, P.V. (2006). "Atmospheric Science: An Introductory
//     Survey." Academic Press, 2nd edition.
//   Finney, M.A. (1998). "FARSITE: Fire Area Simulator." USDA Forest Service
//     Research Paper RMRS-RP-4.
//   Tymstra, C., et al. (2010). "Development and structure of Prometheus."
//     Natural Resources Canada, Information Report NOR-X-417.
// ============================================================================

#ifndef ELEVATION_LAPSE_RATE_H
#define ELEVATION_LAPSE_RATE_H

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_GpuLaunch.H>
#include <cmath>

using namespace amrex;

// ============================================================================
// Lapse Rate Constants
// ============================================================================
namespace LapseRateConst {
    // Standard lapse rates [°C/m]
    constexpr Real LAPSE_DRY = Real(0.0098);      // Dry adiabatic lapse rate
    constexpr Real LAPSE_MOIST = Real(0.0065);    // Moist adiabatic (standard atmosphere)
    constexpr Real LAPSE_DEFAULT = Real(0.0065);  // Default (standard atmosphere)
    
    // Physical constants
    constexpr Real G = Real(9.80665);             // Gravitational acceleration [m/s²]
    constexpr Real M_AIR = Real(0.0289644);       // Molar mass of air [kg/mol]
    constexpr Real R_GAS = Real(8.31447);         // Universal gas constant [J/(mol·K)]
    constexpr Real P_SEA_LEVEL = Real(101325.0);  // Sea level pressure [Pa]
    
    // RH lapse coefficient (empirical)
    constexpr Real RH_LAPSE_COEFF = Real(0.0001); // RH increase per meter [1/m]
}

// ============================================================================
// apply_temperature_lapse
//
// Applies temperature lapse rate correction for elevation.
//
// Formulation:
//   T(z) = T_ref - Γ × (z - z_ref)
//
// Parameters:
//   T_ref      - reference temperature [°C] at reference elevation
//   elevation  - terrain elevation [m]
//   z_ref      - reference elevation [m] (default: 0 = sea level)
//   lapse_rate - environmental lapse rate [°C/m] (default: 0.0065 °C/m)
//
// Returns: lapsed temperature [°C]
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real apply_temperature_lapse(Real T_ref,
                              Real elevation,
                              Real z_ref = Real(0.0),
                              Real lapse_rate = LapseRateConst::LAPSE_DEFAULT)
{
    // Elevation difference
    Real delta_z = elevation - z_ref;
    
    // Apply lapse rate
    Real T_lapsed = T_ref - lapse_rate * delta_z;
    
    return T_lapsed;
}

// ============================================================================
// apply_rh_lapse
//
// Applies relative humidity correction for elevation using exponential model.
//
// Cooler air at higher elevations holds less moisture, increasing RH.
//
// Formulation:
//   RH(z) = RH_ref × exp(k × Δz)
//
// where k ≈ 0.0001 m^-1 (empirical coefficient).
//
// Parameters:
//   RH_ref     - reference relative humidity [%] at reference elevation
//   elevation  - terrain elevation [m]
//   z_ref      - reference elevation [m] (default: 0 = sea level)
//   k          - RH lapse coefficient [1/m] (default: 0.0001)
//
// Returns: lapsed relative humidity [%], clamped to [0, 100]
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real apply_rh_lapse(Real RH_ref,
                    Real elevation,
                    Real z_ref = Real(0.0),
                    Real k = LapseRateConst::RH_LAPSE_COEFF)
{
    // Elevation difference
    Real delta_z = elevation - z_ref;
    
    // Exponential RH increase with elevation
    Real RH_lapsed = RH_ref * std::exp(k * delta_z);
    
    // Clamp to physical range [0, 100]
    return amrex::max(Real(0.0), amrex::min(RH_lapsed, Real(100.0)));
}

// ============================================================================
// compute_barometric_pressure
//
// Computes barometric pressure as a function of elevation.
//
// Barometric formula (standard atmosphere):
//   P(z) = P_0 × (1 - L × z / T_0)^(g × M / (R × L))
//
// where:
//   P_0 - sea level pressure [Pa] (default: 101325 Pa)
//   L   - temperature lapse rate [K/m] (default: 0.0065 K/m)
//   T_0 - sea level temperature [K]
//   g   - gravitational acceleration [m/s²]
//   M   - molar mass of air [kg/mol]
//   R   - universal gas constant [J/(mol·K)]
//
// Parameters:
//   elevation    - terrain elevation [m]
//   T_sea_level  - sea level temperature [°C] (default: 15°C)
//   P_sea_level  - sea level pressure [Pa] (default: 101325 Pa)
//
// Returns: barometric pressure [Pa]
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real compute_barometric_pressure(Real elevation,
                                  Real T_sea_level = Real(15.0),
                                  Real P_sea_level = LapseRateConst::P_SEA_LEVEL)
{
    // Convert temperature to Kelvin
    Real T_0 = T_sea_level + Real(273.15);
    
    // Temperature lapse rate [K/m]
    Real L = LapseRateConst::LAPSE_DEFAULT;
    
    // Barometric formula exponent
    Real exponent = LapseRateConst::G * LapseRateConst::M_AIR / 
                    (LapseRateConst::R_GAS * L);
    
    // Pressure calculation
    // Minimum threshold prevents numerical issues in power calculation (stratospheric limit)
    constexpr Real min_temp_ratio = Real(0.01);
    Real temp_ratio = Real(1.0) - (L * elevation) / T_0;
    Real pressure = P_sea_level * std::pow(amrex::max(min_temp_ratio, temp_ratio), exponent);
    
    return pressure;
}

// ============================================================================
// apply_elevation_lapse_fields
//
// Applies elevation-based lapse rate corrections to temperature and RH fields.
//
// Parameters:
//   temperature_mf - temperature field [°C] (modified in-place)
//   rh_mf          - relative humidity field [%] (modified in-place)
//   elevation_mf   - terrain elevation field [m]
//   T_ref          - reference temperature [°C] at reference elevation
//   RH_ref         - reference RH [%] at reference elevation
//   z_ref          - reference elevation [m] (default: mean elevation)
//   lapse_rate     - temperature lapse rate [°C/m] (default: 0.0065)
// ============================================================================
inline void apply_elevation_lapse_fields(MultiFab& temperature_mf,
                                          MultiFab& rh_mf,
                                          const MultiFab& elevation_mf,
                                          Real T_ref,
                                          Real RH_ref,
                                          Real z_ref,
                                          Real lapse_rate = LapseRateConst::LAPSE_DEFAULT)
{
    for (MFIter mfi(temperature_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        Array4<Real> const& temp = temperature_mf.array(mfi);
        Array4<Real> const& rh = rh_mf.array(mfi);
        Array4<Real const> const& elev = elevation_mf.const_array(mfi);
        
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            // Apply temperature lapse
            temp(i,j,k) = apply_temperature_lapse(T_ref, elev(i,j,k), z_ref, lapse_rate);
            
            // Apply RH lapse
            rh(i,j,k) = apply_rh_lapse(RH_ref, elev(i,j,k), z_ref);
        });
    }
}

// ============================================================================
// compute_air_density_correction
//
// Computes air density correction factor for elevation.
//
// Air density affects oxygen availability and combustion efficiency.
//
// Ideal gas law:
//   ρ(z) / ρ_0 = (P(z) / P_0) × (T_0 / T(z))
//
// Parameters:
//   elevation    - terrain elevation [m]
//   T_sea_level  - sea level temperature [°C] (default: 15°C)
//
// Returns: air density ratio ρ(z) / ρ_0 (dimensionless)
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real compute_air_density_correction(Real elevation, Real T_sea_level = Real(15.0))
{
    // Sea level conditions
    Real P_0 = LapseRateConst::P_SEA_LEVEL;
    Real T_0 = T_sea_level + Real(273.15);  // Convert to Kelvin
    
    // Elevated conditions
    Real P_z = compute_barometric_pressure(elevation, T_sea_level);
    Real T_z = apply_temperature_lapse(T_sea_level, elevation, Real(0.0)) + Real(273.15);
    
    // Density ratio from ideal gas law
    Real rho_ratio = (P_z / P_0) * (T_0 / T_z);
    
    return rho_ratio;
}

#endif // ELEVATION_LAPSE_RATE_H
