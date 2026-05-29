// ============================================================================
// Fine Fuel Moisture Time-Lag Differential Equations
// ============================================================================
// Implements physically-based fuel moisture dynamics using time-lag differential
// equations following Nelson (2000) and Viney (1991). This provides continuous
// moisture evolution driven by ambient temperature, relative humidity, and
// rainfall, accounting for hysteresis and response time.
//
// Physics:
//   - Fuel moisture responds to equilibrium moisture content (EMC) with a time lag
//   - Time lag depends on fuel size class (1-hr, 10-hr, 100-hr, 1000-hr)
//   - Adsorption (wetting) and desorption (drying) follow different EMC curves
//   - Rainfall adds moisture with exponential penetration model
//   - Solar radiation and fuel temperature affect drying rate
//
// Governing differential equation:
//   dM/dt = (M_e - M) / τ + P(t)
//
// where:
//   M    - fuel moisture content [fraction, 0-1]
//   M_e  - equilibrium moisture content [fraction]
//   τ    - time-lag constant [hours]
//   P(t) - precipitation wetting term [1/hour]
//
// Time-lag constants (Van Wagner 1972, Nelson 2000):
//   τ_1hr   = 1 hour    (fine grasses, litter)
//   τ_10hr  = 10 hours  (small twigs, 0.25-1 inch)
//   τ_100hr = 100 hours (medium branches, 1-3 inch)
//   τ_1000hr = 1000 hours (large branches, logs, 3-8 inch)
//
// Equilibrium moisture content (Nelson 2000):
//   Adsorption (wetting):
//     M_e = 0.03229 + 0.2810×H + 0.4093×H² - 1.3560×H³ + 1.6596×H⁴ (H ≤ 0.9)
//   Desorption (drying):
//     M_e = 0.05800 + 0.1985×H + 0.6250×H² - 1.1830×H³ + 1.0570×H⁴ (H ≤ 0.9)
//   where H = RH / 100 (relative humidity as fraction)
//
// References:
//   Nelson, R.M. (2000). "Prediction of diurnal change in 10-h fuel stick
//     moisture content." Canadian Journal of Forest Research, 30, 1071-1087.
//   Viney, N.R. (1991). "A review of fine fuel moisture modelling."
//     International Journal of Wildland Fire, 1(4), 215-234.
//   Van Wagner, C.E. (1972). "Equilibrium moisture contents of some fine
//     forest fuels in eastern Canada." Canadian Forestry Service Report 7.
//   Simard, A.J. (1968). "The moisture content of forest fuels: I. A review
//     of the basic concepts." Forest Fire Research Institute Information
//     Report FF-X-14.
// ============================================================================

#ifndef FUEL_MOISTURE_TIMELAG_DE_H
#define FUEL_MOISTURE_TIMELAG_DE_H

#include "constants.H"
#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_GpuLaunch.H>
#include <cmath>

using namespace amrex;

// Default parameters for fuel moisture time-lag model
namespace FuelMoistureTimelagConst {
    // Time-lag constants [hours]
    constexpr Real TAU_1HR       = 1.0;      // 1-hour fuels
    constexpr Real TAU_10HR      = 10.0;     // 10-hour fuels
    constexpr Real TAU_100HR     = 100.0;    // 100-hour fuels
    constexpr Real TAU_1000HR    = 1000.0;   // 1000-hour fuels
    
    // Minimum/maximum moisture content bounds [fraction]
    constexpr Real M_MIN         = 0.01;     // 1% (very dry)
    constexpr Real M_MAX         = 0.40;     // 40% (saturated fine fuels)
    
    // EMC calculation parameters
    constexpr Real RH_MIN        = 0.01;     // Minimum RH for EMC (1%)
    constexpr Real RH_MAX        = 0.99;     // Maximum RH for EMC (99%)
    
    // Precipitation wetting parameters
    constexpr Real RAIN_WETTING_COEF = 0.01; // Exponential wetting coefficient
    constexpr Real MIN_RAIN_RATE     = 0.1;  // Minimum rain rate [mm/h] for wetting
    
    // Temperature adjustment factor (drying accelerates with higher T)
    constexpr Real TEMP_REFERENCE    = 20.0; // Reference temperature [°C]
    constexpr Real TEMP_FACTOR_K     = 0.015; // Temperature sensitivity [1/°C]
}

// ============================================================================
// compute_emc_adsorption
//
// Computes equilibrium moisture content (EMC) for adsorption (wetting).
//
// Model (Nelson 2000):
//   M_e = a₀ + a₁×H + a₂×H² + a₃×H³ + a₄×H⁴
//
// Coefficients for adsorption (valid for H ≤ 0.9):
//   a₀ = 0.03229, a₁ = 0.2810, a₂ = 0.4093, a₃ = -1.3560, a₄ = 1.6596
//
// Parameters:
//   RH - relative humidity [%] (0-100)
//
// Returns: EMC [fraction] (0-0.35)
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real compute_emc_adsorption(Real RH)
{
    using namespace amrex;
    
    // Convert RH% to fraction and clamp
    Real H = (RH / Real(100.0));
    H = amrex::min(amrex::max(H, FuelMoistureTimelagConst::RH_MIN),
                   FuelMoistureTimelagConst::RH_MAX);
    
    // Nelson (2000) adsorption polynomial
    constexpr Real a0 = Real(0.03229);
    constexpr Real a1 = Real(0.2810);
    constexpr Real a2 = Real(0.4093);
    constexpr Real a3 = Real(-1.3560);
    constexpr Real a4 = Real(1.6596);
    
    const Real H2 = H * H;
    const Real H3 = H2 * H;
    const Real H4 = H3 * H;
    
    Real M_e = a0 + a1*H + a2*H2 + a3*H3 + a4*H4;
    
    // Clamp to physical range
    M_e = amrex::max(Real(0.0), amrex::min(M_e, Real(0.35)));
    
    return M_e;
}

// ============================================================================
// compute_emc_desorption
//
// Computes equilibrium moisture content (EMC) for desorption (drying).
//
// Model (Nelson 2000):
//   M_e = d₀ + d₁×H + d₂×H² + d₃×H³ + d₄×H⁴
//
// Coefficients for desorption (valid for H ≤ 0.9):
//   d₀ = 0.05800, d₁ = 0.1985, d₂ = 0.6250, d₃ = -1.1830, d₄ = 1.0570
//
// Parameters:
//   RH - relative humidity [%] (0-100)
//
// Returns: EMC [fraction] (0-0.35)
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real compute_emc_desorption(Real RH)
{
    using namespace amrex;
    
    // Convert RH% to fraction and clamp
    Real H = (RH / Real(100.0));
    H = amrex::min(amrex::max(H, FuelMoistureTimelagConst::RH_MIN),
                   FuelMoistureTimelagConst::RH_MAX);
    
    // Nelson (2000) desorption polynomial
    constexpr Real d0 = Real(0.05800);
    constexpr Real d1 = Real(0.1985);
    constexpr Real d2 = Real(0.6250);
    constexpr Real d3 = Real(-1.1830);
    constexpr Real d4 = Real(1.0570);
    
    const Real H2 = H * H;
    const Real H3 = H2 * H;
    const Real H4 = H3 * H;
    
    Real M_e = d0 + d1*H + d2*H2 + d3*H3 + d4*H4;
    
    // Clamp to physical range
    M_e = amrex::max(Real(0.0), amrex::min(M_e, Real(0.35)));
    
    return M_e;
}

// ============================================================================
// compute_emc_with_hysteresis
//
// Computes EMC accounting for hysteresis (adsorption vs desorption).
//
// Logic:
//   - If M < M_e_desorption: fuel is drying, use desorption curve
//   - If M > M_e_adsorption: fuel is wetting, use adsorption curve
//   - If between curves: use current state to determine direction
//
// Parameters:
//   RH        - relative humidity [%]
//   M_current - current fuel moisture [fraction]
//
// Returns: EMC [fraction]
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real compute_emc_with_hysteresis(Real RH, Real M_current)
{
    using namespace amrex;
    
    const Real M_e_ads = compute_emc_adsorption(RH);
    const Real M_e_des = compute_emc_desorption(RH);
    
    // Determine wetting or drying based on current moisture
    if (M_current < M_e_des) {
        // Fuel is drying
        return M_e_des;
    } else if (M_current > M_e_ads) {
        // Fuel is wetting
        return M_e_ads;
    } else {
        // Between curves: interpolate or use average
        return (M_e_ads + M_e_des) * Real(0.5);
    }
}

// ============================================================================
// compute_temperature_correction_factor
//
// Computes temperature correction factor for time-lag constant.
//
// Higher temperature → faster drying → shorter effective time lag
//
// Model:
//   f_T = exp(-k × (T - T_ref))
//
// where:
//   T      - ambient temperature [°C]
//   T_ref  - reference temperature [°C] (default: 20°C)
//   k      - temperature sensitivity [1/°C] (default: 0.015)
//
// Parameters:
//   T - ambient temperature [°C]
//
// Returns: temperature correction factor [dimensionless, 0.5-2.0]
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real compute_temperature_correction_factor(Real T)
{
    using namespace amrex;
    
    const Real T_ref = FuelMoistureTimelagConst::TEMP_REFERENCE;
    const Real k     = FuelMoistureTimelagConst::TEMP_FACTOR_K;
    
    Real f_T = std::exp(-k * (T - T_ref));
    
    // Clamp to reasonable range (0.5× to 2× nominal time lag)
    f_T = amrex::min(amrex::max(f_T, Real(0.5)), Real(2.0));
    
    return f_T;
}

// ============================================================================
// compute_precipitation_wetting_rate
//
// Computes precipitation wetting rate P(t) [1/hour].
//
// Model (exponential wetting):
//   P(t) = k_w × rain_rate   for rain_rate > threshold
//        = 0                  for rain_rate ≤ threshold
//
// where:
//   rain_rate - precipitation intensity [mm/h]
//   k_w       - wetting coefficient [h/mm]
//
// Parameters:
//   rain_rate - precipitation intensity [mm/h]
//
// Returns: wetting rate [1/hour]
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real compute_precipitation_wetting_rate(Real rain_rate)
{
    using namespace amrex;
    
    if (rain_rate < FuelMoistureTimelagConst::MIN_RAIN_RATE) {
        return Real(0.0);
    }
    
    const Real k_w = FuelMoistureTimelagConst::RAIN_WETTING_COEF;
    
    return k_w * rain_rate;
}

// ============================================================================
// update_fuel_moisture_timelag_de
//
// Updates fuel moisture content using time-lag differential equation.
//
// Governing equation:
//   dM/dt = (M_e - M) / τ_eff + P(t)
//
// Discretized (forward Euler):
//   M_new = M + dt × [(M_e - M) / τ_eff + P(t)]
//
// where:
//   M       - fuel moisture [fraction]
//   M_e     - equilibrium moisture [fraction]
//   τ_eff   - effective time lag [hours]
//   P(t)    - precipitation wetting rate [1/hour]
//   dt      - time step [hours]
//
// Parameters:
//   M_current - current fuel moisture [fraction]
//   RH        - relative humidity [%]
//   T         - temperature [°C]
//   rain_rate - precipitation rate [mm/h]
//   dt_hours  - time step [hours]
//   tau_hours - base time-lag constant [hours]
//
// Returns: updated fuel moisture [fraction]
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real update_fuel_moisture_timelag_de(
    Real M_current,
    Real RH,
    Real T,
    Real rain_rate,
    Real dt_hours,
    Real tau_hours = FuelMoistureTimelagConst::TAU_1HR)
{
    using namespace amrex;
    
    // Ensure M_current is in valid range
    M_current = amrex::min(amrex::max(M_current, FuelMoistureTimelagConst::M_MIN),
                           FuelMoistureTimelagConst::M_MAX);
    
    // Compute equilibrium moisture with hysteresis
    const Real M_e = compute_emc_with_hysteresis(RH, M_current);
    
    // Compute temperature correction factor
    const Real f_T = compute_temperature_correction_factor(T);
    
    // Effective time-lag constant [hours]
    const Real tau_eff = tau_hours * f_T;
    
    // Compute precipitation wetting rate
    const Real P = compute_precipitation_wetting_rate(rain_rate);
    
    // Time-lag moisture change term
    const Real dM_dt_timelag = (M_e - M_current) / amrex::max(tau_eff, Real(0.1));
    
    // Total moisture change rate
    const Real dM_dt = dM_dt_timelag + P;
    
    // Forward Euler update
    Real M_new = M_current + dt_hours * dM_dt;
    
    // Clamp to physical range
    M_new = amrex::min(amrex::max(M_new, FuelMoistureTimelagConst::M_MIN),
                       FuelMoistureTimelagConst::M_MAX);
    
    return M_new;
}

// ============================================================================
// update_multi_class_fuel_moisture
//
// Updates moisture for multiple fuel size classes simultaneously.
//
// This function applies the time-lag DE to 1-hr, 10-hr, 100-hr, and
// 1000-hr fuel classes, each with their respective time-lag constants.
//
// Parameters:
//   M_1hr, M_10hr, M_100hr, M_1000hr - current moistures [fraction]
//   M_1hr_new, ... (outputs)          - updated moistures [fraction]
//   RH         - relative humidity [%]
//   T          - temperature [°C]
//   rain_rate  - precipitation rate [mm/h]
//   dt_hours   - time step [hours]
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void update_multi_class_fuel_moisture(
    Real  M_1hr,   Real  M_10hr,   Real  M_100hr,   Real  M_1000hr,
    Real& M_1hr_new, Real& M_10hr_new, Real& M_100hr_new, Real& M_1000hr_new,
    Real  RH,
    Real  T,
    Real  rain_rate,
    Real  dt_hours)
{
    using namespace FuelMoistureTimelagConst;
    
    // Update each size class with its respective time-lag constant
    M_1hr_new   = update_fuel_moisture_timelag_de(M_1hr,   RH, T, rain_rate, dt_hours, TAU_1HR);
    M_10hr_new  = update_fuel_moisture_timelag_de(M_10hr,  RH, T, rain_rate, dt_hours, TAU_10HR);
    M_100hr_new = update_fuel_moisture_timelag_de(M_100hr, RH, T, rain_rate, dt_hours, TAU_100HR);
    M_1000hr_new = update_fuel_moisture_timelag_de(M_1000hr, RH, T, rain_rate, dt_hours, TAU_1000HR);
}

#endif // FUEL_MOISTURE_TIMELAG_DE_H
