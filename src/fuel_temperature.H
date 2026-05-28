// ============================================================================
// Fuel Particle Temperature Evolution Model
// ============================================================================
// Tracks the surface temperature of fuel particles ahead of and behind the
// fire front. This determines actual ignition timing (not instantaneous),
// affects spotting ignition probability, and is important for firebrand
// ignition delay.
//
// Physics:
//   Fuel particles are heated by:
//     - Radiative heat flux from flames (Q_rad)
//     - Convective heat flux from hot gases (Q_conv)
//   And cooled by:
//     - Re-radiation to environment (Q_loss)
//
// Energy balance (simple 0-D lumped capacitance model):
//   dT/dt = (Q_rad + Q_conv - Q_loss) / (ρ_p × c_p × V_p)
//
// Simplified form:
//   dT/dt = (Q_net) / C_thermal
//
// where:
//   T         = fuel particle temperature [K]
//   Q_rad     = radiative heat flux from fire [W/m²]
//   Q_conv    = convective heat flux from hot gases [W/m²]
//   Q_loss    = re-radiation loss ≈ ε×σ×(T⁴ - T_a⁴) [W/m²]
//   ρ_p       = particle density [kg/m³]
//   c_p       = specific heat capacity [J/kg/K]
//   V_p       = particle volume per unit area [m³/m²] ≈ fuel depth / packing
//   C_thermal = thermal capacitance [J/K/m²]
//
// Ignition criterion:
//   Ignition occurs when T >= T_ignition (typically ~600 K for wood)
//
// References:
//   Dupuy, J.L., & Maréchal, J. (2011). Slope effect on laboratory fire spread.
//     International Journal of Wildland Fire 20(7):884-892.
//   Pickett, B.M., et al. (2010). Experimental measurements during combustion
//     of moist individual fuel particles. International Journal of Wildland
//     Fire 19(2):153-162.
//   Drysdale, D. (2011). An Introduction to Fire Dynamics, 3rd ed. Wiley.
// ============================================================================

#ifndef FUEL_TEMPERATURE_H
#define FUEL_TEMPERATURE_H

#include "constants.H"
#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_GpuLaunch.H>
#include <cmath>

using namespace amrex;

// Physical constants for fuel temperature model
namespace FuelTempConst {
    constexpr Real T_AMBIENT       = 293.15;  // Ambient temperature [K]
    constexpr Real T_IGNITION      = 600.0;   // Ignition temperature [K] for wood
    constexpr Real RHO_PARTICLE    = 513.0;   // Particle density [kg/m³] (wood)
    constexpr Real C_P_FUEL        = 2500.0;  // Specific heat [J/kg/K]
    constexpr Real EMISSIVITY      = 0.95;    // Surface emissivity
    constexpr Real SIGMA_SB        = 5.67e-8; // Stefan-Boltzmann [W/m²/K⁴]
}

// ============================================================================
// Compute thermal capacitance per unit area
//
// C_thermal [J/K/m²] = ρ_p × c_p × V_char
//
// where V_char is a characteristic fuel volume per unit area [m].
// For a fuel bed: V_char ≈ fuel_depth × (1 - packing_ratio)
//
// Parameters:
//   fuel_depth     – fuel bed depth [m]
//   packing_ratio  – fuel bed packing ratio (dimensionless, ~0.01-0.1)
//   rho_p          – particle density [kg/m³]
//   c_p            – specific heat [J/kg/K]
//
// Returns: thermal capacitance [J/K/m²]
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real compute_thermal_capacitance(Real fuel_depth,
                                 Real packing_ratio,
                                 Real rho_p = FuelTempConst::RHO_PARTICLE,
                                 Real c_p   = FuelTempConst::C_P_FUEL)
{
    // Characteristic volume per unit area [m]
    Real V_char = fuel_depth * amrex::max(Real(1.0) - packing_ratio, Real(0.1));
    
    // Thermal capacitance [J/K/m²]
    return rho_p * c_p * V_char;
}

// ============================================================================
// Compute radiative loss from fuel particle
//
// Q_loss [W/m²] = ε × σ × (T⁴ - T_a⁴)
//
// Parameters:
//   T_fuel      – fuel particle temperature [K]
//   T_ambient   – ambient air temperature [K]
//   emissivity  – surface emissivity (0-1)
//
// Returns: radiative loss [W/m²] (positive = cooling)
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real compute_radiative_loss(Real T_fuel,
                            Real T_ambient = FuelTempConst::T_AMBIENT,
                            Real emissivity = FuelTempConst::EMISSIVITY)
{
    const Real T4_fuel = T_fuel * T_fuel * T_fuel * T_fuel;
    const Real T4_amb  = T_ambient * T_ambient * T_ambient * T_ambient;
    
    return emissivity * FuelTempConst::SIGMA_SB * (T4_fuel - T4_amb);
}

// ============================================================================
// Update fuel temperature field using explicit Euler time integration
//
// Parameters:
//   T_fuel_mf        – input/output: fuel particle temperature [K]
//   Q_rad_mf         – input: radiative heat flux from fire [W/m²]
//   Q_conv_mf        – input: convective heat flux from gases [W/m²]
//   C_thermal_mf     – input: thermal capacitance per unit area [J/K/m²]
//   dt               – time step [s]
//   T_ambient        – ambient temperature [K]
//   T_max            – maximum temperature cap [K] (default: 1500 K)
//
// The temperature is updated as:
//   T^(n+1) = T^n + dt × (Q_rad + Q_conv - Q_loss) / C_thermal
// ============================================================================
inline void update_fuel_temperature(MultiFab& T_fuel_mf,
                                    const MultiFab& Q_rad_mf,
                                    const MultiFab& Q_conv_mf,
                                    const MultiFab& C_thermal_mf,
                                    Real dt,
                                    Real T_ambient = FuelTempConst::T_AMBIENT,
                                    Real T_max     = Real(1500.0))
{
    for (MFIter mfi(T_fuel_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto T_fuel    = T_fuel_mf.array(mfi);
        auto const Q_rad   = Q_rad_mf.const_array(mfi);
        auto const Q_conv  = Q_conv_mf.const_array(mfi);
        auto const C_therm = C_thermal_mf.const_array(mfi);
        const Real dt_local = dt;
        const Real T_amb    = T_ambient;
        const Real T_cap    = T_max;

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const Real T_current = T_fuel(i, j, k);
            const Real C_th      = C_therm(i, j, k);
            
            // Skip if thermal capacitance is too small (numerical stability)
            if (C_th < Real(1.0)) {
                return;
            }

            // Net heat flux [W/m²]
            const Real Q_loss = compute_radiative_loss(T_current, T_amb);
            const Real Q_net  = Q_rad(i, j, k) + Q_conv(i, j, k) - Q_loss;

            // Temperature increment [K]
            const Real dT = dt_local * Q_net / C_th;

            // Update temperature with bounds
            Real T_new = T_current + dT;
            T_new = amrex::max(T_new, T_amb);  // Cannot go below ambient
            T_new = amrex::min(T_new, T_cap);  // Cap at maximum

            T_fuel(i, j, k) = T_new;
        });
    }
}

// ============================================================================
// Check for ignition based on fuel temperature
//
// Sets ignition flag to 1.0 where T_fuel >= T_ignition, 0.0 otherwise.
//
// Parameters:
//   ignition_mf     – output: ignition flag (0 or 1)
//   T_fuel_mf       – input: fuel particle temperature [K]
//   T_ignition      – ignition temperature threshold [K]
// ============================================================================
inline void check_thermal_ignition(MultiFab& ignition_mf,
                                   const MultiFab& T_fuel_mf,
                                   Real T_ignition = FuelTempConst::T_IGNITION)
{
    for (MFIter mfi(ignition_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto ign       = ignition_mf.array(mfi);
        auto const T_f = T_fuel_mf.const_array(mfi);
        const Real T_ig = T_ignition;

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            ign(i, j, k) = (T_f(i, j, k) >= T_ig) ? Real(1.0) : Real(0.0);
        });
    }
}

// ============================================================================
// Initialize fuel temperature field
//
// Sets all cells to ambient temperature.
//
// Parameters:
//   T_fuel_mf     – output: fuel temperature field [K]
//   T_ambient     – ambient temperature [K]
// ============================================================================
inline void init_fuel_temperature(MultiFab& T_fuel_mf,
                                  Real T_ambient = FuelTempConst::T_AMBIENT)
{
    T_fuel_mf.setVal(T_ambient);
}

// ============================================================================
// Compute thermal capacitance field from fuel properties
//
// Parameters:
//   C_thermal_mf    – output: thermal capacitance [J/K/m²]
//   fuel_depth_mf   – input: fuel bed depth [m]
//   packing_mf      – input: packing ratio (dimensionless)
// ============================================================================
inline void compute_thermal_capacitance_field(MultiFab& C_thermal_mf,
                                              const MultiFab& fuel_depth_mf,
                                              const MultiFab& packing_mf)
{
    const Real rho_p = FuelTempConst::RHO_PARTICLE;
    const Real c_p   = FuelTempConst::C_P_FUEL;

    for (MFIter mfi(C_thermal_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto C_therm    = C_thermal_mf.array(mfi);
        auto const depth = fuel_depth_mf.const_array(mfi);
        auto const pack  = packing_mf.const_array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            C_therm(i, j, k) = compute_thermal_capacitance(
                depth(i, j, k), pack(i, j, k), rho_p, c_p
            );
        });
    }
}

#endif // FUEL_TEMPERATURE_H
