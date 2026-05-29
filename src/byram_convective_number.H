// ============================================================================
// Byram's Convective Number (Nc)
// ============================================================================
// Computes Byram's (1959) convective number, which classifies fire behavior as
// wind-driven (Nc < 1) vs. plume-dominated (Nc > 1).
//
// Convective number:
//   Nc = (g × I_B × τ_res) / (ρ_air × C_p × T_a × U³)
//
// where:
//   g       - gravitational acceleration [m/s²]
//   I_B     - Byram fireline intensity [kW/m] = [kW·s/(m·s)]
//   τ_res   - flame residence time [s]
//   ρ_air   - air density [kg/m³]
//   C_p     - specific heat of air [J/(kg·K)]
//   T_a     - ambient temperature [K]
//   U       - wind speed [m/s]
//
// When Nc > 1, fire is plume-dominated (convective); when Nc < 1, wind-dominated.
// Large convective numbers indicate strong vertical plume development and
// increased spotting potential.
//
// References:
//   Byram, G.M. (1959). "Combustion of forest fuels." In: Davis, K.P. (ed.)
//   Forest Fire: Control and Use. McGraw-Hill, New York. pp. 61–89.
//
//   Nelson, R.M. & Adkins, C.W. (1986). "Flame characteristics of wind-driven
//   surface fires." Canadian Journal of Forest Research, 16(6):1293–1300.
// ============================================================================

#ifndef BYRAM_CONVECTIVE_NUMBER_H
#define BYRAM_CONVECTIVE_NUMBER_H

#include "constants.H"
#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_GpuLaunch.H>

using namespace amrex;

// ============================================================================
// Compute Byram's convective number from fireline intensity, residence time,
// and wind speed fields
//
// Parameters:
//   convective_number_mf  – output: Byram convective number [-] (dimensionless)
//   fireline_intensity_mf – input:  Byram fireline intensity [kW/m]
//   residence_time_mf     – input:  flame residence time [s]
//   vel                   – input:  wind velocity field (3 components) [m/s]
//   rho_air               – air density [kg/m³] (default: 1.2)
//   C_p                   – specific heat of air [J/(kg·K)] (default: 1005)
//   T_a                   – ambient temperature [K] (default: 300)
// ============================================================================
inline void compute_byram_convective_number(
    MultiFab& convective_number_mf,
    const MultiFab& fireline_intensity_mf,
    const MultiFab& residence_time_mf,
    const MultiFab& vel,
    Real rho_air = Real(1.2),
    Real C_p     = Real(1005.0),
    Real T_a     = Real(300.0))
{
    const Real g = WildfireConst::GRAVITY;  // 9.81 m/s²

    // Convert kW to W for energy balance: 1 kW = 1000 W = 1000 J/s
    const Real kW_to_W = Real(1000.0);

    // Denominator constant: ρ_air × C_p × T_a
    const Real denom_const = rho_air * C_p * T_a;

    for (MFIter mfi(convective_number_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const I_B    = fireline_intensity_mf.const_array(mfi);
        auto const tau    = residence_time_mf.const_array(mfi);
        auto const u      = vel.const_array(mfi);
        auto nc           = convective_number_mf.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            // Fireline intensity [kW/m] → [W/m]
            const Real I_B_W = I_B(i, j, k) * kW_to_W;

            // Residence time [s]
            const Real tau_s = tau(i, j, k);

            // Wind speed magnitude [m/s]
            const Real u_x = u(i, j, k, 0);
            const Real u_y = u(i, j, k, 1);
            const Real U_mag = std::sqrt(u_x * u_x + u_y * u_y);

            // Avoid division by zero or near-zero wind speeds
            // When U → 0, Nc → ∞ (strongly plume-dominated)
            // Set a floor of 0.1 m/s to prevent numerical overflow
            const Real U_safe = amrex::max(U_mag, Real(0.1));
            const Real U_cubed = U_safe * U_safe * U_safe;

            // Byram convective number (dimensionless)
            //   Nc = (g × I_B × τ_res) / (ρ_air × C_p × T_a × U³)
            const Real numerator = g * I_B_W * tau_s;
            const Real denominator = denom_const * U_cubed;

            nc(i, j, k) = (denominator > Real(0.0)) ? numerator / denominator
                                                      : Real(0.0);
        });
    }
}

#endif // BYRAM_CONVECTIVE_NUMBER_H
