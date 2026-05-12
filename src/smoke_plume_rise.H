// ============================================================================
// Smoke Plume-Rise Model  (Briggs 1965 / 1969)
// ============================================================================
//
// Computes the effective smoke plume-rise height Δh [m] above the fire front
// using the Briggs (1965/1969) buoyancy-dominated plume formula.  The result
// is written per-cell into the `plume_rise_mf` MultiFab and optionally printed
// as a domain maximum in the simulation log.
//
// Theory (Briggs 1969, eq. 3-10 / Slade 1968):
//
//   1. Buoyancy flux per unit fire-line length [m⁴/s³]:
//        F_B = g × I_B / (π × ρ_a × c_p × T_a)
//      where:
//        I_B   = Byram fireline intensity [W/m] (= 1000 × kW/m from fireline_intensity_mf)
//        g     = 9.81 m/s²  (gravitational acceleration)
//        ρ_a   = air density [kg/m³]  (default 1.20)
//        c_p   = specific heat of air [J/(kg·K)]  (default 1005.0)
//        T_a   = ambient temperature [K]  (default 303.15 = 30 °C)
//
//   2. Distance to final rise x_f = 3.5 × x* [m]:
//        x* = 49  × F_B^(5/8)   for F_B < 55 m⁴/s³
//        x* = 120 × F_B^(2/5)   for F_B ≥ 55 m⁴/s³
//
//   3. Final plume rise:
//        Δh = 1.6 × F_B^(1/3) × x_f^(2/3) / u   [m]
//      where u = wind speed [m/s] (clamped to u_min = 0.5 m/s to avoid division by zero).
//
// Only cells with positive fireline intensity (active fire front) receive a
// non-zero plume rise.  All others are set to 0.
//
// References:
//   Briggs, G.A. (1965). A plume rise model compared with observations.
//     JAPCA 15(9):433–438.
//   Briggs, G.A. (1969). Plume rise. USAEC Critical Review Series, TID-25075.
//   Slade, D.H. (ed.) (1968). Meteorology and Atomic Energy.
//     US Atomic Energy Commission, Division of Technical Information.
//
// Input prefix: "smoke_plume."
// ============================================================================

#ifndef SMOKE_PLUME_RISE_H
#define SMOKE_PLUME_RISE_H

#include "parse_inputs.H"
#include "constants.H"

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_GpuLaunch.H>

using namespace amrex;

// ============================================================================
// compute_smoke_plume_rise
//
// Parameters:
//   plume_rise_mf         – output:  single-component MultiFab [m]
//   fireline_intensity_mf – input:   Byram fireline intensity [kW/m]
//   vel                   – input:   wind velocity (used for |u|)
//   sp                    – smoke plume parameters (enable, T_a, rho_a, Cp_a)
// ============================================================================
inline void compute_smoke_plume_rise(
    MultiFab&                                      plume_rise_mf,
    const MultiFab&                                fireline_intensity_mf,
    const MultiFab&                                vel,
    const InputParameters::SmokePlumeParams&       sp)
{
    if (sp.enable == 0) {
        plume_rise_mf.setVal(Real(0.0));
        return;
    }

    // Physical constants
    constexpr Real g = Real(9.81);          // gravitational acceleration [m/s²]
    // Briggs transition threshold and coefficients
    constexpr Real FB_THRESHOLD = Real(55.0); // [m⁴/s³]
    constexpr Real COEFF_SMALL  = Real(49.0); // x* coefficient for F_B < 55
    constexpr Real COEFF_LARGE  = Real(120.0);// x* coefficient for F_B ≥ 55
    constexpr Real EXP_SMALL    = Real(5.0 / 8.0);
    constexpr Real EXP_LARGE    = Real(2.0 / 5.0);
    constexpr Real C_FINAL      = Real(1.6);  // Briggs final-rise coefficient
    constexpr Real F_EXP        = Real(1.0 / 3.0);
    constexpr Real XF_EXP       = Real(2.0 / 3.0);
    constexpr Real U_MIN        = Real(0.5);  // minimum wind speed [m/s]

    // Parameters from struct (trivially copyable scalars only)
    const Real T_a   = Real(sp.T_a);    // ambient temperature [K]
    const Real rho_a = Real(sp.rho_a);  // air density [kg/m³]
    const Real Cp_a  = Real(sp.Cp_a);   // specific heat of air [J/(kg·K)]

    // π × ρ_a × c_p × T_a  (denominator of F_B formula)
    const Real denom = Real(WildfireConst::PI) * rho_a * Cp_a * T_a;

    for (MFIter mfi(plume_rise_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto       pr = plume_rise_mf.array(mfi);
        auto const fi = fireline_intensity_mf.const_array(mfi);
        auto const v  = vel.const_array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const Real I_B_kwm = fi(i, j, k);   // [kW/m]
            if (I_B_kwm <= Real(0.0)) {
                pr(i, j, k) = Real(0.0);
                return;
            }

            // Fireline intensity in W/m
            const Real I_B_wm = I_B_kwm * Real(1000.0);

            // Buoyancy flux per unit fire-line length [m⁴/s³]
            const Real F_B = g * I_B_wm / denom;

            // Distance to final rise x* [m]
            Real x_star;
            if (F_B < FB_THRESHOLD) {
                x_star = COEFF_SMALL * std::pow(F_B, EXP_SMALL);
            } else {
                x_star = COEFF_LARGE * std::pow(F_B, EXP_LARGE);
            }
            // x_f = 3.5 × x*
            const Real x_f = Real(3.5) * x_star;

            // Wind speed (2-D: horizontal components only)
            const Real ux = v(i, j, k, 0);
            const Real uy = v(i, j, k, 1);
            const Real u  = amrex::max(std::sqrt(ux * ux + uy * uy), U_MIN);

            // Briggs final plume rise [m]
            const Real dh = C_FINAL * std::pow(F_B, F_EXP)
                            * std::pow(x_f, XF_EXP)
                            / u;

            pr(i, j, k) = amrex::max(Real(0.0), dh);
        });
    }
}

#endif // SMOKE_PLUME_RISE_H
