// ============================================================================
// Fire Emissions Model  (WRF-Fire / QUIC-Fire style)
// ============================================================================
// Computes per-cell cumulative fire emissions [kg/m²] for three species:
//
//   Component 0 – co2_emissions   [kg CO₂ / m²]
//   Component 1 – co_emissions    [kg CO  / m²]
//   Component 2 – pm25_emissions  [kg PM₂.₅ / m²]
//
// Formula (Seiler & Crutzen 1980; WRF-Fire / QUIC-Fire convention):
//   E_species [kg/m²] = W_fuel [kg/m²] × f_c [-] × EF_species [kg/kg]
//
// where
//   W_fuel  = oven-dry fuel load [kg/m²] (from Rothermel parameters)
//   f_c     = fuel consumption fraction [-]:
//               – from fuel_consumption_mf when bulk model is active (> 0)
//               – default_consumed_frac for burned cells (phi < 0) otherwise
//               – 0 for unburned cells
//   EF_*    = emission factors [kg species / kg dry fuel consumed]:
//               CO₂:   1.570  (Seiler & Crutzen 1980 / WRF-Fire default)
//               CO:    0.102  (WRF-Fire default)
//               PM₂.₅: 0.0162 (WRF-Fire default, from Andreae & Merlet 2001)
//
// The emissions MultiFab represents the total accumulated emissions per unit
// ground area at the current simulation time (kg/m²), not an instantaneous
// rate.  Post-processing tools can divide by cell area (m²) to obtain total
// mass emitted per cell in kg.
//
// Default emission factors are from:
//   Seiler, W. & Crutzen, P.J. (1980). "Estimates of gross and net fluxes of
//     carbon between the biosphere and the atmosphere from biomass burning."
//     Climatic Change, 2(3), 207–247.
//   Andreae, M.O. & Merlet, P. (2001). "Emission of trace gases and aerosols
//     from biomass burning." Global Biogeochemical Cycles, 15(4), 955–966.
//
// These are the same defaults used in:
//   WRF-SFIRE (Mandel et al. 2011, Geosci. Model Dev.)
//   QUIC-Fire (Linn et al. 2020, Fire, 3(1), 8)
// ============================================================================

#ifndef FIRE_EMISSIONS_MODEL_H
#define FIRE_EMISSIONS_MODEL_H

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_GpuLaunch.H>
#include "parse_inputs.H"

using namespace amrex;

// Number of components in the emissions MultiFab
constexpr int EMISSIONS_NCOMP = 3;

// ============================================================================
// compute_fire_emissions
//
// Parameters:
//   emissions_mf          – output: EMISSIONS_NCOMP-component MultiFab
//   phi                   – input:  level-set field (phi < 0 = burned)
//   fuel_consumption_mf   – input:  fuel consumption fraction [0,1]
//                           (non-zero only when bulk consumption model active)
//   rp                    – Rothermel parameters (provides fuel load)
//   ep                    – emissions parameters (EF values, default fraction)
// ============================================================================
static void compute_fire_emissions(
    MultiFab&                                 emissions_mf,
    const MultiFab&                           phi,
    const MultiFab&                           fuel_consumption_mf,
    const RothermelParams&                    rp,
    const InputParameters::EmissionsParams&   ep)
{
    // Fuel load [kg/m²]
    // Use per-class total when multi-class fuel is active; fall back to w0.
    constexpr Real LB_FT2_TO_KG_M2 = Real(4.8824);
    const Real total_per_class = rp.w_d1 + rp.w_d10 + rp.w_d100
                                + rp.w_lh + rp.w_lw;
    const Real w_eff_lb_ft2 = (total_per_class > Real(0.0))
                              ? total_per_class : rp.w0;
    const Real W_fuel = w_eff_lb_ft2 * LB_FT2_TO_KG_M2;  // kg/m²

    // Emission factors [kg species / kg fuel]
    const Real EF_CO2  = Real(ep.EF_CO2);
    const Real EF_CO   = Real(ep.EF_CO);
    const Real EF_PM25 = Real(ep.EF_PM25);
    const Real def_fc  = Real(ep.default_consumed_frac);
    const Real Wf      = W_fuel;

    for (MFIter mfi(emissions_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const p   = phi.const_array(mfi);
        auto const fc  = fuel_consumption_mf.const_array(mfi);
        auto       em  = emissions_mf.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            // Determine fuel consumption fraction for this cell
            const Real fc_val = fc(i, j, k);
            Real f_c;
            if (fc_val > Real(0.0)) {
                // Bulk consumption model result available
                f_c = fc_val;
            } else if (p(i, j, k) < Real(0.0)) {
                // Burned cell without bulk model → use default fraction
                f_c = def_fc;
            } else {
                // Unburned cell
                f_c = Real(0.0);
            }

            const Real burned_load = Wf * f_c;  // kg/m²

            em(i, j, k, 0) = burned_load * EF_CO2;
            em(i, j, k, 1) = burned_load * EF_CO;
            em(i, j, k, 2) = burned_load * EF_PM25;
        });
    }
}

#endif // FIRE_EMISSIONS_MODEL_H
