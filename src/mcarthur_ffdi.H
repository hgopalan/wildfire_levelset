// ============================================================================
// McArthur Forest Fire Danger Index (FFDI) - Mark 5 Meter
// ============================================================================
// Computes the McArthur Mark 5 Forest Fire Danger Index, the Australian
// standard for assessing forest fire danger and potential fire behavior.
//
// The FFDI is an empirical index calibrated for Australian eucalypt forests
// that combines drought factor, temperature, humidity, and wind speed into
// a single fire danger rating.
//
// FFDI formula (McArthur 1967, Noble et al. 1980):
//   FFDI = 2 × exp(-0.45 + 0.987×ln(DF) - 0.0345×H + 0.0338×T + 0.0234×V)
//
// where:
//   DF - Drought Factor (0–10, typically 1–10)
//   H  - Relative humidity [%]
//   T  - Air temperature [°C]
//   V  - Wind speed [km/h]
//
// FFDI interpretation (McArthur 1967):
//   0–5:    Low–Moderate fire danger
//   5–12:   High fire danger
//   12–24:  Very High fire danger
//   24–50:  Severe fire danger
//   50–100: Extreme fire danger
//   100+:   Catastrophic (Code Red in Victoria)
//
// The Drought Factor (DF) is typically computed from the Keetch-Byram Drought
// Index (KBDI), days since rain, or soil moisture deficit. For simplicity,
// this implementation accepts DF as a user-specified parameter (default: 5.0,
// representing moderate drought conditions).
//
// References:
//   McArthur, A.G. (1967). "Fire behaviour in eucalypt forests."
//     Forestry & Timber Bureau, Canberra, Australia. Leaflet No. 107.
//   Noble, I.R., Bary, G.A.V., & Gill, A.M. (1980). "McArthur's fire-danger
//     meters expressed as equations." Australian Journal of Ecology, 5(2):201–203.
//   Luke, R.H. & McArthur, A.G. (1978). "Bushfires in Australia."
//     Australian Government Publishing Service, Canberra.
// ============================================================================

#ifndef MCARTHUR_FFDI_H
#define MCARTHUR_FFDI_H

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_GpuLaunch.H>

using namespace amrex;

// ============================================================================
// Compute McArthur Mark 5 Forest Fire Danger Index
//
// Parameters:
//   ffdi_mf       – output: McArthur FFDI (dimensionless index)
//   temperature_C – air temperature [°C]
//   humidity_pct  – relative humidity [%]
//   wind_kmh      – wind speed [km/h]
//   drought_factor– Drought Factor (0–10, default: 5.0)
//
// The function accepts either spatially-varying or scalar inputs.
// ============================================================================

// Version 1: Spatially-varying temperature, humidity, and wind speed
inline void compute_mcarthur_ffdi(
    MultiFab& ffdi_mf,
    const MultiFab& temperature_mf,  // [°C]
    const MultiFab& humidity_mf,     // [%]
    const MultiFab& vel,             // [m/s] (3 components)
    Real drought_factor = Real(5.0))
{
    // Clamp drought factor to valid range [0, 10]
    const Real DF = amrex::max(Real(0.1), amrex::min(drought_factor, Real(10.0)));

    // Precompute ln(DF)
    const Real ln_DF = std::log(DF);

    // Convert m/s to km/h: 1 m/s = 3.6 km/h
    const Real ms_to_kmh = Real(3.6);

    for (MFIter mfi(ffdi_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const T_C = temperature_mf.const_array(mfi);
        auto const RH  = humidity_mf.const_array(mfi);
        auto const u   = vel.const_array(mfi);
        auto ffdi      = ffdi_mf.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            // Temperature [°C]
            const Real T = T_C(i, j, k);

            // Relative humidity [%]
            // Clamp to valid range [1, 100] to avoid log(0) or negative values
            const Real H = amrex::max(Real(1.0), amrex::min(RH(i, j, k), Real(100.0)));

            // Wind speed [km/h]
            const Real u_x = u(i, j, k, 0);
            const Real u_y = u(i, j, k, 1);
            const Real U_ms = std::sqrt(u_x * u_x + u_y * u_y);
            const Real V = amrex::max(U_ms * ms_to_kmh, Real(0.1));  // minimum 0.1 km/h

            // McArthur FFDI formula (Noble et al. 1980):
            //   FFDI = 2 × exp(-0.45 + 0.987×ln(DF) - 0.0345×H + 0.0338×T + 0.0234×V)
            const Real exponent = -Real(0.45)
                                + Real(0.987)  * ln_DF
                                - Real(0.0345) * H
                                + Real(0.0338) * T
                                + Real(0.0234) * V;

            ffdi(i, j, k) = Real(2.0) * std::exp(exponent);
        });
    }
}

// Version 2: Uniform temperature, humidity, and wind speed (scalars)
inline void compute_mcarthur_ffdi_uniform(
    MultiFab& ffdi_mf,
    Real temperature_C,     // [°C]
    Real humidity_pct,      // [%]
    Real wind_ms,           // [m/s]
    Real drought_factor = Real(5.0))
{
    // Clamp drought factor to valid range [0.1, 10]
    const Real DF = amrex::max(Real(0.1), amrex::min(drought_factor, Real(10.0)));

    // Clamp humidity to valid range [1, 100]
    const Real H = amrex::max(Real(1.0), amrex::min(humidity_pct, Real(100.0)));

    // Convert wind speed from m/s to km/h
    const Real V = amrex::max(wind_ms * Real(3.6), Real(0.1));

    // Precompute FFDI (constant across domain)
    const Real ln_DF = std::log(DF);
    const Real exponent = -Real(0.45)
                        + Real(0.987)  * ln_DF
                        - Real(0.0345) * H
                        + Real(0.0338) * temperature_C
                        + Real(0.0234) * V;
    const Real ffdi_value = Real(2.0) * std::exp(exponent);

    for (MFIter mfi(ffdi_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto ffdi = ffdi_mf.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            ffdi(i, j, k) = ffdi_value;
        });
    }
}

#endif // MCARTHUR_FFDI_H
