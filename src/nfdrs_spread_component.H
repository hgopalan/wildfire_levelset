// ============================================================================
// NFDRS Spread Component (SC)
// ============================================================================
// Implements the Spread Component from the U.S. National Fire Danger Rating
// System (NFDRS). The SC is a numerical rating of the forward rate of spread
// of a fire based on fuel model, fuel moisture, and weather conditions.
//
// The Spread Component is one of the core NFDRS indices used for operational
// fire management planning and resource allocation across the United States.
//
// SC Calculation (NFDRS 1978):
//   SC = f(fuel_model, 1-hr FM, wind_speed, slope)
//
// The SC is calibrated to NFDRS reference fuels (A-U) and ranges from 0
// (no spread potential) to 99+ (extreme spread potential). Values above
// 50 indicate very high fire spread danger.
//
// SC Interpretation:
//     0-10:  Low spread potential
//    10-20:  Moderate spread potential
//    20-40:  High spread potential
//    40-60:  Very high spread potential
//    60+:    Extreme spread potential
//
// Relationship to Other NFDRS Components:
//   - SC feeds into the Burning Index (BI): BI = SC × ERC
//   - Used in fire dispatch planning and initial attack decisions
//   - Combined with other components for overall fire danger rating
//
// NFDRS Fuel Models (simplified):
//   A - Western annual grass
//   B - California chaparral
//   C - Pine-grass savanna
//   D - Southern rough
//   E - Hardwood litter (winter)
//   F - Intermediate brush
//   G - Short needle (heavy dead)
//   H - Short needle (normal dead)
//   I - Heavy slash
//   J - Intermediate slash
//   K - Light slash
//   L - Western perennial grass
//   P - Ponderosa pine
//   Q - Sagebrush grass
//   R - Hardwood litter (normal)
//   S - Tundra
//   T - Sagebrush
//   U - Western pines
//
// References:
//   Deeming, J.E., Burgan, R.E., & Cohen, J.D. (1977). "The National Fire-Danger
//     Rating System—1978." USDA Forest Service General Technical Report INT-39.
//   Bradshaw, L.S., Deeming, J.E., Burgan, R.E., & Cohen, J.D. (1984). "The 1978
//     National Fire-Danger Rating System: Technical documentation."
//     USDA Forest Service General Technical Report INT-169.
//   Andrews, P.L. & Rothermel, R.C. (1982). "Charts for interpreting wildland fire
//     behavior characteristics." USDA Forest Service General Technical Report INT-131.
// ============================================================================

#ifndef NFDRS_SPREAD_COMPONENT_H
#define NFDRS_SPREAD_COMPONENT_H

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_GpuLaunch.H>
#include <cmath>

using namespace amrex;



// ============================================================================
// compute_nfdrs_spread_component
//
// Computes the NFDRS Spread Component using simplified empirical formulation.
//
// Formulation (Deeming et al. 1977, simplified):
//   SC = a × (1 + b × U) × (1 + c × S) × exp(-d × FM / 100)
//
// where:
//   U  - wind speed [mph at 20-ft]
//   S  - slope [%]
//   FM - 1-hour dead fuel moisture [%]
//   a, b, c, d - fuel model coefficients
//
// Parameters:
//   wind_speed_mps  - wind speed [m/s] (converted to mph internally)
//   slope_fraction  - terrain slope (rise/run, dimensionless)
//   fuel_moisture_1hr_pct - 1-hour dead fuel moisture [%]
//   fuel_model_factor - fuel model adjustment factor (default: 1.0)
//
// Returns: Spread Component (0-99+)
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real compute_nfdrs_spread_component(Real wind_speed_mps,
                                     Real slope_fraction,
                                     Real fuel_moisture_1hr_pct,
                                     Real fuel_model_factor = Real(1.0))
{
// Local constexpr constants (safe for device code)
   Real A_COEFF = Real(0.792);
   Real B_COEFF = Real(0.0404);
   Real C_COEFF = Real(2.988);

   Real WIND_MIN_MPH = Real(0.0);
   Real WIND_MAX_MPH = Real(30.0);

   Real FM_MIN_PCT = Real(1.0);
   Real FM_MAX_PCT = Real(50.0);

   Real SLOPE_MIN_PCT = Real(0.0);
   Real SLOPE_MAX_PCT = Real(100.0);

   Real MPS_TO_MPH = Real(2.23694);

    // Convert wind speed to mph
    Real U_mph = wind_speed_mps * MPS_TO_MPH;
    U_mph = amrex::max(WIND_MIN_MPH,
                       amrex::min(U_mph, WIND_MAX_MPH));

    // Convert slope to percent
    Real S_pct = slope_fraction * Real(100.0);
    S_pct = amrex::max(SLOPE_MIN_PCT,
                       amrex::min(S_pct, SLOPE_MAX_PCT));

    // Clamp fuel moisture to valid range
    Real FM = amrex::max(FM_MIN_PCT,
                         amrex::min(fuel_moisture_1hr_pct, FM_MAX_PCT));

    // NFDRS SC formulation (Deeming et al. 1977)
    // SC = base × wind_factor × slope_factor × moisture_factor

    Real base = Real(10.0) * fuel_model_factor;  // Base SC for reference conditions

    // Wind effect (increases SC linearly with wind)
    Real wind_factor = Real(1.0) + B_COEFF * U_mph;

    // Slope effect (increases SC linearly with slope)
    Real slope_factor = Real(1.0) + A_COEFF * S_pct / Real(100.0);

    // Moisture damping (exponential decrease with moisture)
    Real moisture_factor = std::exp(-C_COEFF * FM / Real(100.0));

    // Compute SC
    Real SC = base * wind_factor * slope_factor * moisture_factor;

    // Clamp to reasonable range [0, 99]
    return amrex::max(Real(0.0), amrex::min(SC, Real(99.0)));
}

}

// ============================================================================
// compute_nfdrs_fuel_model_factor
//
// Returns the fuel model adjustment factor for different NFDRS fuel types.
//
// NFDRS fuel models have different base spread characteristics.
//
// Parameters:
//   fuel_model_code - NFDRS fuel model code ('A'-'U')
//
// Returns: fuel model factor (dimensionless multiplier)
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real compute_nfdrs_fuel_model_factor(char fuel_model_code)
{
    // Fuel model spread factors (calibrated to NFDRS reference)
    switch (fuel_model_code) {
        case 'A': return Real(1.5);   // Western annual grass (high spread)
        case 'B': return Real(0.8);   // California chaparral (moderate)
        case 'C': return Real(1.2);   // Pine-grass savanna
        case 'D': return Real(0.9);   // Southern rough
        case 'E': return Real(0.6);   // Hardwood litter (winter)
        case 'F': return Real(0.7);   // Intermediate brush
        case 'G': return Real(1.0);   // Short needle (heavy dead)
        case 'H': return Real(1.1);   // Short needle (normal dead)
        case 'I': return Real(0.5);   // Heavy slash (low spread)
        case 'J': return Real(0.7);   // Intermediate slash
        case 'K': return Real(0.9);   // Light slash
        case 'L': return Real(1.3);   // Western perennial grass
        case 'P': return Real(0.8);   // Ponderosa pine
        case 'Q': return Real(1.1);   // Sagebrush grass
        case 'R': return Real(0.7);   // Hardwood litter (normal)
        case 'S': return Real(0.6);   // Tundra
        case 'T': return Real(0.9);   // Sagebrush
        case 'U': return Real(0.8);   // Western pines
        default:  return Real(1.0);   // Default (model G equivalent)
    }
}

// ============================================================================
// compute_nfdrs_sc_field
//
// Computes NFDRS Spread Component for all cells in a MultiFab.
//
// Parameters:
//   sc_mf           - output Spread Component field (0-99)
//   wind_speed_mf   - wind speed field [m/s]
//   slope_mf        - terrain slope field (rise/run, dimensionless)
//   fuel_moisture_mf - 1-hour dead fuel moisture field [%]
//   fuel_model_factor - fuel model adjustment (default: 1.0 for model G)
// ============================================================================
inline void compute_nfdrs_sc_field(MultiFab& sc_mf,
                                    const MultiFab& wind_speed_mf,
                                    const MultiFab& slope_mf,
                                    const MultiFab& fuel_moisture_mf,
                                    Real fuel_model_factor = Real(1.0))
{
    for (MFIter mfi(sc_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        Array4<Real> const& sc = sc_mf.array(mfi);
        Array4<Real const> const& wind = wind_speed_mf.const_array(mfi);
        Array4<Real const> const& slope = slope_mf.const_array(mfi);
        Array4<Real const> const& fm = fuel_moisture_mf.const_array(mfi);
        
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            sc(i,j,k) = compute_nfdrs_spread_component(wind(i,j,k),
                                                        slope(i,j,k),
                                                        fm(i,j,k),
                                                        fuel_model_factor);
        });
    }
}

// ============================================================================
// compute_nfdrs_burning_index
//
// Computes the NFDRS Burning Index from Spread Component and ERC.
//
// Burning Index (BI) combines spread potential (SC) with fuel energy (ERC):
//   BI = (SC × ERC) / 10
//
// Parameters:
//   SC  - Spread Component (0-99)
//   ERC - Energy Release Component (0-99+)
//
// Returns: Burning Index (0-200+)
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real compute_nfdrs_burning_index(Real SC, Real ERC)
{
    // Burning Index formulation (Deeming et al. 1977)
    Real BI = (SC * ERC) / Real(10.0);
    
    // Clamp to reasonable minimum
    return amrex::max(Real(0.0), BI);
}

#endif // NFDRS_SPREAD_COMPONENT_H
