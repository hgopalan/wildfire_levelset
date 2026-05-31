// ============================================================================
// Haines Index (Lower Atmospheric Severity Index)
// ============================================================================
// Implements the Haines Index, an atmospheric stability index that indicates
// the potential for large plume-dominated wildfires and erratic fire behavior.
//
// The Haines Index combines:
//   1. Atmospheric stability (temperature difference between two levels)
//   2. Moisture content (temperature-dewpoint spread at low level)
//
// Three variants exist for different elevations:
//   - Low elevation (< 3000 ft):   uses 950 mb and 850 mb levels
//   - Mid elevation (3000-7000 ft): uses 850 mb and 700 mb levels
//   - High elevation (> 7000 ft):  uses 700 mb and 500 mb levels
//
// Haines Index values:
//   2-3: Very low potential for large plume-dominated fire
//   4:   Low potential
//   5:   Moderate potential  
//   6:   High potential for extreme fire behavior and rapid spread
//
// The index is widely used operationally for fire weather forecasting and
// to assess the likelihood of blow-up fires and dangerous fire behavior.
//
// References:
//   Haines, D.A. (1988). "A lower atmosphere severity index for wildland fires."
//     National Weather Digest, 13:23-27.
//   Werth, P.A. & Ochoa, R. (1993). "The evaluation of Idaho wildfire growth
//     using the Haines Index." Weather and Forecasting, 8(2):223-234.
//   Mills, G.A. & McCaw, L. (2010). "Atmospheric stability environments and
//     fire weather in Australia." CAWCR Technical Report No. 20.
// ============================================================================

#ifndef HAINES_INDEX_H
#define HAINES_INDEX_H

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_GpuLaunch.H>
#include <cmath>

using namespace amrex;

// ============================================================================
// Haines Index elevation variants
// ============================================================================
enum class HainesVariant {
    LOW,   // Low elevation (< 3000 ft):   950-850 mb
    MID,   // Mid elevation (3000-7000 ft): 850-700 mb
    HIGH   // High elevation (> 7000 ft):  700-500 mb
};

// ============================================================================
// compute_haines_stability_term
//
// Computes the stability component (A term) of the Haines Index.
// This represents the temperature lapse rate between two pressure levels.
//
// For each variant, different thresholds are used:
//   LOW:  ΔT = T_950 - T_850
//         A = 1 if ΔT < 4°C, 2 if 4-7°C, 3 if ≥ 8°C
//   MID:  ΔT = T_850 - T_700
//         A = 1 if ΔT < 6°C, 2 if 6-10°C, 3 if ≥ 11°C
//   HIGH: ΔT = T_700 - T_500
//         A = 1 if ΔT < 18°C, 2 if 18-21°C, 3 if ≥ 22°C
//
// Parameters:
//   T_lower - temperature at lower pressure level [°C]
//   T_upper - temperature at upper pressure level [°C]
//   variant - elevation variant (LOW, MID, or HIGH)
//
// Returns: stability term A (1, 2, or 3)
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int compute_haines_stability_term(Real T_lower, Real T_upper, HainesVariant variant)
{
    // Temperature difference (lapse): T_lower - T_upper
    Real dT = T_lower - T_upper;
    
    int A = 1;  // Default: stable
    
    if (variant == HainesVariant::LOW) {
        // Low elevation: 950-850 mb
        if (dT >= Real(8.0)) {
            A = 3;  // Very unstable
        } else if (dT >= Real(4.0)) {
            A = 2;  // Moderately unstable
        } else {
            A = 1;  // Stable
        }
    } else if (variant == HainesVariant::MID) {
        // Mid elevation: 850-700 mb
        if (dT >= Real(11.0)) {
            A = 3;  // Very unstable
        } else if (dT >= Real(6.0)) {
            A = 2;  // Moderately unstable
        } else {
            A = 1;  // Stable
        }
    } else { // HIGH
        // High elevation: 700-500 mb
        if (dT >= Real(22.0)) {
            A = 3;  // Very unstable
        } else if (dT >= Real(18.0)) {
            A = 2;  // Moderately unstable
        } else {
            A = 1;  // Stable
        }
    }
    
    return A;
}

// ============================================================================
// compute_haines_moisture_term
//
// Computes the moisture component (B term) of the Haines Index.
// This represents the dryness of the lower atmosphere via dewpoint depression.
//
// For each variant, different thresholds and pressure levels are used:
//   LOW:  T-Td at 850 mb
//         B = 1 if < 6°C, 2 if 6-9°C, 3 if ≥ 10°C
//   MID:  T-Td at 850 mb
//         B = 1 if < 6°C, 2 if 6-12°C, 3 if ≥ 13°C
//   HIGH: T-Td at 700 mb
//         B = 1 if < 15°C, 2 if 15-20°C, 3 if ≥ 21°C
//
// Parameters:
//   T  - temperature at measurement level [°C]
//   Td - dewpoint temperature at measurement level [°C]
//   variant - elevation variant (LOW, MID, or HIGH)
//
// Returns: moisture term B (1, 2, or 3)
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int compute_haines_moisture_term(Real T, Real Td, HainesVariant variant)
{
    // Dewpoint depression (dry bulb - dewpoint)
    Real dTd = T - Td;
    
    int B = 1;  // Default: moist
    
    if (variant == HainesVariant::LOW) {
        // Low elevation: T-Td at 850 mb
        if (dTd >= Real(10.0)) {
            B = 3;  // Very dry
        } else if (dTd >= Real(6.0)) {
            B = 2;  // Moderately dry
        } else {
            B = 1;  // Moist
        }
    } else if (variant == HainesVariant::MID) {
        // Mid elevation: T-Td at 850 mb
        if (dTd >= Real(13.0)) {
            B = 3;  // Very dry
        } else if (dTd >= Real(6.0)) {
            B = 2;  // Moderately dry
        } else {
            B = 1;  // Moist
        }
    } else { // HIGH
        // High elevation: T-Td at 700 mb
        if (dTd >= Real(21.0)) {
            B = 3;  // Very dry
        } else if (dTd >= Real(15.0)) {
            B = 2;  // Moderately dry
        } else {
            B = 1;  // Moist
        }
    }
    
    return B;
}

// ============================================================================
// compute_haines_index
//
// Computes the Haines Index from atmospheric temperature and dewpoint profiles.
//
// The Haines Index is the sum of stability and moisture terms:
//   HI = A + B
// where A ∈ {1,2,3} and B ∈ {1,2,3}, giving HI ∈ {2,3,4,5,6}
//
// Parameters (variant-specific):
//   For LOW variant (< 3000 ft):
//     T_950  - temperature at 950 mb [°C]
//     T_850  - temperature at 850 mb [°C]
//     Td_850 - dewpoint at 850 mb [°C]
//   For MID variant (3000-7000 ft):
//     T_850  - temperature at 850 mb [°C]
//     T_700  - temperature at 700 mb [°C]
//     Td_850 - dewpoint at 850 mb [°C]
//   For HIGH variant (> 7000 ft):
//     T_700  - temperature at 700 mb [°C]
//     T_500  - temperature at 500 mb [°C]
//     Td_700 - dewpoint at 700 mb [°C]
//
// Returns: Haines Index (2-6)
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int compute_haines_index(Real T_lower, Real T_upper, Real Td_lower, HainesVariant variant)
{
    int A = compute_haines_stability_term(T_lower, T_upper, variant);
    int B = compute_haines_moisture_term(T_lower, Td_lower, variant);
    
    int HI = A + B;
    
    // Haines Index ranges from 2 to 6
    return HI;
}

// ============================================================================
// Simplified Haines Index using surface observations
//
// For cases where upper-air data is not available, approximate the Haines
// Index using surface temperature and humidity with empirical adjustments.
//
// This uses typical atmospheric lapse rates and humidity profiles to estimate
// conditions aloft from surface observations.
//
// Parameters:
//   T_surface  - surface temperature [°C]
//   RH_surface - surface relative humidity [%]
//   elevation_m - terrain elevation [m above sea level]
//
// Returns: approximate Haines Index (2-6)
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int compute_haines_index_surface(Real T_surface, Real RH_surface, Real elevation_m = Real(0.0))
{
    // Determine variant based on elevation
    HainesVariant variant;
    if (elevation_m < Real(914.4)) {  // < 3000 ft
        variant = HainesVariant::LOW;
    } else if (elevation_m < Real(2133.6)) {  // 3000-7000 ft
        variant = HainesVariant::MID;
    } else {  // > 7000 ft
        variant = HainesVariant::HIGH;
    }
    
    // Estimate upper-level temperatures using standard atmospheric lapse
    // Standard lapse rate: ~6.5°C per 1000m in troposphere
    Real lapse_rate = Real(-6.5);  // °C per 1000m
    
    // Estimate temperatures at standard pressure levels
    // Using approximate height differences for pressure levels
    Real T_lower, T_upper, Td_lower;
    
    if (variant == HainesVariant::LOW) {
        // 950 mb ~ 500m, 850 mb ~ 1500m
        T_lower = T_surface + lapse_rate * (Real(0.5) - elevation_m / Real(1000.0));
        T_upper = T_surface + lapse_rate * (Real(1.5) - elevation_m / Real(1000.0));
        // Dewpoint depression increases with height (typically ~2°C per 1000m)
        Real Td_surface = T_surface - (Real(100.0) - RH_surface) / Real(5.0);
        Td_lower = Td_surface - Real(1.0);  // Approximate dewpoint at 850 mb
    } else if (variant == HainesVariant::MID) {
        // 850 mb ~ 1500m, 700 mb ~ 3000m
        T_lower = T_surface + lapse_rate * (Real(1.5) - elevation_m / Real(1000.0));
        T_upper = T_surface + lapse_rate * (Real(3.0) - elevation_m / Real(1000.0));
        Real Td_surface = T_surface - (Real(100.0) - RH_surface) / Real(5.0);
        Td_lower = Td_surface - Real(1.5);
    } else { // HIGH
        // 700 mb ~ 3000m, 500 mb ~ 5500m
        T_lower = T_surface + lapse_rate * (Real(3.0) - elevation_m / Real(1000.0));
        T_upper = T_surface + lapse_rate * (Real(5.5) - elevation_m / Real(1000.0));
        Real Td_surface = T_surface - (Real(100.0) - RH_surface) / Real(5.0);
        Td_lower = Td_surface - Real(3.0);
    }
    
    return compute_haines_index(T_lower, T_upper, Td_lower, variant);
}

// ============================================================================
// Spatially-varying Haines Index computation over a MultiFab
//
// Computes Haines Index field from temperature and dewpoint profiles.
//
// Parameters:
//   haines_mf    – output: Haines Index field (2-6)
//   T_lower_mf   – input: temperature at lower pressure level [°C]
//   T_upper_mf   – input: temperature at upper pressure level [°C]
//   Td_lower_mf  – input: dewpoint at lower pressure level [°C]
//   variant      – elevation variant (LOW, MID, or HIGH)
// ============================================================================
inline void compute_haines_field(
    MultiFab& haines_mf,
    const MultiFab& T_lower_mf,
    const MultiFab& T_upper_mf,
    const MultiFab& Td_lower_mf,
    HainesVariant variant)
{
    for (MFIter mfi(haines_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto haines = haines_mf.array(mfi);
        auto const T_low = T_lower_mf.const_array(mfi);
        auto const T_up = T_upper_mf.const_array(mfi);
        auto const Td_low = Td_lower_mf.const_array(mfi);
        
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real T_lower = T_low(i, j, k);
            Real T_upper = T_up(i, j, k);
            Real Td_lower = Td_low(i, j, k);
            
            int HI = compute_haines_index(T_lower, T_upper, Td_lower, variant);
            
            haines(i, j, k) = static_cast<Real>(HI);
        });
    }
}

// ============================================================================
// Simplified surface-based Haines Index field computation
//
// Parameters:
//   haines_mf      – output: Haines Index field (2-6)
//   temperature_mf – input: surface temperature [°C]
//   humidity_mf    – input: surface relative humidity [%]
//   elevation_mf   – input: terrain elevation [m] (optional, default 0)
// ============================================================================
inline void compute_haines_field_surface(
    MultiFab& haines_mf,
    const MultiFab& temperature_mf,
    const MultiFab& humidity_mf,
    const MultiFab* elevation_mf = nullptr)
{
    for (MFIter mfi(haines_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto haines = haines_mf.array(mfi);
        auto const T_surf = temperature_mf.const_array(mfi);
        auto const RH_surf = humidity_mf.const_array(mfi);
        
        if (elevation_mf != nullptr) {
            auto const elev = elevation_mf->const_array(mfi);
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                int HI = compute_haines_index_surface(T_surf(i, j, k), RH_surf(i, j, k), elev(i, j, k));
                haines(i, j, k) = static_cast<Real>(HI);
            });
        } else {
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                int HI = compute_haines_index_surface(T_surf(i, j, k), RH_surf(i, j, k));
                haines(i, j, k) = static_cast<Real>(HI);
            });
        }
    }
}

#endif // HAINES_INDEX_H
