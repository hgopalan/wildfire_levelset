// ============================================================================
// Van Wagner's Crown Fire Initiation Model (1977)
// ============================================================================
// This implements Van Wagner's threshold model for predicting the transition
// from surface fire to crown fire based on fuel and fire intensity conditions.
//
// Reference:
// Van Wagner, C.E. (1977). "Conditions for the start and spread of crown fire."
// Canadian Journal of Forest Research, 7(1), 23-34.
//
// Physical Background:
// --------------------
// Crown fires occur when surface fire intensity is sufficient to raise the
// canopy fuels above their ignition temperature. Van Wagner's model predicts
// this transition using:
// 1. Critical surface fire intensity (I_o)
// 2. Canopy base height (CBH)
// 3. Canopy bulk density (CBD)
// 4. Foliar moisture content (FMC)
//
// The model determines if crown fire will initiate and, if so, what the
// active crown fire spread rate will be.
// ============================================================================

#ifndef CROWN_INITIATION_H
#define CROWN_INITIATION_H

#include <AMReX.H>
#include <cmath>
#include "parse_inputs.H"

// ============================================================================
// Compute critical surface fire intensity for crown fire initiation
// ============================================================================
// Van Wagner (1977) equation for the minimum surface fire intensity needed
// to ignite the canopy base.
//
// Parameters:
//   CBH: Canopy base height [m or ft] - height to bottom of live canopy
//   FMC: Foliar moisture content [%] - moisture of canopy foliage (dry basis)
//
// Returns:
//   I_o: Critical surface fire intensity [kW/m or BTU/ft/s]
//
// Formula:
//   I_o = (0.010 * CBH * (460 + 25.9 * FMC))
//
// Physical Interpretation:
//   - Higher canopy base height → higher critical intensity needed
//   - Higher foliar moisture → higher critical intensity needed
//   - Typical values: I_o = 50-500 kW/m for various forest types
//
// Unit Notes:
//   - If CBH in meters and FMC in %, I_o is in kW/m
//   - If CBH in feet and FMC in %, I_o is in BTU/ft/s
//   - For conversion: 1 kW/m ≈ 88.05 BTU/ft/s
// ============================================================================
static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_critical_surface_intensity(
    amrex::Real CBH,
    amrex::Real FMC)
{
    using namespace amrex;
    
    // Ensure positive values
    CBH = amrex::max(CBH, Real(0.1));  // Minimum 0.1 m or ft
    FMC = amrex::max(FMC, Real(50.0)); // Minimum 50% moisture
    
    // Van Wagner's critical intensity equation
    // I_o [kW/m or BTU/ft/s] = 0.010 * CBH * (460 + 25.9 * FMC)
    Real I_o = Real(0.010) * CBH * (Real(460.0) + Real(25.9) * FMC);
    
    return I_o;
}

// ============================================================================
// Check if crown fire will initiate
// ============================================================================
// Compares surface fire intensity to critical threshold.
//
// Parameters:
//   I_surface: Surface fire intensity [kW/m or BTU/ft/s]
//   I_o: Critical surface fire intensity [kW/m or BTU/ft/s]
//
// Returns:
//   true if crown fire initiates (I_surface >= I_o), false otherwise
//
// Physical Interpretation:
//   - If I_surface < I_o: Surface fire only (passive crown fire possible)
//   - If I_surface >= I_o: Active crown fire initiates
// ============================================================================
static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool will_crown_initiate(
    amrex::Real I_surface,
    amrex::Real I_o)
{
    using namespace amrex;
    
    return I_surface >= I_o;
}

// ============================================================================
// Compute active crown fire rate of spread
// ============================================================================
// Van Wagner (1977) model for active crown fire spread rate based on
// canopy bulk density and foliar moisture.
//
// Parameters:
//   CBD: Canopy bulk density [kg/m³] - mass of canopy fuel per volume
//   FMC: Foliar moisture content [%] - moisture of canopy foliage
//
// Returns:
//   R_crown: Crown fire rate of spread [m/min]
//
// Formula:
//   R_crown = 3.0 / CBD   (for FMC ~ 100%)
//
// This is a simplified version. More complex versions account for wind
// and the full relationship with moisture.
//
// Physical Interpretation:
//   - Higher canopy bulk density → slower crown fire spread
//   - CBD = 0.1 kg/m³: R_crown = 30 m/min (very fast, sparse canopy)
//   - CBD = 0.2 kg/m³: R_crown = 15 m/min (moderate)
//   - CBD = 0.3 kg/m³: R_crown = 10 m/min (slower, dense canopy)
//
// Typical CBD values:
//   - Sparse conifer: 0.05-0.10 kg/m³
//   - Moderate conifer: 0.10-0.20 kg/m³
//   - Dense conifer: 0.20-0.30 kg/m³
// ============================================================================
static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_crown_fire_spread_rate(
    amrex::Real CBD,
    amrex::Real FMC)
{
    using namespace amrex;
    
    // Ensure positive CBD
    CBD = amrex::max(CBD, Real(0.01));  // Minimum 0.01 kg/m³
    
    // Simplified Van Wagner formula for crown fire ROS
    // R_crown [m/min] = 3.0 / CBD
    // This is the base formula; can be modified for moisture effects
    Real R_crown = Real(3.0) / CBD;
    
    // Optional: Add moisture effect (reduction factor)
    // Higher moisture reduces spread rate
    Real moisture_factor = Real(1.0) - (FMC - Real(100.0)) / Real(200.0);
    moisture_factor = amrex::max(Real(0.3), amrex::min(Real(1.0), moisture_factor));
    
    R_crown *= moisture_factor;
    
    return R_crown;
}

// ============================================================================
// Compute crown fraction - fraction of fire spread that is crown fire
// ============================================================================
// When crown fire initiates, computes what fraction of the total fire
// spread is from crown fire vs surface fire.
//
// Parameters:
//   R_crown: Crown fire rate of spread [m/min]
//   R_surface: Surface fire rate of spread [m/min]
//   crown_initiated: Whether crown fire has initiated (boolean)
//
// Returns:
//   crown_fraction: Fraction of spread from crown fire (0.0-1.0)
//
// Physical Interpretation:
//   - crown_fraction = 0.0: Pure surface fire
//   - crown_fraction = 0.5: Equal contribution from surface and crown
//   - crown_fraction = 1.0: Pure crown fire (surface fire negligible)
//
// When crown fire is active, total ROS is often dominated by crown component
// because crown fires spread much faster than surface fires.
// ============================================================================
static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_crown_fraction(
    amrex::Real R_crown,
    amrex::Real R_surface,
    bool crown_initiated)
{
    using namespace amrex;
    
    if (!crown_initiated) {
        return Real(0.0);  // No crown fire
    }
    
    // Ensure positive rates
    R_crown = amrex::max(R_crown, Real(0.0));
    R_surface = amrex::max(R_surface, Real(0.0));
    
    // Total rate is weighted combination
    Real R_total = R_crown + R_surface;
    
    if (R_total < Real(1.0e-10)) {
        return Real(0.0);  // Avoid division by zero
    }
    
    // Crown fraction
    Real crown_fraction = R_crown / R_total;
    
    // Ensure in valid range
    crown_fraction = amrex::max(Real(0.0), amrex::min(Real(1.0), crown_fraction));
    
    return crown_fraction;
}

// ============================================================================
// Compute total fire spread rate (surface + crown)
// ============================================================================
// Combines surface and crown fire spread rates when crown fire is active.
//
// Parameters:
//   R_surface: Surface fire rate of spread [m/min or ft/min]
//   R_crown: Crown fire rate of spread [m/min]
//   crown_initiated: Whether crown fire has initiated
//   crown_fraction_weight: How much to weight crown vs surface (0-1)
//
// Returns:
//   R_total: Total fire rate of spread [m/min or ft/min]
//
// Physical Interpretation:
//   - If crown fire not initiated: R_total = R_surface
//   - If crown fire initiated: R_total = R_surface + crown_fraction_weight * R_crown
//
// The crown_fraction_weight allows tuning how much crown fire dominates
// the total spread rate. Common approaches:
//   - Simple addition: crown_fraction_weight = 1.0
//   - Crown dominates: crown_fraction_weight > 1.0 and use max(R_surface, R_crown)
//   - Transition zone: crown_fraction_weight = 0.5-1.0
// ============================================================================
static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_total_spread_rate(
    amrex::Real R_surface,
    amrex::Real R_crown,
    bool crown_initiated,
    amrex::Real crown_fraction_weight = 1.0)
{
    using namespace amrex;
    
    if (!crown_initiated) {
        return R_surface;  // Surface fire only
    }
    
    // Ensure positive rates
    R_surface = amrex::max(R_surface, Real(0.0));
    R_crown = amrex::max(R_crown, Real(0.0));
    
    // When crown fire is active, use the maximum of surface or crown
    // This represents that crown fire typically dominates
    Real R_total = amrex::max(R_surface, R_crown * crown_fraction_weight);
    
    return R_total;
}

// ============================================================================
// Convert Rothermel reaction intensity to surface fire intensity
// ============================================================================
// Converts Rothermel's reaction intensity (I_R) to fireline intensity (I_B)
// which is needed for Van Wagner's crown initiation model.
//
// Parameters:
//   I_R: Reaction intensity [BTU/ft²/min] from Rothermel
//   R_surface: Surface fire rate of spread [ft/min]
//   w_consumed: Weight of fuel consumed per unit area [lb/ft²]
//   h_heat: Heat content of fuel [BTU/lb]
//
// Returns:
//   I_B: Fireline intensity [BTU/ft/s]
//
// Formula (Byram 1959):
//   I_B = h_heat * w_consumed * R_surface / 60
//   where division by 60 converts from per-minute to per-second
//
// Alternatively, from reaction intensity:
//   I_B ≈ I_R * R_surface * residence_time_factor
//
// Physical Interpretation:
//   - Fireline intensity is the power output per unit length of fire front
//   - Higher values indicate more intense fires capable of igniting crowns
//   - Typical range: 50-10000 kW/m (or BTU/ft/s equivalent)
// ============================================================================
static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_fireline_intensity(
    amrex::Real I_R,
    amrex::Real R_surface,
    amrex::Real w_consumed,
    amrex::Real h_heat)
{
    using namespace amrex;
    
    // Byram's fireline intensity formula
    // I_B [BTU/ft/s] = h_heat * w_consumed * R_surface / 60
    Real I_B = h_heat * w_consumed * R_surface / Real(60.0);
    
    return I_B;
}

// ============================================================================
// Convert fireline intensity from BTU/ft/s to kW/m
// ============================================================================
// Van Wagner's model typically uses kW/m, while Rothermel uses BTU/ft/s.
//
// Conversion factor: 1 BTU/ft/s = 3.4613 kW/m
//
// Parameters:
//   I_B_imperial: Fireline intensity [BTU/ft/s]
//
// Returns:
//   I_B_metric: Fireline intensity [kW/m]
// ============================================================================
static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real convert_intensity_to_metric(amrex::Real I_B_imperial)
{
    using namespace amrex;
    
    // Conversion factor: BTU/ft/s to kW/m
    constexpr Real conversion = Real(3.4613);
    
    return I_B_imperial * conversion;
}

// ============================================================================
// Convert fireline intensity from kW/m to BTU/ft/s
// ============================================================================
// Parameters:
//   I_B_metric: Fireline intensity [kW/m]
//
// Returns:
//   I_B_imperial: Fireline intensity [BTU/ft/s]
// ============================================================================
static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real convert_intensity_to_imperial(amrex::Real I_B_metric)
{
    using namespace amrex;
    
    // Conversion factor: kW/m to BTU/ft/s
    constexpr Real conversion = Real(1.0 / 3.4613);
    
    return I_B_metric * conversion;
}

#endif // CROWN_INITIATION_H
