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

#include "constants.H"
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
// Cruz, Alexander & Wakimoto (2005) active crown fire rate of spread
// ============================================================================
// Device-inline version that can be called from within the FARSITE GPU kernel
// and from the levelset crown-override loop in main.cpp.
//
// Formula: R [m/min] = 11.02 × U₁₀^0.90 × CBD^0.19 × exp(−0.17 × MC₁₀)
// Returns ROS in m/s.
//
// Parameters:
//   wind_speed_ms – local 10-m open wind speed [m/s]
//   CBD_kgm3      – canopy bulk density [kg/m³]
//   MC10_pct      – 10-h dead fuel moisture content [%]
// ============================================================================
static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_crown_fire_spread_rate_cruz(
    amrex::Real wind_speed_ms,
    amrex::Real CBD_kgm3,
    amrex::Real MC10_pct) noexcept
{
    using namespace amrex;

    constexpr Real A_LEAD  = WildfireConst::CRUZ_A_LEAD;   // 11.02
    constexpr Real CBD_EXP = WildfireConst::CRUZ_CBD_EXP;  // 0.19
    constexpr Real MC_COEF = WildfireConst::CRUZ_MC_COEFF; // 0.17
    constexpr Real W_EXP   = WildfireConst::CRUZ_WIND_EXP; // 0.90

    // Clamp inputs to physically valid ranges
    const Real CBD  = amrex::max(CBD_kgm3,     Real(0.01));
    const Real MC10 = amrex::max(MC10_pct,     Real(0.0));
    // Convert m/s → km/h; enforce minimum 1 km/h (model not calibrated at calm)
    const Real U_kmh = amrex::max(wind_speed_ms * Real(3.6), Real(1.0));

    // R [m/min] = 11.02 * U^0.90 * CBD^0.19 * exp(-0.17 * MC10)
    const Real R_mmin = A_LEAD
                       * std::pow(U_kmh, W_EXP)
                       * std::pow(CBD,   CBD_EXP)
                       * std::exp(-MC_COEF * MC10);

    // Convert m/min → m/s
    return amrex::max(R_mmin / Real(60.0), Real(0.0));
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
//   crown_fraction_weight: How much to weight crown vs surface (0-2)
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
//   - Simple addition: crown_fraction_weight = 1.0 (default)
//   - Crown dominates: crown_fraction_weight = 1.2-1.5
//   - Transition zone: crown_fraction_weight = 0.5-1.0
//   - Valid range: 0.0-2.0
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
    amrex::Real /*I_R*/,
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
    
    // Conversion factor: BTU/ft/s → kW/m
    constexpr Real conversion = WildfireConst::BTU_FTS_TO_KW_M_VW;
    
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
    
    // Conversion factor: kW/m → BTU/ft/s
    constexpr Real conversion = Real(1.0) / WildfireConst::BTU_FTS_TO_KW_M_VW;
    
    return I_B_metric * conversion;
}

// ============================================================================
// Rothermel (1991) active crown fire rate of spread
// ============================================================================
// Formula: R_crown = 3.34 x R_surface
//
// Reference:
//   Rothermel, R.C. (1991). "Predicting behavior and size of crown fires in
//   the northern Rocky Mountain conifers." USDA Forest Service INT-RP-438.
//
// This is the FARSITE-compatible formula where the active crown fire ROS is
// a fixed multiple of the surface ROS.  Calibrated for northern Rocky Mountain
// conifer forests (MT, ID, WY, CO) - applicable across the continental US.
//
// Parameters:
//   R_surface_ms: surface fire rate of spread [m/s]
//
// Returns:
//   R_crown: active crown fire rate of spread [m/s]
// ============================================================================
static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_rothermel_1991_crown_ros(amrex::Real R_surface_ms) noexcept
{
    return amrex::max(WildfireConst::ROTHERMEL1991_CROWN_MULT * R_surface_ms,
                      amrex::Real(0.0));
}

// ============================================================================
// Van Wagner (1977) passive crown fire interpolation (crowning fraction blend)
// ============================================================================
// When I_B > I_o, a continuous crowning fraction CF blends between surface and
// active crown ROS:
//
//   CF = min((I_B / I_o)^(2/3), 1.0)
//   R_total = R_surface * (1 - CF) + R_crown_active * CF
//
// Reference:
//   Van Wagner, C.E. (1977). "Conditions for the start and spread of crown
//   fire." Canadian Journal of Forest Research, 7(1), 23-34.
//
// Parameters:
//   R_surface_ms:      surface fire ROS [m/s]
//   R_crown_active_ms: active crown fire ROS [m/s]
//   I_B_kwm:           Byram fireline intensity [kW/m]
//   I_o_kwm:           Van Wagner critical surface intensity [kW/m]
//
// Returns:
//   Blended ROS [m/s]; equals R_surface when I_B < I_o,
//   equals R_crown_active when CF = 1 (active crown).
// ============================================================================
static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_van_wagner_passive_blend(
    amrex::Real R_surface_ms,
    amrex::Real R_crown_active_ms,
    amrex::Real I_B_kwm,
    amrex::Real I_o_kwm) noexcept
{
    using namespace amrex;
    if (I_o_kwm < Real(1.0e-6)) return R_surface_ms;
    if (I_B_kwm < I_o_kwm)      return R_surface_ms;

    const Real ratio = I_B_kwm / I_o_kwm;
    const Real CF = amrex::min(std::pow(ratio, Real(2.0 / 3.0)), Real(1.0));

    return R_surface_ms * (Real(1.0) - CF)
         + amrex::max(R_crown_active_ms, R_surface_ms) * CF;
}

#endif // CROWN_INITIATION_H
