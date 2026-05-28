// ============================================================================
// Plume Entrainment Momentum Feedback Model
// ============================================================================
// Enhances the existing heat_flux_model.H with horizontal momentum effects.
// Fire plumes create inflow winds that can strengthen actual ROS by feeding
// the fire with fresh air.
//
// Physics:
//   - Rising plume entrains air from surroundings
//   - Horizontal inflow accelerates toward the fire
//   - Can increase local wind speed and ROS
//   - Important for fire whirl formation and rapid fire growth
//
// Model:
//   dU_plume [m/s] = k_momentum × w_plume × (I / I_ref)^0.5
//   U_eff_total = U_ambient + dU_plume
//
// where:
//   k_momentum = momentum coupling coefficient (typically 0.1-0.3)
//   w_plume    = upward plume velocity [m/s]
//   I          = fireline intensity [kW/m]
//   I_ref      = reference intensity [kW/m] (typically 500-1000)
//
// The upward plume velocity is estimated from buoyancy:
//   w_plume = sqrt(g × Q × H / (ρ × Cp × T_a))
//
// This extends the existing heat_flux_model.H enable_induced capability.
//
// References:
//   Finney, M.A., et al. (2015). Role of buoyant flame dynamics in wildfire
//     spread. PNAS 112(32):9833-9838.
//   Linn, R.R., et al. (2002). Studying wildfire behavior using FIRETEC.
//     International Journal of Wildland Fire 11(4):233-246.
//   Rothermel, R.C. (1983). How to predict the spread and intensity of forest
//     and range fires. USDA Forest Service GTR INT-143.
// ============================================================================

#ifndef PLUME_MOMENTUM_FEEDBACK_H
#define PLUME_MOMENTUM_FEEDBACK_H

#include "constants.H"
#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_GpuLaunch.H>
#include <cmath>

using namespace amrex;

// Default parameters for plume momentum feedback
namespace PlumeMomentumConst {
    constexpr Real K_MOMENTUM_DEFAULT  = 0.2;     // Momentum coupling coefficient
    constexpr Real I_REF_DEFAULT       = 750.0;   // Reference intensity [kW/m]
    constexpr Real RHO_AIR             = 1.2;     // Air density [kg/m³]
    constexpr Real CP_AIR              = 1005.0;  // Specific heat [J/kg/K]
    constexpr Real T_AMBIENT           = 293.15;  // Ambient temperature [K]
    constexpr Real H_REF               = 10.0;    // Reference height [m]
}

// ============================================================================
// Compute upward plume velocity from heat flux
//
// w_plume = sqrt(g × Q × H / (ρ × Cp × T_a))
//
// Parameters:
//   Q_heat    – heat flux [W/m²]
//   H_ref     – reference height [m]
//   rho_air   – air density [kg/m³]
//   Cp_air    – specific heat [J/kg/K]
//   T_ambient – ambient temperature [K]
//
// Returns: upward plume velocity [m/s]
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real compute_plume_velocity(Real Q_heat,
                           Real H_ref     = PlumeMomentumConst::H_REF,
                           Real rho_air   = PlumeMomentumConst::RHO_AIR,
                           Real Cp_air    = PlumeMomentumConst::CP_AIR,
                           Real T_ambient = PlumeMomentumConst::T_AMBIENT)
{
    if (Q_heat <= Real(0.0)) {
        return Real(0.0);
    }
    
    // Buoyancy velocity scale
    // w ~ sqrt(g × Q × H / (ρ × Cp × T))
    const Real g = WildfireConst::GRAVITY;
    const Real denominator = rho_air * Cp_air * T_ambient;
    
    if (denominator <= Real(0.0)) {
        return Real(0.0);
    }
    
    Real w_plume = std::sqrt(g * Q_heat * H_ref / denominator);
    
    // Cap at reasonable maximum (50 m/s)
    return amrex::min(w_plume, Real(50.0));
}

// ============================================================================
// Compute plume-induced horizontal wind increment
//
// dU_plume = k_momentum × w_plume × (I / I_ref)^0.5
//
// Parameters:
//   I_byram     – fireline intensity [kW/m]
//   w_plume     – upward plume velocity [m/s]
//   k_momentum  – momentum coupling coefficient
//   I_ref       – reference intensity [kW/m]
//
// Returns: horizontal wind increment [m/s]
// ============================================================================
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real compute_plume_wind_increment(Real I_byram,
                                  Real w_plume,
                                  Real k_momentum = PlumeMomentumConst::K_MOMENTUM_DEFAULT,
                                  Real I_ref      = PlumeMomentumConst::I_REF_DEFAULT)
{
    if (I_byram <= Real(0.0) || w_plume <= Real(0.0)) {
        return Real(0.0);
    }
    
    // Intensity scaling: (I / I_ref)^0.5
    Real intensity_scale = std::sqrt(I_byram / I_ref);
    
    // Horizontal wind increment
    Real dU = k_momentum * w_plume * intensity_scale;
    
    // Cap at reasonable maximum (10 m/s)
    return amrex::min(dU, Real(10.0));
}

// ============================================================================
// Compute plume-induced wind field from intensity and heat flux
//
// Parameters:
//   dU_plume_mf       – output: plume-induced wind increment [m/s]
//   fireline_intensity_mf – input: Byram intensity [kW/m]
//   Q_heat_mf         – input: heat flux [W/m²]
//   k_momentum        – momentum coupling coefficient
//   I_ref             – reference intensity [kW/m]
// ============================================================================
inline void compute_plume_wind_field(MultiFab& dU_plume_mf,
                                    const MultiFab& fireline_intensity_mf,
                                    const MultiFab& Q_heat_mf,
                                    Real k_momentum = PlumeMomentumConst::K_MOMENTUM_DEFAULT,
                                    Real I_ref      = PlumeMomentumConst::I_REF_DEFAULT)
{
    for (MFIter mfi(dU_plume_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto dU        = dU_plume_mf.array(mfi);
        auto const I_B = fireline_intensity_mf.const_array(mfi);
        auto const Q   = Q_heat_mf.const_array(mfi);
        const Real k_m = k_momentum;
        const Real I_r = I_ref;

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            // Compute upward plume velocity
            Real w_plume = compute_plume_velocity(Q(i, j, k));
            
            // Compute horizontal wind increment
            dU(i, j, k) = compute_plume_wind_increment(
                I_B(i, j, k), w_plume, k_m, I_r
            );
        });
    }
}

// ============================================================================
// Apply plume momentum feedback to wind field (in-place)
//
// Adds the plume-induced wind increment to the existing wind field.
//
// Parameters:
//   U_mf                  – input/output: wind speed field [m/s]
//   fireline_intensity_mf – input: Byram intensity [kW/m]
//   Q_heat_mf             – input: heat flux [W/m²]
//   k_momentum            – momentum coupling coefficient
//   I_ref                 – reference intensity [kW/m]
// ============================================================================
inline void apply_plume_momentum_feedback(MultiFab& U_mf,
                                          const MultiFab& fireline_intensity_mf,
                                          const MultiFab& Q_heat_mf,
                                          Real k_momentum = PlumeMomentumConst::K_MOMENTUM_DEFAULT,
                                          Real I_ref      = PlumeMomentumConst::I_REF_DEFAULT)
{
    for (MFIter mfi(U_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto U         = U_mf.array(mfi);
        auto const I_B = fireline_intensity_mf.const_array(mfi);
        auto const Q   = Q_heat_mf.const_array(mfi);
        const Real k_m = k_momentum;
        const Real I_r = I_ref;

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real w_plume = compute_plume_velocity(Q(i, j, k));
            Real dU = compute_plume_wind_increment(I_B(i, j, k), w_plume, k_m, I_r);
            U(i, j, k) += dU;
        });
    }
}

// ============================================================================
// Compute directional plume wind feedback (2D vectorial)
//
// The plume-induced wind is directed toward the fire (anti-parallel to
// the phi gradient in unburned regions).
//
// Parameters:
//   U_x_mf, U_y_mf        – input/output: wind components [m/s]
//   fireline_intensity_mf – input: Byram intensity [kW/m]
//   Q_heat_mf             – input: heat flux [W/m²]
//   phi_mf                – input: level set field
//   geom                  – geometry for gradient calculation
//   k_momentum            – momentum coupling coefficient
//   I_ref                 – reference intensity [kW/m]
// ============================================================================
inline void apply_directional_plume_feedback(MultiFab& U_x_mf,
                                             MultiFab& U_y_mf,
                                             const MultiFab& fireline_intensity_mf,
                                             const MultiFab& Q_heat_mf,
                                             const MultiFab& phi_mf,
                                             const amrex::Geometry& geom,
                                             Real k_momentum = PlumeMomentumConst::K_MOMENTUM_DEFAULT,
                                             Real I_ref      = PlumeMomentumConst::I_REF_DEFAULT)
{
    const auto dx = geom.CellSizeArray();
    const Real inv_dx = Real(1.0) / dx[0];
    const Real inv_dy = Real(1.0) / dx[1];

    for (MFIter mfi(U_x_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto U_x       = U_x_mf.array(mfi);
        auto U_y       = U_y_mf.array(mfi);
        auto const I_B = fireline_intensity_mf.const_array(mfi);
        auto const Q   = Q_heat_mf.const_array(mfi);
        auto const phi = phi_mf.const_array(mfi);
        const Real k_m = k_momentum;
        const Real I_r = I_ref;

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            // Only apply in unburned region (phi > 0)
            if (phi(i, j, k) <= Real(0.0)) {
                return;
            }
            
            // Compute plume-induced wind magnitude
            Real w_plume = compute_plume_velocity(Q(i, j, k));
            Real dU_mag = compute_plume_wind_increment(I_B(i, j, k), w_plume, k_m, I_r);
            
            if (dU_mag <= Real(0.0)) {
                return;
            }
            
            // Phi gradient (central differences)
            Real dphi_dx = (phi(i+1, j, k) - phi(i-1, j, k)) * Real(0.5) * inv_dx;
            Real dphi_dy = (phi(i, j+1, k) - phi(i, j-1, k)) * Real(0.5) * inv_dy;
            
            Real grad_mag = std::sqrt(dphi_dx * dphi_dx + dphi_dy * dphi_dy);
            
            if (grad_mag > Real(1e-6)) {
                // Normalize and apply direction (toward fire, opposite to gradient)
                Real nx = -dphi_dx / grad_mag;
                Real ny = -dphi_dy / grad_mag;
                
                U_x(i, j, k) += dU_mag * nx;
                U_y(i, j, k) += dU_mag * ny;
            }
        });
    }
}

#endif // PLUME_MOMENTUM_FEEDBACK_H
