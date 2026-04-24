// ============================================================================
// Rothermel (1972) Fire Spread Model - Common Implementation
// ============================================================================
// This file provides the Rothermel fire spread model calculations that are
// shared between the level-set advection and FARSITE ellipse modules.
// ============================================================================

#ifndef ROTHERMEL_MODEL_H
#define ROTHERMEL_MODEL_H

#include "parse_inputs.H"
#include <AMReX.H>

// Structure to hold computed Rothermel parameters for GPU kernel
struct RothermelComputed {
    amrex::Real R0;              // No-wind, no-slope rate of spread [ft/min]
    amrex::Real beta_ratio_E;    // Packing-ratio part of wind factor
    amrex::Real C;               // Wind factor coefficient C (Eq. 47)
    amrex::Real B;               // Wind factor coefficient B (Eq. 48)
    amrex::Real phi_s;           // Slope factor (Eq. 51)
    amrex::Real beta;            // Packing ratio (for terrain-based phi_s calculation)
    amrex::Real wind_conv;       // Unit conversion for wind
    amrex::Real ros_conv;        // Unit conversion for ROS
};

// ============================================================================
// Compute Rothermel (1972) fire spread parameters
// ============================================================================
inline RothermelComputed compute_rothermel_params(const RothermelParams& rp)
{
    using namespace amrex;
    
    RothermelComputed rc;
    
    // Unpack fuel parameters
    const Real w0     = rp.w0;
    const Real sigma  = rp.sigma;
    const Real delta  = rp.delta;
    const Real M_f    = rp.M_f;
    const Real M_x    = rp.M_x;
    const Real h_heat = rp.h_heat;
    const Real S_T    = rp.S_T;
    const Real S_e    = rp.S_e;
    const Real rho_p  = rp.rho_p;

    // Bulk density and packing ratio
    const Real rho_b   = w0 / delta;                                   // Eq. 24
    const Real beta    = rho_b / rho_p;

    // Optimum packing ratio (Eq. 37)
    const Real beta_op = 3.348 * std::pow(sigma, -0.8189);

    // Optimum reaction velocity (Eq. 36)
    const Real sigma15    = std::pow(sigma, 1.5);
    const Real Gamma_max  = sigma15 / (495.0 + 0.0594 * sigma15);

    // A coefficient (Eq. 38)
    const Real A = 133.0 * std::pow(sigma, -0.7913);

    // Actual reaction velocity (Eq. 38)
    const Real beta_ratio = beta / beta_op;
    const Real Gamma_prime = Gamma_max * std::pow(beta_ratio, A)
                             * std::exp(A * (1.0 - beta_ratio));

    // Net fuel load (Eq. 24)
    const Real w_n = w0 * (1.0 - S_T);

    // Moisture damping (Eq. 29)
    const Real rm    = std::min(M_f / M_x, 1.0);
    const Real eta_M = std::max(1.0 - 2.59*rm + 5.11*rm*rm - 3.52*rm*rm*rm, 0.0);

    // Mineral damping (Eq. 30)
    const Real eta_s = 0.174 * std::pow(S_e, -0.19);

    // Reaction intensity [BTU/ft²/min] (Eq. 27)
    const Real I_R = Gamma_prime * w_n * h_heat * eta_M * eta_s;

    // Propagating flux ratio (Eq. 42)
    const Real xi = std::exp((0.792 + 0.681 * std::sqrt(sigma)) * (beta + 0.1))
                    / (192.0 + 0.2595 * sigma);

    // Effective heating number (Eq. 14)
    const Real eps_h = std::exp(-138.0 / sigma);

    // Heat of preignition [BTU/lb] (Eq. 12)
    const Real Q_ig = 250.0 + 1116.0 * M_f;

    // No-wind, no-slope rate of spread [ft/min] (Eq. 1)
    rc.R0 = I_R * xi / (rho_b * eps_h * Q_ig);

    // Wind factor coefficients (Eqs. 47–48)
    rc.C = 7.47  * std::exp(-0.133  * std::pow(sigma, 0.55));
    rc.B = 0.02526 * std::pow(sigma, 0.54);
    const Real E = 0.715  * std::exp(-3.59e-4 * sigma);

    // Packing-ratio part of wind factor (captured for GPU)
    rc.beta_ratio_E = std::pow(beta_ratio, -E);

    // Slope factor (Eq. 51)
    const Real tan_phi_sq = rp.slope_x * rp.slope_x + rp.slope_y * rp.slope_y;
    rc.phi_s = std::max(5.275 * std::pow(beta, -0.3) * tan_phi_sq, 0.0);
    
    // Store beta for terrain-based phi_s calculation
    rc.beta = beta;

    // Unit-conversion scalars for the GPU kernel
    rc.wind_conv = rp.wind_conv;
    rc.ros_conv  = rp.ros_conv;
    
    return rc;
}

#endif // ROTHERMEL_MODEL_H
