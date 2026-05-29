// ============================================================================
// Ember Accumulation and Probabilistic Ignition Model
// ============================================================================
// Tracks the density of landed firebrands at each cell and computes
// probabilistic ignition based on accumulated ember coverage. This extends
// the ember_cascade_flux and Albini spotting models by modeling ember
// persistence and burnout.
//
// Physical Model:
// ───────────────
//  1. Ember density accumulation:
//       ρ_ember(t+dt) = ρ_ember(t) + F_landing × dt - ρ_ember(t) × k_decay × dt
//
//     where:
//       F_landing  – ember landing flux [embers/m²/s] from ember_cascade or Albini
//       k_decay    – decay rate [1/s] representing ember burnout (default: 1/600 = 1/10min)
//
//  2. Probabilistic ignition:
//       P_ignition = 1 - exp(-λ)
//       λ = k_ignition × ρ_ember × dt_check / ρ_threshold
//
//     where:
//       ρ_threshold  – minimum ember density for ignition [embers/m²] (default: 10.0)
//       k_ignition   – ignition rate coefficient (default: 1.0)
//
//  3. Fuel moisture adjustment:
//       When fuel moisture M_f is available, the ignition probability is
//       reduced by a moisture damping factor:
//         f_moisture = max(0, (M_x - M_f) / M_x)
//       P_ignition_eff = P_ignition × f_moisture
//
// Usage:
//   EmberAccumulationState state;
//   update_ember_accumulation(state, ember_flux_mf, dt, k_decay);
//   compute_ember_ignition(phi, state, ember_density_mf, params, dt_check, step);
//
// References:
//   Sardoy, N., et al. (2007). Modeling transport and combustion of firebrands.
//     Combustion and Flame, 150(3):151-169.
//   Koo, E., et al. (2010). Modelling firebrand transport in wildfires.
//     Proceedings of the Combustion Institute, 33(2):2449-2456.
//
// Input prefix: "ember_accumulation."
// ============================================================================

#ifndef EMBER_ACCUMULATION_H
#define EMBER_ACCUMULATION_H

#include "constants.H"
#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuLaunch.H>
#include <cmath>
#include <random>

using namespace amrex;

// ============================================================================
// EmberAccumulationState
//
// Carries per-cell ember density as it evolves with landing flux and decay.
// ============================================================================
struct EmberAccumulationState {
    std::unique_ptr<MultiFab> ember_density_mf;  // ember density [embers/m²]
    bool initialized = false;
};

// ============================================================================
// update_ember_accumulation
//
// Update ember density field with new landing flux and exponential decay.
//
// Parameters:
//   state          – accumulation state (updated in-place)
//   ember_flux_mf  – ember landing flux density [embers/m²/s] (comp 0 from cascade)
//   geom           – AMReX geometry
//   dt             – timestep size [s]
//   k_decay        – decay rate [1/s] (default: 1/600 = 10 min burnout)
// ============================================================================
inline void update_ember_accumulation(
    EmberAccumulationState& state,
    const MultiFab&         ember_flux_mf,
    const Geometry&         geom,
    Real                    dt,
    Real                    k_decay = 1.0/600.0)
{
    // Initialize if needed
    if (!state.initialized || state.ember_density_mf == nullptr) {
        state.ember_density_mf = std::make_unique<MultiFab>(
            ember_flux_mf.boxArray(), ember_flux_mf.DistributionMap(), 1, 0);
        state.ember_density_mf->setVal(0.0);
        state.initialized = true;
    }

    MultiFab& rho = *state.ember_density_mf;

    // Update density: dρ/dt = F_landing - k_decay × ρ
    // Semi-implicit update: ρ_new = (ρ_old + F × dt) / (1 + k_decay × dt)
    for (MFIter mfi(rho); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const& rho_arr = rho.array(mfi);
        auto const& flux_arr = ember_flux_mf.const_array(mfi);

        const Real decay_factor = 1.0 / (1.0 + k_decay * dt);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            const Real flux = flux_arr(i, j, k, 0);  // comp 0 = ember flux
            const Real rho_old = rho_arr(i, j, k, 0);
            
            // Semi-implicit update for stability
            rho_arr(i, j, k, 0) = (rho_old + flux * dt) * decay_factor;
        });
    }
}

// ============================================================================
// compute_ember_ignition
//
// Probabilistic ignition based on accumulated ember density.
//
// Parameters:
//   phi            – level-set field (modified: ignited cells set to -spot_radius)
//   state          – ember accumulation state
//   ember_output_mf – diagnostic output MultiFab (comp 0: density, comp 1: ignition flag)
//   fuel_moisture_mf – optional fuel moisture field (nullptr to disable moisture damping)
//   M_x            – moisture of extinction for damping calculation
//   rho_threshold  – minimum density for ignition [embers/m²]
//   k_ignition     – ignition rate coefficient
//   spot_radius    – radius of ignition zone [m]
//   dt_check       – effective time interval for check [s]
//   step           – current timestep (for RNG seed)
// ============================================================================
inline void compute_ember_ignition(
    MultiFab&                      phi,
    EmberAccumulationState&        state,
    MultiFab&                      ember_output_mf,
    const MultiFab*                fuel_moisture_mf = nullptr,
    Real                           M_x              = 0.30,
    Real                           rho_threshold    = 10.0,
    Real                           k_ignition       = 1.0,
    Real                           spot_radius      = 5.0,
    Real                           dt_check         = 1.0,
    int                            step             = 0)
{
    if (!state.initialized || state.ember_density_mf == nullptr) {
        return;
    }

    const MultiFab& rho = *state.ember_density_mf;

    // Copy density to output (comp 0)
    MultiFab::Copy(ember_output_mf, rho, 0, 0, 1, 0);
    
    // Reset ignition flags (comp 1)
    ember_output_mf.setVal(0.0, 1, 1);

    // Random number generator (one per MPI rank, seeded by step)
    std::mt19937 rng(step + amrex::ParallelDescriptor::MyProc() * 12345);
    std::uniform_real_distribution<Real> uniform(0.0, 1.0);

    for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const& phi_arr = phi.array(mfi);
        auto const& rho_arr = rho.const_array(mfi);
        auto const& output_arr = ember_output_mf.array(mfi);

        // Optional fuel moisture damping
        Real M_f = 0.0;
        bool use_moisture_damping = (fuel_moisture_mf != nullptr);

        for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k) {
            for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
                for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
                    const Real density = rho_arr(i, j, k, 0);
                    
                    // Skip if below threshold or already burning
                    if (density < rho_threshold || phi_arr(i, j, k, 0) < 0.0) {
                        continue;
                    }

                    // Compute ignition probability
                    const Real lambda = k_ignition * density * dt_check / rho_threshold;
                    Real P_ign = 1.0 - std::exp(-lambda);

                    // Apply fuel moisture damping if available
                    if (use_moisture_damping) {
                        M_f = (*fuel_moisture_mf).const_array(mfi)(i, j, k, 0);
                        const Real f_moisture = std::max(0.0, (M_x - M_f) / M_x);
                        P_ign *= f_moisture;
                    }

                    // Stochastic ignition
                    if (uniform(rng) < P_ign) {
                        phi_arr(i, j, k, 0) = -spot_radius;
                        output_arr(i, j, k, 1) = 1.0;  // mark ignition
                    }
                }
            }
        }
    }
}

#endif // EMBER_ACCUMULATION_H
