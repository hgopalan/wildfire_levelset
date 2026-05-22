// ============================================================================
// Fire Acceleration Models
// ============================================================================
// Small fires do not immediately achieve their quasi-steady-state (QSS)
// rate of spread. This module implements two acceleration models:
//
// 1. SIZE-BASED MODEL (Rothermel 1983 / Catchpole et al. 1992)
//    - Simple exponential approach based on current fire size
//    - Applies global scaling factor to ROS field
//    - α_acc = 1 − exp(−r_fire / L_acc)
//
// 2. FARSITE TEMPORAL MODEL (McAlpine & Wakimoto 1991 / VanWagner)
//    - Per-cell temporal tracking of acceleration state
//    - Based on elapsed time since ignition or wind change
//    - R(t) = R_E × (1 − exp(−A × t))
//    - Includes optional wind-onset time-lag capability
//
// References:
//   Rothermel, R.C. (1983). "How to predict the spread and intensity of forest
//       and range fires." USDA Forest Service GTR INT-143.
//   Catchpole, E.A., de Mestre, N.J. & Gill, A.M. (1992). "Intensity of fire
//       at its perimeter." Australian Journal of Ecology, 17(1), 1–4.
//   McAlpine, R.S. & Wakimoto, R.H. (1991). "The acceleration of fire from
//       point source to equilibrium spread." Forest Science, 37(5), 1314–1337.
//   Alexander, M.E., Stocks, B.J. & Lawson, B.D. (1992). "Fire behavior in
//       Black Spruce-lichen woodland." Info. Rep. NOR-X-310. USDA/CFS.
//   Finney, M.A. (1998/2004). "FARSITE: Fire Area Simulator." USDA Forest
//       Service RMRS-RP-4.
//
// Usage:
//   // In main.cpp setup (if using temporal model):
//   if (inputs.acceleration.enable == 1 && inputs.acceleration.use_temporal == 1) {
//       f.accel_state_mf = std::make_unique<MultiFab>(ba, dm, 3, 0);
//       f.accel_state_mf->setVal(Real(0.0));
//   }
//
//   // Before propagation (each timestep):
//   apply_fire_acceleration(R_mf, phi, geom, inputs.acceleration, dt,
//                           f.accel_state_mf.get(), prev_R_equilibrium);
// ============================================================================

#ifndef FIRE_ACCELERATION_H
#define FIRE_ACCELERATION_H

#include "constants.H"
#include "parse_inputs.H"
#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_ParallelDescriptor.H>
#include <cmath>

using namespace amrex;

// ============================================================================
// apply_fire_acceleration_size_based
//
// SIZE-BASED MODEL: Computes fire area from phi and scales R_mf by
// α_acc = 1 − exp(−r_fire / L_acc).
//
// When the fire is large (r_fire >> L_acc), α_acc → 1 and R_mf is unchanged.
// When the fire is small (r_fire << L_acc), α_acc → 0 and R_mf → 0.
//
// Parameters:
//   R_mf    – rate-of-spread field [m/s] (modified in place)
//   phi     – level-set field (phi < 0 = burned)
//   geom    – AMReX geometry (provides cell sizes)
//   L_acc   – acceleration length scale [m]
// ============================================================================
inline void apply_fire_acceleration_size_based(MultiFab&       R_mf,
                                               const MultiFab& phi,
                                               const Geometry& geom,
                                               const Real      L_acc)
{
    const auto dx = geom.CellSize();

    // Count burned cells (CPU loop; same pattern as perimeter counting in main.cpp)
    amrex::Long n_burned = 0;
    for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const p  = phi.const_array(mfi);
        amrex::LoopOnCpu(bx, [&](int i, int j, int k) {
            if (p(i, j, k) < Real(0.0)) ++n_burned;
        });
    }
    ParallelDescriptor::ReduceLongSum(n_burned);

    if (n_burned == 0) return;  // no fire yet — skip

    // Burned area and effective fire radius
    const Real A_fire  = static_cast<Real>(n_burned) * dx[0] * dx[1];
    const Real r_fire  = std::sqrt(A_fire / WildfireConst::PI);
    const Real L_acc_safe = amrex::max(L_acc, Real(1.0));  // floor prevents divide-by-zero
    const Real alpha   = Real(1.0) - std::exp(-r_fire / L_acc_safe);

    // Nothing to do when essentially at full speed (saves a kernel launch)
    if (alpha >= Real(1.0) - Real(1.0e-4)) return;

    amrex::Print() << "  Fire acceleration (size-based): r_fire=" << r_fire
                   << " m  L_acc=" << L_acc << " m  alpha=" << alpha << "\n";

    // Scale R_mf by alpha
    for (MFIter mfi(R_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto R = R_mf.array(mfi);
        const Real a = alpha;
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            R(i, j, k) *= a;
        });
    }
}

// ============================================================================
// apply_fire_acceleration_temporal
//
// FARSITE TEMPORAL MODEL: Per-cell acceleration using VanWagner's equation
// R(t) = R_E × (1 − exp(−A × t)).
//
// Each cell tracks:
//   - ROS_current: current instantaneous ROS [m/s]
//   - elapsed_time: time since entering current state [s]
//   - ROS_equilibrium: target equilibrium ROS [m/s]
//
// When R_mf (quasi-steady Rothermel ROS) changes (e.g., due to wind change),
// the acceleration resets and ramps up exponentially.
//
// Parameters:
//   R_mf            – equilibrium ROS from Rothermel model [m/s] (read; output written to R_out)
//   phi             – level-set field (phi < 0 = burned)
//   geom            – AMReX geometry
//   params          – acceleration parameters
//   dt              – timestep [s]
//   accel_state_mf  – 3-component state (ROS_current, elapsed_time, ROS_equilibrium)
// ============================================================================
inline void apply_fire_acceleration_temporal(MultiFab&       R_mf,
                                             const MultiFab& phi,
                                             const Geometry& geom,
                                             const InputParameters::AccelerationParams& params,
                                             const Real      dt,
                                             MultiFab*       accel_state_mf)
{
    if (accel_state_mf == nullptr) {
        amrex::Abort("apply_fire_acceleration_temporal: accel_state_mf is null");
    }

    const auto dx = geom.CellSize();

    // Compute fire perimeter to determine acceleration constant
    // (point ignition: A_point; line ignition: A_line)
    amrex::Long n_perimeter = 0;
    for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const p  = phi.const_array(mfi);
        amrex::LoopOnCpu(bx, [&](int i, int j, int k) {
            if (p(i, j, k) < Real(0.0)) {
                // Check if this is a perimeter cell (has unburned neighbor)
                bool is_perimeter = false;
                if (p(i+1, j, k) >= Real(0.0)) is_perimeter = true;
                if (p(i-1, j, k) >= Real(0.0)) is_perimeter = true;
                if (p(i, j+1, k) >= Real(0.0)) is_perimeter = true;
                if (p(i, j-1, k) >= Real(0.0)) is_perimeter = true;
                if (is_perimeter) ++n_perimeter;
            }
        });
    }
    ParallelDescriptor::ReduceLongSum(n_perimeter);

    // Estimate perimeter length in meters
    const Real perim_length = static_cast<Real>(n_perimeter) * std::sqrt(dx[0] * dx[0] + dx[1] * dx[1]);

    // Select acceleration constant based on perimeter size
    const Real A = (perim_length > params.perim_limit) ? params.A_line : params.A_point;
    
    // Convert A from [1/min] to [1/s]
    const Real A_per_sec = A / Real(60.0);

    // Update acceleration state for each cell
    for (MFIter mfi(R_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto R_eq    = R_mf.const_array(mfi);          // equilibrium ROS from Rothermel
        auto state   = accel_state_mf->array(mfi);     // [ROS_current, elapsed_time, ROS_equilibrium]
        auto p       = phi.const_array(mfi);           // level set
        auto R_out   = R_mf.array(mfi);                // output (will be modified)

        const Real dt_loc  = dt;
        const Real A_loc   = A_per_sec;
        const Real eps     = Real(1.0e-6);
        const int  enable_wind_lag = params.enable_wind_lag;
        const Real tau_wind = params.tau_wind;

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            // Only accelerate burned cells
            if (p(i, j, k) >= Real(0.0)) {
                state(i, j, k, 0) = Real(0.0);  // ROS_current
                state(i, j, k, 1) = Real(0.0);  // elapsed_time
                state(i, j, k, 2) = Real(0.0);  // ROS_equilibrium
                return;
            }

            const Real R_E = R_eq(i, j, k);  // equilibrium ROS (target)
            Real R_curr    = state(i, j, k, 0);
            Real t_elapsed = state(i, j, k, 1);
            Real R_prev_eq = state(i, j, k, 2);

            // Detect equilibrium ROS change (wind change or fuel change)
            constexpr Real EQUILIBRIUM_CHANGE_THRESHOLD = Real(1.0e-6);
            const Real dR_eq = amrex::Math::abs(R_E - R_prev_eq);
            const bool eq_changed = (dR_eq > EQUILIBRIUM_CHANGE_THRESHOLD * amrex::max(R_E, R_prev_eq, Real(1.0)));

            if (eq_changed) {
                // Equilibrium changed — apply wind-lag or reset
                if (enable_wind_lag == 1 && R_E > R_prev_eq) {
                    // Wind increase: apply exponential lag
                    // New target approaches R_E with time constant tau_wind
                    // R_target(t) = R_prev + (R_E - R_prev) × (1 - exp(-dt/tau_wind))
                    const Real lag_factor = Real(1.0) - std::exp(-dt_loc / tau_wind);
                    const Real R_target = R_prev_eq + (R_E - R_prev_eq) * lag_factor;
                    
                    // Start new acceleration towards lagged target
                    R_prev_eq = R_target;
                    t_elapsed = Real(0.0);
                } else {
                    // Wind decrease or lag disabled: reset immediately to new equilibrium
                    R_prev_eq = R_E;
                    t_elapsed = Real(0.0);
                }
            }

            // VanWagner acceleration equation: R(t) = R_E × (1 - exp(-A × t))
            t_elapsed += dt_loc;
            R_curr = R_prev_eq * (Real(1.0) - std::exp(-A_loc * t_elapsed));

            // Clamp to equilibrium (avoid overshoot due to numerical error)
            if (R_curr > R_prev_eq) R_curr = R_prev_eq;

            // Update state
            state(i, j, k, 0) = R_curr;
            state(i, j, k, 1) = t_elapsed;
            state(i, j, k, 2) = R_prev_eq;

            // Write accelerated ROS to output
            R_out(i, j, k) = R_curr;
        });
    }

    if (ParallelDescriptor::IOProcessor()) {
        amrex::Print() << "  Fire acceleration (temporal): perim=" << perim_length << " m"
                       << "  A=" << A << " 1/min"
                       << "  (using " << (perim_length > params.perim_limit ? "line" : "point") << " constant)\n";
    }
}

// ============================================================================
// apply_fire_acceleration
//
// Main dispatch function: selects size-based or temporal model based on params.
//
// Parameters:
//   R_mf            – rate-of-spread field [m/s] (modified in place)
//   phi             – level-set field (phi < 0 = burned)
//   geom            – AMReX geometry
//   params          – AccelerationParams from InputParameters
//   dt              – timestep [s] (only used for temporal model)
//   accel_state_mf  – per-cell state for temporal model (nullptr for size-based)
// ============================================================================
inline void apply_fire_acceleration(MultiFab&       R_mf,
                                    const MultiFab& phi,
                                    const Geometry& geom,
                                    const InputParameters::AccelerationParams& params,
                                    const Real      dt = Real(0.0),
                                    MultiFab*       accel_state_mf = nullptr)
{
    if (params.enable != 1) return;

    if (params.use_temporal == 1) {
        // FARSITE temporal model
        apply_fire_acceleration_temporal(R_mf, phi, geom, params, dt, accel_state_mf);
    } else {
        // Simple size-based model
        apply_fire_acceleration_size_based(R_mf, phi, geom, params.L_acc);
    }
}

#endif // FIRE_ACCELERATION_H
