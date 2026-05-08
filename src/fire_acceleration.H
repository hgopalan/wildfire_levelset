// ============================================================================
// Fire Acceleration Model (Rothermel 1983 / Catchpole et al. 1992)
// ============================================================================
// Small fires do not immediately achieve their quasi-steady-state (QSS)
// rate of spread.  As the fire grows, convective preheating of unburned fuel
// ahead of the front increases, allowing the spread rate to accelerate from
// near zero (at ignition) toward the QSS value predicted by Rothermel (1972).
//
// This model applies a globally computed acceleration factor to the ROS field
// based on the current fire size.  The factor approaches 1 once the fire has
// grown large enough that preheating is no longer limiting.
//
// Acceleration factor (Catchpole et al. 1992 / Rothermel 1983 style):
//
//   r_fire [m]  = sqrt(A_fire / π)   (effective fire radius from burned area)
//   α_acc       = 1 − exp(−r_fire / L_acc)
//   R_eff       = α_acc × R_QSS
//
// where:
//   A_fire   = burned area [m²]  (= n_burned_cells × dx × dy)
//   L_acc    = acceleration length scale [m]; fire achieves ~63% of QSS at
//              r_fire = L_acc, and ~95% at r_fire = 3 L_acc.
//              Typical values: 30–100 m depending on fuel type.
//
// The model is disabled (α = 1) once r_fire >> L_acc, so it has no effect on
// large fires and does not alter the steady-state behaviour.
//
// References:
//   Rothermel, R.C. (1983). "How to predict the spread and intensity of forest
//       and range fires." USDA Forest Service GTR INT-143.
//   Catchpole, E.A., de Mestre, N.J. & Gill, A.M. (1992). "Intensity of fire
//       at its perimeter." Australian Journal of Ecology, 17(1), 1–4.
//
// Usage (in main.cpp, before propagation):
//   apply_fire_acceleration(R_mf, phi, geom, inputs.acceleration);
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
// apply_fire_acceleration
//
// Computes the current fire area from phi and scales R_mf by the resulting
// acceleration factor α_acc = 1 − exp(−r_fire / L_acc).
//
// When the fire is large (r_fire >> L_acc), α_acc → 1 and R_mf is unchanged.
// When the fire is small (r_fire << L_acc), α_acc → 0 and R_mf → 0.
//
// Parameters:
//   R_mf    – rate-of-spread field [m/s] (modified in place)
//   phi     – level-set field (phi < 0 = burned)
//   geom    – AMReX geometry (provides cell sizes)
//   params  – AccelerationParams from InputParameters
// ============================================================================
inline void apply_fire_acceleration(MultiFab&       R_mf,
                                    const MultiFab& phi,
                                    const Geometry& geom,
                                    const InputParameters::AccelerationParams& params)
{
    if (params.enable != 1) return;

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

    if (n_burned == 0) return;  // no fire yet — skip (no scaling needed)

    // Burned area and effective fire radius
    const Real A_fire  = static_cast<Real>(n_burned) * dx[0] * dx[1];
    const Real r_fire  = std::sqrt(A_fire / WildfireConst::PI);
    const Real L_acc   = amrex::max(static_cast<Real>(params.L_acc), Real(1.0));
    const Real alpha   = Real(1.0) - std::exp(-r_fire / L_acc);

    // Nothing to do when essentially at full speed (saves a kernel launch)
    if (alpha >= Real(1.0) - Real(1.0e-4)) return;

    amrex::Print() << "  Fire acceleration: r_fire=" << r_fire
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

#endif // FIRE_ACCELERATION_H
