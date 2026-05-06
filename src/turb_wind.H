// ============================================================================
// turb_wind.H  –  Turbulent Wind Perturbation Model
// ============================================================================
//
// Provides two turbulent wind models that add stochastic perturbations to the
// ambient wind field at every timestep.  In both cases the unperturbed base
// wind is preserved in a separate MultiFab so the perturbation is always
// computed relative to the clean background field.
//
// Models
// ------
//   "ou_process"
//       Ornstein-Uhlenbeck (OU) temporally-correlated noise (Uhlenbeck &
//       Ornstein 1930).  Two sub-modes:
//
//         L_c = 0  (domain-uniform):
//           A single pair (u', v') is drawn and added uniformly to every cell.
//           All cells experience the same gust event – appropriate when the
//           domain is small relative to the turbulent integral scale.
//
//         L_c > 0  (spatially correlated):
//           Independent per-cell OU states (u'(x,y), v'(x,y)) are maintained
//           in a MultiFab.  The white-noise driving term is smoothed each step
//           with a separable Gaussian kernel of standard deviation
//           sigma_k = L_c / dx [cells] before being injected into the OU
//           update.  This yields spatially correlated perturbations whose
//           e-folding length is approximately L_c.  After smoothing the noise
//           field is rescaled to unit RMS so the stationary standard deviation
//           of each cell's perturbation remains exactly sigma [m/s].
//           Boundary clamping (Neumann condition) is used at subdomain edges,
//           making the method MPI-compatible without FillBoundary calls.
//
//         Discrete-time OU update (exact Euler-Maruyama):
//           alpha       = exp(-theta * dt)
//           sigma_step  = sigma * sqrt(1 - alpha^2)
//           u'(t+dt)    = alpha * u'(t) + sigma_step * eta_u(t)
//           v'(t+dt)    = alpha * v'(t) + sigma_step * eta_v(t)
//         where eta is unit-variance (spatially correlated) white noise.
//         Stationary std dev of each cell = sigma [m/s].
//         Temporal decorrelation time = 1/theta [s].
//
//   "direction_walk"
//       Bounded cumulative random walk of wind direction (domain-uniform).
//       A Gaussian angular increment (sigma_theta [rad/step]) is drawn each
//       timestep and accumulated; the total angle is clamped to ±theta_max.
//       The entire base wind field is rotated by the cumulative angle.
//       Wind speed is preserved exactly; only direction fluctuates.
//
// Parameters (all under the "turb_wind." prefix in the input file)
// ----------------------------------------------------------------
//   model           – "none" / "ou_process" / "direction_walk"  (default: "none")
//   theta           – OU reversion rate [1/s]  (default: 0.1)
//   sigma           – OU stationary std dev [m/s]  (default: 0.5)
//   L_c             – spatial correlation length [m]; 0 = domain-uniform  (default: 0.0)
//   sigma_theta     – direction_walk angular std dev [rad/step]  (default: 0.1)
//   theta_max       – direction_walk max cumulative deviation [rad]  (default: 0.5236 ≈ 30°)
//   random_seed     – RNG seed; 0 = use system clock  (default: 0)
//
// References
// ----------
//   Uhlenbeck, G.E. & Ornstein, L.S. (1930). On the theory of the Brownian
//     motion. Phys. Rev. 36(5):823–841.
//   Finney, M.A. et al. (2011). Role of wind, fuel moisture, and terrain in
//     controlling fire movement. Ecosphere 2(3):art17.
//   Cruz, M.G. et al. (2019). Uncertainty in wildfire behaviour research.
//     Current Forestry Reports 5:155–172.
// ============================================================================

#ifndef TURB_WIND_H
#define TURB_WIND_H

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuLaunch.H>
#include "parse_inputs.H"

#include <algorithm>
#include <cmath>
#include <ctime>
#include <random>
#include <vector>

using namespace amrex;

// ============================================================================
// TurbWindState
// Holds the scalar random state used by domain-uniform OU and direction_walk.
// The per-cell OU state (L_c > 0) is stored in the ou_state MultiFab in main.
// ============================================================================
struct TurbWindState {
    double u_prime;      // domain-uniform OU u-perturbation [m/s]
    double v_prime;      // domain-uniform OU v-perturbation [m/s]
    double angle_cumul;  // accumulated direction angle [rad] (direction_walk)
    std::mt19937 rng;
    std::normal_distribution<double> norm_dist;
};

// ============================================================================
// init_turb_wind_state
// ============================================================================
static void init_turb_wind_state(TurbWindState& state, int seed)
{
    state.u_prime    = 0.0;
    state.v_prime    = 0.0;
    state.angle_cumul = 0.0;
    const unsigned int s = (seed == 0)
                           ? static_cast<unsigned int>(std::time(nullptr))
                           : static_cast<unsigned int>(seed);
    state.rng.seed(s);
    state.norm_dist = std::normal_distribution<double>(0.0, 1.0);
    amrex::Print() << "TurbWindState initialised (seed=" << s << ")\n";
}

// ============================================================================
// apply_gaussian_spatial_kernel  (internal helper)
//
// In-place separable 1-D Gaussian smoothing on a 2-component MultiFab,
// applied as two sequential passes (x then y).  The z-dimension is traversed
// normally; each (i,j) slice at every k receives independent smoothing in
// the horizontal plane, which is appropriate for 2-D surface fire simulations
// embedded in a 3-D domain.
//
// Boundary treatment: Neumann clamping – the index is clamped to the valid
// box edges so no ghost cells or FillBoundary calls are required.  This makes
// the algorithm MPI-safe (each subdomain is processed independently).
//
// After both passes the field is rescaled to unit RMS so the OU sigma
// parameter retains its meaning as the stationary standard deviation.
//
// Parameters:
//   mf      – in/out: 2-component MultiFab (white noise → correlated noise)
//   sigma_k – Gaussian std-dev of the kernel [cells]
// ============================================================================
static void apply_gaussian_spatial_kernel(MultiFab& mf, double sigma_k)
{
    if (sigma_k < 1.0e-6) return;   // nothing to do for zero correlation length

    // Kernel half-radius: 3-sigma rule (captures > 99.7 % of the Gaussian mass)
    const int r = static_cast<int>(std::ceil(3.0 * sigma_k));

    // Pre-compute normalised Gaussian weights for offsets [-r, r]
    const int ksize = 2 * r + 1;
    std::vector<double> kw(ksize);
    double ksum = 0.0;
    for (int d = -r; d <= r; ++d) {
        kw[d + r] = std::exp(-0.5 * d * d / (sigma_k * sigma_k));
        ksum += kw[d + r];
    }
    for (auto& w : kw) w /= ksum;

    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        const int ilo = bx.smallEnd(0), ihi = bx.bigEnd(0);
        const int jlo = bx.smallEnd(1), jhi = bx.bigEnd(1);
        const int klo = bx.smallEnd(2), khi = bx.bigEnd(2);

        const int nx = ihi - ilo + 1;
        const int ny = jhi - jlo + 1;
        const int nz = khi - klo + 1;

        // Scratch buffer for one pass: 2 components × nz × ny × nx
        std::vector<double> tmp(static_cast<std::size_t>(2 * nz * ny * nx), 0.0);

        auto fn = mf.array(mfi);

        // ---- Pass 1: smooth in x-direction --------------------------------
        amrex::LoopOnCpu(bx, [&](int i, int j, int k) noexcept {
            for (int c = 0; c < 2; ++c) {
                double val = 0.0;
                for (int d = -r; d <= r; ++d) {
                    const int ii = std::max(ilo, std::min(ihi, i + d));
                    val += kw[d + r] * static_cast<double>(fn(ii, j, k, c));
                }
                const std::size_t idx =
                    static_cast<std::size_t>(c) * nz * ny * nx
                    + static_cast<std::size_t>(k - klo) * ny * nx
                    + static_cast<std::size_t>(j - jlo) * nx
                    + static_cast<std::size_t>(i - ilo);
                tmp[idx] = val;
            }
        });

        // Write x-smoothed back into mf
        amrex::LoopOnCpu(bx, [&](int i, int j, int k) noexcept {
            for (int c = 0; c < 2; ++c) {
                const std::size_t idx =
                    static_cast<std::size_t>(c) * nz * ny * nx
                    + static_cast<std::size_t>(k - klo) * ny * nx
                    + static_cast<std::size_t>(j - jlo) * nx
                    + static_cast<std::size_t>(i - ilo);
                fn(i, j, k, c) = static_cast<Real>(tmp[idx]);
            }
        });

        // ---- Pass 2: smooth in y-direction --------------------------------
        std::fill(tmp.begin(), tmp.end(), 0.0);

        amrex::LoopOnCpu(bx, [&](int i, int j, int k) noexcept {
            for (int c = 0; c < 2; ++c) {
                double val = 0.0;
                for (int d = -r; d <= r; ++d) {
                    const int jj = std::max(jlo, std::min(jhi, j + d));
                    val += kw[d + r] * static_cast<double>(fn(i, jj, k, c));
                }
                const std::size_t idx =
                    static_cast<std::size_t>(c) * nz * ny * nx
                    + static_cast<std::size_t>(k - klo) * ny * nx
                    + static_cast<std::size_t>(j - jlo) * nx
                    + static_cast<std::size_t>(i - ilo);
                tmp[idx] = val;
            }
        });

        amrex::LoopOnCpu(bx, [&](int i, int j, int k) noexcept {
            for (int c = 0; c < 2; ++c) {
                const std::size_t idx =
                    static_cast<std::size_t>(c) * nz * ny * nx
                    + static_cast<std::size_t>(k - klo) * ny * nx
                    + static_cast<std::size_t>(j - jlo) * nx
                    + static_cast<std::size_t>(i - ilo);
                fn(i, j, k, c) = static_cast<Real>(tmp[idx]);
            }
        });
    }

    // ---- Rescale to unit RMS per component (MPI-safe via mf.norm2) ---------
    // mf.norm2(c) = sqrt(sum_cells(val^2)); divide by sqrt(N) to get RMS.
    const Real N = static_cast<Real>(mf.boxArray().numPts());
    for (int c = 0; c < 2; ++c) {
        const Real rms = (N > Real(0.0)) ? mf.norm2(c) / std::sqrt(N) : Real(1.0);
        if (rms > Real(1.0e-12)) {
            mf.mult(Real(1.0) / rms, c, 1);
        }
    }
}

// ============================================================================
// apply_turb_wind
//
// Updates the turbulent state and writes the perturbed wind into 'vel' by
// copying from 'vel_base' and adding the current perturbation.  Must be
// called after vel_base has been set for this timestep (either constant or
// freshly loaded from a time-dependent wind file).
//
// Parameters:
//   vel       – output: perturbed wind MultiFab (overwritten)
//   vel_base  – input:  unperturbed background wind MultiFab (preserved)
//   ou_state  – in/out: per-cell OU state MultiFab (2 comps); nullptr means
//               domain-uniform mode (L_c = 0) or direction_walk
//   state     – in/out: scalar RNG / domain-uniform OU / walk state
//   dt        – timestep [s]
//   tp        – turbulent wind parameters
//   geom      – domain geometry (used for dx when L_c > 0)
// ============================================================================
static void apply_turb_wind(MultiFab&         vel,
                             const MultiFab&   vel_base,
                             MultiFab*         ou_state,
                             TurbWindState&    state,
                             amrex::Real       dt,
                             const InputParameters::TurbWindParams& tp,
                             const Geometry&   geom)
{
    if (tp.model == "ou_process") {

        const double alpha      = std::exp(-tp.theta * static_cast<double>(dt));
        const double sigma_step = tp.sigma * std::sqrt(1.0 - alpha * alpha);

        if (ou_state == nullptr) {
            // ------------------------------------------------------------------
            // Domain-uniform OU (L_c = 0)
            // ------------------------------------------------------------------
            state.u_prime = alpha * state.u_prime
                            + sigma_step * state.norm_dist(state.rng);
            state.v_prime = alpha * state.v_prime
                            + sigma_step * state.norm_dist(state.rng);

            const Real up = static_cast<Real>(state.u_prime);
            const Real vp = static_cast<Real>(state.v_prime);

            for (MFIter mfi(vel); mfi.isValid(); ++mfi) {
                const Box& bx = mfi.validbox();
                auto const vb = vel_base.const_array(mfi);
                auto       v  = vel.array(mfi);
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                    v(i, j, k, 0) = vb(i, j, k, 0) + up;
                    v(i, j, k, 1) = vb(i, j, k, 1) + vp;
                    v(i, j, k, 2) = vb(i, j, k, 2);
                });
            }
            amrex::Print() << "  OU wind (uniform): du=" << state.u_prime
                           << "  dv=" << state.v_prime << " m/s\n";

        } else {
            // ------------------------------------------------------------------
            // Per-cell spatially-correlated OU (L_c > 0)
            // ------------------------------------------------------------------

            // Step 1: generate per-cell standard normal white noise (CPU)
            MultiFab white_noise(ou_state->boxArray(),
                                 ou_state->DistributionMap(), 2, 0);
            for (MFIter mfi(white_noise); mfi.isValid(); ++mfi) {
                const Box& bx = mfi.validbox();
                auto wn = white_noise.array(mfi);
                amrex::LoopOnCpu(bx, [&](int i, int j, int k) noexcept {
                    wn(i, j, k, 0) = static_cast<Real>(state.norm_dist(state.rng));
                    wn(i, j, k, 1) = static_cast<Real>(state.norm_dist(state.rng));
                });
            }

            // Step 2: apply separable Gaussian spatial kernel
            const double sigma_k = tp.L_c / static_cast<double>(geom.CellSize(0));
            apply_gaussian_spatial_kernel(white_noise, sigma_k);

            // Step 3: update per-cell OU state (CPU, uses smoothed noise)
            const Real a  = static_cast<Real>(alpha);
            const Real ss = static_cast<Real>(sigma_step);
            for (MFIter mfi(*ou_state); mfi.isValid(); ++mfi) {
                const Box& bx = mfi.validbox();
                auto       st = ou_state->array(mfi);
                auto const wn = white_noise.const_array(mfi);
                amrex::LoopOnCpu(bx, [&](int i, int j, int k) noexcept {
                    st(i, j, k, 0) = a * st(i, j, k, 0) + ss * wn(i, j, k, 0);
                    st(i, j, k, 1) = a * st(i, j, k, 1) + ss * wn(i, j, k, 1);
                });
            }

            // Step 4: vel = vel_base + ou_state (GPU-accelerated)
            for (MFIter mfi(vel); mfi.isValid(); ++mfi) {
                const Box& bx = mfi.validbox();
                auto const vb = vel_base.const_array(mfi);
                auto const st = ou_state->const_array(mfi);
                auto       v  = vel.array(mfi);
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                    v(i, j, k, 0) = vb(i, j, k, 0) + st(i, j, k, 0);
                    v(i, j, k, 1) = vb(i, j, k, 1) + st(i, j, k, 1);
                    v(i, j, k, 2) = vb(i, j, k, 2);
                });
            }
        }

    } else if (tp.model == "direction_walk") {
        // ------------------------------------------------------------------
        // Bounded cumulative random walk of wind direction (domain-uniform)
        // ------------------------------------------------------------------
        const double delta = tp.sigma_theta * state.norm_dist(state.rng);
        state.angle_cumul += delta;
        if (state.angle_cumul >  tp.theta_max) state.angle_cumul =  tp.theta_max;
        if (state.angle_cumul < -tp.theta_max) state.angle_cumul = -tp.theta_max;

        const Real cos_a = static_cast<Real>(std::cos(state.angle_cumul));
        const Real sin_a = static_cast<Real>(std::sin(state.angle_cumul));

        for (MFIter mfi(vel); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            auto const vb = vel_base.const_array(mfi);
            auto       v  = vel.array(mfi);
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                const Real ub = vb(i, j, k, 0);
                const Real vb_y = vb(i, j, k, 1);
                v(i, j, k, 0) =  ub * cos_a - vb_y * sin_a;
                v(i, j, k, 1) =  ub * sin_a + vb_y * cos_a;
                v(i, j, k, 2) = vb(i, j, k, 2);
            });
        }
        amrex::Print() << "  Direction-walk angle: "
                       << state.angle_cumul * (180.0 / 3.14159265358979) << " deg\n";
    }
}

#endif // TURB_WIND_H
