// ============================================================================
// turb_wind.H  –  Turbulent Wind Perturbation Model
// ============================================================================
//
// Provides three active turbulent wind models (plus "none") that add stochastic perturbations to
// the ambient wind field at every timestep.  The unperturbed base wind is
// preserved in a separate MultiFab so the perturbation is always computed
// relative to the clean background field.
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
//   "spectral_noise"
//       Random Fourier Feature (RFF) spectral noise with OU temporal
//       evolution (Kraichnan 1970; Rahimi & Recht 2007).  A set of N_modes
//       random wavenumber pairs (kx_n, ky_n) is drawn at initialisation from
//       a 2-D isotropic Gaussian power spectrum with characteristic length
//       scale L_c:
//           kx_n, ky_n  ~  N(0, (1/L_c)^2)   independently
//       Each mode n also has a fixed random phase (phi_u_n, phi_v_n) and a
//       scalar OU amplitude (A_u_n, A_v_n) that evolves according to the
//       standard OU update with the user's theta and unit internal variance.
//       The perturbation field is reconstructed on the GPU at each timestep:
//           u'(x,y) = sigma * sqrt(2/N) * Σ_n A_u_n * cos(kx_n*x + ky_n*y + phi_u_n)
//           v'(x,y) = sigma * sqrt(2/N) * Σ_n A_v_n * cos(kx_n*x + ky_n*y + phi_v_n)
//       This yields:
//         - Temporal decorrelation time = 1/theta [s] (same as ou_process)
//         - Approximate 2-D autocovariance: C(r) ≈ sigma^2 * exp(-r^2 / (2*L_c^2))
//         - Correct energy distribution across wavenumbers (Gaussian PSD)
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
//   model           – "none" / "ou_process" / "spectral_noise" / "direction_walk"
//   theta           – OU reversion rate [1/s]  (default: 0.1)
//   sigma           – OU stationary std dev [m/s]  (default: 0.5)
//   L_c             – spatial correlation length [m]; 0 = domain-uniform  (default: 0.0)
//   N_modes         – number of Fourier modes for spectral_noise  (default: 32)
//   sigma_theta     – direction_walk angular std dev [rad/step]  (default: 0.1)
//   theta_max       – direction_walk max cumulative deviation [rad]  (default: 0.5236 ≈ 30°)
//   random_seed     – RNG seed; 0 = use system clock  (default: 0)
//
// References
// ----------
//   Uhlenbeck, G.E. & Ornstein, L.S. (1930). On the theory of the Brownian
//     motion. Phys. Rev. 36(5):823–841.
//   Kraichnan, R.H. (1970). Diffusion by a random velocity field.
//     Phys. Fluids 13(1):22–31.
//   Rahimi, A. & Recht, B. (2007). Random features for large-scale kernel
//     machines. NIPS 2007.
//   Finney, M.A. et al. (2011). Role of wind, fuel moisture, and terrain in
//     controlling fire movement. Ecosphere 2(3):art17.
//   Cruz, M.G. et al. (2019). Uncertainty in wildfire behaviour research.
//     Current Forestry Reports 5:155–172.
// ============================================================================

#ifndef TURB_WIND_H
#define TURB_WIND_H

#include <AMReX.H>
#include <AMReX_Gpu.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuLaunch.H>
#include "parse_inputs.H"

#include <algorithm>
#include <cmath>
#include <ctime>
#include <random>
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace amrex;

// ============================================================================
// TurbWindState
// Holds all random state required by the turbulent wind models.
//
//   Domain-uniform OU / direction_walk:
//     u_prime, v_prime, angle_cumul – scalar state variables.
//
//   Spectral noise:
//     kx, ky     – fixed wavenumber vectors drawn once at init.
//     phi_u, phi_v – fixed random phases drawn once at init.
//     A_u, A_v   – scalar OU amplitude per mode (evolve each step on CPU).
//     d_*        – device copies of the above for GPU reconstruction.
// ============================================================================
struct TurbWindState {
    // Domain-uniform OU scalars (also used for spectral if L_c = 0)
    double u_prime    = 0.0;
    double v_prime    = 0.0;
    double angle_cumul = 0.0;

    // Random engine (seeded at init; used by all models)
    std::mt19937 rng;
    std::normal_distribution<double> norm_dist{0.0, 1.0};

    // ---------- Spectral noise state ------------------------------------
    int N_modes = 0;

    // Host-side fixed spectral mode data (initialised once, never changed)
    std::vector<amrex::Real> kx;     // wavenumber x-component [rad/m]
    std::vector<amrex::Real> ky;     // wavenumber y-component [rad/m]
    std::vector<amrex::Real> phi_u;  // phase for u' component [rad]
    std::vector<amrex::Real> phi_v;  // phase for v' component [rad]

    // Host-side evolving OU amplitudes (unit internal variance)
    std::vector<amrex::Real> A_u;   // u' OU coefficient per mode
    std::vector<amrex::Real> A_v;   // v' OU coefficient per mode

    // Device-side copies for the GPU reconstruction kernel
    amrex::Gpu::DeviceVector<amrex::Real> d_kx, d_ky;
    amrex::Gpu::DeviceVector<amrex::Real> d_phi_u, d_phi_v;
    amrex::Gpu::DeviceVector<amrex::Real> d_A_u, d_A_v;
};

// ============================================================================
// init_turb_wind_state
//
// Seeds the RNG and (for "spectral_noise") draws the fixed wavenumber vectors
// and random phases from the 2-D isotropic Gaussian spectrum with length
// scale L_c.  Must be called once before the time loop.
// ============================================================================
inline void init_turb_wind_state(TurbWindState& state,
                                  const InputParameters::TurbWindParams& tp)
{
    const unsigned int s = (tp.random_seed == 0)
                           ? static_cast<unsigned int>(std::time(nullptr))
                           : static_cast<unsigned int>(tp.random_seed);
    state.rng.seed(s);
    state.norm_dist = std::normal_distribution<double>(0.0, 1.0);
    state.u_prime    = 0.0;
    state.v_prime    = 0.0;
    state.angle_cumul = 0.0;

    amrex::Print() << "TurbWindState initialised (seed=" << s << ")\n";

    if (tp.model == "spectral_noise") {
        const int N = tp.N_modes;
        state.N_modes = N;

        state.kx.resize(N);   state.ky.resize(N);
        state.phi_u.resize(N); state.phi_v.resize(N);
        state.A_u.assign(N, amrex::Real(0.0));
        state.A_v.assign(N, amrex::Real(0.0));

        // Sample wavenumbers from 2-D isotropic Gaussian PSD:
        //   kx_n, ky_n ~ N(0, (1/L_c)^2)  independently
        // This gives the Gaussian autocorrelation C(r) ~ exp(-r^2 / (2*L_c^2)).
        std::normal_distribution<double> k_dist(0.0, 1.0 / static_cast<double>(tp.L_c));
        std::uniform_real_distribution<double> phi_dist(0.0, 2.0 * M_PI);

        for (int n = 0; n < N; ++n) {
            state.kx[n]    = static_cast<amrex::Real>(k_dist(state.rng));
            state.ky[n]    = static_cast<amrex::Real>(k_dist(state.rng));
            state.phi_u[n] = static_cast<amrex::Real>(phi_dist(state.rng));
            state.phi_v[n] = static_cast<amrex::Real>(phi_dist(state.rng));
        }

        // Allocate and copy fixed data to device
        state.d_kx.resize(N);    state.d_ky.resize(N);
        state.d_phi_u.resize(N); state.d_phi_v.resize(N);
        state.d_A_u.resize(N);   state.d_A_v.resize(N);

        amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                         state.kx.begin(), state.kx.end(), state.d_kx.begin());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                         state.ky.begin(), state.ky.end(), state.d_ky.begin());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                         state.phi_u.begin(), state.phi_u.end(), state.d_phi_u.begin());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                         state.phi_v.begin(), state.phi_v.end(), state.d_phi_v.begin());
        // Amplitude device vectors initialised to zero via Gpu::copy below
        amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                         state.A_u.begin(), state.A_u.end(), state.d_A_u.begin());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                         state.A_v.begin(), state.A_v.end(), state.d_A_v.begin());

        amrex::Print() << "  Spectral noise: N_modes=" << N
                       << "  L_c=" << tp.L_c << " m\n";
    }
}

// ============================================================================
// apply_gaussian_spatial_kernel  (internal helper for "ou_process" + L_c > 0)
//
// In-place separable 1-D Gaussian smoothing on a 2-component MultiFab,
// applied as two sequential passes (x then y).  The z-dimension is traversed
// normally; each (i,j) slice at every k receives independent smoothing in
// the horizontal plane, appropriate for 2-D surface fire simulations embedded
// in a 3-D domain.
//
// Boundary treatment: Neumann clamping – the index is clamped to the valid
// box edges so no ghost cells or FillBoundary calls are required, making the
// algorithm MPI-safe (each subdomain processed independently).
//
// After both passes the field is rescaled to unit RMS so the OU sigma
// parameter retains its meaning as the stationary standard deviation.
// ============================================================================
inline void apply_gaussian_spatial_kernel(MultiFab& mf, double sigma_k)
{
    if (sigma_k < 1.0e-6) return;

    // Kernel half-radius: 3-sigma rule
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
#if (AMREX_SPACEDIM >= 3)
        const int klo = bx.smallEnd(2), khi = bx.bigEnd(2);
#else
        const int klo = 0, khi = 0;
#endif

        const int nx = ihi - ilo + 1;
        const int ny = jhi - jlo + 1;
        const int nz = khi - klo + 1;

        std::vector<double> tmp(static_cast<std::size_t>(2 * nz * ny * nx), 0.0);

        auto fn = mf.array(mfi);

        // ---- Pass 1: smooth in x ------------------------------------------
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

        // ---- Pass 2: smooth in y ------------------------------------------
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

    // ---- Rescale to unit RMS per component --------------------------------
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
// called after vel_base has been set for this timestep.
//
// Parameters:
//   vel       – output: perturbed wind MultiFab (overwritten)
//   vel_base  – input:  unperturbed background wind MultiFab (preserved)
//   ou_state  – in/out: per-cell OU state MultiFab (2 comps); nullptr means
//               domain-uniform mode (L_c = 0), direction_walk, or spectral_noise
//   state     – in/out: scalar + spectral mode state
//   dt        – timestep [s]
//   tp        – turbulent wind parameters
//   geom      – domain geometry (cell sizes and problo)
// ============================================================================
inline void apply_turb_wind(MultiFab&         vel,
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
            // ---- Domain-uniform OU (L_c = 0) --------------------------------
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
            // ---- Per-cell spatially-correlated OU (L_c > 0) -----------------
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

            const double sigma_k = tp.L_c / static_cast<double>(geom.CellSize(0));
            apply_gaussian_spatial_kernel(white_noise, sigma_k);

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

    } else if (tp.model == "spectral_noise") {
        // ---- Random Fourier Feature spectral noise with OU amplitudes -------
        //
        // Step 1: update scalar OU amplitude per mode on CPU.
        //   Internal variance = 1; user's sigma is applied as a scale factor
        //   in the GPU reconstruction kernel.
        const double alpha      = std::exp(-tp.theta * static_cast<double>(dt));
        const double sigma_step = std::sqrt(1.0 - alpha * alpha);  // unit internal sigma
        const int N = state.N_modes;

        for (int n = 0; n < N; ++n) {
            state.A_u[n] = static_cast<amrex::Real>(
                alpha * static_cast<double>(state.A_u[n])
                + sigma_step * state.norm_dist(state.rng));
            state.A_v[n] = static_cast<amrex::Real>(
                alpha * static_cast<double>(state.A_v[n])
                + sigma_step * state.norm_dist(state.rng));
        }

        // Step 2: copy updated amplitudes to device
        amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                         state.A_u.begin(), state.A_u.end(), state.d_A_u.begin());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                         state.A_v.begin(), state.A_v.end(), state.d_A_v.begin());

        // Step 3: reconstruct field on GPU
        //   u'(x,y) = sigma * sqrt(2/N) * Σ_n A_u_n * cos(kx_n*x + ky_n*y + phi_u_n)
        const auto dx     = geom.CellSizeArray();
        const auto problo = geom.ProbLoArray();
        const Real scale  = static_cast<Real>(
            tp.sigma * std::sqrt(2.0 / static_cast<double>(N)));

        const Real* kx_ptr    = state.d_kx.data();
        const Real* ky_ptr    = state.d_ky.data();
        const Real* phi_u_ptr = state.d_phi_u.data();
        const Real* phi_v_ptr = state.d_phi_v.data();
        const Real* A_u_ptr   = state.d_A_u.data();
        const Real* A_v_ptr   = state.d_A_v.data();

        for (MFIter mfi(vel); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            auto const vb = vel_base.const_array(mfi);
            auto       v  = vel.array(mfi);

            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                const Real x = problo[0] + (i + Real(0.5)) * dx[0];
                const Real y = problo[1] + (j + Real(0.5)) * dx[1];
                Real up = Real(0.0);
                Real vp = Real(0.0);
                for (int n = 0; n < N; ++n) {
                    up += A_u_ptr[n] * std::cos(kx_ptr[n] * x + ky_ptr[n] * y + phi_u_ptr[n]);
                    vp += A_v_ptr[n] * std::cos(kx_ptr[n] * x + ky_ptr[n] * y + phi_v_ptr[n]);
                }
                v(i, j, k, 0) = vb(i, j, k, 0) + scale * up;
                v(i, j, k, 1) = vb(i, j, k, 1) + scale * vp;
                v(i, j, k, 2) = vb(i, j, k, 2);
            });
        }

    } else if (tp.model == "direction_walk") {
        // ---- Bounded cumulative random walk of wind direction ---------------
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
                const Real ub  = vb(i, j, k, 0);
                const Real vb_y = vb(i, j, k, 1);
                v(i, j, k, 0) =  ub * cos_a - vb_y * sin_a;
                v(i, j, k, 1) =  ub * sin_a + vb_y * cos_a;
                v(i, j, k, 2) = vb(i, j, k, 2);
            });
        }
        amrex::Print() << "  Direction-walk angle: "
                       << state.angle_cumul * (180.0 / M_PI) << " deg\n";
    }
}

#endif // TURB_WIND_H

