// ============================================================================
// Flux-Based Ember Cascade Model (plume-rise driven)
// ============================================================================
//
// Models dense firebrand showers produced by rising smoke plumes as a
// continuous ember-landing-flux field rather than tracking individual
// firebrands (as the Albini 1983 trajectory model does).  The approach is
// appropriate when fire intensity is high enough to produce large ember
// cascades from the convective plume, where the single-brand Monte-Carlo
// treatment becomes computationally or physically inadequate.
//
// Physical model
// ──────────────
// For every active fire-front cell (phi ≤ phi_thresh, I_B > I_B_min):
//
//  1. Plume lofting height (Albini 1983):
//       H_z = 12.2 × I_B^(1/3)   [I_B in kW/m, H_z in m]
//
//  2. Ember source flux per unit fire-line length [embers/m/s]:
//       F_src = k_flux × (I_B / I_B_ref)^flux_exp
//
//  3. Flight time and mean wind-advection displacement:
//       t_f = H_z / terminal_velocity  [s]
//       dx  = u × t_f,   dy = v × t_f  [m]
//     where (u, v) is the horizontal wind at the source cell.
//     When use_3d_wind = 1 (or require_3d_wind = 1) the height-averaged wind
//     from the 3-D massconsistent_amr plt file is used instead of the 2-D
//     ambient velocity MultiFab.
//
//  4. Gaussian spread radius (turbulent diffusion from plume):
//       sigma = sigma_base + k_sigma × H_z  [m]
//     This captures both near-source turbulent diffusion and the effective
//     widening of the ember cloud with plume rise height.
//
//  5. Landing flux density contribution at cell (x_land, y_land) [embers/m²/s]:
//       F_land(x, y) += F_src × (dx_cell / (2π σ²))
//                       × exp(-(Δx² + Δy²) / (2 σ²))
//     where Δx = x_land − x_src − dx,  Δy = y_land − y_src − dy.
//     Only cells within n_sigma_cutoff × sigma of the mean landing position
//     receive contributions (Gaussian truncation for efficiency).
//
//  6. Stochastic Poisson ignition per check interval:
//       λ = F_land × dt_check / N_min_density
//       P_ign = 1 − exp(−λ)   ≈ min(1, λ)  for small λ
//     If a uniform random number on [0,1] is less than P_ign, a spot fire
//     of radius spot_radius is placed at that cell.
//
// 3-D wind assertion
// ──────────────────
// When require_3d_wind = 1 the model asserts that a valid PltWindData struct
// (from an AMReX massconsistent_amr plt file) is supplied.  If the plt file
// is missing or could not be read, the run aborts with a diagnostic message.
// When use_3d_wind = 1 but require_3d_wind = 0, the model falls back to the
// 2-D ambient velocity MultiFab if no 3-D wind data are available.
//
// Diagnostic MultiFab (ember_cascade_mf, 2 components)
// ─────────────────────────────────────────────────────
//   comp 0: accumulated ember landing flux density [embers/m²/s]
//   comp 1: ignition flag (1.0 = spot fire placed this check, 0.0 otherwise)
//
// References
// ──────────
//   Albini, F.A. (1983). "Potential spotting distance from wind-driven
//       surface fires." USDA Forest Service Research Paper INT-309.
//   Sardoy, N., Consalvi, J.L., Porterie, B., & Fernandez-Pello, A.C. (2007).
//       "Modeling transport and combustion of firebrands from burning trees."
//       Combustion and Flame 150(3):151–169.
//   Ellis, P.F.M. (2000). "The aerodynamic and combustion characteristics of
//       eucalypt bark – a firebrand study." PhD thesis, ANU.
//
// Input prefix: "ember_cascade."
// ============================================================================

#ifndef EMBER_CASCADE_FLUX_H
#define EMBER_CASCADE_FLUX_H

#include "constants.H"
#include "parse_inputs.H"
#include "plt_wind_reader.H"

#include <AMReX.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>

#include <cmath>
#include <random>
#include <vector>

using namespace amrex;

// ============================================================================
// compute_ember_cascade_flux
//
// Parameters:
//   phi                – level-set field (modified in-place: ignited cells
//                        are set to −spot_radius to start a new fire circle)
//   ember_cascade_mf   – output diagnostic MultiFab (2 comps; zeroed each call)
//                          comp 0: ember landing flux [embers/m²/s]
//                          comp 1: ignition flag (1=new spot fire placed)
//   fireline_intensity_mf – Byram fireline intensity [kW/m] at each cell
//   vel                – 2-D ambient wind velocity MultiFab (u=comp0, v=comp1)
//   geom               – AMReX geometry
//   ec                 – EmberCascadeParams from InputParameters
//   dt_check           – effective time interval for this check [s]
//                        (= step_dt × check_interval, used in Poisson λ)
//   step               – current timestep (used to seed RNG)
//   plt_wind           – optional 3-D wind from massconsistent_amr (may be
//                        nullptr; required when ec.require_3d_wind == 1)
// ============================================================================
inline void compute_ember_cascade_flux(
    MultiFab&                                          phi,
    MultiFab&                                          ember_cascade_mf,
    const MultiFab&                                    fireline_intensity_mf,
    const MultiFab&                                    vel,
    const Geometry&                                    geom,
    const InputParameters::EmberCascadeParams&         ec,
    Real                                               dt_check,
    int                                                step,
    const PltWindData*                                 plt_wind = nullptr)
{
    if (ec.enable == 0) return;

    // ---- 3-D wind assertion ----
    // When require_3d_wind = 1 the caller must supply a valid PltWindData.
    // Abort early with a clear message rather than silently falling back.
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        !(ec.require_3d_wind == 1 && (plt_wind == nullptr || !plt_wind->valid)),
        "ember_cascade.require_3d_wind = 1 but no valid 3-D plt wind data were "
        "provided.  Set ember_cascade.plt_wind_file to a massconsistent_amr "
        "plotfile directory.");

    const bool use_3d = (ec.use_3d_wind == 1 || ec.require_3d_wind == 1)
                        && (plt_wind != nullptr && plt_wind->valid);

    // ---- Host copies (ember cascade integration is serial on the CPU) ----
    amrex::Gpu::streamSynchronize();

    amrex::Arena* pinned = amrex::The_Pinned_Arena();
    amrex::MFInfo host_info;
    host_info.SetArena(pinned);

    amrex::MultiFab phi_h(phi.boxArray(), phi.DistributionMap(),
                          phi.nComp(), phi.nGrow(), host_info);
    amrex::MultiFab::Copy(phi_h, phi, 0, 0, phi.nComp(), phi.nGrow());

    amrex::MultiFab fi_h(fireline_intensity_mf.boxArray(),
                         fireline_intensity_mf.DistributionMap(),
                         1, 0, host_info);
    amrex::MultiFab::Copy(fi_h, fireline_intensity_mf, 0, 0, 1, 0);

    amrex::MultiFab vel_h(vel.boxArray(), vel.DistributionMap(),
                          vel.nComp(), vel.nGrow(), host_info);
    amrex::MultiFab::Copy(vel_h, vel, 0, 0, vel.nComp(), vel.nGrow());

    amrex::Gpu::streamSynchronize();

    const auto  dx       = geom.CellSize();
    const Real* prob_lo  = geom.ProbLo();
    const int   nx       = geom.Domain().length(0);
    const int   ny       = geom.Domain().length(1);
    const Real  cell_area = dx[0] * dx[1];

    // ---- Build flat wind arrays for transport ----
    std::vector<Real> u_arr(nx * ny, Real(0.0));
    std::vector<Real> v_arr(nx * ny, Real(0.0));

    if (use_3d) {
        // Bilinearly interpolate height-averaged wind from 3-D plt grid onto
        // the fire-model 2-D cell grid (same logic as in albini_spotting.H).
        const int   pw_nx = plt_wind->nx_2d;
        const int   pw_ny = plt_wind->ny_2d;
        const Real  pw_x0 = plt_wind->lo2d[0];
        const Real  pw_y0 = plt_wind->lo2d[1];
        const Real  pw_x1 = plt_wind->hi2d[0];
        const Real  pw_y1 = plt_wind->hi2d[1];
        const Real  pw_dx = (pw_nx > 0) ? (pw_x1 - pw_x0) / pw_nx : Real(1.0);
        const Real  pw_dy = (pw_ny > 0) ? (pw_y1 - pw_y0) / pw_ny : Real(1.0);

        const auto& u2d_h = plt_wind->u2d;
        const auto& v2d_h = plt_wind->v2d;

        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                Real xc = prob_lo[0] + (i + Real(0.5)) * dx[0];
                Real yc = prob_lo[1] + (j + Real(0.5)) * dx[1];

                Real xi = (xc - pw_x0) / pw_dx - Real(0.5);
                Real yj = (yc - pw_y0) / pw_dy - Real(0.5);

                int pi  = std::max(0, std::min(static_cast<int>(std::floor(xi)), pw_nx - 2));
                int pj  = std::max(0, std::min(static_cast<int>(std::floor(yj)), pw_ny - 2));
                Real fx = std::max(Real(0.0), std::min(xi - static_cast<Real>(pi), Real(1.0)));
                Real fy = std::max(Real(0.0), std::min(yj - static_cast<Real>(pj), Real(1.0)));

                auto pidx = [pw_nx](int ii, int jj){ return jj * pw_nx + ii; };
                u_arr[j * nx + i] =
                    (Real(1.0)-fx)*(Real(1.0)-fy)*u2d_h[pidx(pi,   pj  )]
                    + fx           *(Real(1.0)-fy)*u2d_h[pidx(pi+1, pj  )]
                    + (Real(1.0)-fx)*fy           *u2d_h[pidx(pi,   pj+1)]
                    + fx            *fy           *u2d_h[pidx(pi+1, pj+1)];
                v_arr[j * nx + i] =
                    (Real(1.0)-fx)*(Real(1.0)-fy)*v2d_h[pidx(pi,   pj  )]
                    + fx           *(Real(1.0)-fy)*v2d_h[pidx(pi+1, pj  )]
                    + (Real(1.0)-fx)*fy           *v2d_h[pidx(pi,   pj+1)]
                    + fx            *fy           *v2d_h[pidx(pi+1, pj+1)];
            }
        }
    } else {
        // Use the 2-D ambient velocity MultiFab
        for (MFIter mfi(vel_h); mfi.isValid(); ++mfi) {
            const Box& bx    = mfi.validbox();
            auto const v_arr4 = vel_h.const_array(mfi);
            const auto lo    = lbound(bx);
            const auto hi    = ubound(bx);
            for (int j = lo.y; j <= hi.y; ++j) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    u_arr[j * nx + i] = v_arr4(i, j, 0, 0);
                    v_arr[j * nx + i] = v_arr4(i, j, 0, 1);
                }
            }
        }
#ifdef AMREX_USE_MPI
        MPI_Allreduce(MPI_IN_PLACE, u_arr.data(), nx * ny,
                      ParallelDescriptor::Mpi_typemap<Real>::type(),
                      MPI_SUM, ParallelDescriptor::Communicator());
        MPI_Allreduce(MPI_IN_PLACE, v_arr.data(), nx * ny,
                      ParallelDescriptor::Mpi_typemap<Real>::type(),
                      MPI_SUM, ParallelDescriptor::Communicator());
#endif
    }

    // ---- Accumulate landing-flux density field ----
    // Flat host array: flux_field[j * nx + i] = F_land [embers/m²/s]
    std::vector<Real> flux_field(nx * ny, Real(0.0));

    const Real phi_thresh    = Real(0.1);
    const Real I_B_min       = Real(ec.I_B_min);
    const Real v_t           = Real(ec.terminal_velocity);
    const Real k_flux        = Real(ec.k_flux);
    const Real I_B_ref       = Real(ec.I_B_ref);
    const Real flux_exp      = Real(ec.flux_exp);
    const Real sigma_base    = Real(ec.sigma_base);
    const Real k_sigma       = Real(ec.k_sigma);
    const Real n_sigma_cut   = Real(ec.n_sigma_cutoff);
    const Real two_pi        = Real(2.0) * Real(WildfireConst::PI);

    for (MFIter mfi(phi_h); mfi.isValid(); ++mfi) {
        const Box& bx    = mfi.validbox();
        auto const p_arr = phi_h.const_array(mfi);
        auto const fi_arr= fi_h.const_array(mfi);
        const auto lo    = lbound(bx);
        const auto hi    = ubound(bx);

        for (int js = lo.y; js <= hi.y; ++js) {
            for (int is = lo.x; is <= hi.x; ++is) {

                // Only fire-front cells (phi near zero) with sufficient intensity
                if (p_arr(is, js, 0) > phi_thresh) continue;
                const Real I_B = fi_arr(is, js, 0);
                if (I_B < I_B_min) continue;

                // Source physical position (cell centre)
                const Real xs = prob_lo[0] + (is + Real(0.5)) * dx[0];
                const Real ys = prob_lo[1] + (js + Real(0.5)) * dx[1];

                // Plume lofting height (Albini 1983)
                const Real H_z = Real(12.2) * std::pow(I_B, Real(1.0) / Real(3.0));

                // Ember source flux [embers/m/s]
                const Real F_src = k_flux
                    * std::pow(I_B / I_B_ref, flux_exp);

                // Flight time and mean transport displacement
                const Real t_f  = H_z / v_t;
                const Real u_s  = u_arr[js * nx + is];
                const Real v_s  = v_arr[js * nx + is];
                const Real dx_mean = u_s * t_f;
                const Real dy_mean = v_s * t_f;

                // Gaussian spread (turbulent diffusion + plume width)
                const Real sigma = sigma_base + k_sigma * H_z;
                const Real inv_2s2 = Real(1.0) / (Real(2.0) * sigma * sigma);
                const Real gauss_norm = F_src / (two_pi * sigma * sigma);

                // Mean landing cell centre
                const Real xc_mean = xs + dx_mean;
                const Real yc_mean = ys + dy_mean;

                // Bounding box in cell indices for the Gaussian truncation
                const Real cut_dist = n_sigma_cut * sigma;
                const int i0 = std::max(0, static_cast<int>(
                    std::floor((xc_mean - cut_dist - prob_lo[0]) / dx[0])));
                const int i1 = std::min(nx - 1, static_cast<int>(
                    std::floor((xc_mean + cut_dist - prob_lo[0]) / dx[0])));
                const int j0 = std::max(0, static_cast<int>(
                    std::floor((yc_mean - cut_dist - prob_lo[1]) / dx[1])));
                const int j1 = std::min(ny - 1, static_cast<int>(
                    std::floor((yc_mean + cut_dist - prob_lo[1]) / dx[1])));

                // Skip if mean landing is entirely outside the domain
                if (i0 > i1 || j0 > j1) continue;

                // Accumulate Gaussian contributions
                for (int jl = j0; jl <= j1; ++jl) {
                    const Real yc_l = prob_lo[1] + (jl + Real(0.5)) * dx[1];
                    const Real dy2  = yc_l - yc_mean;
                    for (int il = i0; il <= i1; ++il) {
                        const Real xc_l = prob_lo[0] + (il + Real(0.5)) * dx[0];
                        const Real dx2  = xc_l - xc_mean;
                        const Real r2   = dx2 * dx2 + dy2 * dy2;
                        // gauss_norm = F_src / (2π σ²) is the Gaussian flux density
                        // [embers/m²/s] at the landing point.  Multiply by cell_area
                        // [m²] to obtain the total landing rate [embers/s] in this cell.
                        flux_field[jl * nx + il] +=
                            gauss_norm * std::exp(-r2 * inv_2s2) * cell_area;
                    }
                }
            }
        }
    }

#ifdef AMREX_USE_MPI
    // Reduce flux field across all MPI ranks so every rank has the complete
    // accumulated landing flux before applying ignitions.
    MPI_Allreduce(MPI_IN_PLACE, flux_field.data(), nx * ny,
                  ParallelDescriptor::Mpi_typemap<Real>::type(),
                  MPI_SUM, ParallelDescriptor::Communicator());
#endif

    // ---- Write landing flux to ember_cascade_mf comp 0 ----
    // Zero both components first; comp 0 gets the flux field, comp 1 the
    // per-cell ignition flag (written below).
    ember_cascade_mf.setVal(Real(0.0));

    amrex::MultiFab ecmf_h(ember_cascade_mf.boxArray(),
                           ember_cascade_mf.DistributionMap(),
                           ember_cascade_mf.nComp(), 0, host_info);
    ecmf_h.setVal(Real(0.0));

    for (MFIter mfi(ecmf_h); mfi.isValid(); ++mfi) {
        const Box& bx  = mfi.validbox();
        auto ec_arr    = ecmf_h.array(mfi);
        const auto lo  = lbound(bx);
        const auto hi  = ubound(bx);
        for (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                ec_arr(i, j, 0, 0) = flux_field[j * nx + i];
            }
        }
    }

    // ---- Stochastic Poisson ignition ----
    // RNG seed: identical on all ranks for deterministic results.
    unsigned rng_seed;
    if (ec.random_seed == 0) {
#ifdef AMREX_USE_MPI
        if (ParallelDescriptor::IOProcessor()) rng_seed = std::random_device{}();
        MPI_Bcast(&rng_seed, 1,
                  ParallelDescriptor::Mpi_typemap<unsigned>::type(),
                  ParallelDescriptor::IOProcessorNumber(),
                  ParallelDescriptor::Communicator());
#else
        rng_seed = std::random_device{}();
#endif
    } else {
        rng_seed = static_cast<unsigned>(ec.random_seed + step);
    }
    std::mt19937 gen(rng_seed);
    std::uniform_real_distribution<Real> uniform(Real(0.0), Real(1.0));

    const Real N_min = Real(ec.N_min_density);
    const Real spot_r = Real(ec.spot_radius);
    const int radius_cells = static_cast<int>(std::ceil(spot_r / dx[0])) + 1;

    int n_ignitions = 0;

    for (MFIter mfi(phi_h); mfi.isValid(); ++mfi) {
        const Box& bx   = mfi.validbox();
        auto p_arr      = phi_h.array(mfi);
        auto ec_arr     = ecmf_h.array(mfi);
        const auto lo   = lbound(bx);
        const auto hi   = ubound(bx);

        for (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {

                const Real flux = flux_field[j * nx + i];
                if (flux <= Real(0.0)) continue;

                // Poisson arrival probability
                const Real lambda = flux * dt_check / N_min;
                const Real P_ign  = Real(1.0) - std::exp(-lambda);

                if (uniform(gen) >= P_ign) continue;

                // Mark ignition flag
                ec_arr(i, j, 0, 1) = Real(1.0);
                ++n_ignitions;

                // Apply spot-fire circle in phi
                const Real xc_ign = prob_lo[0] + (i + Real(0.5)) * dx[0];
                const Real yc_ign = prob_lo[1] + (j + Real(0.5)) * dx[1];

                const int i0 = std::max(lo.x, i - radius_cells);
                const int i1 = std::min(hi.x, i + radius_cells);
                const int j0 = std::max(lo.y, j - radius_cells);
                const int j1 = std::min(hi.y, j + radius_cells);

                for (int jj = j0; jj <= j1; ++jj) {
                    for (int ii = i0; ii <= i1; ++ii) {
                        const Real xc2 = prob_lo[0] + (ii + Real(0.5)) * dx[0];
                        const Real yc2 = prob_lo[1] + (jj + Real(0.5)) * dx[1];
                        const Real dd  = std::sqrt((xc2 - xc_ign) * (xc2 - xc_ign)
                                                 + (yc2 - yc_ign) * (yc2 - yc_ign));
                        if (dd < spot_r) {
                            p_arr(ii, jj, 0) =
                                amrex::min(p_arr(ii, jj, 0), dd - spot_r);
                        }
                    }
                }
            }
        }
    }

    amrex::Print() << "Ember cascade flux: " << n_ignitions
                   << " spot-fire ignition(s) at step " << step
                   << (use_3d ? " [3-D wind]" : " [2-D wind]") << "\n";

    // ---- Copy results from host back to device ----
    amrex::MultiFab::Copy(phi,              phi_h,   0, 0, phi.nComp(),              phi.nGrow());
    amrex::MultiFab::Copy(ember_cascade_mf, ecmf_h,  0, 0, ember_cascade_mf.nComp(), 0);
    amrex::Gpu::streamSynchronize();
}

#endif // EMBER_CASCADE_FLUX_H
