#include <AMReX.H>
#include <AMReX_Array4.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_GPULaunch.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_BoxArray.H>
#include "numerical_schemes.H"
#include <cmath>
#include <vector>
#include <string>

using namespace amrex;

// ============================================================================
// UPDATED: Build RHS with WENO5-Z using Eq. (2) with wind magnitude R and Laplacian
// ============================================================================
static void build_rhs_weno5z(const MultiFab& phi,
                             const MultiFab& vel,
                             MultiFab& rhs,
                             const Geometry& geom)
{
    const auto dx = geom.CellSize();
    GpuArray<Real,AMREX_SPACEDIM> gdx{dx[0], dx[1], dx[2]};

    // Artificial viscosity ε per Eq. (2); use a single uniform value across the domain.
    // Paper default ~0.4 is common operationally. [1](https://nrel-my.sharepoint.com/personal/hgopalan_nrel_gov/Documents/Microsoft%20Copilot%20Chat%20Files/gmd-2024-124.pdf)
    constexpr Real eps_visc = 0.4;

    // Helper: discrete Laplacian Δφ (second-order)
    auto laplacian_phi = [=] AMREX_GPU_HOST_DEVICE (Array4<const Real> const& p,
                                                    int i, int j, int k) noexcept -> Real {
        Real ddx = (p(i+1,j,k) - 2.0*p(i,j,k) + p(i-1,j,k)) / (gdx[0]*gdx[0]);
        Real ddy = (p(i,j+1,k) - 2.0*p(i,j,k) + p(i,j-1,k)) / (gdx[1]*gdx[1]);
#if (AMREX_SPACEDIM==3)
        Real ddz = (p(i,j,k+1) - 2.0*p(i,j,k) + p(i,j,k-1)) / (gdx[2]*gdx[2]);
        return ddx + ddy + ddz;
#else
        return ddx + ddy;
#endif
    };

    // RHS per Eq. (2): ∂φ/∂t = - R ( |∇φ| - ε Δφ ), with R = ‖u‖ (wind magnitude). [1](https://nrel-my.sharepoint.com/personal/hgopalan_nrel_gov/Documents/Microsoft%20Copilot%20Chat%20Files/gmd-2024-124.pdf)
    for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
        const Box& bx  = mfi.validbox();
        auto const p   = phi.const_array(mfi);
        auto const v   = vel.const_array(mfi);
        auto       fr  = rhs.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            // |∇φ| using your existing Godunov/WENO5 routine (structure preserved)
            Real grad_norm = godunov_norm_grad_phi(p, i, j, k, gdx);

            // Laplacian Δφ
            Real lap = laplacian_phi(p, i, j, k);

            // Wind speed magnitude R
            Real ux = v(i,j,k,0);
            Real uy = v(i,j,k,1);
#if (AMREX_SPACEDIM==3)
            Real uz = v(i,j,k,2);
            Real R  = std::sqrt(ux*ux + uy*uy + uz*uz);
#else
            Real R  = std::sqrt(ux*ux + uy*uy);
#endif

            // Final RHS (keeps your WENO5 & RK3 integrators intact)
            fr(i,j,k) = - R * (grad_norm - eps_visc * lap);
        });
    }
}
// ======================= Outflow boundary fill ======================
// Fills ghost cells using linear extrapolation (outflow treatment) for phi.
// Linear extrapolation is applied unconditionally for all boundaries.
static void fill_boundary_extrap(MultiFab& mf, const Geometry& geom)
{
    // Fill interior ghost cells for MPI exchange between ranks.
    // Since is_periodic is {0,0,0}, this does not copy across periodic boundaries.
    mf.FillBoundary(geom.periodicity());

    const auto& domain = geom.Domain();

    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        const Box& gbx = mfi.fabbox(); // includes ghost cells
        auto const arr = mf.array(mfi);

        // Get domain bounds
        int dom_ilo = domain.smallEnd(0);
        int dom_ihi = domain.bigEnd(0);
        int dom_jlo = domain.smallEnd(1);
        int dom_jhi = domain.bigEnd(1);
        int dom_klo = domain.smallEnd(2);
        int dom_khi = domain.bigEnd(2);

        // Get box bounds
        int bx_ilo = bx.smallEnd(0);
        int bx_ihi = bx.bigEnd(0);
        int bx_jlo = bx.smallEnd(1);
        int bx_jhi = bx.bigEnd(1);
        int bx_klo = bx.smallEnd(2);
        int bx_khi = bx.bigEnd(2);

        int gbx_ilo = gbx.smallEnd(0);
        int gbx_ihi = gbx.bigEnd(0);
        int gbx_jlo = gbx.smallEnd(1);
        int gbx_jhi = gbx.bigEnd(1);
        int gbx_klo = gbx.smallEnd(2);
        int gbx_khi = gbx.bigEnd(2);

        int ncomp = mf.nComp();

        // X-lo boundary: linear extrapolation (outflow)
        if (bx_ilo == dom_ilo) {
            ParallelFor(Box(IntVect(gbx_ilo, gbx_jlo, gbx_klo),
                           IntVect(dom_ilo-1, gbx_jhi, gbx_khi)),
                       ncomp,
                       [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept {
                           arr(i,j,k,n) = arr(dom_ilo,j,k,n)
                               + Real(dom_ilo-i)*(arr(dom_ilo,j,k,n) - arr(dom_ilo+1,j,k,n));
                       });
        }

        // X-hi boundary: linear extrapolation (outflow)
        if (bx_ihi == dom_ihi) {
            ParallelFor(Box(IntVect(dom_ihi+1, gbx_jlo, gbx_klo),
                           IntVect(gbx_ihi, gbx_jhi, gbx_khi)),
                       ncomp,
                       [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept {
                           arr(i,j,k,n) = arr(dom_ihi,j,k,n)
                               + Real(i-dom_ihi)*(arr(dom_ihi,j,k,n) - arr(dom_ihi-1,j,k,n));
                       });
        }

        // Y-lo boundary: linear extrapolation (outflow)
        if (bx_jlo == dom_jlo) {
            ParallelFor(Box(IntVect(gbx_ilo, gbx_jlo, gbx_klo),
                           IntVect(gbx_ihi, dom_jlo-1, gbx_khi)),
                       ncomp,
                       [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept {
                           arr(i,j,k,n) = arr(i,dom_jlo,k,n)
                               + Real(dom_jlo-j)*(arr(i,dom_jlo,k,n) - arr(i,dom_jlo+1,k,n));
                       });
        }

        // Y-hi boundary: linear extrapolation (outflow)
        if (bx_jhi == dom_jhi) {
            ParallelFor(Box(IntVect(gbx_ilo, dom_jhi+1, gbx_klo),
                           IntVect(gbx_ihi, gbx_jhi, gbx_khi)),
                       ncomp,
                       [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept {
                           arr(i,j,k,n) = arr(i,dom_jhi,k,n)
                               + Real(j-dom_jhi)*(arr(i,dom_jhi,k,n) - arr(i,dom_jhi-1,k,n));
                       });
        }

        // Z-lo boundary: linear extrapolation (outflow)
        if (bx_klo == dom_klo) {
            ParallelFor(Box(IntVect(gbx_ilo, gbx_jlo, gbx_klo),
                           IntVect(gbx_ihi, gbx_jhi, dom_klo-1)),
                       ncomp,
                       [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept {
                           arr(i,j,k,n) = arr(i,j,dom_klo,n)
                               + Real(dom_klo-k)*(arr(i,j,dom_klo,n) - arr(i,j,dom_klo+1,n));
                       });
        }

        // Z-hi boundary: linear extrapolation (outflow)
        if (bx_khi == dom_khi) {
            ParallelFor(Box(IntVect(gbx_ilo, gbx_jlo, dom_khi+1),
                           IntVect(gbx_ihi, gbx_jhi, gbx_khi)),
                       ncomp,
                       [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept {
                           arr(i,j,k,n) = arr(i,j,dom_khi,n)
                               + Real(k-dom_khi)*(arr(i,j,dom_khi,n) - arr(i,j,dom_khi-1,n));
                       });
        }
    }
}

// ======================= WENO5-Z advection + RK3 ==================
static void advect_levelset_weno5z_rk3(MultiFab& phi,
                                       const MultiFab& vel,
                                       const Geometry& geom,
                                       Real dt)
{
    const auto& ba = phi.boxArray();
    const auto& dm = phi.DistributionMap();

    // Stage buffers with 3 ghosts (needed for WENO stencils)
    MultiFab phi1(ba, dm, 1, 3);
    MultiFab phi2(ba, dm, 1, 3);
    MultiFab rhs(ba, dm, 1, 0);


    // ---------- Stage 1 ----------
    build_rhs_weno5z(phi, vel, rhs, geom);
    for (MFIter mfi(phi, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& tbx = mfi.tilebox();
        auto const p  = phi.const_array(mfi);
        auto const r  = rhs.const_array(mfi);
        auto const q1 = phi1.array(mfi);
        ParallelFor(tbx, [=] AMREX_GPU_DEVICE (int i,int j,int k) noexcept {
            q1(i,j,k) = p(i,j,k) + dt * r(i,j,k);
        });
    }
    fill_boundary_extrap(phi1, geom);

    // ---------- Stage 2 ----------
    build_rhs_weno5z(phi1, vel, rhs, geom);
    for (MFIter mfi(phi, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& tbx = mfi.tilebox();
        auto const p   = phi.const_array(mfi);
        auto const q1c = phi1.const_array(mfi);
        auto const r   = rhs.const_array(mfi);
        auto const q2  = phi2.array(mfi);
        ParallelFor(tbx, [=] AMREX_GPU_DEVICE (int i,int j,int k) noexcept {
            q2(i,j,k) = 0.75 * p(i,j,k) + 0.25 * ( q1c(i,j,k) + dt * r(i,j,k) );
        });
    }
    fill_boundary_extrap(phi2, geom);

    // ---------- Stage 3 ----------
    build_rhs_weno5z(phi2, vel, rhs, geom);
    for (MFIter mfi(phi, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& tbx = mfi.tilebox();
        auto const p   = phi.const_array(mfi);
        auto const q2c = phi2.const_array(mfi);
        auto const r   = rhs.const_array(mfi);
        auto const pout= phi.array(mfi);
        ParallelFor(tbx, [=] AMREX_GPU_DEVICE (int i,int j,int k) noexcept {
           // pout(i,j,k) = (1.0/3.0) * p(i,j,k) + (2.0/3.0) * ( q2c(i,j,k) + dt * r(i,j,k) );

        Real newphi = (1.0/3.0) * p(i,j,k) + (2.0/3.0) * ( q2c(i,j,k) + dt * r(i,j,k) );

        // Clip non-positive values:
        if (newphi < 0.0) newphi = 0.0;

        pout(i,j,k) = newphi;

        });
    }
}

// ================================================================
// Smoothed sign function
// ================================================================
static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real smoothed_sign(Real phi0, Real eps)
{
    return phi0 / std::sqrt(phi0*phi0 + eps*eps);
}


// ================================================================
// Corrected reinitialization routine
// ================================================================
static void reinitialize_phi(MultiFab& phi,
                             const Geometry& geom,
                             int n_iters,
                             Real dtau)
{
    const auto dx = geom.CellSize();
    GpuArray<Real,AMREX_SPACEDIM> gdx{dx[0], dx[1], dx[2]};

    // Epsilon proportional to smallest dx
    Real eps = std::min(dx[0], std::min(dx[1], dx[2])) * 1.0;

    // phi0 must have 1 ghost cell for ∇phi operations
    MultiFab phi0(phi.boxArray(), phi.DistributionMap(), 1, 1);
    MultiFab::Copy(phi0, phi, 0, 0, 1, 1);
    fill_boundary_extrap(phi0, geom);

    // rhs needs same ghost depth as phi to allow growntilebox(1)
    MultiFab rhs(phi.boxArray(), phi.DistributionMap(), 1, 1);

    for (int it = 0; it < n_iters; ++it) {

        // Ensure ghost cells of phi are fresh each iteration using linear extrapolation
        fill_boundary_extrap(phi, geom);

        for (MFIter mfi(phi, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.growntilebox(1); // SAFE for ±1 stencil

            auto const p  = phi.const_array(mfi);
            auto const p0 = phi0.const_array(mfi);
            auto const r  = rhs.array(mfi);

            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i,int j,int k) noexcept {
                Real s = smoothed_sign(p0(i,j,k), eps);
                Real g = godunov_norm_grad_phi(p, i, j, k, gdx);
                r(i,j,k) = - s * (g - 1.0);
            });
        }

        // Update phi
        for (MFIter mfi(phi, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.growntilebox(1);

            auto const p  = phi.array(mfi);
            auto const r  = rhs.const_array(mfi);

            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i,int j,int k) noexcept {
                p(i,j,k) += dtau * r(i,j,k);
            });
        }
    }
}

// ======================= Initialize signed distance (sphere) ======================
static void init_phi_sphere(MultiFab& phi, const Geometry& geom,
                            Real cx, Real cy, Real cz, Real radius)
{
    const auto dx  = geom.CellSize();
    const auto plo = geom.ProbLo();

    for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const p = phi.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real x = plo[0] + (i + 0.5) * dx[0];
            Real y = plo[1] + (j + 0.5) * dx[1];
            Real z = plo[2] + (k + 0.5) * dx[2];
            Real d = std::sqrt((x - cx)*(x - cx) +
                               (y - cy)*(y - cy) +
                               (z - cz)*(z - cz)) - radius;
            p(i,j,k) = d; // negative inside the sphere
        });
    }
}

// ======================= Initialize velocity (constant) ==========================
static void init_velocity_constant(MultiFab& vel, Real ux, Real uy, Real uz)
{
    for (MFIter mfi(vel); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const v = vel.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            v(i,j,k,0) = ux;
            v(i,j,k,1) = uy;
            v(i,j,k,2) = uz;
        });
    }
}

// ======================= Compute dt from CFL and max|u| =========================
static Real compute_dt(const MultiFab& vel, const Geometry& geom, Real cfl)
{
    // Compute |u| into a temp field, then take infinity norm
    MultiFab umag(vel.boxArray(), vel.DistributionMap(), 1, 0);

    for (MFIter mfi(umag, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& tbx = mfi.tilebox();
        auto const v = vel.const_array(mfi);
        auto const u = umag.array(mfi);

        ParallelFor(tbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real mag2 = v(i,j,k,0)*v(i,j,k,0)
                      + v(i,j,k,1)*v(i,j,k,1)
                      + v(i,j,k,2)*v(i,j,k,2);
            u(i,j,k) = std::sqrt(mag2);
        });
    }

    Real umax = umag.norminf(0, 0, true); // comp 0, nghost 0, global

    const auto dx = geom.CellSize();
    Real dmin = dx[0];
    dmin = std::min(dmin, dx[1]);
    dmin = std::min(dmin, dx[2]);

    if (umax <= 1.0e-14) {
        return 1.0; // stationary flow: generous dt
    }
    return cfl * dmin / umax;
}

// ======================= Main ================================================
int main(int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        ParmParse pp;

        // ---------------- Inputs: grid & domain ----------------
        // You can specify either n_cell (cubic) or n_cell_x/y/z
        int n_cell      = 64;   pp.query("n_cell", n_cell);
        int n_cell_x    = n_cell; pp.query("n_cell_x", n_cell_x);
        int n_cell_y    = n_cell; pp.query("n_cell_y", n_cell_y);
        int n_cell_z    = n_cell; pp.query("n_cell_z", n_cell_z);

        int max_grid    = 32;   pp.query("max_grid_size", max_grid);

        Real plo_x = 0.0, plo_y = 0.0, plo_z = 0.0;
        Real phi_x = 1.0, phi_y = 1.0, phi_z = 1.0;
        pp.query("prob_lo_x", plo_x);
        pp.query("prob_lo_y", plo_y);
        pp.query("prob_lo_z", plo_z);
        pp.query("prob_hi_x", phi_x);
        pp.query("prob_hi_y", phi_y);
        pp.query("prob_hi_z", phi_z);

        // ---------------- Inputs: time & output ----------------
        int reinit_int = 20;
        int reinit_iters = 20;
        Real reinit_dtau = 0.5;

        pp.query("reinit_int", reinit_int);
        pp.query("reinit_iters", reinit_iters);
        pp.query("reinit_dtau", reinit_dtau);
        int nsteps    = 300;    pp.query("nsteps", nsteps);
        Real cfl      = 0.5;    pp.query("cfl", cfl);
        int plot_int  = 50;     pp.query("plot_int", plot_int);

        // ---------------- Inputs: velocity ---------------------
        Real ux = 0.25, uy = 0.0, uz = 0.0;
        pp.query("u_x", ux);
        pp.query("u_y", uy);
        pp.query("u_z", uz);


        // ---------------- Inputs: source selection ----------------
        std::string source_type = "sphere"; // "line" or "sphere"
        pp.query("source_type", source_type);

        // ---------------- Inputs: sphere -----------------------
        Real cx = 0.5, cy = 0.5, cz = 0.5, radius = 0.25;
        pp.query("sphere_center_x", cx);
        pp.query("sphere_center_y", cy);
        pp.query("sphere_center_z", cz);
        pp.query("sphere_radius",   radius);

        // ---------------- Geometry setup -----------------------
        IntVect dom_lo(0, 0, 0);
        IntVect dom_hi(n_cell_x-1, n_cell_y-1, n_cell_z-1);
        Box domain(dom_lo, dom_hi);

        RealBox rb({plo_x, plo_y, plo_z}, {phi_x, phi_y, phi_z});
        Array<int,AMREX_SPACEDIM> is_periodic{0, 0, 0};
        Geometry geom(domain, &rb, 0, is_periodic.data());

        // ---------------- Grids & distribution -----------------
        BoxArray ba(domain);
        ba.maxSize(max_grid);
        DistributionMapping dm(ba);

        // ---------------- Fields: phi (1 comp), vel (3 comps) --
        const int ng_phi = 3; // 3 ghost cells for stencil operations (WENO5-Z flux divergence uses up to ±3 cells)
        MultiFab phi(ba, dm, 1, ng_phi);
        MultiFab vel(ba, dm, 3, 1);

        // ---------------- Initialize ---------------------------
        //init_phi_sphere(phi, geom, cx, cy, cz, radius);

        // ---------------- Initialize ---------------------------
        if (source_type == "sphere") {

            // Fallback: sphere SDF (old behavior)
            init_phi_sphere(phi, geom, cx, cy, cz, radius);
        }

        init_velocity_constant(vel, ux, uy, uz);

        // ---------------- dt from CFL --------------------------
        Real dt = compute_dt(vel, geom, cfl);
        amrex::Print() << "Computed dt = " << dt << "\n";

        // ---------------- Write initial plotfile ---------------
        {
            Vector<std::string> names = {"phi"};
            WriteSingleLevelPlotfile("plt0000", phi, names, geom, 0.0, 0);
        }

        // ---------------- Time stepping ------------------------
        for (int step = 1; step <= nsteps; ++step) {
            advect_levelset_weno5z_rk3 (phi, vel, geom, dt);
            Real philomax = phi.min(0);
            Real phihimax = phi.max(0);

            amrex::Print() << "Step " << step
                        << " : phi_min = " << philomax
                        << " , phi_max = " << phihimax << "\n";

            if (reinit_int > 0 && (step % reinit_int == 0)) {
                amrex::Print() << "Reinitializing at step " << step << "\n";
                reinitialize_phi(phi, geom, reinit_iters, reinit_dtau);
            }

            if (plot_int > 0 && (step % plot_int == 0)) {
                char buf[64];
                std::snprintf(buf, sizeof(buf), "plt%04d", step);
                Vector<std::string> names = {"phi"};
                WriteSingleLevelPlotfile(buf, phi, names, geom, step*dt, step);
                amrex::Print() << "Wrote " << buf << "\n";
            }
        }

        // ---------------- Final write --------------------------
        {
            char buf[64];
            std::snprintf(buf, sizeof(buf), "plt%04d", nsteps);
            Vector<std::string> names = {"phi"};
            WriteSingleLevelPlotfile(buf, phi, names, geom, nsteps*dt, nsteps);
            amrex::Print() << "Wrote final " << buf << "\n";
        }
    }
    amrex::Finalize();
    return 0;
}