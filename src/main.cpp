#include <AMReX.H>
#include <AMReX_Array4.H>
#include <AMReX_Gpu.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelFor.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_BoxArray.H>

#include <cmath>
#include <algorithm>
#include <string>

using namespace amrex;

// ============================================================================
// Helpers
// ============================================================================

static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real safe_norm2(Real ax, Real ay, Real eps = 1e-12)
{
    return std::sqrt(ax*ax + ay*ay + eps);
}

static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool in_front_band(Real phi0, Real dx, Real dy)
{
    // Small band around zero level-set where we can apply a "front" viscosity if desired.
    Real band = 3.0 * std::max(dx, dy);
    return std::abs(phi0) <= band;
}

// ============================================================================
// Initialize signed-distance phi of a sphere
// phi(x) = distance_to_center - radius
// ============================================================================
static void init_phi_sphere(MultiFab& phi, const Geometry& geom,
                            Real cx, Real cy, Real cz, Real radius)
{
    const auto dx  = geom.CellSize();
    const auto plo = geom.ProbLo();

    for (MFIter mfi(phi, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto p = phi.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            Real x = plo[0] + (i + 0.5_rt) * dx[0];
            Real y = plo[1] + (j + 0.5_rt) * dx[1];
#if (AMREX_SPACEDIM==3)
            Real z = plo[2] + (k + 0.5_rt) * dx[2];
#else
            Real z = cz;
#endif
            Real dist = std::sqrt( (x-cx)*(x-cx) + (y-cy)*(y-cy) + (z-cz)*(z-cz) );
            p(i,j,k) = dist - radius;
        });
    }
}

// ============================================================================
// Initialize velocity (constant); keep vel layout as-is: components (ux,uy,uz)
// ============================================================================
static void init_velocity_constant(MultiFab& vel, Real ux, Real uy, Real uz)
{
    for (MFIter mfi(vel, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto v = vel.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            v(i,j,k,0) = ux;
            v(i,j,k,1) = uy;
#if (AMREX_SPACEDIM==3)
            v(i,j,k,2) = uz;
#else
            v(i,j,k,2) = 0.0_rt;
#endif
        });
    }
}

// ============================================================================
// Compute CFL dt from velocity magnitude
// dt = cfl * min(dx) / max(|u|)
// ============================================================================
static Real compute_dt(const MultiFab& vel, const Geometry& geom, Real cfl)
{
    Real umax = 1e-12;
    for (MFIter mfi(vel, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto v = vel.array(mfi);

        // Reduce on-GPU
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            Real ux = v(i,j,k,0);
            Real uy = v(i,j,k,1);
#if (AMREX_SPACEDIM==3)
            Real uz = v(i,j,k,2);
#else
            Real uz = 0.0_rt;
#endif
            Real um = std::sqrt(ux*ux + uy*uy + uz*uz);
            // atomic max
            amrex::Gpu::AtomicMax(&umax, um);
        });
    }
    // Copy back umax from device
    Real umax_h = umax;
#ifdef AMREX_USE_GPU
    amrex::Gpu::streamSynchronize();
#endif
    const auto dx = geom.CellSize();
    Real mindx = dx[0];
    for (int d=1; d<AMREX_SPACEDIM; ++d) mindx = std::min(mindx, dx[d]);

    if (umax_h < 1e-12) {
        // stationary flow; choose conservative dt
        return cfl * mindx;
    } else {
        return cfl * mindx / umax_h;
    }
}

// ============================================================================
// Build RHS per Eq. (2) in the paper:
// ∂φ/∂t = - R |∇φ| + R ε (Δx φ_xx + Δy φ_yy)
// Here R is the normal-front speed. Placeholder: R = max(0, u · n).
// You can later replace R with Rothermel ROS projected onto n.
// [1](https://nrel-my.sharepoint.com/personal/hgopalan_nrel_gov/Documents/Microsoft%20Copilot%20Chat%20Files/gmd-2024-124.pdf)
// ============================================================================
static void build_rhs_levelset_eq2(const MultiFab& phi,
                                   const MultiFab& vel,
                                   MultiFab& rhs,
                                   const Geometry& geom,
                                   Real epsilon_front = 0.4,   // near front
                                   Real epsilon_far   = 0.4)   // elsewhere
{
    const auto dx = geom.CellSize();
    GpuArray<Real,AMREX_SPACEDIM> gdx{dx[0], dx[1], (AMREX_SPACEDIM>2 ? dx[2] : 1.0_rt)};

    for (MFIter mfi(phi, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx  = mfi.validbox();
        auto const p   = phi.const_array(mfi);
        auto const v   = vel.const_array(mfi);
        auto const fr  = rhs.array(mfi);

        // need interior points for derivatives; shrink by one
        const Box ibx = amrex::grow(bx, -1);

        amrex::ParallelFor(ibx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            // First derivatives (central)
            Real dphidx = (p(i+1,j,k) - p(i-1,j,k)) / (2.0_rt * gdx[0]);
            Real dphidy = (p(i,j+1,k) - p(i,j-1,k)) / (2.0_rt * gdx[1]);
            Real gradmag = safe_norm2(dphidx, dphidy);

            // Unit normal to the front (horizontal)
            Real nx = dphidx / gradmag;
            Real ny = dphidy / gradmag;

            // Second derivatives
            Real d2phidx2 = (p(i+1,j,k) - 2.0_rt*p(i,j,k) + p(i-1,j,k)) / (gdx[0]*gdx[0]);
            Real d2phidy2 = (p(i,j+1,k) - 2.0_rt*p(i,j,k) + p(i,j-1,k)) / (gdx[1]*gdx[1]);

            // Paper's scaled Laplacian: Δφ = Δx φ_xx + Δy φ_yy
            Real lap_scaled = gdx[0]*d2phidx2 + gdx[1]*d2phidy2;

            // Local artificial viscosity (front vs. far)
            Real eps_loc = in_front_band(p(i,j,k), gdx[0], gdx[1]) ? epsilon_front : epsilon_far;

            // Velocity kept as-is: components (ux, uy, uz)
            Real ux = v(i,j,k,0);
            Real uy = v(i,j,k,1);

            // Normal-front speed placeholder: project ambient wind onto front normal
            Real un = ux*nx + uy*ny;
            Real R  = amrex::max(0.0_rt, un);

            // RHS of Eq. (2)
            fr(i,j,k) = - R * gradmag + R * eps_loc * lap_scaled;
        });

        // simple boundary: zero RHS on 1-cell rim
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            if (i==bx.smallEnd(0) || i==bx.bigEnd(0) ||
                j==bx.smallEnd(1) || j==bx.bigEnd(1)) {
                fr(i,j,k) = 0.0_rt;
            }
        });
    }
}

// ============================================================================
// RK3 (Shu–Osher) advect step using the Eq. (2) RHS
// ============================================================================
static void advect_levelset_rk3(MultiFab& phi,
                                const MultiFab& vel,
                                const Geometry& geom,
                                Real dt,
                                Real epsilon_front = 0.4,
                                Real epsilon_far   = 0.4)
{
    const auto& ba = phi.boxArray();
    const auto& dm = phi.DistributionMap();

    MultiFab rhs(ba, dm, 1, 1);
    MultiFab phi1(ba, dm, 1, 1);
    MultiFab phi2(ba, dm, 1, 1);
    MultiFab phi_new(ba, dm, 1, 1);

    // ---- Stage 1 ----
    build_rhs_levelset_eq2(phi, vel, rhs, geom, epsilon_front, epsilon_far);

    for (MFIter mfi(phi, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto p  = phi.array(mfi);
        auto r  = rhs.array(mfi);
        auto p1 = phi1.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            p1(i,j,k) = p(i,j,k) + dt * r(i,j,k);
        });
    }

    // ---- Stage 2 ----
    build_rhs_levelset_eq2(phi1, vel, rhs, geom, epsilon_front, epsilon_far);

    for (MFIter mfi(phi, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto p   = phi.array(mfi);
        auto p1  = phi1.array(mfi);
        auto r   = rhs.array(mfi);
        auto p2  = phi2.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            Real tmp = p1(i,j,k) + dt * r(i,j,k);
            p2(i,j,k) = 0.75_rt * p(i,j,k) + 0.25_rt * tmp;
        });
    }

    // ---- Stage 3 ----
    build_rhs_levelset_eq2(phi2, vel, rhs, geom, epsilon_front, epsilon_far);

    for (MFIter mfi(phi, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto p   = phi.array(mfi);
        auto p2  = phi2.array(mfi);
        auto r   = rhs.array(mfi);
        auto pn  = phi_new.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            Real tmp = p2(i,j,k) + dt * r(i,j,k);
            pn(i,j,k) = (1.0_rt/3.0_rt) * p(i,j,k) + (2.0_rt/3.0_rt) * tmp;
        });
    }

    // overwrite phi with phi_new
    for (MFIter mfi(phi, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto p  = phi.array(mfi);
        auto pn = phi_new.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            p(i,j,k) = pn(i,j,k);
        });
    }
}

// ============================================================================
// (Optional) simple reinitialization stub -- keep φ as signed-distance-ish.
// The paper solves reinit Eq. (3) each step; you can wire that in similarly.
// [1](https://nrel-my.sharepoint.com/personal/hgopalan_nrel_gov/Documents/Microsoft%20Copilot%20Chat%20Files/gmd-2024-124.pdf)
// ============================================================================
static void reinitialize_phi(MultiFab& /*phi*/,
                             const Geometry& /*geom*/,
                             int /*n_iters*/,
                             Real /*dtau*/)
{
    // Stub: no-op by default. Implement Sussman reinitialization if desired.  [1](https://nrel-my.sharepoint.com/personal/hgopalan_nrel_gov/Documents/Microsoft%20Copilot%20Chat%20Files/gmd-2024-124.pdf)
}

// ============================================================================
// Main
// ============================================================================
int main(int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    {
        // ---------------- Parameters ----------------
        ParmParse pp;

        // domain cells
        int nx = 128, ny = 128, nz = (AMREX_SPACEDIM==3 ? 64 : 1);
        pp.query("nx", nx);
        pp.query("ny", ny);
        pp.query("nz", nz);

        // physical box
        Real prob_lo[3] = {0.0_rt, 0.0_rt, 0.0_rt};
        Real prob_hi[3] = {1.0_rt, 1.0_rt, (AMREX_SPACEDIM==3 ? 1.0_rt : 0.0_rt)};
        pp.queryarr("prob_lo", prob_lo, 0, 3);
        pp.queryarr("prob_hi", prob_hi, 0, 3);

        // sphere init
        Real cx = 0.5_rt, cy = 0.5_rt, cz = (AMREX_SPACEDIM==3 ? 0.5_rt : 0.0_rt);
        Real radius = 0.15_rt;
        pp.query("cx", cx);
        pp.query("cy", cy);
        pp.query("cz", cz);
        pp.query("radius", radius);

        // velocity
        Real ux = 1.0_rt, uy = 0.0_rt, uz = 0.0_rt;
        pp.query("ux", ux);
        pp.query("uy", uy);
        pp.query("uz", uz);

        // time integration
        int  nsteps = 100;
        Real cfl    = 0.5_rt;
        Real dt     = -1.0_rt; // if <0, compute from CFL
        pp.query("nsteps", nsteps);
        pp.query("cfl", cfl);
        pp.query("dt", dt);

        // artificial viscosity values (front/far)
        Real epsilon_front = 0.4_rt;
        Real epsilon_far   = 0.4_rt;
        pp.query("epsilon_front", epsilon_front);
        pp.query("epsilon_far", epsilon_far);

        int plot_int = 50;
        pp.query("plot_int", plot_int);

        // ---------------- Geometry & grids ----------------
        IntVect dom_lo(0,0,0);
        IntVect dom_hi(nx-1, ny-1, (AMREX_SPACEDIM==3 ? nz-1 : 0));
        Box domain(dom_lo, dom_hi);

        RealBox real_box({prob_lo[0], prob_lo[1], prob_lo[2]},
                         {prob_hi[0], prob_hi[1], prob_hi[2]});

#ifdef AMREX_USE_GPU
        // prefer tiling on GPU
        amrex::Gpu::setLaunchRegion(1);
#endif

        int coord_sys = 0; // Cartesian
        Vector<int> is_periodic(AMREX_SPACEDIM, 0); // non-periodic by default
        Geometry geom(domain, &real_box, coord_sys, is_periodic.data());

        BoxArray ba(domain);
        int max_grid_size = 64;
        pp.query("max_grid_size", max_grid_size);
        ba.maxSize(max_grid_size);

        DistributionMapping dm(ba);

        // ---------------- State fields ----------------
        MultiFab phi(ba, dm, 1, 1);           // level set
        MultiFab vel(ba, dm, 3, 0);           // velocity (ux,uy,uz)

        // ---------------- Initialization ----------------
        init_phi_sphere(phi, geom, cx, cy, cz, radius);
        init_velocity_constant(vel, ux, uy, uz);

        // Compute dt if unspecified
        if (dt <= 0.0_rt) {
            dt = compute_dt(vel, geom, cfl);
        }

        // ---------------- Time loop ----------------
        for (int step = 0; step < nsteps; ++step)
        {
            advect_levelset_rk3(phi, vel, geom, dt, epsilon_front, epsilon_far);

            // Optional: reinitialize to maintain signed-distance properties
            // reinitialize_phi(phi, geom, /*n_iters=*/1, /*dtau=*/0.5_rt);

            if (plot_int > 0 && (step % plot_int == 0))
            {
                std::string plt = amrex::Concatenate("plt", step, 5);
                Vector<std::string> varnames = {"phi"};
                amrex::WriteSingleLevelPlotfile(plt, phi, varnames, geom, step, 0.0_rt);
            }
        }

        // final write
        {
            std::string plt = amrex::Concatenate("plt", nsteps, 5);
            Vector<std::string> varnames = {"phi"};
            amrex::WriteSingleLevelPlotfile(plt, phi, varnames, geom, nsteps, nsteps*dt);
        }
    }
    amrex::Finalize();
    return 0;
}