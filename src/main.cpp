#include <AMReX.H>
#include <AMReX_Array4.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_GPULaunch.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_BoxArray.H>
#include <cmath>
#include <vector>
#include <string>

using namespace amrex;

// ======================= WENO5-Z helpers ===========================
static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real weno5z_left(Real vm2, Real vm1, Real v0, Real vp1, Real vp2)
{
    // Candidate polys at i+1/2 (left-biased)
    Real p0 = ( 2.0*vm2 - 7.0*vm1 + 11.0*v0 ) / 6.0;
    Real p1 = ( -1.0*vm1 + 5.0*v0 + 2.0*vp1 ) / 6.0;
    Real p2 = ( 2.0*v0 + 5.0*vp1 - 1.0*vp2 ) / 6.0;

    // Smoothness indicators β_k
    Real b0 = (13.0/12.0) * (vm2 - 2.0*vm1 + v0) * (vm2 - 2.0*vm1 + v0)
            + 0.25 * (vm2 - 4.0*vm1 + 3.0*v0) * (vm2 - 4.0*vm1 + 3.0*v0);

    Real b1 = (13.0/12.0) * (vm1 - 2.0*v0 + vp1) * (vm1 - 2.0*v0 + vp1)
            + 0.25 * (vm1 - vp1) * (vm1 - vp1);

    Real b2 = (13.0/12.0) * (v0 - 2.0*vp1 + vp2) * (v0 - 2.0*vp1 + vp2)
            + 0.25 * (3.0*v0 - 4.0*vp1 + vp2) * (3.0*v0 - 4.0*vp1 + vp2);

    // WENO-Z: tau5 = |β0 - β2|
    Real tau5 = amrex::Math::abs(b0 - b2);
    const Real eps = 1e-12;        // small number to avoid /0
    const Real d0 = 0.1, d1 = 0.6, d2 = 0.3; // linear weights

    Real a0 = d0 * (1.0 + tau5 / (b0 + eps));
    Real a1 = d1 * (1.0 + tau5 / (b1 + eps));
    Real a2 = d2 * (1.0 + tau5 / (b2 + eps));

    Real sum = a0 + a1 + a2;
    a0 /= sum; a1 /= sum; a2 /= sum;

    return a0*p0 + a1*p1 + a2*p2;
}

// ================================================================
// Godunov |grad phi| operator
// ================================================================
static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real godunov_norm_grad_phi(const Array4<Real const>& phi,
                           int i, int j, int k,
                           const GpuArray<Real,AMREX_SPACEDIM>& dx)
{
    // Forward/backward differences
    Real px_pos = (phi(i+1,j,k) - phi(i,j,k)) / dx[0];
    Real px_neg = (phi(i,j,k) - phi(i-1,j,k)) / dx[0];

    Real py_pos = (phi(i,j+1,k) - phi(i,j,k)) / dx[1];
    Real py_neg = (phi(i,j,k) - phi(i,j-1,k)) / dx[1];

    Real pz_pos = (phi(i,j,k+1) - phi(i,j,k)) / dx[2];
    Real pz_neg = (phi(i,j,k) - phi(i,j,k-1)) / dx[2];

    // Godunov choices for reinitialization
    Real gx = (px_pos > 0.0 ? px_pos*px_pos : 0.0)
            + (px_neg < 0.0 ? px_neg*px_neg : 0.0);

    Real gy = (py_pos > 0.0 ? py_pos*py_pos : 0.0)
            + (py_neg < 0.0 ? py_neg*py_neg : 0.0);

    Real gz = (pz_pos > 0.0 ? pz_pos*pz_pos : 0.0)
            + (pz_neg < 0.0 ? pz_neg*pz_neg : 0.0);

    return std::sqrt(gx + gy + gz);
}

// Right-biased value at i+1/2 from the "right" side.
// Implemented via mirrored call to the left-biased routine.
static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real weno5z_right(Real vm1, Real v0, Real vp1, Real vp2, Real vp3)
{
    // Mirror: (vp3, vp2, vp1, v0, vm1) acts like "left" in reversed indexing.
    return weno5z_left(vp3, vp2, vp1, v0, vm1);
}

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
// ======================= First-order extrapolation boundary fill ======================
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

        // X-lo boundary
        if (bx_ilo == dom_ilo) {
            ParallelFor(Box(IntVect(gbx_ilo, gbx_jlo, gbx_klo),
                           IntVect(dom_ilo-1, gbx_jhi, gbx_khi)),
                       ncomp,
                       [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept {
                           arr(i,j,k,n) = arr(dom_ilo,j,k,n);
                       });
        }

        // X-hi boundary
        if (bx_ihi == dom_ihi) {
            ParallelFor(Box(IntVect(dom_ihi+1, gbx_jlo, gbx_klo),
                           IntVect(gbx_ihi, gbx_jhi, gbx_khi)),
                       ncomp,
                       [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept {
                           arr(i,j,k,n) = arr(dom_ihi,j,k,n);
                       });
        }

        // Y-lo boundary
        if (bx_jlo == dom_jlo) {
            ParallelFor(Box(IntVect(gbx_ilo, gbx_jlo, gbx_klo),
                           IntVect(gbx_ihi, dom_jlo-1, gbx_khi)),
                       ncomp,
                       [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept {
                           arr(i,j,k,n) = arr(i,dom_jlo,k,n);
                       });
        }

        // Y-hi boundary
        if (bx_jhi == dom_jhi) {
            ParallelFor(Box(IntVect(gbx_ilo, dom_jhi+1, gbx_klo),
                           IntVect(gbx_ihi, gbx_jhi, gbx_khi)),
                       ncomp,
                       [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept {
                           arr(i,j,k,n) = arr(i,dom_jhi,k,n);
                       });
        }

        // Z-lo boundary
        if (bx_klo == dom_klo) {
            ParallelFor(Box(IntVect(gbx_ilo, gbx_jlo, gbx_klo),
                           IntVect(gbx_ihi, gbx_jhi, dom_klo-1)),
                       ncomp,
                       [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept {
                           arr(i,j,k,n) = arr(i,j,dom_klo,n);
                       });
        }

        // Z-hi boundary
        if (bx_khi == dom_khi) {
            ParallelFor(Box(IntVect(gbx_ilo, gbx_jlo, dom_khi+1),
                           IntVect(gbx_ihi, gbx_jhi, gbx_khi)),
                       ncomp,
                       [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept {
                           arr(i,j,k,n) = arr(i,j,dom_khi,n);
                       });
        }
    }
}

// ======================= Build RHS with WENO5-Z ====================
static void build_rhs_weno5z_old(const MultiFab& phi,
                             const MultiFab& vel,
                             MultiFab& rhs,
                             const Geometry& geom)
{
    const auto dx = geom.CellSize();
    GpuArray<Real,AMREX_SPACEDIM> gdx{dx[0], dx[1], dx[2]};

    rhs.setVal(0.0);

    // Fill ghost cells of phi for stencils using first-order extrapolation
    fill_boundary_extrap(const_cast<MultiFab&>(phi), geom);

    for (MFIter mfi(rhs, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& tbx = mfi.tilebox();
        auto const p = phi.const_array(mfi);
        auto const v = vel.const_array(mfi);
        auto const f = rhs.array(mfi);

        ParallelFor(tbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            // Velocity components (constant per cell in your setup)
            Real ux = v(i,j,k,0);
            Real uy = v(i,j,k,1);
            Real uz = v(i,j,k,2);

            // Flux splitting
            Real upx = std::max(ux, 0.0);
            Real umx = std::min(ux, 0.0);

            Real upy = std::max(uy, 0.0);
            Real umy = std::min(uy, 0.0);

            Real upz = std::max(uz, 0.0);
            Real umz = std::min(uz, 0.0);

            // ---------- X-direction: F_{i+1/2} and F_{i-1/2} ----------
            Real phiL_p = weno5z_left( p(i-2,j,k), p(i-1,j,k), p(i,j,k), p(i+1,j,k), p(i+2,j,k) );
            Real phiR_p = weno5z_right( p(i-1,j,k), p(i,j,k), p(i+1,j,k), p(i+2,j,k), p(i+3,j,k) );
            Real Fx_p   = upx * phiL_p + umx * phiR_p;

            Real phiL_m = weno5z_left( p(i-3,j,k), p(i-2,j,k), p(i-1,j,k), p(i,j,k), p(i+1,j,k) );
            Real phiR_m = weno5z_right( p(i-2,j,k), p(i-1,j,k), p(i,j,k), p(i+1,j,k), p(i+2,j,k) );
            Real Fx_m   = upx * phiL_m + umx * phiR_m;

            Real dFdx = (Fx_p - Fx_m) / gdx[0];

            // ---------- Y-direction ----------
            Real phiL_p_y = weno5z_left( p(i,j-2,k), p(i,j-1,k), p(i,j,k), p(i,j+1,k), p(i,j+2,k) );
            Real phiR_p_y = weno5z_right( p(i,j-1,k), p(i,j,k), p(i,j+1,k), p(i,j+2,k), p(i,j+3,k) );
            Real Fy_p     = upy * phiL_p_y + umy * phiR_p_y;

            Real phiL_m_y = weno5z_left( p(i,j-3,k), p(i,j-2,k), p(i,j-1,k), p(i,j,k), p(i,j+1,k) );
            Real phiR_m_y = weno5z_right( p(i,j-2,k), p(i,j-1,k), p(i,j,k), p(i,j+1,k), p(i,j+2,k) );
            Real Fy_m     = upy * phiL_m_y + umy * phiR_m_y;

            Real dFdy = (Fy_p - Fy_m) / gdx[1];

            // ---------- Z-direction ----------
            Real phiL_p_z = weno5z_left( p(i,j,k-2), p(i,j,k-1), p(i,j,k), p(i,j,k+1), p(i,j,k+2) );
            Real phiR_p_z = weno5z_right( p(i,j,k-1), p(i,j,k), p(i,j,k+1), p(i,j,k+2), p(i,j,k+3) );
            Real Fz_p     = upz * phiL_p_z + umz * phiR_p_z;

            Real phiL_m_z = weno5z_left( p(i,j,k-3), p(i,j,k-2), p(i,j,k-1), p(i,j,k), p(i,j,k+1) );
            Real phiR_m_z = weno5z_right( p(i,j,k-2), p(i,j,k-1), p(i,j,k), p(i,j,k+1), p(i,j,k+2) );
            Real Fz_m     = upz * phiL_m_z + umz * phiR_m_z;

            Real dFdz = (Fz_p - Fz_m) / gdx[2];

            // φ_t = - (∂F/∂x + ∂F/∂y + ∂F/∂z)
            f(i,j,k) = - (dFdx + dFdy + dFdz);
        });
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



// Distance from point P to line segment A->B (clamped projection)
static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real distance_point_to_segment(Real px, Real py, Real pz,
                               Real ax, Real ay, Real az,
                               Real bx, Real by, Real bz)
{
    Real vx = bx - ax;
    Real vy = by - ay;
    Real vz = bz - az;

    Real wx = px - ax;
    Real wy = py - ay;
    Real wz = pz - az;

    Real vv = vx*vx + vy*vy + vz*vz;

    // Handle degenerate (A==B): distance to point A
    if (vv <= 1.0e-30) {
        Real dx = wx, dy = wy, dz = wz;
        return std::sqrt(dx*dx + dy*dy + dz*dz);
    }

    Real t = (wx*vx + wy*vy + wz*vz) / vv;
    t = std::max(0.0, std::min(1.0, t));

    Real cx = ax + t * vx;
    Real cy = ay + t * vy;
    Real cz = az + t * vz;

    Real dx = px - cx;
    Real dy = py - cy;
    Real dz = pz - cz;

    return std::sqrt(dx*dx + dy*dy + dz*dz);
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

        // Ensure ghost cells of phi are fresh each iteration using first-order extrapolation
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

// ======================= Upwind advection RHS: u · ∇phi =========================
static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real upwind_u_dot_grad_phi(Array4<Real const> const& phi,
                           Array4<Real const> const& vel,
                           int i, int j, int k,
                           const GpuArray<Real,AMREX_SPACEDIM>& dx)
{
    Real ux = vel(i,j,k,0);
    Real uy = vel(i,j,k,1);
    Real uz = vel(i,j,k,2);

    Real dphidx = (ux > 0.0)
                  ? (phi(i,j,k) - phi(i-1,j,k)) / dx[0]
                  : (phi(i+1,j,k) - phi(i,j,k)) / dx[0];

    Real dphidy = (uy > 0.0)
                  ? (phi(i,j,k) - phi(i,j-1,k)) / dx[1]
                  : (phi(i,j+1,k) - phi(i,j,k)) / dx[1];

    Real dphidz = (uz > 0.0)
                  ? (phi(i,j,k) - phi(i,j,k-1)) / dx[2]
                  : (phi(i,j,k+1) - phi(i,j,k)) / dx[2];

    return ux * dphidx + uy * dphidy + uz * dphidz;
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

// ======================= Single-step upwind advection ===========================
static void advect_levelset_upwind(MultiFab& phi,
                                   const MultiFab& vel,
                                   const Geometry& geom,
                                   Real dt)
{
    const auto dx = geom.CellSize();
    GpuArray<Real,AMREX_SPACEDIM> gdx{dx[0], dx[1], dx[2]};

    // One ghost cell for upwind stencils; fill using first-order extrapolation
    fill_boundary_extrap(phi, geom);

    for (MFIter mfi(phi, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& tbx = mfi.tilebox();
        auto const p = phi.array(mfi);

        ParallelFor(tbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real val = p(i,j,k);
            p(i,j,k) = (val < 0.0) ? 0.0 : val;    // clamp to non-negative
        });
    }

    MultiFab rhs(phi.boxArray(), phi.DistributionMap(), 1, 0);

    for (MFIter mfi(rhs, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& tbx = mfi.tilebox();

        auto const p  = phi.const_array(mfi);
        auto const v  = vel.const_array(mfi);
        auto const f  = rhs.array(mfi);

        ParallelFor(tbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            f(i,j,k) = - upwind_u_dot_grad_phi(p, v, i, j, k, gdx);
        });
    }

    for (MFIter mfi(phi, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& tbx = mfi.tilebox();
        auto const p = phi.array(mfi);
        auto const f = rhs.const_array(mfi);

        ParallelFor(tbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            p(i,j,k) += dt * f(i,j,k);
        });
    }
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
        const int ng_phi = 3; // one ghost cell for upwind stencil
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
            //advect_levelset_upwind(phi, vel, geom, dt);
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