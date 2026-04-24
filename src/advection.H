// ============================================================================
// UPDATED: Build RHS with WENO5-Z using Eq. (2) with Rothermel (1972) ROS
// ============================================================================
#include "parse_inputs.H"
#include "rothermel_model.H"

static void build_rhs_weno5z(const MultiFab& phi,
                             const MultiFab& vel,
                             MultiFab& rhs,
                             const Geometry& geom,
                             const RothermelParams& rp)
{
    const auto dx = geom.CellSize();
#if (AMREX_SPACEDIM == 3)
    GpuArray<Real,AMREX_SPACEDIM> gdx{dx[0], dx[1], dx[2]};
#else
    GpuArray<Real,AMREX_SPACEDIM> gdx{dx[0], dx[1]};
#endif

    // Artificial viscosity ε per Eq. (2); use a single uniform value across the domain.
    // Paper default ~0.4 is common operationally. 
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

    // Compute Rothermel (1972) fire spread parameters
    const RothermelComputed rc = compute_rothermel_params(rp);
    const Real R0 = rc.R0;
    const Real C = rc.C;
    const Real B = rc.B;
    const Real beta_ratio_E = rc.beta_ratio_E;
    const Real phi_s = rc.phi_s;
    const Real wind_conv = rc.wind_conv;
    const Real ros_conv = rc.ros_conv;

    // RHS per Eq. (2): ∂φ/∂t = - R ( |∇φ| - ε Δφ ), with R from Rothermel (1972). 
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

            // Horizontal wind components
            Real ux = v(i,j,k,0);
            Real uy = v(i,j,k,1);

            // Wind speed in ft/min (Rothermel units)
            Real U_ftmin = std::sqrt(ux*ux + uy*uy) * wind_conv;

            // Wind factor φ_w (Rothermel Eq. 47)
            Real phi_w = (U_ftmin > 0.0) ? C * std::pow(U_ftmin, B) * beta_ratio_E : 0.0;

            // Rate of spread [ft/min] → simulation units
            Real R = std::max(R0 * (1.0 + phi_w + phi_s), 0.0) * ros_conv;

            // Final RHS (keeps your WENO5 & RK3 integrators intact)
            fr(i,j,k) = -R * (grad_norm - eps_visc * lap);
        });
    }
}


// ======================= WENO5-Z advection + RK3 ==================
static void advect_levelset_weno5z_rk3(MultiFab& phi,
                                       const MultiFab& vel,
                                       const Geometry& geom,
                                       Real dt,
                                       const RothermelParams& rp)
{
    const auto& ba = phi.boxArray();
    const auto& dm = phi.DistributionMap();

    // Stage buffers with 3 ghosts (needed for WENO stencils)
    MultiFab phi1(ba, dm, 1, 3);
    MultiFab phi2(ba, dm, 1, 3);
    MultiFab rhs(ba, dm, 1, 0);


    // ---------- Stage 1 ----------
    build_rhs_weno5z(phi, vel, rhs, geom, rp);
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
    build_rhs_weno5z(phi1, vel, rhs, geom, rp);
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
    build_rhs_weno5z(phi2, vel, rhs, geom, rp);
    for (MFIter mfi(phi, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& tbx = mfi.tilebox();
        auto const p   = phi.const_array(mfi);
        auto const q2c = phi2.const_array(mfi);
        auto const r   = rhs.const_array(mfi);
        auto const pout= phi.array(mfi);
        ParallelFor(tbx, [=] AMREX_GPU_DEVICE (int i,int j,int k) noexcept {
           // pout(i,j,k) = (1.0/3.0) * p(i,j,k) + (2.0/3.0) * ( q2c(i,j,k) + dt * r(i,j,k) );

        Real newphi = (1.0/3.0) * p(i,j,k) + (2.0/3.0) * ( q2c(i,j,k) + dt * r(i,j,k) );
        pout(i,j,k) = newphi;

        });
    }
}
