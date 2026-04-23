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