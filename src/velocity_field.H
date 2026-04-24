// ======================= Initialize velocity (constant) ==========================
static void init_velocity_constant(MultiFab& vel, const Geometry& geom, Real ux, Real uy, Real uz)
{
#if (AMREX_SPACEDIM == 3)
    (void)geom;  // Suppress unused parameter warning in case we don't need plo/dx
#else
    (void)geom;
    (void)uz;  // uz not used in 2D
#endif
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