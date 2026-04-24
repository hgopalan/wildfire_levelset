// ======================= Initialize velocity (constant) ==========================
static void init_velocity_constant(MultiFab& vel, const Geometry& geom, Real ux, Real uy, Real uz)
{
    (void)geom;  // Suppress unused parameter warning
#if (AMREX_SPACEDIM == 2)
    (void)uz;  // uz not used in 2D
#endif
    for (MFIter mfi(vel); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const v = vel.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            v(i,j,k,0) = ux;
            v(i,j,k,1) = uy;
#if (AMREX_SPACEDIM == 3)
            v(i,j,k,2) = uz;
#else
            v(i,j,k,2) = Real(0.0);  // In 2D, set z-component to 0
#endif
        });
    }
}