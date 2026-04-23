// ======================= Initialize velocity (constant) ==========================
static void init_velocity_constant(MultiFab& vel, const Geometry& geom, Real ux, Real uy, Real uz)
{
    const auto dx  = geom.CellSize();
    const auto plo = geom.ProbLo();
    for (MFIter mfi(vel); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const v = vel.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const Real x = plo[0] + (i + 0.5) * dx[0];
            const Real y = plo[1] + (j + 0.5) * dx[1];
            const Real z = plo[2] + (k + 0.5) * dx[2];
            v(i,j,k,0) = ux;
            v(i,j,k,1) = uy;
            v(i,j,k,2) = uz;
        });
    }
}