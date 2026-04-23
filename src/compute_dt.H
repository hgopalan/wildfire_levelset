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