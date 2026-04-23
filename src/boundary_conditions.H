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