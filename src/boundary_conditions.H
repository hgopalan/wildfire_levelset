// ======================= Outflow boundary fill ======================
// Fills ghost cells using linear extrapolation (outflow treatment) for phi.
//
// Uses a hierarchical three-pass strategy to ensure deterministic, race-free
// behavior under arbitrary MPI decompositions:
//
//   Pass 1 – X faces only:  fill X-ghost cells, reading from domain interior
//             only (j, k clamped to the intersection of the domain range and
//             this FAB's ghost box).  No reads from ghost cells; no overlap.
//
//   Pass 2 – Y faces (+ X edges):  fill Y-ghost cells over the full X range
//             of this FAB (including X-ghost columns filled in Pass 1).
//             K is clamped to the domain range so Z-corners are not yet read.
//
//   Pass 3 – Z faces (+ edges/corners):  fill Z-ghost cells over the full
//             X,Y range of this FAB (all columns filled by Passes 1-2).
//             All source cells are guaranteed valid at this point.
//
// Motivation: the original code filled the complete ghost slab for each face
// (including edges/corners) in a single ParallelFor per face.  Edge/corner
// ghost cells were written by multiple kernels AND the source data for those
// cells (ghost rows in the transverse directions) had not yet been populated,
// creating circular reads from uninitialised memory.  Different MPI
// decompositions exposed the uninitialised reads differently, producing
// nondeterministic crashes.
static void fill_boundary_extrap(MultiFab& mf, const Geometry& geom)
{
    // Sanity-check: at least one ghost cell layer must be present.
    // (AMREX_ASSERT is active unless AMREX_ASSERT_OFF is defined.)
    AMREX_ASSERT(mf.nGrow() >= 1);

    // Fill interior ghost cells via MPI exchange between ranks.
    // is_periodic is {0,0,0}, so this only copies between neighbouring ranks
    // (no periodic wrap-around).
    mf.FillBoundary(geom.periodicity());

    const auto& domain = geom.Domain();

    const int dom_ilo = domain.smallEnd(0);
    const int dom_ihi = domain.bigEnd(0);
    const int dom_jlo = domain.smallEnd(1);
    const int dom_jhi = domain.bigEnd(1);
    const int dom_klo = domain.smallEnd(2);
    const int dom_khi = domain.bigEnd(2);

    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const Box& bx  = mfi.validbox();
        const Box& gbx = mfi.fabbox(); // includes ghost cells
        auto const arr = mf.array(mfi);

        const int bx_ilo = bx.smallEnd(0);
        const int bx_ihi = bx.bigEnd(0);
        const int bx_jlo = bx.smallEnd(1);
        const int bx_jhi = bx.bigEnd(1);
        const int bx_klo = bx.smallEnd(2);
        const int bx_khi = bx.bigEnd(2);

        const int gbx_ilo = gbx.smallEnd(0);
        const int gbx_ihi = gbx.bigEnd(0);
        const int gbx_jlo = gbx.smallEnd(1);
        const int gbx_jhi = gbx.bigEnd(1);
        const int gbx_klo = gbx.smallEnd(2);
        const int gbx_khi = gbx.bigEnd(2);

        const int ncomp = mf.nComp();

        // Clamp ranges: intersect domain interior with this FAB's ghost box.
        // Prevents out-of-bounds reads/writes on FABs that do not span the
        // full domain extent (the common case for multi-rank decompositions).
        const int cl_jlo = std::max(gbx_jlo, dom_jlo);
        const int cl_jhi = std::min(gbx_jhi, dom_jhi);
        const int cl_klo = std::max(gbx_klo, dom_klo);
        const int cl_khi = std::min(gbx_khi, dom_khi);

        // ------------------------------------------------------------
        // Pass 1: X-face ghost cells only.
        // Y,Z restricted to the clamped domain range so all source reads
        // come from valid domain-interior cells.
        // ------------------------------------------------------------
        if (bx_ilo == dom_ilo) {
            ParallelFor(Box(IntVect(gbx_ilo, cl_jlo, cl_klo),
                           IntVect(dom_ilo-1, cl_jhi, cl_khi)),
                       ncomp,
                       [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept {
                           arr(i,j,k,n) = arr(dom_ilo,j,k,n)
                               + Real(dom_ilo-i)*(arr(dom_ilo,j,k,n) - arr(dom_ilo+1,j,k,n));
                       });
        }

        if (bx_ihi == dom_ihi) {
            ParallelFor(Box(IntVect(dom_ihi+1, cl_jlo, cl_klo),
                           IntVect(gbx_ihi, cl_jhi, cl_khi)),
                       ncomp,
                       [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept {
                           arr(i,j,k,n) = arr(dom_ihi,j,k,n)
                               + Real(i-dom_ihi)*(arr(dom_ihi,j,k,n) - arr(dom_ihi-1,j,k,n));
                       });
        }

        // ------------------------------------------------------------
        // Pass 2: Y-face ghost cells (includes X-ghost edge columns).
        // Full X range (gbx_ilo..gbx_ihi) is safe: X-ghost cells at
        // domain Y,Z levels were populated in Pass 1.
        // K clamped to domain range so Z-corners are not yet read.
        // ------------------------------------------------------------
        if (bx_jlo == dom_jlo) {
            ParallelFor(Box(IntVect(gbx_ilo, gbx_jlo, cl_klo),
                           IntVect(gbx_ihi, dom_jlo-1, cl_khi)),
                       ncomp,
                       [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept {
                           arr(i,j,k,n) = arr(i,dom_jlo,k,n)
                               + Real(dom_jlo-j)*(arr(i,dom_jlo,k,n) - arr(i,dom_jlo+1,k,n));
                       });
        }

        if (bx_jhi == dom_jhi) {
            ParallelFor(Box(IntVect(gbx_ilo, dom_jhi+1, cl_klo),
                           IntVect(gbx_ihi, gbx_jhi, cl_khi)),
                       ncomp,
                       [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept {
                           arr(i,j,k,n) = arr(i,dom_jhi,k,n)
                               + Real(j-dom_jhi)*(arr(i,dom_jhi,k,n) - arr(i,dom_jhi-1,k,n));
                       });
        }

        // ------------------------------------------------------------
        // Pass 3: Z-face ghost cells (includes all edge/corner columns).
        // Full X,Y range (gbx_*) is safe: all transverse ghost cells at
        // Z-face levels were populated by Passes 1-2.
        // ------------------------------------------------------------
        if (bx_klo == dom_klo) {
            ParallelFor(Box(IntVect(gbx_ilo, gbx_jlo, gbx_klo),
                           IntVect(gbx_ihi, gbx_jhi, dom_klo-1)),
                       ncomp,
                       [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept {
                           arr(i,j,k,n) = arr(i,j,dom_klo,n)
                               + Real(dom_klo-k)*(arr(i,j,dom_klo,n) - arr(i,j,dom_klo+1,n));
                       });
        }

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