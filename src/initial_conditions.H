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
                               (z - cz)* (z - cz)
                               ) - radius;
            p(i,j,k) = std::max(-1.0,std::min(d,1.0)); // negative inside the sphere
        });
    }
}

// ======================= Initialize signed distance (box) ======================
static void init_phi_box(MultiFab& phi, const Geometry& geom,
                         Real cx, Real cy, Real cz, 
                         Real sx, Real sy, Real sz) 
{
    const auto dx  = geom.CellSize();
    const auto plo = geom.ProbLo();
    const auto hx = 0.5 * sx;
    const auto hy = 0.5 * sy;
    const auto hz = 0.5 * sz;
    for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const p = phi.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            // Get cell center coordinates
            const Real x = plo[0] + (i + 0.5) * dx[0];
            const Real y = plo[1] + (j + 0.5) * dx[1];
            const Real z = plo[2] + (k + 0.5) * dx[2];

            // Compute distance to each face relative to the center
            Real dx_val = std::abs(x - cx) - hx;
            Real dy_val = std::abs(y - cy) - hy;
            Real dz_val = std::abs(z - cz) - hz;

            // The signed distance to a box is the max of the component distances
            // This yields d < 0 inside, d = 0 on the surface, d > 0 outside
            Real d = std::max({dx_val, dy_val, dz_val});

            p(i,j,k) = d; 
        });
    }
}