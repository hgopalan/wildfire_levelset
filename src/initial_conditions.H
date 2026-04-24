// ======================= Initialize signed distance (sphere) ======================
static void init_phi_sphere(MultiFab& phi, const Geometry& geom,
                            Real cx, Real cy, Real cz, Real radius)
{
#if (AMREX_SPACEDIM == 2)
    (void)cz;  // Suppress unused parameter warning in 2D
#endif
    const auto dx  = geom.CellSize();
    const auto plo = geom.ProbLo();

    for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const p = phi.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real x = plo[0] + (i + 0.5) * dx[0];
            Real y = plo[1] + (j + 0.5) * dx[1];
#if (AMREX_SPACEDIM == 3)
            Real z = plo[2] + (k + 0.5) * dx[2];
            Real d = std::sqrt((x - cx)*(x - cx) +
                               (y - cy)*(y - cy) +
                               (z - cz)* (z - cz)
                               ) - radius;
#else
            (void)k;  // Suppress unused parameter warning
            Real d = std::sqrt((x - cx)*(x - cx) +
                               (y - cy)*(y - cy)
                               ) - radius;
#endif
            p(i,j,k) = d; // negative inside the sphere
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
#if (AMREX_SPACEDIM == 3)
    const auto hz = 0.5 * sz;
#else
    (void)cz;  // Suppress unused parameter warning
    (void)sz;
#endif
    for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const p = phi.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            // Get cell center coordinates
            const Real x = plo[0] + (i + 0.5) * dx[0];
            const Real y = plo[1] + (j + 0.5) * dx[1];

            // Compute distance to each face relative to the center
            Real dx_val = std::abs(x - cx) - hx;
            Real dy_val = std::abs(y - cy) - hy;

#if (AMREX_SPACEDIM == 3)
            const Real z = plo[2] + (k + 0.5) * dx[2];
            Real dz_val = std::abs(z - cz) - hz;
            // The signed distance to a box is the max of the component distances
            // This yields d < 0 inside, d = 0 on the surface, d > 0 outside
            Real d = std::max({dx_val, dy_val, dz_val});
#else
            (void)k;  // Suppress unused parameter warning
            Real d = std::max(dx_val, dy_val);
#endif

            p(i,j,k) = d; 
        });
    }
}

// ======================= Initialize indicator function (sphere) ======================
// For FARSITE mode: phi = 1 inside the sphere, phi = 0 outside
static void init_phi_sphere_indicator(MultiFab& phi, const Geometry& geom,
                                       Real cx, Real cy, Real cz, Real radius)
{
#if (AMREX_SPACEDIM == 2)
    (void)cz;  // Suppress unused parameter warning in 2D
#endif
    const auto dx  = geom.CellSize();
    const auto plo = geom.ProbLo();

    for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const p = phi.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real x = plo[0] + (i + 0.5) * dx[0];
            Real y = plo[1] + (j + 0.5) * dx[1];
#if (AMREX_SPACEDIM == 3)
            Real z = plo[2] + (k + 0.5) * dx[2];
            Real d = std::sqrt((x - cx)*(x - cx) +
                               (y - cy)*(y - cy) +
                               (z - cz)*(z - cz)
                               );
#else
            (void)k;  // Suppress unused parameter warning
            Real d = std::sqrt((x - cx)*(x - cx) +
                               (y - cy)*(y - cy)
                               );
#endif
            // Indicator function: 1 inside, 0 outside
            p(i,j,k) = (d < radius) ? 1.0 : 0.0;
        });
    }
}

// ======================= Initialize indicator function (box) ======================
// For FARSITE mode: phi = 1 inside the box, phi = 0 outside
static void init_phi_box_indicator(MultiFab& phi, const Geometry& geom,
                                    Real cx, Real cy, Real cz, 
                                    Real sx, Real sy, Real sz) 
{
    const auto dx  = geom.CellSize();
    const auto plo = geom.ProbLo();
    const auto hx = 0.5 * sx;
    const auto hy = 0.5 * sy;
#if (AMREX_SPACEDIM == 3)
    const auto hz = 0.5 * sz;
#else
    (void)cz;  // Suppress unused parameter warning
    (void)sz;
#endif
    for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const p = phi.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            // Get cell center coordinates
            const Real x = plo[0] + (i + 0.5) * dx[0];
            const Real y = plo[1] + (j + 0.5) * dx[1];

            // Check if point is inside the box
            bool inside_x = std::abs(x - cx) < hx;
            bool inside_y = std::abs(y - cy) < hy;

#if (AMREX_SPACEDIM == 3)
            const Real z = plo[2] + (k + 0.5) * dx[2];
            bool inside_z = std::abs(z - cz) < hz;
            bool inside = inside_x && inside_y && inside_z;
#else
            (void)k;  // Suppress unused parameter warning
            bool inside = inside_x && inside_y;
#endif
            
            // Indicator function: 1 inside, 0 outside
            p(i,j,k) = inside ? 1.0 : 0.0;
        });
    }
}