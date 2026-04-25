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
                         Real xmin, Real ymin, Real zmin, 
                         Real xmax, Real ymax, Real zmax) 
{
    const auto dx  = geom.CellSize();
    const auto plo = geom.ProbLo();

#if (AMREX_SPACEDIM == 2)
    (void)zmin;  // Suppress unused parameter warning
    (void)zmax;
#endif
    for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const p = phi.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            // Get cell center coordinates
            const Real x = plo[0] + (i + 0.5) * dx[0];
            const Real y = plo[1] + (j + 0.5) * dx[1];

#if (AMREX_SPACEDIM == 3)
            const Real z = plo[2] + (k + 0.5) * dx[2];
#else
            (void)k;  // Suppress unused parameter warning
#endif
        if (x > xmin && y > ymin && x< xmax && y < ymax) {
            p(i,j,k) = -10; 
        }
        else
        {
            p(i,j,k)=0.0;
        }
        });
    }
}

// ======================= Initialize indicator function (sphere) ======================
// For FARSITE mode: phi = -10 inside the sphere, phi = 0 outside
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
            // Indicator function: -10 inside, 0 outside
            p(i,j,k) = (d < radius) ? -10.0 : 0.0;
        });
    }
}

// ======================= Initialize indicator function (box) ======================
// For FARSITE mode: phi = -10 inside the box, phi = 0 outside
static void init_phi_box_indicator(MultiFab& phi, const Geometry& geom,
                         Real xmin, Real ymin, Real zmin, 
                         Real xmax, Real ymax, Real zmax) 
{
    const auto dx  = geom.CellSize();
    const auto plo = geom.ProbLo();

#if (AMREX_SPACEDIM == 2)
    (void)zmin;  // Suppress unused parameter warning
    (void)zmax;
#endif
    for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const p = phi.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            // Get cell center coordinates
            const Real x = plo[0] + (i + 0.5) * dx[0];
            const Real y = plo[1] + (j + 0.5) * dx[1];

#if (AMREX_SPACEDIM == 3)
            const Real z = plo[2] + (k + 0.5) * dx[2];
#else
            (void)k;  // Suppress unused parameter warning
#endif
        if (x > xmin && y > ymin && x< xmax && y < ymax) {
            p(i,j,k) = -10; 
        }
        else
        {
            p(i,j,k)=0.0;
        }
        });
    }
}