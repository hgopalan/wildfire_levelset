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

// ======================= Initialize signed distance (ellipse) ======================
static void init_phi_ellipse(MultiFab& phi, const Geometry& geom,
                             Real cx, Real cy, Real cz, Real rx, Real ry, Real rz)
{
#if (AMREX_SPACEDIM == 2)
    (void)cz;  // Suppress unused parameter warning in 2D
    (void)rz;
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
            // For an ellipse, we compute the approximate signed distance
            // The exact SDF is complex, so we use an approximation
            Real dx_norm = (x - cx) / rx;
            Real dy_norm = (y - cy) / ry;
            Real dz_norm = (z - cz) / rz;
            Real r = std::sqrt(dx_norm*dx_norm + dy_norm*dy_norm + dz_norm*dz_norm);
            
            // Approximate SDF for ellipsoid
            // This gives the correct sign and is approximately a distance
            Real d = r - 1.0;
            // Scale back to get approximate distance
            // Add small epsilon (1.0e-16) to denominator to prevent division by zero
            Real scale = std::sqrt((dx_norm*dx_norm*rx*rx + 
                                   dy_norm*dy_norm*ry*ry + 
                                   dz_norm*dz_norm*rz*rz) / 
                                  (dx_norm*dx_norm + dy_norm*dy_norm + dz_norm*dz_norm + 1.0e-16));
            d *= scale;
#else
            (void)k;  // Suppress unused parameter warning
            // For 2D ellipse
            Real dx_norm = (x - cx) / rx;
            Real dy_norm = (y - cy) / ry;
            Real r = std::sqrt(dx_norm*dx_norm + dy_norm*dy_norm);
            
            // Approximate SDF for ellipse
            Real d = r - 1.0;
            // Scale back to get approximate distance
            // Add small epsilon (1.0e-16) to denominator to prevent division by zero
            Real scale = std::sqrt((dx_norm*dx_norm*rx*rx + dy_norm*dy_norm*ry*ry) / 
                                  (dx_norm*dx_norm + dy_norm*dy_norm + 1.0e-16));
            d *= scale;
#endif
            p(i,j,k) = d; // negative inside the ellipse
        });
    }
}

// ======================= Initialize indicator function (ellipse) ======================
// For FARSITE mode: phi = -10 inside the ellipse, phi = 0 outside
static void init_phi_ellipse_indicator(MultiFab& phi, const Geometry& geom,
                                        Real cx, Real cy, Real cz, Real rx, Real ry, Real rz)
{
#if (AMREX_SPACEDIM == 2)
    (void)cz;  // Suppress unused parameter warning in 2D
    (void)rz;
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
            Real dx_norm = (x - cx) / rx;
            Real dy_norm = (y - cy) / ry;
            Real dz_norm = (z - cz) / rz;
            Real r = std::sqrt(dx_norm*dx_norm + dy_norm*dy_norm + dz_norm*dz_norm);
#else
            (void)k;  // Suppress unused parameter warning
            Real dx_norm = (x - cx) / rx;
            Real dy_norm = (y - cy) / ry;
            Real r = std::sqrt(dx_norm*dx_norm + dy_norm*dy_norm);
#endif
            // Indicator function: -10 inside, 0 outside
            p(i,j,k) = (r < 1.0) ? -10.0 : 0.0;
        });
    }
}

// ======================= Initialize from EB signed distance function ======================
// This function initializes phi using an implicit function that could represent
// an embedded boundary geometry
static void init_phi_from_eb_implicit(MultiFab& phi, const Geometry& geom,
                                       const std::string& eb_type,
                                       Real param1, Real param2, Real param3,
                                       Real param4, Real param5, Real param6)
{
    const auto dx  = geom.CellSize();
    const auto plo = geom.ProbLo();

    for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const p = phi.array(mfi);

        if (eb_type == "plane") {
            // Plane: param1*x + param2*y + param3*z + param4 = 0
            // Normalized form where sqrt(param1^2 + param2^2 + param3^2) = 1
            Real norm = std::sqrt(param1*param1 + param2*param2 + param3*param3);
            Real nx = param1 / norm;
            Real ny = param2 / norm;
            Real nz = param3 / norm;
            Real d0 = param4 / norm;
            
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                Real x = plo[0] + (i + 0.5) * dx[0];
                Real y = plo[1] + (j + 0.5) * dx[1];
#if (AMREX_SPACEDIM == 3)
                Real z = plo[2] + (k + 0.5) * dx[2];
                p(i,j,k) = nx*x + ny*y + nz*z + d0;
#else
                (void)k;
                p(i,j,k) = nx*x + ny*y + d0;
#endif
            });
        }
        else if (eb_type == "cylinder") {
            // Cylinder: center (param1, param2, param3), radius param4, axis along z
            Real cx = param1;
            Real cy = param2;
            Real cz = param3;
            Real radius = param4;
            (void)cz; // unused in 2D
            
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                Real x = plo[0] + (i + 0.5) * dx[0];
                Real y = plo[1] + (j + 0.5) * dx[1];
                (void)k;
                Real r = std::sqrt((x - cx)*(x - cx) + (y - cy)*(y - cy));
                p(i,j,k) = r - radius;
            });
        }
        else if (eb_type == "sphere") {
            // Sphere: center (param1, param2, param3), radius param4
            Real cx = param1;
            Real cy = param2;
            Real cz = param3;
            Real radius = param4;
            
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                Real x = plo[0] + (i + 0.5) * dx[0];
                Real y = plo[1] + (j + 0.5) * dx[1];
#if (AMREX_SPACEDIM == 3)
                Real z = plo[2] + (k + 0.5) * dx[2];
                Real r = std::sqrt((x - cx)*(x - cx) + (y - cy)*(y - cy) + (z - cz)*(z - cz));
#else
                (void)k;
                Real r = std::sqrt((x - cx)*(x - cx) + (y - cy)*(y - cy));
#endif
                p(i,j,k) = r - radius;
            });
        }
        else if (eb_type == "ellipsoid") {
            // Ellipsoid: center (param1, param2, param3), semi-axes (param4, param5, param6)
            Real cx = param1;
            Real cy = param2;
            Real cz = param3;
            Real rx = param4;
            Real ry = param5;
            Real rz = param6;
            
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                Real x = plo[0] + (i + 0.5) * dx[0];
                Real y = plo[1] + (j + 0.5) * dx[1];
#if (AMREX_SPACEDIM == 3)
                Real z = plo[2] + (k + 0.5) * dx[2];
                Real dx_norm = (x - cx) / rx;
                Real dy_norm = (y - cy) / ry;
                Real dz_norm = (z - cz) / rz;
                Real r = std::sqrt(dx_norm*dx_norm + dy_norm*dy_norm + dz_norm*dz_norm);
                Real d = r - 1.0;
                // Add small epsilon (1.0e-16) to denominator to prevent division by zero
                Real scale = std::sqrt((dx_norm*dx_norm*rx*rx + 
                                       dy_norm*dy_norm*ry*ry + 
                                       dz_norm*dz_norm*rz*rz) / 
                                      (dx_norm*dx_norm + dy_norm*dy_norm + dz_norm*dz_norm + 1.0e-16));
                p(i,j,k) = d * scale;
#else
                (void)k;
                Real dx_norm = (x - cx) / rx;
                Real dy_norm = (y - cy) / ry;
                Real r = std::sqrt(dx_norm*dx_norm + dy_norm*dy_norm);
                Real d = r - 1.0;
                // Add small epsilon (1.0e-16) to denominator to prevent division by zero
                Real scale = std::sqrt((dx_norm*dx_norm*rx*rx + dy_norm*dy_norm*ry*ry) / 
                                      (dx_norm*dx_norm + dy_norm*dy_norm + 1.0e-16));
                p(i,j,k) = d * scale;
#endif
            });
        }
        else {
            amrex::Abort("Unknown EB type: " + eb_type);
        }
    }
}