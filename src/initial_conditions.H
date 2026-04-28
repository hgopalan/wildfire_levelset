// ======================= Initialize from CSV fire points (Gaussian blobs) =========
// Each row in the CSV provides the coordinates of an ignition point (X Y in 2D,
// X Y Z in 3D).  For every grid cell a Gaussian contribution
//   phi_contrib = -exp( -r^2 / (2*sigma^2) )
// is computed for each ignition point; phi at a cell is set to the minimum (most
// negative) value across all points.  Values are in [-1, 0): negative near fire
// centres and approaching 0 far away, which is consistent with the signed-distance
// convention used by the level-set solver (negative = inside fire).
//
// sigma   > 0  use provided value (physical units)
// sigma  <= 0  auto-select as 3 * min(dx, dy [, dz])
static void init_phi_from_fire_points_csv(MultiFab& phi, const Geometry& geom,
                                          const std::string& filename,
                                          amrex::Real sigma)
{
    // ---- Read fire points on rank 0 and broadcast ----
    std::vector<amrex::Real> pts_x, pts_y, pts_z;

    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ifstream fin(filename);
        if (!fin.is_open()) {
            amrex::Abort("init_phi_from_fire_points_csv: cannot open file: " + filename);
        }
        std::string line;
        while (std::getline(fin, line)) {
            // Skip blank lines and '#' comment lines (short-circuit ensures
            // line[0] is only accessed when the line is non-empty)
            if (line.empty() || line[0] == '#') continue;
            // Replace commas with spaces so both CSV and space-delimited work
            for (char& c : line) { if (c == ',') c = ' '; }
            std::istringstream iss(line);
            amrex::Real x = 0.0, y = 0.0, z = 0.0;
            if (!(iss >> x >> y)) continue;  // skip malformed lines
            iss >> z;                         // optional Z (stays 0 if not present)
            pts_x.push_back(x);
            pts_y.push_back(y);
            pts_z.push_back(z);
        }
        amrex::Print() << "init_phi_from_fire_points_csv: read " << pts_x.size()
                       << " ignition point(s) from '" << filename << "'\n";
    }

    // Broadcast the point arrays to all MPI ranks
    {
        int n = static_cast<int>(pts_x.size());
        amrex::ParallelDescriptor::Bcast(&n, 1, amrex::ParallelDescriptor::IOProcessorNumber());
        pts_x.resize(n); pts_y.resize(n); pts_z.resize(n);
        amrex::ParallelDescriptor::Bcast(pts_x.data(), n, amrex::ParallelDescriptor::IOProcessorNumber());
        amrex::ParallelDescriptor::Bcast(pts_y.data(), n, amrex::ParallelDescriptor::IOProcessorNumber());
        amrex::ParallelDescriptor::Bcast(pts_z.data(), n, amrex::ParallelDescriptor::IOProcessorNumber());
    }

    if (pts_x.empty()) {
        amrex::Abort("init_phi_from_fire_points_csv: no valid ignition points found in '" + filename + "'");
    }

    const auto dx  = geom.CellSize();
    const auto plo = geom.ProbLo();

    // Auto-compute sigma if not specified (3 cells)
    if (sigma <= 0.0) {
#if (AMREX_SPACEDIM == 3)
        sigma = 3.0 * std::min({dx[0], dx[1], dx[2]});
#else
        sigma = 3.0 * std::min(dx[0], dx[1]);
#endif
        amrex::Print() << "init_phi_from_fire_points_csv: auto sigma = " << sigma << "\n";
    }

    const amrex::Real inv2sig2 = 1.0 / (2.0 * sigma * sigma);
    const int npts = static_cast<int>(pts_x.size());

    // Raw pointers captured by [=] in the ParallelFor lambda.
    // The project currently builds with AMReX_GPU_BACKEND=NONE so these host
    // pointers are safe to access inside the (CPU) ParallelFor / AMREX_GPU_DEVICE
    // lambda.  When GPU support is enabled the vectors must be replaced with
    // amrex::Gpu::DeviceVector and the data copied to device before launching.
    const amrex::Real* px = pts_x.data();
    const amrex::Real* py = pts_y.data();
    const amrex::Real* pz = pts_z.data();

    for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const p = phi.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const amrex::Real x = plo[0] + (i + 0.5) * dx[0];
            const amrex::Real y = plo[1] + (j + 0.5) * dx[1];
#if (AMREX_SPACEDIM == 3)
            const amrex::Real z = plo[2] + (k + 0.5) * dx[2];
#else
            (void)k;
            (void)pz;  // Z coordinates unused in 2D
#endif
            amrex::Real phi_min = 0.0;
            for (int n = 0; n < npts; ++n) {
#if (AMREX_SPACEDIM == 3)
                amrex::Real r2 = (x - px[n])*(x - px[n])
                               + (y - py[n])*(y - py[n])
                               + (z - pz[n])*(z - pz[n]);
#else
                amrex::Real r2 = (x - px[n])*(x - px[n])
                               + (y - py[n])*(y - py[n]);
#endif
                amrex::Real contrib = -std::exp(-r2 * inv2sig2);
                if (contrib < phi_min) phi_min = contrib;
            }
            p(i,j,k) = phi_min;
        });
    }
}

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
#if (AMREX_SPACEDIM == 3)
            Real nz = param3 / norm;
#endif
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
            Real radius = param4;
            
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                Real x = plo[0] + (i + 0.5) * dx[0];
                Real y = plo[1] + (j + 0.5) * dx[1];
#if (AMREX_SPACEDIM == 3)
                Real cz = param3;
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
            Real rx = param4;
            Real ry = param5;
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                Real x = plo[0] + (i + 0.5) * dx[0];
                Real y = plo[1] + (j + 0.5) * dx[1];
#if (AMREX_SPACEDIM == 3)
                Real cz = param3;
                Real rz = param6;   
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
                (void) param3; // unused in 2D
                (void) param6; // unused in 2D
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