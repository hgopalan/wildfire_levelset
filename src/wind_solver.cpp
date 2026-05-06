// ==========================================================================
// wind_solver.cpp
// Terrain-following mass-consistent 3-D wind solver (QUIC-URB style)
//
// Reads a terrain file (X Y Z), a reference wind vector at a specified
// height, and surface roughness.  Constructs a log-law initial wind field
// over the terrain and enforces mass consistency by solving the anisotropic
// Poisson equation
//
//   -(α_h² ∂²λ/∂x² + α_h² ∂²λ/∂y² + α_v² ∂²λ/∂z²) = -(∇·u₀)
//
// via AMReX MLMG (MLABecLaplacian) on a single-level (level-0) 3-D
// Cartesian grid.  The corrected divergence-free wind field is written as
// an AMReX plotfile and, optionally, as a terrain-aligned CSV slice.
//
// Terrain-following initialisation:
//   For each horizontal column (i, j) the local terrain elevation z_s is
//   interpolated from the data file with inverse-distance weighting (IDW).
//   The height above ground level (AGL) of cell centre (i, j, k) is
//       z_agl = (k + 0.5) * dz  [uniform vertical grid in AGL space]
//   so the lowest model level always lies dz/2 above the terrain surface.
//   The log-law profile is evaluated at z_agl:
//       u(z_agl) = u* / κ  * ln((z_agl + z0) / z0)
//       u*       = κ * |U_ref| / ln((z_ref + z0) / z0)
//   Cells at z_agl <= 0 (sub-surface) are zeroed and masked in the solve.
//
// Mass-consistent correction (Lagrange multiplier method, Sherman 1978):
//   Minimise E = ∫[(u−u₀)²/α_h² + (v−v₀)²/α_h² + (w−w₀)²/α_v²] dV
//   subject to ∇·u = 0.
//   Euler-Lagrange conditions give:
//       u = u₀ − α_h² ∂λ/∂x
//       v = v₀ − α_h² ∂λ/∂y
//       w = w₀ − α_v² ∂λ/∂z
//   Substituting into ∇·u = 0 yields the anisotropic Poisson equation
//   solved above.
//
// Terrain-aligned extraction:
//   After the corrected wind field is computed, an optional 2-D slice can be
//   extracted and written as a CSV.  The slice follows the terrain surface at
//   a constant height above ground level (AGL), so each row records the local
//   terrain elevation and the resulting physical extraction height:
//       z_physical(i,j) = z_terrain(i,j) + z_agl_requested
//   Specify the extraction level with ONE of:
//       extract_agl  = <height_m>   # AGL height [m] (snapped to nearest cell)
//       extract_k    = <k_index>    # explicit k-index (0 = lowest cell)
//   The output CSV (extract_file) has columns:
//       x, y, z_terrain, z_physical, z_agl, u, v, w, speed
//
// Usage:  wind_solver inputs.i   (or  wind_solver key=value ...)
//
// Key parameters (with defaults):
//   terrain_file  = terrain.csv   # X Y Z, whitespace- or comma-separated
//   U_ref         = 10.0          # reference wind x-component [m/s]
//   V_ref         = 0.0           # reference wind y-component [m/s]
//   z_ref         = 10.0          # reference height above local terrain [m]
//   z0            = 0.1           # aerodynamic roughness length [m]
//   dx            = 30.0          # grid spacing x [m]
//   dy            = 30.0          # grid spacing y [m]
//   dz            = 30.0          # grid spacing z (AGL) [m]
//   domain_height = 300.0         # vertical extent above terrain base [m]
//   alpha_h       = 1.0           # horizontal Lagrange anisotropy factor
//   alpha_v       = 1.0           # vertical   Lagrange anisotropy factor
//   mlmg_verbose  = 1             # MLMG verbosity (0 = silent, 4 = max)
//   tol_rel       = 1.e-8         # MLMG relative tolerance
//   max_grid_size = 32            # maximum AMReX box size (per dimension)
//   plot_file     = plt_wind      # output plotfile prefix
//   extract_agl   = -1.0          # terrain-aligned extraction AGL [m] (<0 = off)
//   extract_k     = -1            # explicit k-index extraction (<0 = off)
//   extract_file  = wind_extract.csv  # terrain-aligned CSV output filename
// ==========================================================================

#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MLMG.H>
#include <AMReX_LO_BCTYPES.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_VisMF.H>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include <stdexcept>

using namespace amrex;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

// Read an X Y Z terrain file (whitespace or comma separated; '#' comments).
static void read_terrain_file(const std::string& filename,
                               std::vector<Real>& xd,
                               std::vector<Real>& yd,
                               std::vector<Real>& zd)
{
    std::ifstream f(filename);
    if (!f.is_open())
        amrex::Abort("wind_solver: cannot open terrain file: " + filename);

    std::string line;
    while (std::getline(f, line)) {
        // strip comments
        auto pos = line.find('#');
        if (pos != std::string::npos) line = line.substr(0, pos);
        // replace commas with spaces
        std::replace(line.begin(), line.end(), ',', ' ');
        std::istringstream ss(line);
        Real x, y, z;
        if (ss >> x >> y >> z) {
            xd.push_back(x);
            yd.push_back(y);
            zd.push_back(z);
        }
    }
    if (xd.empty())
        amrex::Abort("wind_solver: no data read from terrain file: " + filename);

    amrex::Print() << "wind_solver: read " << xd.size()
                   << " terrain points from " << filename << "\n";
}

// IDW interpolation: terrain height at query point (xq, yq)
// Uses k nearest data points with inverse-square-distance weights.
static Real idw_terrain(Real xq, Real yq,
                        const std::vector<Real>& x,
                        const std::vector<Real>& y,
                        const std::vector<Real>& z,
                        int k = 6)
{
    int n = static_cast<int>(x.size());
    k = std::min(k, n);

    // Squared distances to all data points
    std::vector<std::pair<Real, int>> d2(n);
    for (int i = 0; i < n; ++i) {
        Real dx = x[i] - xq;
        Real dy = y[i] - yq;
        d2[i] = {dx * dx + dy * dy, i};
    }
    // Partial sort: first k elements are the k nearest
    std::partial_sort(d2.begin(), d2.begin() + k, d2.end());

    Real wsum = 0.0, zval = 0.0;
    for (int i = 0; i < k; ++i) {
        if (d2[i].first < Real(1.0e-12)) return z[d2[i].second]; // exact hit
        Real w = Real(1.0) / d2[i].first;  // inverse-square-distance weight
        wsum += w;
        zval += w * z[d2[i].second];
    }
    return zval / wsum;
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------
int main(int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        // ----------------------------------------------------------------
        // 1. Parse user inputs
        // ----------------------------------------------------------------
        ParmParse pp;

        std::string terrain_file = "terrain.csv";
        pp.query("terrain_file", terrain_file);

        Real U_ref = 10.0;  // x-component of reference wind [m/s]
        Real V_ref =  0.0;  // y-component of reference wind [m/s]
        Real z_ref = 10.0;  // reference height above local terrain [m]
        Real z0    =  0.1;  // aerodynamic roughness length [m]
        pp.query("U_ref", U_ref);
        pp.query("V_ref", V_ref);
        pp.query("z_ref", z_ref);
        pp.query("z0",    z0);

        Real dx_req = 30.0;
        Real dy_req = 30.0;
        Real dz_req = 30.0;
        pp.query("dx", dx_req);
        pp.query("dy", dy_req);
        pp.query("dz", dz_req);

        Real domain_height = 300.0;  // [m] vertical domain extent
        pp.query("domain_height", domain_height);

        Real alpha_h = 1.0;  // horizontal Lagrange anisotropy coeff
        Real alpha_v = 1.0;  // vertical   Lagrange anisotropy coeff
        pp.query("alpha_h", alpha_h);
        pp.query("alpha_v", alpha_v);

        int  mlmg_verbose = 1;
        Real tol_rel      = 1.e-8;
        int  max_grid_size = 32;
        std::string plot_file = "plt_wind";
        pp.query("mlmg_verbose",  mlmg_verbose);
        pp.query("tol_rel",       tol_rel);
        pp.query("max_grid_size", max_grid_size);
        pp.query("plot_file",     plot_file);

        // Terrain-aligned extraction parameters
        // extract_agl  : sample at this height above local terrain [m]; snapped to
        //                the nearest cell-centre level.  Takes priority over extract_k.
        // extract_k    : sample at this k-index (0 = lowest model level).
        // Either < 0 disables that mode.  If both are < 0, no extraction is written.
        Real extract_agl = -1.0;
        int  extract_k   = -1;
        std::string extract_file = "wind_extract.csv";
        pp.query("extract_agl",  extract_agl);
        pp.query("extract_k",    extract_k);
        pp.query("extract_file", extract_file);

        // ----------------------------------------------------------------
        // 2. Read terrain file and determine horizontal domain bounds
        // ----------------------------------------------------------------
        std::vector<Real> x_terr, y_terr, z_terr;
        read_terrain_file(terrain_file, x_terr, y_terr, z_terr);

        Real x_lo = *std::min_element(x_terr.begin(), x_terr.end());
        Real x_hi = *std::max_element(x_terr.begin(), x_terr.end());
        Real y_lo = *std::min_element(y_terr.begin(), y_terr.end());
        Real y_hi = *std::max_element(y_terr.begin(), y_terr.end());

        amrex::Print() << "wind_solver: terrain x [" << x_lo << ", " << x_hi << "] m\n";
        amrex::Print() << "wind_solver: terrain y [" << y_lo << ", " << y_hi << "] m\n";

        // ----------------------------------------------------------------
        // 3. Determine grid dimensions from requested spacing
        // ----------------------------------------------------------------
        int nx = std::max(1, static_cast<int>(std::round((x_hi - x_lo) / dx_req)));
        int ny = std::max(1, static_cast<int>(std::round((y_hi - y_lo) / dy_req)));
        int nz = std::max(1, static_cast<int>(std::round(domain_height / dz_req)));

        // Actual cell sizes (may differ slightly from requested if domain
        // size is not an exact multiple of dx_req / dy_req / dz_req).
        Real dx = (x_hi - x_lo) / nx;
        Real dy = (y_hi - y_lo) / ny;
        Real dz = domain_height / nz;

        amrex::Print() << "wind_solver: grid " << nx << " x " << ny << " x " << nz
                       << "  (dx=" << dx << " m, dy=" << dy << " m, dz=" << dz << " m)\n";
        amrex::Print() << "wind_solver: vertical domain [0, " << domain_height << "] m\n";

        // ----------------------------------------------------------------
        // 4. Build AMReX geometry (Cartesian, non-periodic)
        // ----------------------------------------------------------------
        IntVect dom_lo(0, 0, 0);
        IntVect dom_hi(nx - 1, ny - 1, nz - 1);
        Box domain(dom_lo, dom_hi);

        RealBox rb({x_lo, y_lo, 0.0}, {x_hi, y_hi, domain_height});
        Array<int, AMREX_SPACEDIM> is_periodic{0, 0, 0};
        Geometry geom(domain, &rb, CoordSys::cartesian, is_periodic.data());

        BoxArray ba(domain);
        ba.maxSize(max_grid_size);
        DistributionMapping dm(ba);

        // ----------------------------------------------------------------
        // 5. Precompute per-column terrain height via IDW (host side)
        // ----------------------------------------------------------------
        // terrain_h[j*nx + i] = interpolated elevation at column (i,j) [m]
        std::vector<Real> terrain_h(static_cast<std::size_t>(nx) * ny);

        for (int j = 0; j < ny; ++j) {
            Real yc = y_lo + (j + 0.5) * dy;
            for (int i = 0; i < nx; ++i) {
                Real xc = x_lo + (i + 0.5) * dx;
                terrain_h[static_cast<std::size_t>(j) * nx + i] =
                    idw_terrain(xc, yc, x_terr, y_terr, z_terr);
            }
        }

        // Copy to device for use in GPU kernels
        Gpu::DeviceVector<Real> d_terr(terrain_h.size());
        amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                         terrain_h.begin(), terrain_h.end(), d_terr.begin());
        Real const* d_terr_ptr = d_terr.data();

        // Summary statistics
        Real zs_min = *std::min_element(terrain_h.begin(), terrain_h.end());
        Real zs_max = *std::max_element(terrain_h.begin(), terrain_h.end());
        amrex::Print() << "wind_solver: terrain elevation [" << zs_min
                       << ", " << zs_max << "] m\n";

        // ----------------------------------------------------------------
        // 6. Allocate MultiFabs
        //    vel0  – initial log-law wind (u0, v0, w0)       [3 comps, ng=1]
        //    lam   – Lagrange multiplier λ                   [1 comp,  ng=1]
        //    rhs   – Poisson RHS = -(∇·u0)                  [1 comp,  ng=0]
        // ----------------------------------------------------------------
        MultiFab vel0(ba, dm, 3, 1);
        MultiFab lam (ba, dm, 1, 1);
        MultiFab rhs (ba, dm, 1, 0);

        vel0.setVal(0.0);
        lam .setVal(0.0);
        rhs .setVal(0.0);

        // ----------------------------------------------------------------
        // 7. Fill initial log-law wind field
        // ----------------------------------------------------------------
        Real speed_ref = std::sqrt(U_ref * U_ref + V_ref * V_ref);
        const Real kappa = 0.41;  // von Karman constant

        // Compute friction velocity from reference speed and height
        // u* = κ * |U_ref| / ln((z_ref + z0) / z0)
        Real ustar = (speed_ref > Real(1.0e-10))
                   ? kappa * speed_ref / std::log((z_ref + z0) / z0)
                   : Real(0.0);

        Real ux_hat = (speed_ref > Real(1.0e-10)) ? U_ref / speed_ref : Real(1.0);
        Real uy_hat = (speed_ref > Real(1.0e-10)) ? V_ref / speed_ref : Real(0.0);

        // Capture parameters for GPU lambda
        const Real ustar_cap = ustar;
        const Real kappa_cap = kappa;
        const Real z0_cap    = z0;
        const Real ux_h      = ux_hat;
        const Real uy_h      = uy_hat;
        const Real dz_cap    = dz;
        const int  nx_cap    = nx;

        for (MFIter mfi(vel0); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            auto vel = vel0.array(mfi);

            amrex::ParallelFor(bx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real z_agl = (k + Real(0.5)) * dz_cap; // height above local terrain

                if (z_agl <= Real(0.0)) {
                    vel(i, j, k, 0) = Real(0.0);
                    vel(i, j, k, 1) = Real(0.0);
                    vel(i, j, k, 2) = Real(0.0);
                } else {
                    Real speed = (ustar_cap / kappa_cap)
                               * std::log((z_agl + z0_cap) / z0_cap);
                    vel(i, j, k, 0) = speed * ux_h;
                    vel(i, j, k, 1) = speed * uy_h;
                    vel(i, j, k, 2) = Real(0.0);
                }
            });
        }

        // Fill interior (inter-box) ghost cells via MPI exchange
        vel0.FillBoundary(geom.periodicity());

        // ----------------------------------------------------------------
        // 8. Compute divergence of initial wind  →  RHS = -(∇·u0)
        //    One-sided differences at physical domain boundaries;
        //    centred differences in the interior.
        //    Terrain (sub-surface) cells: rhs = 0 (not enforced).
        // ----------------------------------------------------------------
        const IntVect glo = domain.smallEnd();
        const IntVect ghi = domain.bigEnd();
        const int ilo = glo[0], ihi = ghi[0];
        const int jlo = glo[1], jhi = ghi[1];
        const int klo = glo[2], khi = ghi[2];
        const Real inv2dx = Real(0.5) / dx;
        const Real inv2dy = Real(0.5) / dy;
        const Real inv2dz = Real(0.5) / dz;
        const Real inv1dx = Real(1.0) / dx;
        const Real inv1dy = Real(1.0) / dy;
        const Real inv1dz = Real(1.0) / dz;

        for (MFIter mfi(rhs); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            const auto vel = vel0.const_array(mfi);
            auto rh = rhs.array(mfi);

            amrex::ParallelFor(bx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                // Height above terrain for this cell
                Real z_agl = (k + Real(0.5)) * dz_cap;
                if (z_agl <= Real(0.0)) { rh(i, j, k) = Real(0.0); return; }

                // du/dx
                Real du;
                if (i == ilo)
                    du = (vel(i+1,j,k,0) - vel(i,j,k,0)) * inv1dx;
                else if (i == ihi)
                    du = (vel(i,j,k,0) - vel(i-1,j,k,0)) * inv1dx;
                else
                    du = (vel(i+1,j,k,0) - vel(i-1,j,k,0)) * inv2dx;

                // dv/dy
                Real dv;
                if (j == jlo)
                    dv = (vel(i,j+1,k,1) - vel(i,j,k,1)) * inv1dy;
                else if (j == jhi)
                    dv = (vel(i,j,k,1) - vel(i,j-1,k,1)) * inv1dy;
                else
                    dv = (vel(i,j+1,k,1) - vel(i,j-1,k,1)) * inv2dy;

                // dw/dz
                Real dw;
                if (k == klo)
                    dw = (vel(i,j,k+1,2) - vel(i,j,k,2)) * inv1dz;
                else if (k == khi)
                    dw = (vel(i,j,k,2) - vel(i,j,k-1,2)) * inv1dz;
                else
                    dw = (vel(i,j,k+1,2) - vel(i,j,k-1,2)) * inv2dz;

                rh(i, j, k) = -(du + dv + dw);   // rhs = -div(u0)
            });
        }

        // ----------------------------------------------------------------
        // 9. Set up MLABecLaplacian and MLMG for the Poisson solve
        //
        //   Operator:  -(α_h² ∂²λ/∂x² + α_h² ∂²λ/∂y² + α_v² ∂²λ/∂z²) = rhs
        //
        //   Domain BCs:
        //     x-faces (lo, hi): Dirichlet λ = 0  (inflow / outflow)
        //     y-faces (lo, hi): Neumann ∂λ/∂y = 0 (lateral symmetry)
        //     z-faces (lo, hi): Neumann ∂λ/∂z = 0 (ground, top)
        // ----------------------------------------------------------------
        LPInfo info;
        info.setAgglomeration(true);
        info.setConsolidation(true);

        MLABecLaplacian mlabec({geom}, {ba}, {dm}, info);
        mlabec.setMaxOrder(2);

        // Boundary condition types
        Array<LinOpBCType, AMREX_SPACEDIM> lo_bc, hi_bc;
        lo_bc[0] = LinOpBCType::Dirichlet;
        hi_bc[0] = LinOpBCType::Dirichlet;
        lo_bc[1] = LinOpBCType::Neumann;
        hi_bc[1] = LinOpBCType::Neumann;
        lo_bc[2] = LinOpBCType::Neumann;
        hi_bc[2] = LinOpBCType::Neumann;
        mlabec.setDomainBC(lo_bc, hi_bc);

        // Scalars: α_a = 0 (no identity term), β_b = 1 (full diffusion)
        mlabec.setScalars(0.0, 1.0);

        // A coefficients (not used since α_a = 0, but must be set)
        MultiFab acoef(ba, dm, 1, 0);
        acoef.setVal(0.0);
        mlabec.setACoeffs(0, acoef);

        // B coefficients (face-centred, anisotropic)
        //   b_x = b_y = alpha_h², b_z = alpha_v²
        const Real bh = alpha_h * alpha_h;
        const Real bv = alpha_v * alpha_v;
        Array<MultiFab, AMREX_SPACEDIM> bcoef;
        bcoef[0].define(convert(ba, IntVect(1, 0, 0)), dm, 1, 0);
        bcoef[1].define(convert(ba, IntVect(0, 1, 0)), dm, 1, 0);
        bcoef[2].define(convert(ba, IntVect(0, 0, 1)), dm, 1, 0);
        bcoef[0].setVal(bh);
        bcoef[1].setVal(bh);
        bcoef[2].setVal(bv);
        mlabec.setBCoeffs(0, GetArrOfConstPtrs(bcoef));

        // Level BC: homogeneous (λ = 0 on Dirichlet faces)
        mlabec.setLevelBC(0, nullptr);

        // ----------------------------------------------------------------
        // 10. Solve with MLMG
        // ----------------------------------------------------------------
        MLMG mlmg(mlabec);
        mlmg.setMaxIter(200);
        mlmg.setMaxFmgIter(20);
        mlmg.setVerbose(mlmg_verbose);
        mlmg.setBottomVerbose(0);

        lam.setVal(0.0);  // initial guess

        amrex::Print() << "wind_solver: starting MLMG Poisson solve...\n";
        mlmg.solve({&lam}, {&rhs}, tol_rel, Real(0.0));
        amrex::Print() << "wind_solver: MLMG solve complete.\n";

        // Fill interior ghost cells of λ (needed for gradient computation)
        lam.FillBoundary(geom.periodicity());

        // ----------------------------------------------------------------
        // 11. Correct velocity field:  u = u0 - α_h² ∂λ/∂x  etc.
        //     One-sided gradient at physical domain boundaries.
        //     Terrain cells are reset to zero.
        // ----------------------------------------------------------------
        MultiFab vel_c(ba, dm, 3, 0);

        for (MFIter mfi(vel_c); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            const auto v0  = vel0.const_array(mfi);
            const auto la  = lam.const_array(mfi);
            auto       vc  = vel_c.array(mfi);

            amrex::ParallelFor(bx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real z_agl = (k + Real(0.5)) * dz_cap;
                if (z_agl <= Real(0.0)) {
                    vc(i, j, k, 0) = Real(0.0);
                    vc(i, j, k, 1) = Real(0.0);
                    vc(i, j, k, 2) = Real(0.0);
                    return;
                }

                // ∂λ/∂x
                Real dlx;
                if (i == ilo)
                    dlx = (la(i+1,j,k) - la(i,j,k)) * inv1dx;
                else if (i == ihi)
                    dlx = (la(i,j,k) - la(i-1,j,k)) * inv1dx;
                else
                    dlx = (la(i+1,j,k) - la(i-1,j,k)) * inv2dx;

                // ∂λ/∂y
                Real dly;
                if (j == jlo)
                    dly = (la(i,j+1,k) - la(i,j,k)) * inv1dy;
                else if (j == jhi)
                    dly = (la(i,j,k) - la(i,j-1,k)) * inv1dy;
                else
                    dly = (la(i,j+1,k) - la(i,j-1,k)) * inv2dy;

                // ∂λ/∂z
                Real dlz;
                if (k == klo)
                    dlz = (la(i,j,k+1) - la(i,j,k)) * inv1dz;
                else if (k == khi)
                    dlz = (la(i,j,k) - la(i,j,k-1)) * inv1dz;
                else
                    dlz = (la(i,j,k+1) - la(i,j,k-1)) * inv2dz;

                vc(i, j, k, 0) = v0(i, j, k, 0) - bh * dlx;
                vc(i, j, k, 1) = v0(i, j, k, 1) - bh * dly;
                vc(i, j, k, 2) = v0(i, j, k, 2) - bv * dlz;
            });
        }

        // ----------------------------------------------------------------
        // 12. Compute diagnostics: divergence before and after correction
        // ----------------------------------------------------------------
        MultiFab div_before(ba, dm, 1, 0);
        MultiFab div_after (ba, dm, 1, 0);

        // Need ghost cells for vel_c (for div_after stencil)
        MultiFab vel_c_g(ba, dm, 3, 1);
        MultiFab::Copy(vel_c_g, vel_c, 0, 0, 3, 0);
        vel_c_g.FillBoundary(geom.periodicity());

        for (MFIter mfi(div_before); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            const auto v0b = vel0.const_array(mfi);
            const auto vcg = vel_c_g.const_array(mfi);
            auto db = div_before.array(mfi);
            auto da = div_after .array(mfi);

            amrex::ParallelFor(bx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real z_agl = (k + Real(0.5)) * dz_cap;

                // --- divergence before ---
                Real du_b, dv_b, dw_b;
                if (i == ilo) du_b = (v0b(i+1,j,k,0)-v0b(i,j,k,0))*inv1dx;
                else if (i == ihi) du_b = (v0b(i,j,k,0)-v0b(i-1,j,k,0))*inv1dx;
                else du_b = (v0b(i+1,j,k,0)-v0b(i-1,j,k,0))*inv2dx;

                if (j == jlo) dv_b = (v0b(i,j+1,k,1)-v0b(i,j,k,1))*inv1dy;
                else if (j == jhi) dv_b = (v0b(i,j,k,1)-v0b(i,j-1,k,1))*inv1dy;
                else dv_b = (v0b(i,j+1,k,1)-v0b(i,j-1,k,1))*inv2dy;

                if (k == klo) dw_b = (v0b(i,j,k+1,2)-v0b(i,j,k,2))*inv1dz;
                else if (k == khi) dw_b = (v0b(i,j,k,2)-v0b(i,j,k-1,2))*inv1dz;
                else dw_b = (v0b(i,j,k+1,2)-v0b(i,j,k-1,2))*inv2dz;

                db(i,j,k) = (z_agl <= Real(0.0)) ? Real(0.0) : (du_b+dv_b+dw_b);

                // --- divergence after ---
                Real du_a, dv_a, dw_a;
                if (i == ilo) du_a = (vcg(i+1,j,k,0)-vcg(i,j,k,0))*inv1dx;
                else if (i == ihi) du_a = (vcg(i,j,k,0)-vcg(i-1,j,k,0))*inv1dx;
                else du_a = (vcg(i+1,j,k,0)-vcg(i-1,j,k,0))*inv2dx;

                if (j == jlo) dv_a = (vcg(i,j+1,k,1)-vcg(i,j,k,1))*inv1dy;
                else if (j == jhi) dv_a = (vcg(i,j,k,1)-vcg(i,j-1,k,1))*inv1dy;
                else dv_a = (vcg(i,j+1,k,1)-vcg(i,j-1,k,1))*inv2dy;

                if (k == klo) dw_a = (vcg(i,j,k+1,2)-vcg(i,j,k,2))*inv1dz;
                else if (k == khi) dw_a = (vcg(i,j,k,2)-vcg(i,j,k-1,2))*inv1dz;
                else dw_a = (vcg(i,j,k+1,2)-vcg(i,j,k-1,2))*inv2dz;

                da(i,j,k) = (z_agl <= Real(0.0)) ? Real(0.0) : (du_a+dv_a+dw_a);
            });
        }

        Real div_b_max = div_before.norm0();
        Real div_a_max = div_after .norm0();
        amrex::Print() << "wind_solver: max |div(u)| before correction = "
                       << div_b_max << " s⁻¹\n";
        amrex::Print() << "wind_solver: max |div(u)| after  correction = "
                       << div_a_max << " s⁻¹\n";

        // ----------------------------------------------------------------
        // 13. Assemble output MultiFab and write plotfile
        //
        //    Components:
        //      0  u             corrected x-wind [m/s]
        //      1  v             corrected y-wind [m/s]
        //      2  w             corrected z-wind [m/s]
        //      3  vel_magnitude |U| [m/s]
        //      4  u0            initial (log-law) x-wind [m/s]
        //      5  v0            initial (log-law) y-wind [m/s]
        //      6  w0            initial (log-law) z-wind [m/s]
        //      7  lambda        Lagrange multiplier [m²/s]
        //      8  div_before    ∇·u₀ before correction [s⁻¹]
        //      9  div_after     ∇·u  after  correction [s⁻¹]
        //     10  terrain_z     terrain elevation at column [m]
        // ----------------------------------------------------------------
        const int nout = 11;
        MultiFab output(ba, dm, nout, 0);

        for (MFIter mfi(output); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            const auto vc   = vel_c.const_array(mfi);
            const auto v0a  = vel0.const_array(mfi);
            const auto la   = lam.const_array(mfi);
            const auto dib  = div_before.const_array(mfi);
            const auto dia  = div_after.const_array(mfi);
            auto out = output.array(mfi);

            amrex::ParallelFor(bx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real u = vc(i,j,k,0), v = vc(i,j,k,1), w = vc(i,j,k,2);
                out(i,j,k, 0) = u;
                out(i,j,k, 1) = v;
                out(i,j,k, 2) = w;
                out(i,j,k, 3) = std::sqrt(u*u + v*v + w*w);
                out(i,j,k, 4) = v0a(i,j,k,0);
                out(i,j,k, 5) = v0a(i,j,k,1);
                out(i,j,k, 6) = v0a(i,j,k,2);
                out(i,j,k, 7) = la(i,j,k);
                out(i,j,k, 8) = dib(i,j,k);
                out(i,j,k, 9) = dia(i,j,k);
                out(i,j,k,10) = d_terr_ptr[j * nx_cap + i];
            });
        }

        Vector<std::string> var_names = {
            "u", "v", "w", "vel_magnitude",
            "u0", "v0", "w0",
            "lambda",
            "div_before", "div_after",
            "terrain_z"
        };

        WriteSingleLevelPlotfile(plot_file, output, var_names, geom, 0.0, 0);
        amrex::Print() << "wind_solver: plotfile written to " << plot_file << "\n";

        // ----------------------------------------------------------------
        // 14. Optional terrain-aligned extraction
        //
        //  Determine the extraction k-index:
        //   • If extract_agl >= 0: snap to the cell whose centre lies closest
        //     to the requested AGL height,
        //         k_ext = clamp( floor(extract_agl / dz), 0, nz-1 )
        //     so that cell k_ext straddles [k*dz, (k+1)*dz) with centre at
        //     (k+0.5)*dz.
        //   • Else if extract_k >= 0: use that index directly (clamped).
        //   • Otherwise skip.
        //
        //  For each horizontal column (i, j) the extracted point has:
        //     z_terrain  = interpolated terrain elevation [m]
        //     z_agl      = (k_ext + 0.5) * dz            [m above terrain]
        //     z_physical = z_terrain + z_agl              [m above sea-level]
        //  This defines a surface that follows the terrain at a constant AGL
        //  offset — the "terrain-aligned" representation.
        //
        //  Output CSV columns:
        //     x, y, z_terrain, z_physical, z_agl, u, v, w, speed
        // ----------------------------------------------------------------
        const bool do_extract = (extract_agl >= Real(0.0)) || (extract_k >= 0);

        if (do_extract) {
            // Determine k_ext
            int k_ext = -1;
            Real z_agl_ext = Real(0.0);

            if (extract_agl >= Real(0.0)) {
                // Snap requested AGL to the nearest cell-centre level
                k_ext = static_cast<int>(std::floor(extract_agl / dz));
                k_ext = std::max(0, std::min(nz - 1, k_ext));
                z_agl_ext = (k_ext + Real(0.5)) * dz;
                amrex::Print() << "wind_solver: terrain-aligned extraction at AGL = "
                               << extract_agl << " m  →  k = " << k_ext
                               << "  (cell-centre AGL = " << z_agl_ext << " m)\n";
            } else {
                k_ext = std::max(0, std::min(nz - 1, extract_k));
                z_agl_ext = (k_ext + Real(0.5)) * dz;
                amrex::Print() << "wind_solver: terrain-aligned extraction at k = "
                               << k_ext << "  (cell-centre AGL = " << z_agl_ext << " m)\n";
            }

            // Ensure all GPU work is complete before host-side data access
            amrex::Gpu::streamSynchronize();

            // Collect (x, y, z_terrain, z_physical, z_agl, u, v, w, speed)
            // per column for the k_ext level.
            // Each MPI rank collects its own portion; all ranks write
            // sequentially to produce a complete file.
            struct ExtPt {
                Real x, y, z_terrain, z_phys, z_agl_val;
                Real u, v, w, speed;
                int gi, gj;   // global cell indices for sort-order
            };
            std::vector<ExtPt> local_pts;
            local_pts.reserve(static_cast<std::size_t>(nx) * ny / 4 + 1);

            for (MFIter mfi(vel_c, false /*no tiling*/); mfi.isValid(); ++mfi) {
                const Box& bx = mfi.validbox();
                // Skip boxes that do not contain the extraction level
                if (k_ext < bx.smallEnd(2) || k_ext > bx.bigEnd(2)) continue;

                // On CPU builds const_array() returns host-accessible data.
                // On GPU builds a Gpu::streamSynchronize() above ensures the
                // data is up to date; array() here still accesses device memory,
                // so copy the slice to a host FArrayBox first.
#ifdef AMREX_USE_GPU
                Box slice_bx(IntVect(bx.smallEnd(0), bx.smallEnd(1), k_ext),
                              IntVect(bx.bigEnd(0),   bx.bigEnd(1),   k_ext));
                FArrayBox slice_fab(slice_bx, 3, The_Pinned_Arena());
                slice_fab.copy<RunOn::Device>(vel_c[mfi], slice_bx);
                amrex::Gpu::streamSynchronize();
                auto const& vc = slice_fab.const_array();
#else
                auto const& vc = vel_c.const_array(mfi);
#endif

                for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
                    for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
                        Real zs  = terrain_h[static_cast<std::size_t>(j) * nx + i];
                        Real xc  = x_lo + (i + Real(0.5)) * dx;
                        Real yc  = y_lo + (j + Real(0.5)) * dy;
                        Real u_  = vc(i, j, k_ext, 0);
                        Real v_  = vc(i, j, k_ext, 1);
                        Real w_  = vc(i, j, k_ext, 2);
                        Real spd = std::sqrt(u_*u_ + v_*v_ + w_*w_);
                        local_pts.push_back({xc, yc, zs, zs + z_agl_ext,
                                             z_agl_ext, u_, v_, w_, spd,
                                             i, j});
                    }
                }
            }

            // Sort local portion by (j, i) for reproducible output ordering
            std::sort(local_pts.begin(), local_pts.end(),
                      [](const ExtPt& a, const ExtPt& b) {
                          return (a.gj != b.gj) ? (a.gj < b.gj) : (a.gi < b.gi);
                      });

            // Sequential write: rank 0 creates the file with the header;
            // higher ranks append their portion in rank order.
            const int nranks = amrex::ParallelDescriptor::NProcs();
            const int myrank = amrex::ParallelDescriptor::MyProc();

            auto write_pts = [&](bool write_header) {
                std::ofstream outf(extract_file,
                                   write_header ? std::ios::out
                                                : std::ios::app);
                outf << std::scientific << std::setprecision(6);
                if (write_header) {
                    outf << "x,y,z_terrain,z_physical,z_agl,u,v,w,speed\n";
                }
                for (const auto& p : local_pts) {
                    outf << p.x       << ","
                         << p.y       << ","
                         << p.z_terrain << ","
                         << p.z_phys  << ","
                         << p.z_agl_val << ","
                         << p.u       << ","
                         << p.v       << ","
                         << p.w       << ","
                         << p.speed   << "\n";
                }
            };

            if (myrank == 0) {
                write_pts(true /*header*/);
            }
            for (int r = 1; r < nranks; ++r) {
                amrex::ParallelDescriptor::Barrier();
                if (myrank == r) {
                    write_pts(false /*no header — append*/);
                }
            }
            amrex::ParallelDescriptor::Barrier();

            amrex::Print() << "wind_solver: terrain-aligned extraction written to "
                           << extract_file << "  (" << (nx * ny) << " points)\n";
        }

        amrex::Print() << "wind_solver: done.\n";
    }
    amrex::Finalize();
    return 0;
}
