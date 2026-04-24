#include <AMReX.H>
#include <AMReX_Array4.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_GPULaunch.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_BoxArray.H>
#include <cmath>
#include <vector>
#include <string>

using namespace amrex;

// ================================================================
// Fire Spread Model Constants and Functions
// ================================================================

// Anderson (1983) L/W ratio calculation based on wind speed
// Used in FARSITE model for elliptical fire shape
static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real anderson_LW_ratio(Real wind_speed_mph)
{
    // Anderson (1983): L/W = 0.936 * exp(0.2566 * U) + 0.461 * exp(-0.1548 * U) - 0.397
    // where U is wind speed at midflame height in mph
    // For U = 0, L/W = 1.0 (circular fire)
    // This relationship is used in FARSITE
    
    if (wind_speed_mph < 0.01) {
        return Real(1.0); // Nearly circular for no wind
    }
    
    Real LW = Real(0.936) * std::exp(Real(0.2566) * wind_speed_mph) 
            + Real(0.461) * std::exp(Real(-0.1548) * wind_speed_mph) 
            - Real(0.397);
    
    // Ensure physically reasonable bounds
    return amrex::max(Real(1.0), amrex::min(LW, Real(8.0)));
}

// Rothermel slope correction factor
// Accounts for terrain slope effects on fire spread
static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real rothermel_slope_factor(Real slope_degrees)
{
    // Rothermel (1972): Slope factor φ_s
    // φ_s = 5.275 * β^(-0.3) * (tan(θ))^2
    // Simplified version using tan^2(slope) approximation
    // For typical fuel packing ratios β ≈ 0.005-0.01, the coefficient is ~10-12
    
    Real slope_rad = slope_degrees * Real(M_PI) / Real(180.0);
    Real tan_slope = std::tan(slope_rad);
    
    // Using simplified coefficient of 5.275 for typical fuels
    // In practice, this should be adjusted based on fuel bed properties
    Real phi_s = Real(5.275) * tan_slope * tan_slope;
    
    return phi_s;
}

// FARSITE combined wind and slope factor
// Accounts for both wind and terrain effects on spread rate
static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real farsite_combined_factor(Real wind_speed_mph, Real slope_degrees, 
                             Real wind_dir_rad, Real slope_aspect_rad)
{
    // FARSITE combines wind and slope effects vectorially
    // The effective wind speed accounts for terrain-induced draft
    
    // Wind factor (simplified)
    Real phi_w = wind_speed_mph > 0.01 ? 
                 std::exp(Real(0.05) * wind_speed_mph) : Real(1.0);
    
    // Slope factor from Rothermel
    Real phi_s = rothermel_slope_factor(slope_degrees);
    
    // Vectorial combination based on wind direction and slope aspect
    // If wind and slope are aligned, effects add; if opposed, they partially cancel
    Real dir_diff = wind_dir_rad - slope_aspect_rad;
    Real alignment = std::cos(dir_diff);
    
    // Combined factor (simplified formulation)
    // In full FARSITE, this involves more complex vector addition
    Real combined = Real(1.0) + phi_w + phi_s * amrex::max(Real(0.0), alignment);
    
    return amrex::max(Real(1.0), combined);
}

// ================================================================
// WENO3 one-sided reconstruction helper
// Returns the WENO3 reconstructed derivative (left or right) in one direction.
// f0, f1, f2 are three consecutive cell values; the reconstruction targets
// the interface between f1 and f2 (right-biased) or f0 and f1 (left-biased).
//
// For a right-biased (positive upwind) reconstruction:
//   stencil 0: uses f0, f1  -> q0 = f1 - f0  (divided by dx, caller does that)
//   stencil 1: uses f1, f2  -> q1 = f2 - f1
// Ideal weights: C0 = 1/3, C1 = 2/3
// ================================================================
static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real weno3_reconstruct(Real fm1, Real f0, Real fp1)
{
    // Smoothness indicators for the two 2-point stencils
    constexpr Real weno_eps = Real(1.0e-6);
    Real beta0 = (f0  - fm1) * (f0  - fm1);
    Real beta1 = (fp1 - f0)  * (fp1 - f0);

    // Ideal weights
    constexpr Real C0 = Real(1.0/3.0);
    constexpr Real C1 = Real(2.0/3.0);

    // Un-normalised weights
    Real alpha0 = C0 / ((weno_eps + beta0) * (weno_eps + beta0));
    Real alpha1 = C1 / ((weno_eps + beta1) * (weno_eps + beta1));
    Real alpha_sum = alpha0 + alpha1;

    // Normalised weights
    Real w0 = alpha0 / alpha_sum;
    Real w1 = alpha1 / alpha_sum;

    // Candidate stencil reconstructions (finite differences, not divided by dx)
    Real q0 = f0  - fm1;   // stencil 0: left pair
    Real q1 = fp1 - f0;    // stencil 1: right pair

    return w0 * q0 + w1 * q1;
}

// ================================================================
// Godunov |grad phi| operator using WENO3 reconstruction
// ================================================================
static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real godunov_norm_grad_phi(const Array4<Real const>& phi,
                           int i, int j, int k,
                           const GpuArray<Real,AMREX_SPACEDIM>& dx)
{
    // Grid-aware epsilon for regularization (prevents near-zero gradient issues)
    Real eps = Real(1.0e-8) * amrex::min(dx[0], amrex::min(dx[1], dx[2]));

    // WENO3 reconstructed upwind derivatives in each direction.
    // For each direction, two one-sided estimates are computed:
    //   px_right: Godunov "positive" (forward-biased) derivative, using stencil (i-1, i, i+1)
    //   px_left:  Godunov "negative" (backward-biased) derivative, using stencil (i-2, i-1, i)
    // The Godunov upwind selection then picks the appropriate one-sided value.

    // X direction
    Real px_right = weno3_reconstruct(phi(i-1,j,k), phi(i,j,k), phi(i+1,j,k)) / dx[0];
    Real px_left  = weno3_reconstruct(phi(i-2,j,k), phi(i-1,j,k), phi(i,j,k)) / dx[0];

    // Y direction
    Real py_right = weno3_reconstruct(phi(i,j-1,k), phi(i,j,k), phi(i,j+1,k)) / dx[1];
    Real py_left  = weno3_reconstruct(phi(i,j-2,k), phi(i,j-1,k), phi(i,j,k)) / dx[1];

    // Z direction
    Real pz_right = weno3_reconstruct(phi(i,j,k-1), phi(i,j,k), phi(i,j,k+1)) / dx[2];
    Real pz_left  = weno3_reconstruct(phi(i,j,k-2), phi(i,j,k-1), phi(i,j,k)) / dx[2];

    // Godunov upwind selection: keep only contributing sides
    Real px_p = amrex::max(px_right, Real(0.0));
    Real px_n = amrex::min(px_left,  Real(0.0));
    Real py_p = amrex::max(py_right, Real(0.0));
    Real py_n = amrex::min(py_left,  Real(0.0));
    Real pz_p = amrex::max(pz_right, Real(0.0));
    Real pz_n = amrex::min(pz_left,  Real(0.0));

    Real gx = px_p*px_p + px_n*px_n;
    Real gy = py_p*py_p + py_n*py_n;
    Real gz = pz_p*pz_p + pz_n*pz_n;

    // Add epsilon regularization to smooth the gradient magnitude
    return std::sqrt(gx + gy + gz + eps*eps);
}

// ============================================================================
// UPDATED: Build RHS with WENO5-Z using fire spread models
// Includes Anderson L/W ratio and terrain effects (Rothermel/FARSITE)
// ============================================================================
static void build_rhs_weno5z(const MultiFab& phi,
                             const MultiFab& vel,
                             const MultiFab& terrain,
                             MultiFab& rhs,
                             const Geometry& geom,
                             bool use_terrain_effects,
                             bool use_farsite_model)
{
    const auto dx = geom.CellSize();
    GpuArray<Real,AMREX_SPACEDIM> gdx{dx[0], dx[1], dx[2]};

    // Artificial viscosity ε per Eq. (2); use a single uniform value across the domain.
    // Paper default ~0.4 is common operationally.
    constexpr Real eps_visc = 0.4;

    // Conversion factor: m/s to mph for Anderson L/W ratio
    constexpr Real ms_to_mph = 2.23694;

    // Helper: discrete Laplacian Δφ (second-order)
    auto laplacian_phi = [=] AMREX_GPU_HOST_DEVICE (Array4<const Real> const& p,
                                                    int i, int j, int k) noexcept -> Real {
        Real ddx = (p(i+1,j,k) - 2.0*p(i,j,k) + p(i-1,j,k)) / (gdx[0]*gdx[0]);
        Real ddy = (p(i,j+1,k) - 2.0*p(i,j,k) + p(i,j-1,k)) / (gdx[1]*gdx[1]);
#if (AMREX_SPACEDIM==3)
        Real ddz = (p(i,j,k+1) - 2.0*p(i,j,k) + p(i,j,k-1)) / (gdx[2]*gdx[2]);
        return ddx + ddy + ddz;
#else
        return ddx + ddy;
#endif
    };

    // RHS with fire spread model corrections
    for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
        const Box& bx  = mfi.validbox();
        auto const p   = phi.const_array(mfi);
        auto const v   = vel.const_array(mfi);
        auto const t   = terrain.const_array(mfi);
        auto       fr  = rhs.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            // |∇φ| using existing Godunov/WENO5 routine
            Real grad_norm = godunov_norm_grad_phi(p, i, j, k, gdx);

            // Laplacian Δφ
            Real lap = laplacian_phi(p, i, j, k);

            // Wind velocity components
            Real ux = v(i,j,k,0);
            Real uy = v(i,j,k,1);
#if (AMREX_SPACEDIM==3)
            Real uz = v(i,j,k,2);
            Real wind_mag = std::sqrt(ux*ux + uy*uy + uz*uz);
#else
            Real uz = 0.0;
            Real wind_mag = std::sqrt(ux*ux + uy*uy);
#endif

            // Base rate of spread (ROS)
            Real R = wind_mag;

            // Apply terrain and wind effects if enabled
            if (use_terrain_effects) {
                // Terrain data: t(i,j,k,0) = slope (degrees), t(i,j,k,1) = aspect (radians)
                Real slope_deg = t(i,j,k,0);
                Real aspect_rad = t(i,j,k,1);

                // Wind direction (radians)
                Real wind_dir = std::atan2(uy, ux);

                if (use_farsite_model) {
                    // FARSITE model: Anderson L/W ratio + combined wind/slope effects
                    Real wind_mph = wind_mag * ms_to_mph;
                    Real LW_ratio = anderson_LW_ratio(wind_mph);

                    // Combined wind and slope factor
                    Real combined_factor = farsite_combined_factor(wind_mph, slope_deg, 
                                                                   wind_dir, aspect_rad);

                    // Modify ROS with L/W ratio and combined effects
                    // The L/W ratio affects the elliptical spread pattern
                    // Higher L/W means faster spread in wind direction
                    R = wind_mag * combined_factor * (Real(1.0) + (LW_ratio - Real(1.0)) * Real(0.3));

                } else {
                    // Rothermel model: slope correction only
                    Real slope_factor = rothermel_slope_factor(slope_deg);
                    
                    // Apply slope correction to base ROS
                    // Slope factor increases spread uphill, decreases downhill
                    Real dir_diff = wind_dir - aspect_rad;
                    Real uphill_component = std::cos(dir_diff);
                    
                    R = wind_mag * (Real(1.0) + slope_factor * amrex::max(Real(0.0), uphill_component));
                }
            }

            // Final RHS with fire spread model corrections
            fr(i,j,k) = - R * (grad_norm - eps_visc * lap);
        });
    }
}
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

// ======================= WENO5-Z advection + RK3 ==================
static void advect_levelset_weno5z_rk3(MultiFab& phi,
                                       const MultiFab& vel,
                                       const MultiFab& terrain,
                                       const Geometry& geom,
                                       Real dt,
                                       bool use_terrain_effects,
                                       bool use_farsite_model)
{
    const auto& ba = phi.boxArray();
    const auto& dm = phi.DistributionMap();

    // Stage buffers with 3 ghosts (needed for WENO stencils)
    MultiFab phi1(ba, dm, 1, 3);
    MultiFab phi2(ba, dm, 1, 3);
    MultiFab rhs(ba, dm, 1, 0);


    // ---------- Stage 1 ----------
    build_rhs_weno5z(phi, vel, terrain, rhs, geom, use_terrain_effects, use_farsite_model);
    for (MFIter mfi(phi, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& tbx = mfi.tilebox();
        auto const p  = phi.const_array(mfi);
        auto const r  = rhs.const_array(mfi);
        auto const q1 = phi1.array(mfi);
        ParallelFor(tbx, [=] AMREX_GPU_DEVICE (int i,int j,int k) noexcept {
            q1(i,j,k) = p(i,j,k) + dt * r(i,j,k);
        });
    }
    fill_boundary_extrap(phi1, geom);

    // ---------- Stage 2 ----------
    build_rhs_weno5z(phi1, vel, terrain, rhs, geom, use_terrain_effects, use_farsite_model);
    for (MFIter mfi(phi, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& tbx = mfi.tilebox();
        auto const p   = phi.const_array(mfi);
        auto const q1c = phi1.const_array(mfi);
        auto const r   = rhs.const_array(mfi);
        auto const q2  = phi2.array(mfi);
        ParallelFor(tbx, [=] AMREX_GPU_DEVICE (int i,int j,int k) noexcept {
            q2(i,j,k) = 0.75 * p(i,j,k) + 0.25 * ( q1c(i,j,k) + dt * r(i,j,k) );
        });
    }
    fill_boundary_extrap(phi2, geom);

    // ---------- Stage 3 ----------
    build_rhs_weno5z(phi2, vel, terrain, rhs, geom, use_terrain_effects, use_farsite_model);
    for (MFIter mfi(phi, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& tbx = mfi.tilebox();
        auto const p   = phi.const_array(mfi);
        auto const q2c = phi2.const_array(mfi);
        auto const r   = rhs.const_array(mfi);
        auto const pout= phi.array(mfi);
        ParallelFor(tbx, [=] AMREX_GPU_DEVICE (int i,int j,int k) noexcept {
           // pout(i,j,k) = (1.0/3.0) * p(i,j,k) + (2.0/3.0) * ( q2c(i,j,k) + dt * r(i,j,k) );

        Real newphi = (1.0/3.0) * p(i,j,k) + (2.0/3.0) * ( q2c(i,j,k) + dt * r(i,j,k) );

        // Clip non-positive values:
        if (newphi < 0.0) newphi = 0.0;

        pout(i,j,k) = newphi;

        });
    }
}

// ================================================================
// Smoothed sign function
// ================================================================
static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real smoothed_sign(Real phi0, Real eps)
{
    return phi0 / std::sqrt(phi0*phi0 + eps*eps);
}


// ================================================================
// Corrected reinitialization routine
// ================================================================
static void reinitialize_phi(MultiFab& phi,
                             const Geometry& geom,
                             int n_iters,
                             Real dtau)
{
    const auto dx = geom.CellSize();
    GpuArray<Real,AMREX_SPACEDIM> gdx{dx[0], dx[1], dx[2]};

    // Epsilon proportional to smallest dx
    Real eps = std::min(dx[0], std::min(dx[1], dx[2])) * 1.0;

    // phi0 must have 1 ghost cell for ∇phi operations
    MultiFab phi0(phi.boxArray(), phi.DistributionMap(), 1, 1);
    MultiFab::Copy(phi0, phi, 0, 0, 1, 1);
    fill_boundary_extrap(phi0, geom);

    // rhs needs same ghost depth as phi to allow growntilebox(1)
    MultiFab rhs(phi.boxArray(), phi.DistributionMap(), 1, 1);

    for (int it = 0; it < n_iters; ++it) {

        // Ensure ghost cells of phi are fresh each iteration using linear extrapolation
        fill_boundary_extrap(phi, geom);

        for (MFIter mfi(phi, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.growntilebox(1); // SAFE for ±1 stencil

            auto const p  = phi.const_array(mfi);
            auto const p0 = phi0.const_array(mfi);
            auto const r  = rhs.array(mfi);

            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i,int j,int k) noexcept {
                Real s = smoothed_sign(p0(i,j,k), eps);
                Real g = godunov_norm_grad_phi(p, i, j, k, gdx);
                r(i,j,k) = - s * (g - 1.0);
            });
        }

        // Update phi
        for (MFIter mfi(phi, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.growntilebox(1);

            auto const p  = phi.array(mfi);
            auto const r  = rhs.const_array(mfi);

            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i,int j,int k) noexcept {
                p(i,j,k) += dtau * r(i,j,k);
            });
        }
    }
}

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
                               (z - cz)*(z - cz)) - radius;
            p(i,j,k) = d; // negative inside the sphere
        });
    }
}

// ======================= Initialize velocity (constant) ==========================
static void init_velocity_constant(MultiFab& vel, Real ux, Real uy, Real uz)
{
    for (MFIter mfi(vel); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const v = vel.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            v(i,j,k,0) = ux;
            v(i,j,k,1) = uy;
            v(i,j,k,2) = uz;
        });
    }
}

// ======================= Initialize terrain data ==================================
// Terrain data has 2 components: 
//   comp 0: slope (degrees)
//   comp 1: aspect (radians, 0 = East, π/2 = North)
static void init_terrain(MultiFab& terrain, const Geometry& geom,
                        Real slope_degrees, Real aspect_degrees)
{
    Real aspect_rad = aspect_degrees * M_PI / 180.0;
    
    for (MFIter mfi(terrain); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const t = terrain.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            t(i,j,k,0) = slope_degrees;  // slope in degrees
            t(i,j,k,1) = aspect_rad;     // aspect in radians
        });
    }
}

// ======================= Initialize terrain from elevation data ===================
// For more complex terrain, compute slope and aspect from elevation gradients
static void init_terrain_from_elevation(MultiFab& terrain, const MultiFab& elevation,
                                        const Geometry& geom)
{
    const auto dx = geom.CellSize();
    GpuArray<Real,AMREX_SPACEDIM> gdx{dx[0], dx[1], dx[2]};
    
    for (MFIter mfi(terrain); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const t = terrain.array(mfi);
        auto const e = elevation.const_array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            // Compute elevation gradients (centered differences)
            Real dz_dx = (e(i+1,j,k) - e(i-1,j,k)) / (2.0 * gdx[0]);
            Real dz_dy = (e(i,j+1,k) - e(i,j-1,k)) / (2.0 * gdx[1]);
            
            // Slope (degrees)
            Real slope_rad = std::atan(std::sqrt(dz_dx*dz_dx + dz_dy*dz_dy));
            Real slope_deg = slope_rad * 180.0 / M_PI;
            
            // Aspect (radians): direction of steepest uphill slope
            // 0 = East, π/2 = North, π = West, 3π/2 = South
            Real aspect_rad = std::atan2(dz_dy, dz_dx);
            
            t(i,j,k,0) = slope_deg;
            t(i,j,k,1) = aspect_rad;
        });
    }
}

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

// ======================= Main ================================================
int main(int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        ParmParse pp;

        // ---------------- Inputs: grid & domain ----------------
        // You can specify either n_cell (cubic) or n_cell_x/y/z
        int n_cell      = 64;   pp.query("n_cell", n_cell);
        int n_cell_x    = n_cell; pp.query("n_cell_x", n_cell_x);
        int n_cell_y    = n_cell; pp.query("n_cell_y", n_cell_y);
        int n_cell_z    = n_cell; pp.query("n_cell_z", n_cell_z);

        int max_grid    = 32;   pp.query("max_grid_size", max_grid);

        Real plo_x = 0.0, plo_y = 0.0, plo_z = 0.0;
        Real phi_x = 1.0, phi_y = 1.0, phi_z = 1.0;
        pp.query("prob_lo_x", plo_x);
        pp.query("prob_lo_y", plo_y);
        pp.query("prob_lo_z", plo_z);
        pp.query("prob_hi_x", phi_x);
        pp.query("prob_hi_y", phi_y);
        pp.query("prob_hi_z", phi_z);

        // ---------------- Inputs: time & output ----------------
        int reinit_int = 20;
        int reinit_iters = 20;
        Real reinit_dtau = 0.5;

        pp.query("reinit_int", reinit_int);
        pp.query("reinit_iters", reinit_iters);
        pp.query("reinit_dtau", reinit_dtau);
        int nsteps    = 300;    pp.query("nsteps", nsteps);
        Real cfl      = 0.5;    pp.query("cfl", cfl);
        int plot_int  = 50;     pp.query("plot_int", plot_int);

        // ---------------- Inputs: velocity ---------------------
        Real ux = 0.25, uy = 0.0, uz = 0.0;
        pp.query("u_x", ux);
        pp.query("u_y", uy);
        pp.query("u_z", uz);

        // ---------------- Inputs: terrain ----------------------
        bool use_terrain_effects = false;
        bool use_farsite_model = false;
        Real terrain_slope = 0.0;      // slope in degrees
        Real terrain_aspect = 0.0;     // aspect in degrees (0=East, 90=North)
        
        pp.query("use_terrain_effects", use_terrain_effects);
        pp.query("use_farsite_model", use_farsite_model);
        pp.query("terrain_slope", terrain_slope);
        pp.query("terrain_aspect", terrain_aspect);

        // ---------------- Inputs: source selection ----------------
        std::string source_type = "sphere"; // "line" or "sphere"
        pp.query("source_type", source_type);

        // ---------------- Inputs: sphere -----------------------
        Real cx = 0.5, cy = 0.5, cz = 0.5, radius = 0.25;
        pp.query("sphere_center_x", cx);
        pp.query("sphere_center_y", cy);
        pp.query("sphere_center_z", cz);
        pp.query("sphere_radius",   radius);

        // ---------------- Geometry setup -----------------------
        IntVect dom_lo(0, 0, 0);
        IntVect dom_hi(n_cell_x-1, n_cell_y-1, n_cell_z-1);
        Box domain(dom_lo, dom_hi);

        RealBox rb({plo_x, plo_y, plo_z}, {phi_x, phi_y, phi_z});
        Array<int,AMREX_SPACEDIM> is_periodic{0, 0, 0};
        Geometry geom(domain, &rb, 0, is_periodic.data());

        // ---------------- Grids & distribution -----------------
        BoxArray ba(domain);
        ba.maxSize(max_grid);
        DistributionMapping dm(ba);

        // ---------------- Fields: phi (1 comp), vel (3 comps), terrain (2 comps) --
        const int ng_phi = 3; // 3 ghost cells for stencil operations (WENO5-Z flux divergence uses up to ±3 cells)
        MultiFab phi(ba, dm, 1, ng_phi);
        MultiFab vel(ba, dm, 3, 1);
        MultiFab terrain(ba, dm, 2, 1); // 2 components: slope (degrees), aspect (radians)

        // ---------------- Initialize ---------------------------
        //init_phi_sphere(phi, geom, cx, cy, cz, radius);

        // ---------------- Initialize ---------------------------
        if (source_type == "sphere") {

            // Fallback: sphere SDF (old behavior)
            init_phi_sphere(phi, geom, cx, cy, cz, radius);
        }

        init_velocity_constant(vel, ux, uy, uz);
        init_terrain(terrain, geom, terrain_slope, terrain_aspect);
        
        // Print fire spread model configuration
        amrex::Print() << "Fire spread model configuration:\n";
        if (use_terrain_effects) {
            if (use_farsite_model) {
                amrex::Print() << "  Using FARSITE model with Anderson L/W ratio\n";
                amrex::Print() << "  Terrain slope: " << terrain_slope << " degrees\n";
                amrex::Print() << "  Terrain aspect: " << terrain_aspect << " degrees\n";
            } else {
                amrex::Print() << "  Using Rothermel model with terrain effects\n";
                amrex::Print() << "  Terrain slope: " << terrain_slope << " degrees\n";
                amrex::Print() << "  Terrain aspect: " << terrain_aspect << " degrees\n";
            }
        } else {
            amrex::Print() << "  Using basic model (no terrain effects)\n";
        }

        // ---------------- dt from CFL --------------------------
        Real dt = compute_dt(vel, geom, cfl);
        amrex::Print() << "Computed dt = " << dt << "\n";

        // ---------------- Write initial plotfile ---------------
        {
            Vector<std::string> names = {"phi"};
            WriteSingleLevelPlotfile("plt0000", phi, names, geom, 0.0, 0);
        }

        // ---------------- Time stepping ------------------------
        for (int step = 1; step <= nsteps; ++step) {
            advect_levelset_weno5z_rk3(phi, vel, terrain, geom, dt, 
                                      use_terrain_effects, use_farsite_model);
            Real philomax = phi.min(0);
            Real phihimax = phi.max(0);

            amrex::Print() << "Step " << step
                        << " : phi_min = " << philomax
                        << " , phi_max = " << phihimax << "\n";

            if (reinit_int > 0 && (step % reinit_int == 0)) {
                amrex::Print() << "Reinitializing at step " << step << "\n";
                reinitialize_phi(phi, geom, reinit_iters, reinit_dtau);
            }

            if (plot_int > 0 && (step % plot_int == 0)) {
                char buf[64];
                std::snprintf(buf, sizeof(buf), "plt%04d", step);
                Vector<std::string> names = {"phi"};
                WriteSingleLevelPlotfile(buf, phi, names, geom, step*dt, step);
                amrex::Print() << "Wrote " << buf << "\n";
            }
        }

        // ---------------- Final write --------------------------
        {
            char buf[64];
            std::snprintf(buf, sizeof(buf), "plt%04d", nsteps);
            Vector<std::string> names = {"phi"};
            WriteSingleLevelPlotfile(buf, phi, names, geom, nsteps*dt, nsteps);
            amrex::Print() << "Wrote final " << buf << "\n";
        }
    }
    amrex::Finalize();
    return 0;
}