// ============================================================================
// UPDATED: Build RHS with WENO5-Z using Eq. (2) with Rothermel (1972) ROS
// ============================================================================
#include "constants.H"
#include "parse_inputs.H"
#include "rothermel_model.H"

inline void build_rhs_weno5z(const MultiFab& phi,
                             const MultiFab& vel,
                             MultiFab& rhs,
                             const Geometry& geom,
                             const RothermelParams& rp,
                             const MultiFab* slopes = nullptr,
                             const MultiFab* R_precomputed = nullptr)
{
    const auto dx = geom.CellSize();
#if (AMREX_SPACEDIM == 3)
    GpuArray<Real,AMREX_SPACEDIM> gdx{dx[0], dx[1], dx[2]};
#else
    GpuArray<Real,AMREX_SPACEDIM> gdx{dx[0], dx[1]};
#endif

    // Artificial viscosity ε per Eq. (2); use a single uniform value across the domain.
    // Paper default ~0.4 is common operationally.
    constexpr Real eps_visc = WildfireConst::ADVECTION_VISC;

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

    // Compute Rothermel (1972) fire spread parameters
    const RothermelComputed rc = compute_rothermel_params(rp);
    const Real R0 = rc.R0;
    const Real C = rc.C;
    const Real B = rc.B;
    const Real beta_ratio_E = rc.beta_ratio_E;
    const Real phi_s_const = rc.phi_s;  // Constant phi_s for backward compatibility
    const Real beta = rc.beta;  // Packing ratio for terrain-based phi_s calculation
    const Real wind_conv = rc.wind_conv;
    const Real ros_conv = rc.ros_conv;
    // Slope/wind cross-term
    const int  use_sw_cross = rc.use_slope_wind_cross;
    const Real k_sw_cross   = rc.k_slope_wind_cross;
    // FARSITE-style vectorial slope/wind combination
    const int  use_sw_vectors = rc.use_slope_wind_vectors;

    // RHS per Eq. (2): ∂φ/∂t = - R ( |∇φ| - ε Δφ ), with R from Rothermel (1972)
    // or from a pre-computed R MultiFab (e.g. Balbi model).
    for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
        const Box& bx  = mfi.validbox();
        auto const p   = phi.const_array(mfi);
        auto const v   = vel.const_array(mfi);
        auto       fr  = rhs.array(mfi);
        
        // Get slopes array if terrain file was provided
        const bool use_terrain = (slopes != nullptr);
        Array4<const Real> s;
        if (use_terrain) {
            s = slopes->const_array(mfi);
        }

        // Pre-computed R field (e.g. from Balbi model)
        const bool use_R_precomp = (R_precomputed != nullptr);
        Array4<const Real> R_pre;
        if (use_R_precomp) {
            R_pre = R_precomputed->const_array(mfi);
        }

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            // Terrain slope components at this cell (zero on flat terrain)
            Real sx = use_terrain ? s(i,j,k,0) : Real(0.0);
            Real sy = use_terrain ? s(i,j,k,1) : Real(0.0);

            // |∇φ| using Godunov/WENO3 with terrain-corrected arc-length spacings
            Real grad_norm = godunov_norm_grad_phi(p, i, j, k, gdx, sx, sy);

            // Laplacian Δφ
            Real lap = laplacian_phi(p, i, j, k);

            // Rate of spread: use pre-computed field or compute Rothermel inline
            Real R;
            if (use_R_precomp) {
                // Pre-computed ROS [m/s] (from Balbi or other external model)
                R = amrex::max(R_pre(i, j, k), Real(0.0));
            } else {
                // Horizontal wind components
                Real ux = v(i,j,k,0);
                Real uy = v(i,j,k,1);

                // Wind speed in ft/min (Rothermel units)
                Real U_ftmin = std::sqrt(ux*ux + uy*uy) * wind_conv;

                // Wind factor φ_w (Rothermel Eq. 47)
                Real phi_w = (U_ftmin > 0.0) ? C * std::pow(U_ftmin, B) * beta_ratio_E : 0.0;

                // Slope factor φ_s: use terrain-based if available, otherwise constant
                Real phi_s;
                if (use_terrain) {
                    // Compute phi_s from terrain slopes at this cell
                    // Using full Rothermel Eq. 51: φ_s = 5.275 * β^(-0.3) * tan²(φ)
                    Real slope_x = s(i,j,k,0);
                    Real slope_y = s(i,j,k,1);
                    Real tan_phi_sq = slope_x * slope_x + slope_y * slope_y;
                    phi_s = std::max(Real(5.275) * std::pow(beta, Real(-0.3)) * tan_phi_sq, Real(0.0));
                } else {
                    // Use constant phi_s from RothermelComputed
                    phi_s = phi_s_const;
                }

                // Combine φ_w and φ_s using the selected mode.
                // Mode A – FARSITE vectorial (use_slope_wind_vectors = 1):
                //   φ_combined = |φ_w·wind_hat + φ_s·slope_hat|
                // Mode B – Scalar cross-term (use_slope_wind_cross = 1):
                //   φ_combined = φ_w + φ_s + k × φ_w × φ_s
                // Mode C – Default scalar additive:
                //   φ_combined = φ_w + φ_s
                Real phi_combined;
                if (use_sw_vectors == 1 && use_terrain) {
                    Real slope_x = s(i,j,k,0);
                    Real slope_y = s(i,j,k,1);
                    Real tan_phi = std::sqrt(slope_x*slope_x + slope_y*slope_y);
                    Real U_mag   = std::sqrt(ux*ux + uy*uy);
                    if (phi_w > Real(0.0) && U_mag > Real(0.0) && tan_phi > Real(0.0)) {
                        Real wind_hat_x  = ux / U_mag;
                        Real wind_hat_y  = uy / U_mag;
                        Real slope_hat_x = slope_x / tan_phi;
                        Real slope_hat_y = slope_y / tan_phi;
                        Real comb_x = phi_w * wind_hat_x + phi_s * slope_hat_x;
                        Real comb_y = phi_w * wind_hat_y + phi_s * slope_hat_y;
                        phi_combined = std::sqrt(comb_x*comb_x + comb_y*comb_y);
                    } else {
                        phi_combined = phi_w + phi_s;
                    }
                } else {
                    Real cross = (use_sw_cross == 1) ? k_sw_cross * phi_w * phi_s : Real(0.0);
                    phi_combined = phi_w + phi_s + cross;
                }

                R = std::max(R0 * (Real(1.0) + phi_combined), Real(0.0)) * ros_conv;
            }

            // Final RHS (keeps your WENO5 & RK3 integrators intact)
            fr(i,j,k) = -R * (grad_norm - eps_visc * lap);
        });
    }
}


// ======================= WENO5-Z advection + RK3 ==================
inline void advect_levelset_weno5z_rk3(MultiFab& phi,
                                       const MultiFab& vel,
                                       const Geometry& geom,
                                       Real dt,
                                       const RothermelParams& rp,
                                       const MultiFab* slopes = nullptr,
                                       const MultiFab* R_precomputed = nullptr)
{
    const auto& ba = phi.boxArray();
    const auto& dm = phi.DistributionMap();

    // Stage buffers with 3 ghosts (needed for WENO stencils)
    MultiFab phi1(ba, dm, 1, 3);
    MultiFab phi2(ba, dm, 1, 3);
    MultiFab rhs(ba, dm, 1, 0);


    // ---------- Stage 1 ----------
    build_rhs_weno5z(phi, vel, rhs, geom, rp, slopes, R_precomputed);
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
    build_rhs_weno5z(phi1, vel, rhs, geom, rp, slopes, R_precomputed);
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
    build_rhs_weno5z(phi2, vel, rhs, geom, rp, slopes, R_precomputed);
    for (MFIter mfi(phi, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& tbx = mfi.tilebox();
        auto const p   = phi.const_array(mfi);
        auto const q2c = phi2.const_array(mfi);
        auto const r   = rhs.const_array(mfi);
        auto const pout= phi.array(mfi);
        ParallelFor(tbx, [=] AMREX_GPU_DEVICE (int i,int j,int k) noexcept {
           // pout(i,j,k) = (1.0/3.0) * p(i,j,k) + (2.0/3.0) * ( q2c(i,j,k) + dt * r(i,j,k) );

        Real newphi = (1.0/3.0) * p(i,j,k) + (2.0/3.0) * ( q2c(i,j,k) + dt * r(i,j,k) );
        pout(i,j,k) = newphi;

        });
    }
}
