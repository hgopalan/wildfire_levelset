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
    const Real beta0 = (f0  - fm1) * (f0  - fm1);
    const Real beta1 = (fp1 - f0)  * (fp1 - f0);

    // Ideal weights
    constexpr Real C0 = Real(1.0/3.0);
    constexpr Real C1 = Real(2.0/3.0);

    // Un-normalised weights
    const Real alpha0 = C0 / ((weno_eps + beta0) * (weno_eps + beta0));
    const Real alpha1 = C1 / ((weno_eps + beta1) * (weno_eps + beta1));
    const Real alpha_sum = alpha0 + alpha1;

    // Normalised weights
    const Real w0 = alpha0 / alpha_sum;
    const Real w1 = alpha1 / alpha_sum;

    // Candidate stencil reconstructions (finite differences, not divided by dx)
    const Real q0 = f0  - fm1;   // stencil 0: left pair
    const Real q1 = fp1 - f0;    // stencil 1: right pair

    return w0 * q0 + w1 * q1;
}

// ================================================================
// Godunov |grad phi| operator — WENO3 with WENO1 fallback when phi is unbounded
// ================================================================
static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real godunov_norm_grad_phi(const Array4<Real const>& phi,
                           int i, int j, int k,
                           const GpuArray<Real,AMREX_SPACEDIM>& dx)
{
    // Grid-aware epsilon for regularization
    Real eps = Real(1.0e-8) * amrex::min(dx[0], amrex::min(dx[1], dx[2]));

    // ----------------------------------------------------------------
    // Detect unbounded phi: check all stencil values (WENO3 needs ±2 cells).
    // If any value is NaN/Inf or beyond the threshold, drop to WENO1.
    // ----------------------------------------------------------------
    constexpr Real phi_bound = Real(1.1);

    auto is_bad = [&](Real v) -> bool {
        return !std::isfinite(v) || std::abs(v) > phi_bound;
    };

    bool fallback =
        is_bad(phi(i-2,j,k)) || is_bad(phi(i-1,j,k)) ||
        is_bad(phi(i,  j,k)) ||
        is_bad(phi(i+1,j,k)) || is_bad(phi(i+2,j,k)) ||
        is_bad(phi(i,j-2,k)) || is_bad(phi(i,j-1,k)) ||
        is_bad(phi(i,j+1,k)) || is_bad(phi(i,j+2,k));
#if (AMREX_SPACEDIM==3)
    fallback = fallback ||
        is_bad(phi(i,j,k-2)) || is_bad(phi(i,j,k-1)) ||
        is_bad(phi(i,j,k+1)) || is_bad(phi(i,j,k+2));
#endif

    Real px_right, px_left;
    Real py_right, py_left;
#if (AMREX_SPACEDIM==3)
    Real pz_right, pz_left;
#endif

    if (fallback) {
        // ---- WENO1: plain first-order one-sided differences ----
        px_right = (phi(i+1,j,k) - phi(i,  j,k)) / dx[0];   // D+x
        px_left  = (phi(i,  j,k) - phi(i-1,j,k)) / dx[0];   // D-x
        py_right = (phi(i,j+1,k) - phi(i,j,  k)) / dx[1];   // D+y
        py_left  = (phi(i,j,  k) - phi(i,j-1,k)) / dx[1];   // D-y
#if (AMREX_SPACEDIM==3)
        pz_right = (phi(i,j,k+1) - phi(i,j,k  )) / dx[2];   // D+z
        pz_left  = (phi(i,j,k  ) - phi(i,j,k-1)) / dx[2];   // D-z
#endif
    } else {
        // ---- WENO3: blended high-order reconstruction ----
        px_right = weno3_reconstruct(phi(i-1,j,k), phi(i,j,k), phi(i+1,j,k)) / dx[0];
        px_left  = weno3_reconstruct(phi(i-2,j,k), phi(i-1,j,k), phi(i,j,k)) / dx[0];
        py_right = weno3_reconstruct(phi(i,j-1,k), phi(i,j,k), phi(i,j+1,k)) / dx[1];
        py_left  = weno3_reconstruct(phi(i,j-2,k), phi(i,j-1,k), phi(i,j,k)) / dx[1];
#if (AMREX_SPACEDIM==3)
        pz_right = weno3_reconstruct(phi(i,j,k-1), phi(i,j,k), phi(i,j,k+1)) / dx[2];
        pz_left  = weno3_reconstruct(phi(i,j,k-2), phi(i,j,k-1), phi(i,j,k)) / dx[2];
#endif
    }

    // Godunov upwind selection
    Real px_p = amrex::max(px_right, Real(0.0));
    Real px_n = amrex::min(px_left,  Real(0.0));
    Real py_p = amrex::max(py_right, Real(0.0));
    Real py_n = amrex::min(py_left,  Real(0.0));
#if (AMREX_SPACEDIM==3)
    Real pz_p = amrex::max(pz_right, Real(0.0));
    Real pz_n = amrex::min(pz_left,  Real(0.0));
#endif

    Real gx = px_p*px_p + px_n*px_n;
    Real gy = py_p*py_p + py_n*py_n;
#if (AMREX_SPACEDIM==3)
    Real gz = pz_p*pz_p + pz_n*pz_n;
    return std::sqrt(gx + gy + gz + eps*eps);
#else
    return std::sqrt(gx + gy + eps*eps);
#endif
}

