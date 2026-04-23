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
// Godunov |grad phi| operator using WENO3 reconstruction
// ================================================================
static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real godunov_norm_grad_phi(const amrex::Array4<Real const>& phi,
                           int i, int j, int k,
                           const amrex::GpuArray<Real,AMREX_SPACEDIM>& dx)
{
    // Grid-aware epsilon for regularization (prevents near-zero gradient issues)
    const Real eps = Real(1.0e-8) * amrex::min(dx[0], amrex::min(dx[1], dx[2]));

    // WENO3 reconstructed upwind derivatives in each direction.
    // For each direction, two one-sided estimates are computed:
    //   px_right: Godunov "positive" (forward-biased) derivative, using stencil (i-1, i, i+1)
    //   px_left:  Godunov "negative" (backward-biased) derivative, using stencil (i-2, i-1, i)
    // The Godunov upwind selection then picks the appropriate one-sided value.

    // X direction
    const Real px_right = weno3_reconstruct(phi(i-1,j,k), phi(i,j,k), phi(i+1,j,k)) / dx[0];
    const Real px_left  = weno3_reconstruct(phi(i-2,j,k), phi(i-1,j,k), phi(i,j,k)) / dx[0];

    // Y direction
    const Real py_right = weno3_reconstruct(phi(i,j-1,k), phi(i,j,k), phi(i,j+1,k)) / dx[1];
    const Real py_left  = weno3_reconstruct(phi(i,j-2,k), phi(i,j-1,k), phi(i,j,k)) / dx[1];

    // Z direction
    const Real pz_right = weno3_reconstruct(phi(i,j,k-1), phi(i,j,k), phi(i,j,k+1)) / dx[2];
    const Real pz_left  = weno3_reconstruct(phi(i,j,k-2), phi(i,j,k-1), phi(i,j,k)) / dx[2];

    // Godunov upwind selection: keep only contributing sides
    const Real px_p = amrex::max(px_right, Real(0.0));
    const Real px_n = amrex::min(px_left,  Real(0.0));
    const Real py_p = amrex::max(py_right, Real(0.0));
    const Real py_n = amrex::min(py_left,  Real(0.0));
    const Real pz_p = amrex::max(pz_right, Real(0.0));
    const Real pz_n = amrex::min(pz_left,  Real(0.0));

    const Real gx = px_p*px_p + px_n*px_n;
    const Real gy = py_p*py_p + py_n*py_n;
    const Real gz = pz_p*pz_p + pz_n*pz_n;

    // Add epsilon regularization to smooth the gradient magnitude
    return std::sqrt(gx + gy + gz + eps*eps);
}