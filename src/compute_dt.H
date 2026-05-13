// ======================= Compute dt from CFL and max R =========================
static Real compute_dt(const MultiFab& R_mf, const Geometry& geom, Real cfl)
{
    // Take infinity norm of Rothermel wind speed R
    Real Rmax = R_mf.norminf(0, 0, true); // comp 0, nghost 0, global

    const auto dx = geom.CellSize();
    Real dmin = dx[0];
    dmin = std::min(dmin, dx[1]);
#if (AMREX_SPACEDIM == 3)
    dmin = std::min(dmin, dx[2]);
#endif

    if (Rmax <= 1.0e-14) {
        return 1.0; // stationary flow: generous dt
    }
    return cfl * dmin / Rmax;
}

// ======================= FARSITE CFL dt (ellipse-scaled) =======================
// FARSITE propagates an elliptical wavelet whose major semi-axis grows at
// ROS * LWR per unit time.  Stability requires that no ellipse grows larger
// than one grid cell per step, giving:
//   dt = cfl * dmin / (Rmax * lwr)
// This mirrors the original FARSITE timestep criterion (Finney 1998).
static Real compute_dt(const MultiFab& R_mf, const Geometry& geom, Real cfl, Real lwr)
{
    Real Rmax = R_mf.norminf(0, 0, true); // comp 0, nghost 0, global

    const auto dx = geom.CellSize();
    Real dmin = dx[0];
    dmin = std::min(dmin, dx[1]);
#if (AMREX_SPACEDIM == 3)
    dmin = std::min(dmin, dx[2]);
#endif

    if (Rmax <= 1.0e-14) {
        return 1.0; // stationary flow: generous dt
    }
    // lwr must be >= 1: the ellipse major axis cannot be shorter than the minor axis
    const Real lwr_safe = std::max(lwr, Real(1.0));
    return cfl * dmin / (Rmax * lwr_safe);
}