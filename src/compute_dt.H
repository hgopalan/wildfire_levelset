// ======================= Compute dt from CFL and max R =========================
static Real compute_dt(const MultiFab& R_mf, const Geometry& geom, Real cfl)
{
    // Take infinity norm of Rothermel wind speed R
    Real Rmax = R_mf.norminf(0, 0, true); // comp 0, nghost 0, global

    const auto dx = geom.CellSize();
    Real dmin = dx[0];
    dmin = std::min(dmin, dx[1]);
    dmin = std::min(dmin, dx[2]);

    if (Rmax <= 1.0e-14) {
        return 1.0; // stationary flow: generous dt
    }
    return cfl * dmin / Rmax;
}