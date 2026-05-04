// ============================================================================
// Compute Balbi (2009) rate of spread and store in a MultiFab
// ============================================================================
// Mirrors compute_rothermel_R.H but uses the Balbi physical model.
// Outputs R in simulation length/time units (m/s by default).
//
// Optional spatially-varying fuel support:
//   fuel_model_mf  – per-cell integer fuel codes (from landscape file)
//   fuel_table     – GPU-accessible BalbiComputed array indexed by code
//   fuel_table_size – number of entries
//
// Codes outside [0, fuel_table_size) fall back to the global default.
// ============================================================================

#ifndef COMPUTE_BALBI_R_H
#define COMPUTE_BALBI_R_H

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_Gpu.H>
#include "parse_inputs.H"
#include "balbi_model.H"

using namespace amrex;

static void compute_balbi_R(MultiFab& R_mf,
                             const MultiFab& vel,
                             const Geometry& /*geom*/,
                             const RothermelParams& rp,
                             const InputParameters::BalbiParams& bp,
                             const MultiFab* slopes = nullptr,
                             const MultiFab* fuel_model_mf = nullptr,
                             const BalbiComputed* fuel_table = nullptr,
                             int fuel_table_size = 0)
{
    // Pre-compute global default Balbi coefficients
    const BalbiComputed bc_default = compute_balbi_params(rp, bp);
    const Real A_coeff_default = bc_default.A_coeff;
    const Real v_b_default     = bc_default.v_b;

    const bool use_fuel_table = (fuel_model_mf != nullptr &&
                                  fuel_table    != nullptr &&
                                  fuel_table_size > 0);

    // Domain-mean wind for clamping extreme cells (same strategy as Rothermel)
    const Long ncells_vel = vel.boxArray().numPts();
    const Real ux_global_mean = (ncells_vel > 0) ? vel.sum(0) / Real(ncells_vel) : Real(0.0);
    const Real uy_global_mean = (ncells_vel > 0) ? vel.sum(1) / Real(ncells_vel) : Real(0.0);

    for (MFIter mfi(R_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const v  = vel.const_array(mfi);
        auto       R  = R_mf.array(mfi);

        const bool use_terrain = (slopes != nullptr);
        Array4<const Real> s;
        if (use_terrain) s = slopes->const_array(mfi);

        Array4<const Real> fm;
        if (use_fuel_table) fm = fuel_model_mf->const_array(mfi);

        const BalbiComputed* tbl     = fuel_table;
        const int            tbl_sz  = fuel_table_size;

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            // Wind components – clamp extreme wind speeds
            Real ux = v(i,j,k,0);
            Real uy = v(i,j,k,1);
            constexpr Real MAX_WIND_SQ = 25.0 * 25.0;  // 25 m/s threshold
            if (ux*ux + uy*uy > MAX_WIND_SQ) {
                ux = ux_global_mean;
                uy = uy_global_mean;
            }

            // Select per-cell or global Balbi coefficients
            Real cell_A  = A_coeff_default;
            Real cell_vb = v_b_default;
            if (use_fuel_table) {
                int code = static_cast<int>(fm(i, j, k));
                if (code >= 0 && code < tbl_sz) {
                    cell_A  = tbl[code].A_coeff;
                    cell_vb = tbl[code].v_b;
                }
            }

            // Wind speed magnitude [m/s]
            const Real wind_speed = std::sqrt(ux*ux + uy*uy);

            // Flame tilt angle from wind: tan(α_wind) = wind_speed / v_b
            Real tan_alpha = wind_speed / cell_vb;

            // Add slope contribution: tan(α_slope) = sqrt(sx² + sy²)
            if (use_terrain) {
                const Real sx = s(i,j,k,0);
                const Real sy = s(i,j,k,1);
                tan_alpha += std::sqrt(sx*sx + sy*sy);
            }

            // Geometric factors from flame tilt
            //   sin α = tan α / sqrt(1 + tan²α)
            //   cos α = 1     / sqrt(1 + tan²α)
            const Real denom   = std::sqrt(Real(1.0) + tan_alpha * tan_alpha);
            const Real sin_a   = tan_alpha / denom;
            const Real cos_a   = Real(1.0)  / denom;

            // Balbi (2009) ROS [m/s]:  R = A · (1 + sin α − cos α)
            R(i,j,k) = amrex::max(cell_A * (Real(1.0) + sin_a - cos_a), Real(0.0));
        });
    }
}

#endif // COMPUTE_BALBI_R_H
