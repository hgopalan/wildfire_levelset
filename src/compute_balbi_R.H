// ============================================================================
// Compute Balbi (2009) rate of spread and store in a MultiFab
// ============================================================================
// Mirrors compute_rothermel_R.H but uses the Balbi physical model.
// Outputs R in simulation length/time units (m/s by default).
//
// Optional spatially-varying fuel support:
//   fuel_model_mf  -- per-cell integer fuel codes (from landscape file)
//   fuel_table     -- GPU-accessible BalbiComputed array indexed by code
//   fuel_table_size -- number of entries
//
// Optional heat flux augmentation:
//   heat_flux_mf   -- per-cell heat flux [W/m2]; when provided, the per-cell
//                    buoyancy velocity is augmented by the fire-induced
//                    buoyancy component:
//                      v_b_Q  = k_upward * sqrt(g * Q * H / (rho_air * Cp_air * T_a))
//                      v_b_eff = sqrt(v_b_fuel^2 + v_b_Q^2)
//   hp             -- heat flux model parameters (rho_air, Cp_air, T_a, ref_height)
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

// Physical constant for heat flux buoyancy
static constexpr Real BALBI_HF_G = Real(9.81);

inline void compute_balbi_R(MultiFab& R_mf,
                             const MultiFab& vel,
                             const Geometry& /*geom*/,
                             const RothermelParams& rp,
                             const InputParameters::BalbiParams& bp,
                             const MultiFab* slopes = nullptr,
                             const MultiFab* fuel_model_mf = nullptr,
                             const BalbiComputed* fuel_table = nullptr,
                             int fuel_table_size = 0,
                             const MultiFab* heat_flux_mf = nullptr,
                             const InputParameters::HeatFluxParams* hp = nullptr)
{
    // Pre-compute global default Balbi coefficients
    const BalbiComputed bc_default = compute_balbi_params(rp, bp);
    const Real A_coeff_default = bc_default.A_coeff;
    const Real v_b_default     = bc_default.v_b;

    const bool use_fuel_table = (fuel_model_mf != nullptr &&
                                  fuel_table    != nullptr &&
                                  fuel_table_size > 0);

    // Heat flux augmentation of buoyancy velocity
    const bool use_heat_flux = (heat_flux_mf != nullptr && hp != nullptr &&
                                (hp->enable_upward == 1 || hp->enable_induced == 1));
    const Real hf_rho_air = use_heat_flux ? hp->rho_air    : Real(1.2);
    const Real hf_Cp_air  = use_heat_flux ? hp->Cp_air     : Real(1005.0);
    const Real hf_T_a     = use_heat_flux ? hp->T_a        : Real(300.0);
    const Real hf_refH    = use_heat_flux ? hp->ref_height : Real(10.0);
    const Real hf_k_up    = use_heat_flux ? hp->k_upward   : Real(1.0);

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

        Array4<const Real> Q_arr;
        if (use_heat_flux) Q_arr = heat_flux_mf->const_array(mfi);

        const BalbiComputed* tbl    = fuel_table;
        const int            tbl_sz = fuel_table_size;

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            // Wind components -- clamp extreme wind speeds
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

            // Heat flux augmentation of buoyancy velocity
            //   v_b_Q  = k_upward * sqrt(g * Q * H / (rho * Cp * T_a))
            //   v_b_eff = sqrt(v_b_fuel^2 + v_b_Q^2)
            if (use_heat_flux) {
                const Real Q_cell = amrex::max(Q_arr(i, j, k), Real(0.0));
                if (Q_cell > Real(0.0)) {
                    const Real denom = hf_rho_air * hf_Cp_air * hf_T_a;
                    if (denom > Real(0.0)) {
                        const Real v_b_Q = hf_k_up *
                                           std::sqrt(BALBI_HF_G * Q_cell * hf_refH / denom);
                        cell_vb = std::sqrt(cell_vb * cell_vb + v_b_Q * v_b_Q);
                    }
                }
                cell_vb = amrex::max(cell_vb, BALBI_V_B_MIN);
            }

            // Wind speed magnitude [m/s]
            const Real wind_speed = std::sqrt(ux*ux + uy*uy);

            // Flame tilt angle from wind: tan(alpha_wind) = wind_speed / v_b
            Real tan_alpha = wind_speed / cell_vb;

            // Add slope contribution: tan(alpha_slope) = sqrt(sx^2 + sy^2)
            if (use_terrain) {
                const Real sx = s(i,j,k,0);
                const Real sy = s(i,j,k,1);
                tan_alpha += std::sqrt(sx*sx + sy*sy);
            }

            // Geometric factors from flame tilt
            //   sin a = tan a / sqrt(1 + tan^2 a)
            //   cos a = 1     / sqrt(1 + tan^2 a)
            const Real denom   = std::sqrt(Real(1.0) + tan_alpha * tan_alpha);
            const Real sin_a   = tan_alpha / denom;
            const Real cos_a   = Real(1.0)  / denom;

            // Balbi (2009) ROS [m/s]:  R = A * (1 + sin a - cos a)
            R(i,j,k) = amrex::max(cell_A * (Real(1.0) + sin_a - cos_a), Real(0.0));
        });
    }
}

#endif // COMPUTE_BALBI_R_H
