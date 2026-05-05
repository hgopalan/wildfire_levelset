// ============================================================================
// Heat Flux MultiFab and Fire-Induced Wind Model
// ============================================================================
// Provides heat-flux-driven wind corrections based on WindNinja-style
// parameterisations of upward convective velocity and horizontal inflow
// induced by the fire heat release.
//
// The heat flux Q [W/m²] can be provided as:
//   (a) a uniform value (heat_flux.value > 0), or
//   (b) a spatially-varying ASCII file (heat_flux.file), X Y Q columns,
//       interpolated to the grid using the same Wendland C² IDW kernel
//       used for terrain and wind files.
//
// Fire-induced wind terms (WindNinja framework):
// -----------------------------------------------
//
//  Upward velocity (enable_upward = 1):
//    The fire plume drives a vertical velocity proportional to the
//    buoyancy of the hot gas column above the fuel bed:
//
//      w_up = k_upward * sqrt(g * Q * ref_height / (rho_air * Cp_air * T_a))
//
//    In the 2D projection this upward motion entrains horizontal inflow.
//    The horizontal inflow component is added to vel_eff as:
//
//      U_inflow = k_upward * w_up, directed anti-parallel to ambient wind
//                 (toward the fire centre)
//
//  Induced horizontal inflow (enable_induced = 1):
//    The fire-plume entrainment creates a horizontal inflow velocity
//    directed toward the fire (anti-parallel to the phi gradient in the
//    unburned region):
//
//      U_ind = k_induced * w_up
//
//    Direction: toward decreasing phi (phi gradient points outward).
//    Only applied outside the fire perimeter (phi > 0).
//
// Balbi v_b augmentation:
// ------------------------
//  When heat flux is active, the Balbi buoyancy velocity is augmented by
//  the fire-induced buoyancy:
//
//      v_b_Q  = sqrt(g * Q * ref_height / (rho_air * Cp_air * T_a))
//      v_b_eff = sqrt(v_b_fuel^2 + v_b_Q^2)
//
//  A larger v_b makes the flame more vertical (less forward tilt), reducing
//  radiant heat to unburned fuel and thereby reducing ROS – physically
//  consistent with strongly buoyant plumes above high-intensity fires.
//
// References:
//   Rothermel (1983). USDA Forest Service GTR INT-143.
//   Finney, M.A. et al. (2010). FARSITE: Fire Area Simulator.
//   WindNinja documentation: https://weather.firelab.org/windninja/
// ============================================================================

#ifndef HEAT_FLUX_MODEL_H
#define HEAT_FLUX_MODEL_H

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuLaunch.H>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "parse_inputs.H"
#include "spatial_grid.H"

using namespace amrex;

// Gravitational acceleration [m/s²] (local to this translation unit)
static constexpr Real HF_G = Real(9.81);

// ============================================================================
// init_heat_flux_from_value
//
// Fills Q_mf uniformly with the given value [W/m²].
// ============================================================================
static void init_heat_flux_from_value(MultiFab& Q_mf, Real value)
{
    Q_mf.setVal(std::max(value, Real(0.0)));
}

// ============================================================================
// init_heat_flux_from_file
//
// Reads a spatially-varying heat flux field from an ASCII file with columns
//   X  Y  Q
// and maps it to Q_mf using the same Wendland C² IDW kernel used throughout
// the codebase (terrain files, wind files, etc.).
//
// Lines starting with '#' or blank lines are ignored.
// ============================================================================
#if (AMREX_SPACEDIM == 2)
static void init_heat_flux_from_file(MultiFab& Q_mf,
                                     const Geometry& geom,
                                     const std::string& filename)
{
    std::vector<Real> x_data, y_data, q_data;

    std::ifstream infile(filename);
    if (!infile.is_open()) {
        amrex::Abort("Cannot open heat flux file: " + filename);
    }

    std::string line;
    while (std::getline(infile, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream iss(line);
        Real x, y, q;
        if (iss >> x >> y >> q) {
            x_data.push_back(x);
            y_data.push_back(y);
            q_data.push_back(std::max(q, Real(0.0)));
        }
    }
    infile.close();

    if (x_data.empty()) {
        amrex::Abort("No data read from heat flux file: " + filename);
    }

    amrex::Print() << "Read " << x_data.size()
                   << " heat-flux data points from " << filename << "\n";

    int n_points = static_cast<int>(x_data.size());

    // Characteristic spacing
    Real x_min_d = *std::min_element(x_data.begin(), x_data.end());
    Real x_max_d = *std::max_element(x_data.begin(), x_data.end());
    Real y_min_d = *std::min_element(y_data.begin(), y_data.end());
    Real y_max_d = *std::max_element(y_data.begin(), y_data.end());
    Real area    = (x_max_d - x_min_d) * (y_max_d - y_min_d);
    Real h = (area > Real(0.0) && n_points > 1)
                 ? std::sqrt(area / static_cast<Real>(n_points))
                 : Real(1.0);

    // Build cell-list spatial index
    Gpu::DeviceVector<Real> d_sx, d_sy, d_sq;
    Gpu::DeviceVector<int>  d_bstart, d_bcount;
    CellListParams cl;

    std::vector<int> sidx = build_cell_list(x_data, y_data, h,
                                            d_sx, d_sy,
                                            d_bstart, d_bcount, cl);

    std::vector<Real> q_sorted = reorder_values(q_data, sidx);
    d_sq.resize(n_points);
    Gpu::copy(Gpu::hostToDevice, q_sorted.begin(), q_sorted.end(), d_sq.begin());

    cl.sx     = d_sx.data();
    cl.sy     = d_sy.data();
    cl.bstart = d_bstart.data();
    cl.bcount = d_bcount.data();

    const Real* q_ptr = d_sq.data();

    for (MFIter mfi(Q_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto Q = Q_mf.array(mfi);

        const auto& geom_ref = geom;
        const auto dx = geom_ref.CellSizeArray();
        const auto plo = geom_ref.ProbLoArray();

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real cx = plo[0] + (i + Real(0.5)) * dx[0];
            Real cy = plo[1] + (j + Real(0.5)) * dx[1];

            Real val = cell_list_interp(cl, q_ptr, cx, cy);
            Q(i, j, k) = amrex::max(val, Real(0.0));
        });
    }
}
#endif // AMREX_SPACEDIM == 2

// ============================================================================
// apply_heatflux_wind
//
// Adds heat-flux-driven wind corrections to vel_eff.
// vel_eff should be pre-populated with the ambient or terrain-corrected
// wind field before this function is called.
//
// Parameters:
//   vel_eff   – in/out: effective wind field [m/s] (3 components)
//   vel       – input:  ambient wind [m/s]
//   Q_mf      – input:  heat flux field [W/m²]
//   phi_mf    – optional: level-set field for inflow direction (nullptr = off)
//   geom      – grid geometry (used for phi gradient, may be nullptr if phi_mf is nullptr)
//   hp        – heat flux model parameters
//
// The function handles both 2D and 3D builds:
//   3D: w_up is added directly to the z-component of vel_eff.
//   2D: upward motion is projected onto the horizontal plane as an inflow
//       directed upstream (anti-parallel to the ambient wind direction).
// ============================================================================
static void apply_heatflux_wind(MultiFab&              vel_eff,
                                const MultiFab&        vel,
                                const MultiFab&        Q_mf,
                                const MultiFab*        phi_mf,
                                const InputParameters::HeatFluxParams& hp)
{
    if (hp.enable_upward == 0 && hp.enable_induced == 0) return;

    // Parameters captured for GPU lambda
    const Real rho_air   = hp.rho_air;
    const Real Cp_air    = hp.Cp_air;
    const Real T_a       = hp.T_a;
    const Real ref_height = hp.ref_height;
    const Real k_upward   = hp.k_upward;
    const Real k_induced  = hp.k_induced;
    const int  do_upward  = hp.enable_upward;
    const int  do_induced = hp.enable_induced;

    const bool has_phi = (phi_mf != nullptr);

    // Domain-mean wind for extreme-cell clamping
    const Long ncells_vel = vel.boxArray().numPts();
    const Real ux_mean = (ncells_vel > 0) ? vel.sum(0) / Real(ncells_vel) : Real(0.0);
    const Real uy_mean = (ncells_vel > 0) ? vel.sum(1) / Real(ncells_vel) : Real(0.0);

    for (MFIter mfi(vel_eff); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const v   = vel.const_array(mfi);
        auto const Q   = Q_mf.const_array(mfi);
        auto       ve  = vel_eff.array(mfi);

        Array4<const Real> phi_arr;
        if (has_phi) phi_arr = phi_mf->const_array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            // Ambient wind, clamping extreme values
            Real ux = v(i, j, k, 0);
            Real uy = v(i, j, k, 1);
            constexpr Real MAX_WIND_SQ = Real(25.0) * Real(25.0);
            if (ux * ux + uy * uy > MAX_WIND_SQ) {
                ux = ux_mean;
                uy = uy_mean;
            }

            // Heat flux at this cell [W/m²]
            const Real Q_cell = amrex::max(Q(i, j, k), Real(0.0));

            // Upward velocity from fire buoyancy [m/s]
            // w_up = k_upward * sqrt(g * Q * H / (rho * Cp * T_a))
            Real w_up = Real(0.0);
            if (do_upward == 1 && Q_cell > Real(0.0)) {
                const Real denom = rho_air * Cp_air * T_a;
                if (denom > Real(0.0)) {
                    w_up = k_upward * std::sqrt(HF_G * Q_cell * ref_height / denom);
                }
            }

#if (AMREX_SPACEDIM == 3)
            // In 3D: add w_up directly to the z-component
            if (do_upward == 1) {
                ve(i, j, k, 2) += w_up;
            }
#else
            // In 2D: project upward plume as upstream inflow (anti-parallel to wind)
            // directed toward the fire: opposite to ambient wind direction
            if (do_upward == 1 && w_up > Real(0.0)) {
                const Real wind_spd = std::sqrt(ux * ux + uy * uy);
                if (wind_spd > Real(1.0e-8)) {
                    // Inflow is opposite to ambient wind direction
                    Real ix = -ux / wind_spd;
                    Real iy = -uy / wind_spd;
                    ve(i, j, k, 0) += w_up * ix;
                    ve(i, j, k, 1) += w_up * iy;
                }
            }
#endif

            // Induced horizontal inflow toward fire
            // Direction: anti-parallel to grad(phi) ≈ inward toward phi < 0
            if (do_induced == 1 && w_up > Real(0.0) && has_phi) {
                // Only apply outside the fire perimeter (phi >= 0)
                if (phi_arr(i, j, k) >= Real(0.0)) {
                    const Real U_ind = k_induced * w_up;
                    // grad(phi) points outward; inflow = anti-parallel to grad(phi)
                    // Use centered differences where possible, fallback to one-sided at boundaries
                    const Real gpx = phi_arr(amrex::min(i+1, bx.bigEnd(0)), j, k)
                                   - phi_arr(amrex::max(i-1, bx.smallEnd(0)), j, k);
                    const Real gpy = phi_arr(i, amrex::min(j+1, bx.bigEnd(1)), k)
                                   - phi_arr(i, amrex::max(j-1, bx.smallEnd(1)), k);
                    const Real gp_mag = std::sqrt(gpx * gpx + gpy * gpy);
                    if (gp_mag > Real(1.0e-8)) {
                        // inflow directed anti-parallel to grad(phi)
                        ve(i, j, k, 0) -= U_ind * gpx / gp_mag;
                        ve(i, j, k, 1) -= U_ind * gpy / gp_mag;
                    }
                }
            }
        });
    }
}

#endif // HEAT_FLUX_MODEL_H
