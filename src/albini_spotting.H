// ============================================================================
// Albini (1983) Firebrand Spotting Model with 2-D Trajectory Integration
// ============================================================================
// Implements firebrand spotting in two stages:
//   1. Lofting height – uses Albini's (1983) thermal-plume formula to compute
//      the maximum height H_z reached by a firebrand above the fire front,
//      driven by Byram's fire line intensity.
//   2. 2-D trajectory – treats the lofted firebrand as a particle that descends
//      at a constant terminal velocity v_t while being carried horizontally by
//      the 2-D wind field.  The horizontal position is integrated with a simple
//      forward-Euler scheme interpolating the cell-centered velocity MultiFab.
//
// References:
//   Albini (1983). "Potential spotting distance from wind-driven surface fires."
//       Research Paper INT-309, USDA Forest Service.
//   Byram (1959). "Combustion of forest fuels." In Forest Fire: Control and Use.
//
// Key formulas (SI units throughout):
//   Byram's fire line intensity:
//       I_B [kW/m] = R [ft/min] * w0 [lb/ft²] * h_heat [BTU/lb] / 60 * 3.459
//   Albini lofting height:
//       H_z [m] = 12.2 * I_B [kW/m]^(1/3)
//   Flight time:
//       t_f [s] = H_z / v_t
//   2-D trajectory (Euler, n_traj_steps sub-steps):
//       x(t + dt) = x(t) + u(x,y) * dt
//       y(t + dt) = y(t) + v(x,y) * dt
//
// Serial execution is assumed (same constraint as firebrand_spotting.H).
// ============================================================================

#ifndef ALBINI_SPOTTING_H
#define ALBINI_SPOTTING_H

#include "parse_inputs.H"
#include "rothermel_model.H"
#include <AMReX.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <random>
#include <vector>
#include <cmath>

using namespace amrex;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// BTU/ft/s -> kW/m conversion factor
static constexpr Real BTU_FTS_TO_KW_M = Real(3.459);

// ============================================================================
// Compute Byram's fire line intensity [kW/m] from Rothermel parameters
//
//   I_B = R [ft/min] * w0 [lb/ft²] * h_heat [BTU/lb] / 60  [BTU/ft/s]
//       * 3.459  [kW/m per BTU/ft/s]
//
// R_ros is the actual rate of spread in simulation length/time units (m/s).
// rp.ros_conv [m/s per ft/min] converts simulation ROS to ft/min.
// ============================================================================
static AMREX_FORCE_INLINE
Real albini_byram_intensity(Real R_ros_sim, const RothermelParams& rp)
{
    Real R_ftmin = (rp.ros_conv > Real(0.0)) ? (R_ros_sim / rp.ros_conv) : Real(0.0);
    Real I_B_BTU = R_ftmin * rp.w0 * rp.h_heat / Real(60.0);
    return I_B_BTU * BTU_FTS_TO_KW_M;
}

// ============================================================================
// Compute Albini (1983) maximum lofting height [m]
//
//   H_z = 12.2 * I_B^(1/3)   [I_B in kW/m, H_z in m]
//
// Returns 0 if I_B is below the minimum threshold.
// ============================================================================
static AMREX_FORCE_INLINE
Real albini_lofting_height(Real I_B_kW, Real I_B_min)
{
    if (I_B_kW < I_B_min) {
        return Real(0.0);
    }
    return Real(12.2) * std::pow(I_B_kW, Real(1.0) / Real(3.0));
}

// ============================================================================
// Bilinear interpolation of the wind field at an arbitrary (x, y) position.
//
// The velocity MultiFab vel stores cell-centered u (component 0) and
// v (component 1) on a 2-D grid.  Wind data are copied into flat arrays
// u_arr / v_arr (size nx * ny) and indexed as j*nx + i.
//
// Returns true if (x,y) is inside the domain; false if it has left the domain.
// ============================================================================
static bool albini_interp_wind(
    Real x, Real y,
    const std::vector<Real>& u_arr,
    const std::vector<Real>& v_arr,
    int nx, int ny,
    const Real* prob_lo,
    const Real* prob_hi,
    const Real* dx,
    Real& u_out, Real& v_out)
{
    // Check domain bounds
    if (x < prob_lo[0] || x > prob_hi[0] || y < prob_lo[1] || y > prob_hi[1]) {
        u_out = Real(0.0);
        v_out = Real(0.0);
        return false;
    }

    // Fractional cell index (cell (i,j) is centred at prob_lo + (i+0.5)*dx)
    Real xi = (x - prob_lo[0]) / dx[0] - Real(0.5);
    Real yj = (y - prob_lo[1]) / dx[1] - Real(0.5);

    int i = static_cast<int>(std::floor(xi));
    int j = static_cast<int>(std::floor(yj));
    Real fx = xi - static_cast<Real>(i);
    Real fy = yj - static_cast<Real>(j);

    // Clamp so that i+1 and j+1 remain valid
    i  = std::max(0, std::min(i, nx - 2));
    j  = std::max(0, std::min(j, ny - 2));
    fx = std::max(Real(0.0), std::min(fx, Real(1.0)));
    fy = std::max(Real(0.0), std::min(fy, Real(1.0)));

    auto idx = [nx](int ii, int jj) { return jj * nx + ii; };

    u_out = (Real(1.0) - fx) * (Real(1.0) - fy) * u_arr[idx(i,   j  )]
          + fx               * (Real(1.0) - fy) * u_arr[idx(i+1, j  )]
          + (Real(1.0) - fx) * fy               * u_arr[idx(i,   j+1)]
          + fx               * fy               * u_arr[idx(i+1, j+1)];

    v_out = (Real(1.0) - fx) * (Real(1.0) - fy) * v_arr[idx(i,   j  )]
          + fx               * (Real(1.0) - fy) * v_arr[idx(i+1, j  )]
          + (Real(1.0) - fx) * fy               * v_arr[idx(i,   j+1)]
          + fx               * fy               * v_arr[idx(i+1, j+1)];

    return true;
}

// ============================================================================
// Main Albini spotting function
//
// For every fire-front cell (phi <= phi_thresh):
//   1. Compute fire line intensity I_B via Rothermel/Byram.
//   2. Compute lofting height H_z via Albini (1983).
//   3. Stochastically decide whether a firebrand is launched (probability
//      proportional to I_B, capped at asp.P_base).
//   4. Integrate the 2-D horizontal trajectory forward-Euler until the
//      particle has descended from H_z to the ground (flight time = H_z /
//      asp.terminal_velocity).
//   5. Place a new spot-fire circle in phi at the landing location.
//
// albini_data (4 components):
//   0 – lofting height H_z [m] at fire-front source cells
//   1 – number of firebrands launched from each source cell
//   2 – maximum landing distance [m] from each source cell
//   3 – active flag (1.0) at cells that received a spot ignition
//
// NOTE: This function requires serial execution (no MPI / OMP / GPU),
//       identical to the constraint in firebrand_spotting.H.
// ============================================================================
static void compute_albini_spotting(
    MultiFab& phi,
    MultiFab& albini_data,
    const MultiFab& vel,
    const MultiFab& R_mf,        // Rothermel rate-of-spread field [m/s]
    const Geometry& geom,
    const RothermelParams& rp,
    const InputParameters::AlbiniSpottingParams& asp,
    int step)
{
    if (asp.enable == 0) {
        return;
    }

#if defined(AMREX_USE_MPI) || defined(AMREX_USE_OMP) || defined(AMREX_USE_GPU)
    amrex::Abort("Albini spotting requires serial execution. "
                 "Rebuild with AMReX_MPI=OFF, AMReX_OMP=OFF, AMReX_GPU_BACKEND=NONE.");
#endif

    const Real phi_thresh = Real(0.1);
    const auto dx        = geom.CellSize();
    const Real* prob_lo  = geom.ProbLo();
    const Real* prob_hi  = geom.ProbHi();

    const int nx = geom.Domain().length(0);
    const int ny = geom.Domain().length(1);

    // ---- Reset albini_data ----
    albini_data.setVal(Real(0.0));

    // ---- Copy wind field into flat arrays for trajectory integration ----
    std::vector<Real> u_arr(nx * ny, Real(0.0));
    std::vector<Real> v_arr(nx * ny, Real(0.0));

    for (MFIter mfi(vel); mfi.isValid(); ++mfi) {
        const Box& bx     = mfi.validbox();
        auto const v_data = vel.const_array(mfi);
        const auto lo = lbound(bx);
        const auto hi = ubound(bx);
        for (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                u_arr[j * nx + i] = v_data(i, j, 0, 0);
                v_arr[j * nx + i] = v_data(i, j, 0, 1);
            }
        }
    }

    // ---- Identify fire-front cells and compute lofting heights ----
    struct FrontCell {
        Real x, y;                // physical position of cell centre
        int  i, j;                // grid indices
        Real H_z;                 // Albini lofting height [m]
        Real I_B;                 // Byram fire line intensity [kW/m]
    };
    std::vector<FrontCell> front_cells;

    for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
        const Box& bx   = mfi.validbox();
        auto const p    = phi.const_array(mfi);
        auto const r    = R_mf.const_array(mfi);
        auto       adata = albini_data.array(mfi);

        const auto lo = lbound(bx);
        const auto hi = ubound(bx);

        for (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                if (p(i, j, 0) > phi_thresh) {
                    continue; // not near fire front
                }
                Real R_sim = r(i, j, 0); // ROS in simulation units (m/s)
                if (R_sim <= Real(0.0)) {
                    continue;
                }
                Real I_B = albini_byram_intensity(R_sim, rp);
                Real H_z = albini_lofting_height(I_B, asp.I_B_min);
                if (H_z <= Real(0.0)) {
                    continue;
                }

                // Store lofting height diagnostic
                adata(i, j, 0, 0) = H_z;

                Real x = prob_lo[0] + (i + Real(0.5)) * dx[0];
                Real y = prob_lo[1] + (j + Real(0.5)) * dx[1];
                front_cells.push_back({x, y, i, j, H_z, I_B});
            }
        }
    }

    // ---- Stochastic firebrand generation and trajectory integration ----
    std::mt19937 gen(asp.random_seed == 0
                     ? std::random_device{}()
                     : static_cast<unsigned>(asp.random_seed + step));
    std::uniform_real_distribution<Real> uniform(Real(0.0), Real(1.0));

    // Stores landing locations (x, y) and source index
    struct LandingSpot {
        Real x, y;
        int  src_i, src_j;
        int  ic, jc;      // grid indices of the landing cell
        Real distance;
    };
    std::vector<LandingSpot> spots;

    for (const auto& fc : front_cells) {
        // Intensity-weighted launch probability (capped at P_base)
        Real p_launch = asp.P_base *
            std::min(Real(1.0), fc.I_B / std::max(asp.I_B_min, Real(1.0e-10)));

        if (uniform(gen) > p_launch) {
            continue; // no firebrand from this cell
        }

        // Flight time based on lofting height and terminal velocity
        Real t_flight = fc.H_z / asp.terminal_velocity;
        Real dt_traj  = t_flight / static_cast<Real>(asp.n_traj_steps);

        // Integrate 2-D trajectory
        Real px = fc.x;
        Real py = fc.y;
        bool in_domain = true;

        for (int s = 0; s < asp.n_traj_steps; ++s) {
            Real u_w, v_w;
            in_domain = albini_interp_wind(px, py, u_arr, v_arr, nx, ny,
                                           prob_lo, prob_hi, dx,
                                           u_w, v_w);
            if (!in_domain) {
                break;
            }
            px += u_w * dt_traj;
            py += v_w * dt_traj;
        }

        if (!in_domain) {
            continue; // firebrand left the domain
        }

        Real dist = std::sqrt((px - fc.x) * (px - fc.x) +
                              (py - fc.y) * (py - fc.y));
        // Require the firebrand to have travelled at least one cell width
        if (dist < dx[0]) {
            continue;
        }

        int ic = static_cast<int>(std::floor((px - prob_lo[0]) / dx[0]));
        int jc = static_cast<int>(std::floor((py - prob_lo[1]) / dx[1]));
        spots.push_back({px, py, fc.i, fc.j, ic, jc, dist});
    }

    amrex::Print() << "Albini spotting: generated " << spots.size()
                   << " firebrands at step " << step << "\n";

    if (spots.empty()) {
        return;
    }

    // ---- Apply spot ignitions to phi ----
    int radius_cells = static_cast<int>(std::ceil(asp.spot_radius / dx[0])) + 1;

    for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto p = phi.array(mfi);
        const auto lo = lbound(bx);
        const auto hi = ubound(bx);

        for (const auto& sp : spots) {
            int ic = sp.ic;
            int jc = sp.jc;

            int i0 = std::max(ic - radius_cells, lo.x);
            int i1 = std::min(ic + radius_cells, hi.x);
            int j0 = std::max(jc - radius_cells, lo.y);
            int j1 = std::min(jc + radius_cells, hi.y);

            for (int j = j0; j <= j1; ++j) {
                for (int i = i0; i <= i1; ++i) {
                    Real cx = prob_lo[0] + (i + Real(0.5)) * dx[0];
                    Real cy = prob_lo[1] + (j + Real(0.5)) * dx[1];
                    Real d  = std::sqrt((cx - sp.x) * (cx - sp.x) +
                                        (cy - sp.y) * (cy - sp.y));
                    if (d < asp.spot_radius) {
                        p(i, j, 0) = amrex::min(p(i, j, 0), d - asp.spot_radius);
                    }
                }
            }
        }
    }

    // ---- Update albini_data diagnostics ----
    for (MFIter mfi(albini_data); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto adata = albini_data.array(mfi);
        const auto lo = lbound(bx);
        const auto hi = ubound(bx);

        for (const auto& sp : spots) {
            // Update source-cell diagnostics
            if (sp.src_i >= lo.x && sp.src_i <= hi.x &&
                sp.src_j >= lo.y && sp.src_j <= hi.y) {
                adata(sp.src_i, sp.src_j, 0, 1) += Real(1.0);            // count
                adata(sp.src_i, sp.src_j, 0, 2) =
                    amrex::max(adata(sp.src_i, sp.src_j, 0, 2), sp.distance); // max dist
            }

            // Mark landing cell as active
            int ic = sp.ic;
            int jc = sp.jc;
            if (ic >= lo.x && ic <= hi.x && jc >= lo.y && jc <= hi.y) {
                adata(ic, jc, 0, 3) = Real(1.0);
            }
        }
    }
}

#endif // ALBINI_SPOTTING_H
