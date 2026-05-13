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
//      When 3-D wind data from massconsistent_amr is provided (via
//      albini_spotting.use_3d_wind = 1 and albini_spotting.plt_wind_file),
//      the height-averaged horizontal wind from that dataset replaces the 2-D
//      MultiFab wind for trajectory integration.
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
// Serial execution is assumed for trajectory integration (same constraint as
// firebrand_spotting.H).  In GPU builds the function synchronizes the device,
// processes trajectory data on the CPU host, and copies results back to device.
//
// GPU pathway for 3-D wind interpolation (AMREX_USE_GPU):
//   When use_3d_wind = 1 and a GPU build is active, the bilinear interpolation
//   of the height-averaged plt wind onto the fire-model 2-D cell grid is
//   offloaded to the device via amrex::ParallelFor, using the device-side
//   PltWindData::d_u2d / d_v2d arrays.  The resulting wind field is copied back
//   to the host u_arr/v_arr for the serial trajectory integration step.
//   In CPU-only builds the same interpolation is performed as a plain nested
//   double loop over host arrays (CPU pathway, kept intact).
// ============================================================================

#ifndef ALBINI_SPOTTING_H
#define ALBINI_SPOTTING_H

#include "constants.H"
#include "parse_inputs.H"
#include "plt_wind_reader.H"
#include "rothermel_model.H"
#include "scott_spotting_table.H"
#include <AMReX.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <random>
#include <vector>
#include <cmath>

using namespace amrex;

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
    return I_B_BTU * WildfireConst::BTU_FTS_TO_KW_M;
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
//   5. Cap the trajectory travel distance against the Scott/Albini (1979)
//      per-fuel-model maximum spotting distance table (when fuel_model_mf
//      and fuel_system are provided).
//   6. Place a new spot-fire circle in phi at the landing location.
//
// albini_data (4 components):
//   0 – lofting height H_z [m] at fire-front source cells
//   1 – number of firebrands launched from each source cell
//   2 – maximum landing distance [m] from each source cell
//   3 – active flag (1.0) at cells that received a spot ignition
//
// When asp.use_3d_wind == 1 and plt_wind is valid, the column-averaged
// horizontal wind from the 3-D plt file is used for trajectory integration
// instead of the 2-D velocity MultiFab.
//
// GPU safety: All stochastic trajectory computation runs on the CPU host.
// In GPU builds (AMREX_USE_GPU) the function synchronizes the device stream,
// creates host-pinned copies of phi, vel, R_mf, and albini_data, then
// performs trajectory integration on the host before copying results back.
//
// GPU pathway for 3-D wind (AMREX_USE_GPU + use_3d_wind = 1):
//   The bilinear interpolation of the height-averaged plt wind onto the
//   fire-model grid is offloaded to a GPU kernel (amrex::ParallelFor) that
//   reads PltWindData::d_u2d / d_v2d (device arrays populated by
//   read_plt_wind_file).  The resulting wind field is then copied back to
//   host for the serial trajectory step.  In CPU-only builds the same
//   interpolation is a plain double loop over host u2d/v2d arrays.
//
// The function is safe to call from any build configuration (CPU/GPU).
// ============================================================================
inline void compute_albini_spotting(
    MultiFab& phi,
    MultiFab& albini_data,
    const MultiFab& vel,
    const MultiFab& R_mf,        // Rothermel rate-of-spread field [m/s]
    const Geometry& geom,
    const RothermelParams& rp,
    const InputParameters::AlbiniSpottingParams& asp,
    int step,
    const MultiFab* fuel_model_mf   = nullptr,  // per-cell fuel model code (optional)
    const std::string* fuel_system  = nullptr,   // "13" or "40" (optional)
    const PltWindData* plt_wind     = nullptr,   // 3-D wind from massconsistent_amr (optional)
    const MultiFab* residual_fuel_mf = nullptr,  // per-cell residual fuel fraction (optional)
    Real reentry_fuel_threshold      = Real(0.05)) // min fuel fraction to allow ignition
{
    if (asp.enable == 0) {
        return;
    }

    // ---- GPU safety: synchronize device work, then operate on host copies ----
    // Albini trajectory integration is inherently serial (stochastic branching,
    // dynamic landing zones) and runs on the CPU host.  In GPU builds the
    // MultiFab data lives on the device; we create explicit host-pinned copies
    // (using The_Pinned_Arena, which is host-accessible in all build configs)
    // before the CPU loops and write results back to device afterward.
    // In CPU-only builds The_Pinned_Arena allocates regular heap memory and
    // the copies are simple memcpy operations.
    amrex::Gpu::streamSynchronize();  // flush any pending GPU kernels

    // Host-pinned arena: accessible from CPU in all builds (GPU and CPU-only)
    amrex::Arena* pinned = amrex::The_Pinned_Arena();
    amrex::MFInfo host_info;
    host_info.SetArena(pinned);

    amrex::MultiFab phi_h(phi.boxArray(), phi.DistributionMap(),
                          phi.nComp(), phi.nGrow(), host_info);
    amrex::MultiFab::Copy(phi_h, phi, 0, 0, phi.nComp(), phi.nGrow());

    amrex::MultiFab vel_h(vel.boxArray(), vel.DistributionMap(),
                          vel.nComp(), vel.nGrow(), host_info);
    amrex::MultiFab::Copy(vel_h, vel, 0, 0, vel.nComp(), vel.nGrow());

    amrex::MultiFab R_h(R_mf.boxArray(), R_mf.DistributionMap(),
                        R_mf.nComp(), R_mf.nGrow(), host_info);
    amrex::MultiFab::Copy(R_h, R_mf, 0, 0, R_mf.nComp(), R_mf.nGrow());

    amrex::MultiFab adata_h(albini_data.boxArray(), albini_data.DistributionMap(),
                            albini_data.nComp(), albini_data.nGrow(), host_info);
    adata_h.setVal(Real(0.0));

    std::unique_ptr<amrex::MultiFab> fm_h_ptr;
    if (fuel_model_mf != nullptr) {
        fm_h_ptr = std::make_unique<amrex::MultiFab>(
            fuel_model_mf->boxArray(), fuel_model_mf->DistributionMap(),
            fuel_model_mf->nComp(), 0, host_info);
        amrex::MultiFab::Copy(*fm_h_ptr, *fuel_model_mf, 0, 0,
                              fuel_model_mf->nComp(), 0);
    }
    amrex::Gpu::streamSynchronize();  // wait for device→host copies

    const Real phi_thresh = Real(0.1);
    const auto dx        = geom.CellSize();
    const Real* prob_lo  = geom.ProbLo();
    const Real* prob_hi  = geom.ProbHi();

    const int nx = geom.Domain().length(0);
    const int ny = geom.Domain().length(1);

    // ---- Build flat wind arrays for trajectory integration ----
    // When use_3d_wind = 1 and plt_wind is valid, use the height-averaged
    // horizontal wind from the 3-D plt file.  Otherwise use the 2-D vel MultiFab.
    std::vector<Real> u_arr(nx * ny, Real(0.0));
    std::vector<Real> v_arr(nx * ny, Real(0.0));

    const bool use_3d = (asp.use_3d_wind == 1 && plt_wind != nullptr && plt_wind->valid);

    if (use_3d) {
        // 3-D wind: bilinearly interpolate the height-averaged u2d/v2d grid
        // from PltWindData onto the fire-model 2-D cell grid.
        const int   pw_nx = plt_wind->nx_2d;
        const int   pw_ny = plt_wind->ny_2d;
        const Real  pw_x0 = plt_wind->lo2d[0];
        const Real  pw_y0 = plt_wind->lo2d[1];
        const Real  pw_x1 = plt_wind->hi2d[0];
        const Real  pw_y1 = plt_wind->hi2d[1];
        const Real  pw_dx = (pw_nx > 0) ? (pw_x1 - pw_x0) / pw_nx : Real(1.0);
        const Real  pw_dy = (pw_ny > 0) ? (pw_y1 - pw_y0) / pw_ny : Real(1.0);

        amrex::Print() << "Albini spotting: using 3-D plt wind ("
                       << pw_nx << "x" << pw_ny << " 2-D averaged) for trajectory";

#ifdef AMREX_USE_GPU
        // ---- GPU pathway ----
        // Use device-side d_u2d/d_v2d to build the interpolated wind field
        // in parallel on the GPU, then copy the result to host u_arr/v_arr
        // for the serial stochastic trajectory integration step.
        amrex::Print() << " [GPU kernel]\n";
        if (pw_nx >= 2 && pw_ny >= 2) {
            amrex::Gpu::DeviceVector<Real> u_dev(nx * ny);
            amrex::Gpu::DeviceVector<Real> v_dev(nx * ny);

            const Real* d_u2d_ptr = plt_wind->d_u2d.data();
            const Real* d_v2d_ptr = plt_wind->d_v2d.data();
            Real* u_dev_ptr = u_dev.data();
            Real* v_dev_ptr = v_dev.data();

            // Capture all scalar values needed inside the kernel by value
            // (host pointer dereference is not allowed inside device lambdas).
            const int  g_nx    = nx;
            const int  g_pw_nx = pw_nx;
            const int  g_pw_ny = pw_ny;
            const Real g_pw_x0 = pw_x0;
            const Real g_pw_y0 = pw_y0;
            const Real g_pw_dx = pw_dx;
            const Real g_pw_dy = pw_dy;
            const Real g_lo_x  = prob_lo[0];
            const Real g_lo_y  = prob_lo[1];
            const Real g_dx0   = dx[0];
            const Real g_dx1   = dx[1];

            amrex::ParallelFor(nx * ny,
                [=] AMREX_GPU_DEVICE (int n) noexcept {
                    const int i  = n % g_nx;
                    const int j  = n / g_nx;

                    // Cell-centre physical coordinates in the fire-model domain
                    const Real xc = g_lo_x + (i + Real(0.5)) * g_dx0;
                    const Real yc = g_lo_y + (j + Real(0.5)) * g_dx1;

                    // Fractional index on the plt wind 2-D grid
                    const Real xi = (xc - g_pw_x0) / g_pw_dx - Real(0.5);
                    const Real yj = (yc - g_pw_y0) / g_pw_dy - Real(0.5);

                    const int pi = amrex::max(0, amrex::min(
                        static_cast<int>(amrex::Math::floor(xi)), g_pw_nx - 2));
                    const int pj = amrex::max(0, amrex::min(
                        static_cast<int>(amrex::Math::floor(yj)), g_pw_ny - 2));

                    const Real fx = amrex::max(Real(0.0),
                        amrex::min(xi - static_cast<Real>(pi), Real(1.0)));
                    const Real fy = amrex::max(Real(0.0),
                        amrex::min(yj - static_cast<Real>(pj), Real(1.0)));

                    const int idx00 =  pj      * g_pw_nx + pi;
                    const int idx10 =  pj      * g_pw_nx + pi + 1;
                    const int idx01 = (pj + 1) * g_pw_nx + pi;
                    const int idx11 = (pj + 1) * g_pw_nx + pi + 1;

                    u_dev_ptr[n] =
                        (Real(1.0)-fx)*(Real(1.0)-fy)*d_u2d_ptr[idx00]
                        + fx           *(Real(1.0)-fy)*d_u2d_ptr[idx10]
                        + (Real(1.0)-fx)*fy            *d_u2d_ptr[idx01]
                        + fx            *fy            *d_u2d_ptr[idx11];

                    v_dev_ptr[n] =
                        (Real(1.0)-fx)*(Real(1.0)-fy)*d_v2d_ptr[idx00]
                        + fx           *(Real(1.0)-fy)*d_v2d_ptr[idx10]
                        + (Real(1.0)-fx)*fy            *d_v2d_ptr[idx01]
                        + fx            *fy            *d_v2d_ptr[idx11];
                });
            amrex::Gpu::streamSynchronize();

            // Copy GPU-interpolated wind field back to host u_arr/v_arr so
            // the serial trajectory integration loop can use it.
            amrex::Gpu::copy(amrex::Gpu::deviceToHost,
                             u_dev.begin(), u_dev.end(), u_arr.begin());
            amrex::Gpu::copy(amrex::Gpu::deviceToHost,
                             v_dev.begin(), v_dev.end(), v_arr.begin());
            amrex::Gpu::streamSynchronize();
        }
#else
        // ---- CPU pathway (kept intact for non-GPU builds) ----
        // Bilinearly interpolate host-side u2d/v2d onto the fire-model grid.
        amrex::Print() << "\n";
        {
            const auto& u2d_h = plt_wind->u2d;
            const auto& v2d_h = plt_wind->v2d;

            for (int j = 0; j < ny; ++j) {
                for (int i = 0; i < nx; ++i) {
                    // Cell-centre physical coordinates in the fire-model domain
                    Real xc = prob_lo[0] + (i + Real(0.5)) * dx[0];
                    Real yc = prob_lo[1] + (j + Real(0.5)) * dx[1];

                    // Bilinear interpolation on the plt wind 2-D grid
                    Real xi = (xc - pw_x0) / pw_dx - Real(0.5);
                    Real yj = (yc - pw_y0) / pw_dy - Real(0.5);
                    int pi = std::max(0, std::min(static_cast<int>(std::floor(xi)), pw_nx - 2));
                    int pj = std::max(0, std::min(static_cast<int>(std::floor(yj)), pw_ny - 2));
                    Real fx = std::max(Real(0.0), std::min(xi - static_cast<Real>(pi), Real(1.0)));
                    Real fy = std::max(Real(0.0), std::min(yj - static_cast<Real>(pj), Real(1.0)));

                    auto pidx = [pw_nx](int ii, int jj){ return jj * pw_nx + ii; };
                    Real u_interp =
                        (Real(1.0)-fx)*(Real(1.0)-fy)*u2d_h[pidx(pi,  pj  )]
                       + fx           *(Real(1.0)-fy)*u2d_h[pidx(pi+1,pj  )]
                       + (Real(1.0)-fx)*fy           *u2d_h[pidx(pi,  pj+1)]
                       + fx            *fy           *u2d_h[pidx(pi+1,pj+1)];
                    Real v_interp =
                        (Real(1.0)-fx)*(Real(1.0)-fy)*v2d_h[pidx(pi,  pj  )]
                       + fx           *(Real(1.0)-fy)*v2d_h[pidx(pi+1,pj  )]
                       + (Real(1.0)-fx)*fy           *v2d_h[pidx(pi,  pj+1)]
                       + fx            *fy           *v2d_h[pidx(pi+1,pj+1)];

                    u_arr[j * nx + i] = u_interp;
                    v_arr[j * nx + i] = v_interp;
                }
            }
        }
#endif
    } else {
        // Default: copy from the 2-D velocity host copy
        for (MFIter mfi(vel_h); mfi.isValid(); ++mfi) {
            const Box& bx     = mfi.validbox();
            auto const v_data = vel_h.const_array(mfi);
            const auto lo = lbound(bx);
            const auto hi = ubound(bx);
            for (int j = lo.y; j <= hi.y; ++j) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    u_arr[j * nx + i] = v_data(i, j, 0, 0);
                    v_arr[j * nx + i] = v_data(i, j, 0, 1);
                }
            }
        }

#ifdef AMREX_USE_MPI
        // Each MFIter rank iterates over disjoint valid boxes (guaranteed by AMReX's
        // DistributionMapping), so summing across ranks assembles the full wind field
        // without double-counting.
        MPI_Allreduce(MPI_IN_PLACE, u_arr.data(), nx * ny,
                      ParallelDescriptor::Mpi_typemap<Real>::type(),
                      MPI_SUM, ParallelDescriptor::Communicator());
        MPI_Allreduce(MPI_IN_PLACE, v_arr.data(), nx * ny,
                      ParallelDescriptor::Mpi_typemap<Real>::type(),
                      MPI_SUM, ParallelDescriptor::Communicator());
#endif
    }

    // ---- Identify fire-front cells and compute lofting heights ----
    struct FrontCell {
        Real x, y;                // physical position of cell centre
        int  i, j;                // grid indices
        Real H_z;                 // Albini lofting height [m]
        Real I_B;                 // Byram fire line intensity [kW/m]
        int  fuel_model;          // FBFM fuel model code at this cell (0 = unknown)
    };
    std::vector<FrontCell> front_cells;

    const bool has_fm = (fm_h_ptr != nullptr);

    for (MFIter mfi(phi_h); mfi.isValid(); ++mfi) {
        const Box& bx   = mfi.validbox();
        auto const p    = phi_h.const_array(mfi);
        auto const r    = R_h.const_array(mfi);
        auto       ad   = adata_h.array(mfi);

        Array4<const Real> fm_arr;
        if (has_fm) fm_arr = fm_h_ptr->const_array(mfi);

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
                ad(i, j, 0, 0) = H_z;

                int fm_code = has_fm ? static_cast<int>(fm_arr(i, j, 0)) : 0;
                Real x = prob_lo[0] + (i + Real(0.5)) * dx[0];
                Real y = prob_lo[1] + (j + Real(0.5)) * dx[1];
                front_cells.push_back({x, y, i, j, H_z, I_B, fm_code});
            }
        }
    }

#ifdef AMREX_USE_MPI
    // Gather fire-front cells from all MPI ranks so every rank integrates
    // the same trajectories and applies the same spot ignitions.
    {
        int local_n = static_cast<int>(front_cells.size());
        int nprocs = ParallelDescriptor::NProcs();
        std::vector<int> counts(nprocs), displs(nprocs);
        MPI_Allgather(&local_n, 1, MPI_INT,
                      counts.data(), 1, MPI_INT,
                      ParallelDescriptor::Communicator());
        int total = 0;
        for (int r = 0; r < nprocs; ++r) { displs[r] = total; total += counts[r]; }

        // Serialise struct fields into flat arrays
        std::vector<Real> fc_x(local_n), fc_y(local_n), fc_Hz(local_n), fc_IB(local_n);
        std::vector<int>  fc_i(local_n), fc_j(local_n), fc_fm(local_n);
        for (int n = 0; n < local_n; ++n) {
            fc_x[n]  = front_cells[n].x;
            fc_y[n]  = front_cells[n].y;
            fc_Hz[n] = front_cells[n].H_z;
            fc_IB[n] = front_cells[n].I_B;
            fc_i[n]  = front_cells[n].i;
            fc_j[n]  = front_cells[n].j;
            fc_fm[n] = front_cells[n].fuel_model;
        }

        auto allgatherv_real = [&](std::vector<Real>& v) {
            std::vector<Real> all(total);
            MPI_Allgatherv(v.data(), local_n,
                           ParallelDescriptor::Mpi_typemap<Real>::type(),
                           all.data(), counts.data(), displs.data(),
                           ParallelDescriptor::Mpi_typemap<Real>::type(),
                           ParallelDescriptor::Communicator());
            v = std::move(all);
        };
        auto allgatherv_int = [&](std::vector<int>& v) {
            std::vector<int> all(total);
            MPI_Allgatherv(v.data(), local_n, MPI_INT,
                           all.data(), counts.data(), displs.data(), MPI_INT,
                           ParallelDescriptor::Communicator());
            v = std::move(all);
        };

        allgatherv_real(fc_x);  allgatherv_real(fc_y);
        allgatherv_real(fc_Hz); allgatherv_real(fc_IB);
        allgatherv_int(fc_i);   allgatherv_int(fc_j);
        allgatherv_int(fc_fm);

        front_cells.clear();
        front_cells.reserve(total);
        for (int n = 0; n < total; ++n) {
            front_cells.push_back({fc_x[n], fc_y[n], fc_i[n], fc_j[n],
                                   fc_Hz[n], fc_IB[n], fc_fm[n]});
        }
    }
#endif

    // ---- Stochastic firebrand generation and trajectory integration ----
    // Use a seed that is identical on all MPI ranks for deterministic results.
    unsigned rng_seed;
    if (asp.random_seed == 0) {
#ifdef AMREX_USE_MPI
        if (ParallelDescriptor::IOProcessor()) {
            rng_seed = std::random_device{}();
        }
        MPI_Bcast(&rng_seed, 1,
                  ParallelDescriptor::Mpi_typemap<unsigned>::type(),
                  ParallelDescriptor::IOProcessorNumber(),
                  ParallelDescriptor::Communicator());
#else
        rng_seed = std::random_device{}();
#endif
    } else {
        rng_seed = static_cast<unsigned>(asp.random_seed + step);
    }
    std::mt19937 gen(rng_seed);
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

        // ---- Scott/Albini (1979) per-fuel-model maximum spotting distance cap ----
        // When a fuel model MultiFab is available, look up the table maximum for
        // the fuel model at the source cell and discard the firebrand if the
        // integrated distance exceeds it.
        if (has_fm && fuel_system != nullptr && fc.fuel_model > 0) {
            float d_max = (*fuel_system == "40")
                          ? get_max_spot_dist_fbfm40(fc.fuel_model)
                          : get_max_spot_dist_fbfm13(fc.fuel_model);
            if (d_max > 0.0f && dist > static_cast<Real>(d_max)) {
                // Clamp landing position to the cap distance along the same direction
                Real dir_x = px - fc.x;
                Real dir_y = py - fc.y;
                Real inv_dist = Real(1.0) / dist;
                px = fc.x + dir_x * inv_dist * static_cast<Real>(d_max);
                py = fc.y + dir_y * inv_dist * static_cast<Real>(d_max);
                dist = static_cast<Real>(d_max);
            }
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

    // ---- Catching (ignition) probability filter ----
    // Each successfully landed firebrand ignites the receiving fuel only with
    // probability asp.P_catch (Anderson 1970 / Rothermel 1972 concept).
    // P_catch = 1.0 (default) → always ignite; original behaviour preserved.
    if (asp.P_catch < Real(1.0)) {
        const Real P_c = std::max(Real(0.0), std::min(Real(1.0), asp.P_catch));
        std::uniform_real_distribution<Real> catch_roll(Real(0.0), Real(1.0));
        spots.erase(
            std::remove_if(spots.begin(), spots.end(),
                           [&](const LandingSpot&) { return catch_roll(gen) >= P_c; }),
            spots.end());
        amrex::Print() << "  After P_catch=" << P_c << " filter: "
                       << spots.size() << " Albini brands caught fire\n";
    }

    // ---- Post-fire re-entry fuel adjustment ----
    // When residual_fuel_mf is provided, remove or reduce Albini spot fires
    // that land in previously-burned cells with insufficient residual fuel.
    if (residual_fuel_mf != nullptr && !spots.empty()) {
        // Make a host copy of the residual fuel MultiFab
        amrex::MultiFab rf_h(residual_fuel_mf->boxArray(),
                             residual_fuel_mf->DistributionMap(),
                             1, 0, host_info);
        amrex::MultiFab::Copy(rf_h, *residual_fuel_mf, 0, 0, 1, 0);
        amrex::Gpu::streamSynchronize();

        std::uniform_real_distribution<Real> fuel_roll(Real(0.0), Real(1.0));
        const size_t prev_size = spots.size();
        spots.erase(
            std::remove_if(spots.begin(), spots.end(),
                [&](const LandingSpot& sp) -> bool {
                    // Look up residual fuel at landing cell
                    Real f_res = Real(1.0);
                    for (MFIter mfi(rf_h); mfi.isValid(); ++mfi) {
                        const Box& bx = mfi.validbox();
                        if (bx.contains(IntVect(sp.ic, sp.jc))) {
                            f_res = rf_h.array(mfi)(sp.ic, sp.jc, 0);
                            break;
                        }
                    }
                    // Remove if below absolute threshold
                    if (f_res < reentry_fuel_threshold) return true;
                    // Scale P_ignite by residual fuel fraction
                    return (fuel_roll(gen) >= f_res);
                }),
            spots.end());
        if (spots.size() < prev_size)
            amrex::Print() << "  Re-entry fuel filter: " << spots.size()
                           << " of " << prev_size << " Albini brands in fueled area\n";
    }

    if (spots.empty()) {
        return;
    }

    // ---- Apply spot ignitions to phi ----
    int radius_cells = static_cast<int>(std::ceil(asp.spot_radius / dx[0])) + 1;

    for (MFIter mfi(phi_h); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto p = phi_h.array(mfi);
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
    for (MFIter mfi(adata_h); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto ad = adata_h.array(mfi);
        const auto lo = lbound(bx);
        const auto hi = ubound(bx);

        for (const auto& sp : spots) {
            // Update source-cell diagnostics
            if (sp.src_i >= lo.x && sp.src_i <= hi.x &&
                sp.src_j >= lo.y && sp.src_j <= hi.y) {
                ad(sp.src_i, sp.src_j, 0, 1) += Real(1.0);            // count
                ad(sp.src_i, sp.src_j, 0, 2) =
                    amrex::max(ad(sp.src_i, sp.src_j, 0, 2), sp.distance); // max dist
            }

            // Mark landing cell as active
            int ic = sp.ic;
            int jc = sp.jc;
            if (ic >= lo.x && ic <= hi.x && jc >= lo.y && jc <= hi.y) {
                ad(ic, jc, 0, 3) = Real(1.0);
            }
        }
    }

    // ---- Copy results from host back to device (no-op in CPU-only builds) ----
    amrex::MultiFab::Copy(phi,         phi_h,   0, 0, phi.nComp(),         phi.nGrow());
    amrex::MultiFab::Copy(albini_data, adata_h, 0, 0, albini_data.nComp(), 0);
    amrex::Gpu::streamSynchronize();
}

#endif // ALBINI_SPOTTING_H
