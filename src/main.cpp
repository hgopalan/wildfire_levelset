#include <AMReX.H>
#include <AMReX_Array4.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_BoxArray.H>
#include <cmath>
#include <vector>
#include <string>
#include <memory>

using namespace amrex;
#include "numerical_schemes.H"
#include "initial_conditions.H"
#include "compute_dt.H"
#include "boundary_conditions.H"
#include "reinitialization.H"
#include "advection.H"
#include "velocity_field.H"
#include "parse_inputs.H"
#include "regrid_negative_phi.H"



// ======================= Main ================================================
int main(int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        // --- new: parse all inputs in one place
        InputParameters inputs;
        parse_inputs(inputs);

        // ---------------- Geometry setup -----------------------
#if (AMREX_SPACEDIM == 3)
        IntVect dom_lo(0, 0, 0);
        IntVect dom_hi(inputs.n_cell_x-1, inputs.n_cell_y-1, inputs.n_cell_z-1);
        Box domain(dom_lo, dom_hi);

        RealBox rb({inputs.plo_x, inputs.plo_y, inputs.plo_z},
                   {inputs.phi_x, inputs.phi_y, inputs.phi_z});
        
        Array<int,AMREX_SPACEDIM> is_periodic{0, 0, 0};
#else
        IntVect dom_lo(0, 0);
        IntVect dom_hi(inputs.n_cell_x-1, inputs.n_cell_y-1);
        Box domain(dom_lo, dom_hi);

        RealBox rb({inputs.plo_x, inputs.plo_y},
                   {inputs.phi_x, inputs.phi_y});
        
        Array<int,AMREX_SPACEDIM> is_periodic{0, 0};
#endif
        Geometry geom(domain, &rb, 0, is_periodic.data());


        // ---------------- Grids & distribution -----------------
        BoxArray ba(domain);
        ba.maxSize(inputs.max_grid);
        DistributionMapping dm(ba);

        // ---------------- Fields: phi (1 comp), vel (3 comps) --
        const int ng_phi = 3; // 3 ghost cells for stencil operations (WENO5-Z flux divergence uses up to ±3 cells)
        MultiFab phi(ba, dm, 1, ng_phi);
        MultiFab vel(ba, dm, 3, 1);


        // ---------------- Initialize ---------------------------
        if (inputs.source_type == "sphere") {
            init_phi_sphere(phi, geom, inputs.cx, inputs.cy, inputs.cz, inputs.radius);
        }
        else if(inputs.source_type == "box") {
            init_phi_box(phi, geom, inputs.cx, inputs.cy, inputs.cz, inputs.bx, inputs.by, inputs.bz);
        }
        else {
            amrex::Abort("Invalid source_type: " + inputs.source_type);
        }
        fill_boundary_extrap(phi, geom);
        
        // Initialize velocity field
#if (AMREX_SPACEDIM == 2)
        if (!inputs.velocity_file.empty()) {
            init_velocity_from_file(vel, geom, inputs.velocity_file);
        } else {
            init_velocity_constant(vel, geom, inputs.ux, inputs.uy, inputs.uz);
        }
#else
        init_velocity_constant(vel, geom, inputs.ux, inputs.uy, inputs.uz);
#endif

        // ---------------- dt from CFL --------------------------
        Real dt = compute_dt(vel, geom, inputs.cfl);
        amrex::Print() << "Computed dt = " << dt << "\n";
        Real time = 0.0;
        bool has_fine_level = false;
        BoxArray fine_ba;
        DistributionMapping fine_dm;
        std::unique_ptr<Geometry> fine_geom;
        std::unique_ptr<MultiFab> fine_phi;
        std::unique_ptr<MultiFab> fine_vel;


        // ---------------- Write initial plotfile ---------------
        {
            Vector<std::string> names = {"phi", "velx", "vely"
#if (AMREX_SPACEDIM == 3)
                , "velz"
#endif
            };
            MultiFab plotmf(ba, dm, 1 + AMREX_SPACEDIM, 0);
            MultiFab::Copy(plotmf, phi, 0, 0, 1, 0);
            MultiFab::Copy(plotmf, vel, 0, 1, AMREX_SPACEDIM, 0);
            WriteSingleLevelPlotfile("plt0000", plotmf, names, geom, 0.0, 0);
        }

        // ---------------- Time stepping ------------------------
        for (int step = 1; step <= inputs.nsteps; ++step) {
            fill_boundary_extrap(phi, geom);
            if (has_fine_level) {
                // Fill fine ghost cells from coarse data first, then apply boundary conditions
                fill_fine_ghost_from_coarse(*fine_phi, phi, *fine_geom, geom, inputs.amr_refine_ratio);
                fill_boundary_extrap(*fine_phi, *fine_geom);
            }
            const Real dt_step = dt;
            advect_levelset_weno5z_rk3 (phi, vel, geom, dt_step, inputs.rothermel);
            if (has_fine_level) {
                advect_levelset_weno5z_rk3(*fine_phi, *fine_vel, *fine_geom, dt_step, inputs.rothermel);
                synchronize_coarse_from_fine(phi, vel, *fine_phi, *fine_vel, inputs.amr_refine_ratio);
            }

            // Grid tagging is disabled in 2D mode
#if (AMREX_SPACEDIM == 3)
            if (inputs.amr_enable_negative_phi_refine == 1 &&
                inputs.amr_regrid_int > 0 &&
                (step % inputs.amr_regrid_int) == 0)
            {
                bool refined = regrid_negative_phi(phi, vel, geom,
                                                   ng_phi,
                                                   inputs.amr_refine_ratio,
                                                   inputs.amr_tag_phi_threshold,
                                                   inputs.amr_max_refinements,
                                                   inputs.max_grid,
                                                   has_fine_level,
                                                   fine_ba,
                                                   fine_dm,
                                                   fine_geom,
                                                   fine_phi,
                                                   fine_vel);
                amrex::Print() << "Regrid at step " << step << (refined ? " (refined)" : " (no refinement)") << "\n";
                if (refined) {
                    fill_boundary_extrap(phi, geom);
                    if (has_fine_level) {
                        fill_fine_ghost_from_coarse(*fine_phi, phi, *fine_geom, geom, inputs.amr_refine_ratio);
                        fill_boundary_extrap(*fine_phi, *fine_geom);
                    }
                }
                amrex::Print() << "After regrid: has_fine_level = " << has_fine_level << "\n";
            }
#endif  // AMREX_SPACEDIM == 3

            dt = compute_dt(vel, geom, inputs.cfl);
            if (has_fine_level) {
                dt = std::min(dt, compute_dt(*fine_vel, *fine_geom, inputs.cfl));
            }
            Real phi_min = phi.min(0);
            Real phi_max = phi.max(0);

            amrex::Print() << "Step " << step
                        << " : phi_min = " << phi_min
                        << " , phi_max = " << phi_max << "\n";

        if (inputs.reinit_int > 0 && (step % inputs.reinit_int == 0)) {
            amrex::Print() << "Reinitializing at step " << step << "\n";

            // --- Coarse level: dtau and iters from coarse dx ---
            {
                const auto dx = geom.CellSize();
#if (AMREX_SPACEDIM == 3)
                Real dx_min   = std::min({dx[0], dx[1], dx[2]});
#else
                Real dx_min   = std::min(dx[0], dx[1]);
#endif
                Real dtau     = 0.1 * dx_min;               // CFL = 0.5 for forward-Euler + Godunov
                int  niters   = static_cast<int>(std::ceil(ng_phi / 0.5)); // = 6 for ng_phi=3
                reinitialize_phi(phi, geom, niters, dtau);
                amrex::Print() << "Reinitialized coarse level phi\n" << " with dtau = " << dtau << " and niters = " << niters << "\n";
            }
            
            // --- Fine level: dtau and iters from fine dx (finer spacing → smaller dtau) ---
            if (has_fine_level) {
                const auto dx_f = fine_geom->CellSize();
#if (AMREX_SPACEDIM == 3)
                Real dx_min_f   = std::min({dx_f[0], dx_f[1], dx_f[2]});
#else
                Real dx_min_f   = std::min(dx_f[0], dx_f[1]);
#endif
                Real dtau_f     = 0.5 * dx_min_f;           // same CFL, but scales with fine dx
                int  niters_f   = static_cast<int>(std::ceil(ng_phi / 0.5)); // same band depth in cells
                reinitialize_phi(*fine_phi, *fine_geom, niters_f, dtau_f);
                synchronize_coarse_phi_from_fine(phi, *fine_phi, inputs.amr_refine_ratio);
            }
        }

            time += dt_step;

            if (inputs.plot_int > 0 && (step % inputs.plot_int == 0)) {
                char buf[64];
                std::snprintf(buf, sizeof(buf), "plt%04d", step);
                Vector<std::string> names = {"phi", "velx", "vely"
#if (AMREX_SPACEDIM == 3)
                    , "velz"
#endif
                };
                if (has_fine_level) {
                    MultiFab coarse_plotmf(ba, dm, 1 + AMREX_SPACEDIM, 0);
                    MultiFab::Copy(coarse_plotmf, phi, 0, 0, 1, 0);
                    MultiFab::Copy(coarse_plotmf, vel, 0, 1, AMREX_SPACEDIM, 0);
                    MultiFab fine_plotmf(fine_ba, fine_dm, 1 + AMREX_SPACEDIM, 0);
                    MultiFab::Copy(fine_plotmf, *fine_phi, 0, 0, 1, 0);
                    MultiFab::Copy(fine_plotmf, *fine_vel, 0, 1, AMREX_SPACEDIM, 0);
                    Vector<const MultiFab*> plot_data = {&coarse_plotmf, &fine_plotmf};
                    Vector<Geometry> geoms = {geom, *fine_geom};
                    Vector<int> isteps = {step, step};
                    Vector<IntVect> ref_ratio = {IntVect(AMREX_D_DECL(inputs.amr_refine_ratio,
                                                                       inputs.amr_refine_ratio,
                                                                       inputs.amr_refine_ratio))};
                    WriteMultiLevelPlotfile(buf, 2, plot_data, names, geoms, time, isteps, ref_ratio);
                } else {
                    MultiFab plotmf(ba, dm, 1 + AMREX_SPACEDIM, 0);
                    MultiFab::Copy(plotmf, phi, 0, 0, 1, 0);
                    MultiFab::Copy(plotmf, vel, 0, 1, AMREX_SPACEDIM, 0);
                    WriteSingleLevelPlotfile(buf, plotmf, names, geom, time, step);
                }
                amrex::Print() << "Wrote " << buf << "\n";
            }
        }

        // ---------------- Final write --------------------------
        {
            char buf[64];
            std::snprintf(buf, sizeof(buf), "plt%04d", inputs.nsteps);
            Vector<std::string> names = {"phi", "velx", "vely"
#if (AMREX_SPACEDIM == 3)
                , "velz"
#endif
            };
            if (has_fine_level) {
                MultiFab coarse_plotmf(ba, dm, 1 + AMREX_SPACEDIM, 0);
                MultiFab::Copy(coarse_plotmf, phi, 0, 0, 1, 0);
                MultiFab::Copy(coarse_plotmf, vel, 0, 1, AMREX_SPACEDIM, 0);
                MultiFab fine_plotmf(fine_ba, fine_dm, 1 + AMREX_SPACEDIM, 0);
                MultiFab::Copy(fine_plotmf, *fine_phi, 0, 0, 1, 0);
                MultiFab::Copy(fine_plotmf, *fine_vel, 0, 1, AMREX_SPACEDIM, 0);
                Vector<const MultiFab*> plot_data = {&coarse_plotmf, &fine_plotmf};
                Vector<Geometry> geoms = {geom, *fine_geom};
                Vector<int> isteps = {inputs.nsteps, inputs.nsteps};
                Vector<IntVect> ref_ratio = {IntVect(AMREX_D_DECL(inputs.amr_refine_ratio,
                                                                   inputs.amr_refine_ratio,
                                                                   inputs.amr_refine_ratio))};
                WriteMultiLevelPlotfile(buf, 2, plot_data, names, geoms,
                                        time, isteps, ref_ratio);
            } else {
                MultiFab plotmf(ba, dm, 1 + AMREX_SPACEDIM, 0);
                MultiFab::Copy(plotmf, phi, 0, 0, 1, 0);
                MultiFab::Copy(plotmf, vel, 0, 1, AMREX_SPACEDIM, 0);
                WriteSingleLevelPlotfile(buf, plotmf, names, geom, time, inputs.nsteps);
            }
            amrex::Print() << "Wrote final " << buf << "\n";
        }
    }
    amrex::Finalize();
    return 0;
}
