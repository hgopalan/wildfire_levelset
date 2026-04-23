#include <AMReX.H>
#include <AMReX_Array4.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_GPULaunch.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_BoxArray.H>
#include <cmath>
#include <vector>
#include <string>

using namespace amrex;
#include "numerical_schemes.H"
#include "initial_conditions.H"
#include "compute_dt.H"
#include "boundary_conditions.H"
#include "reinitialization.H"
#include "advection.H"
#include "velocity_field.H"



// ======================= Main ================================================
int main(int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        ParmParse pp;

        // ---------------- Inputs: grid & domain ----------------
        // You can specify either n_cell (cubic) or n_cell_x/y/z
        int n_cell      = 64;   pp.query("n_cell", n_cell);
        int n_cell_x    = n_cell; pp.query("n_cell_x", n_cell_x);
        int n_cell_y    = n_cell; pp.query("n_cell_y", n_cell_y);
        int n_cell_z    = n_cell; pp.query("n_cell_z", n_cell_z);

        int max_grid    = 32;   pp.query("max_grid_size", max_grid);

        Real plo_x = 0.0, plo_y = 0.0, plo_z = 0.0;
        Real phi_x = 1.0, phi_y = 1.0, phi_z = 1.0;
        pp.query("prob_lo_x", plo_x);
        pp.query("prob_lo_y", plo_y);
        pp.query("prob_lo_z", plo_z);
        pp.query("prob_hi_x", phi_x);
        pp.query("prob_hi_y", phi_y);
        pp.query("prob_hi_z", phi_z);

        // ---------------- Inputs: time & output ----------------
        int reinit_int = 20;
        int reinit_iters = 20;
        Real reinit_dtau = 0.5;

        pp.query("reinit_int", reinit_int);
        pp.query("reinit_iters", reinit_iters);
        pp.query("reinit_dtau", reinit_dtau);
        int nsteps    = 300;    pp.query("nsteps", nsteps);
        Real cfl      = 0.5;    pp.query("cfl", cfl);
        int plot_int  = 50;     pp.query("plot_int", plot_int);

        // ---------------- Inputs: velocity ---------------------
        Real ux = 0.25, uy = 0.0, uz = 0.0;
        pp.query("u_x", ux);
        pp.query("u_y", uy);
        pp.query("u_z", uz);


        // ---------------- Inputs: source selection ----------------
        std::string source_type = "sphere"; // "line" or "sphere"
        pp.query("source_type", source_type);

        // ---------------- Inputs: sphere -----------------------
        Real cx = 0.5, cy = 0.5, cz = 0.5, radius = 0.25;
        pp.query("sphere_center_x", cx);
        pp.query("sphere_center_y", cy);
        pp.query("sphere_center_z", cz);
        pp.query("sphere_radius",   radius);

        // ---------------- Geometry setup -----------------------
        IntVect dom_lo(0, 0, 0);
        IntVect dom_hi(n_cell_x-1, n_cell_y-1, n_cell_z-1);
        Box domain(dom_lo, dom_hi);

        RealBox rb({plo_x, plo_y, plo_z}, {phi_x, phi_y, phi_z});
        Array<int,AMREX_SPACEDIM> is_periodic{0, 0, 0};
        Geometry geom(domain, &rb, 0, is_periodic.data());

        // ---------------- Grids & distribution -----------------
        BoxArray ba(domain);
        ba.maxSize(max_grid);
        DistributionMapping dm(ba);

        // ---------------- Fields: phi (1 comp), vel (3 comps) --
        const int ng_phi = 3; // 3 ghost cells for stencil operations (WENO5-Z flux divergence uses up to ±3 cells)
        MultiFab phi(ba, dm, 1, ng_phi);
        MultiFab vel(ba, dm, 3, 1);

        // ---------------- Initialize ---------------------------
        //init_phi_sphere(phi, geom, cx, cy, cz, radius);

        // ---------------- Initialize ---------------------------
        if (source_type == "sphere") {

            // Fallback: sphere SDF (old behavior)
            init_phi_sphere(phi, geom, cx, cy, cz, radius);
        }

        init_velocity_constant(vel, ux, uy, uz);

        // ---------------- dt from CFL --------------------------
        Real dt = compute_dt(vel, geom, cfl);
        amrex::Print() << "Computed dt = " << dt << "\n";

        // ---------------- Write initial plotfile ---------------
        {
            Vector<std::string> names = {"phi"};
            WriteSingleLevelPlotfile("plt0000", phi, names, geom, 0.0, 0);
        }

        // ---------------- Time stepping ------------------------
        for (int step = 1; step <= nsteps; ++step) {
            advect_levelset_weno5z_rk3 (phi, vel, geom, dt);
            Real philomax = phi.min(0);
            Real phihimax = phi.max(0);

            amrex::Print() << "Step " << step
                        << " : phi_min = " << philomax
                        << " , phi_max = " << phihimax << "\n";

            if (reinit_int > 0 && (step % reinit_int == 0)) {
                amrex::Print() << "Reinitializing at step " << step << "\n";
                reinitialize_phi(phi, geom, reinit_iters, reinit_dtau);
            }

            if (plot_int > 0 && (step % plot_int == 0)) {
                char buf[64];
                std::snprintf(buf, sizeof(buf), "plt%04d", step);
                Vector<std::string> names = {"phi"};
                WriteSingleLevelPlotfile(buf, phi, names, geom, step*dt, step);
                amrex::Print() << "Wrote " << buf << "\n";
            }
        }

        // ---------------- Final write --------------------------
        {
            char buf[64];
            std::snprintf(buf, sizeof(buf), "plt%04d", nsteps);
            Vector<std::string> names = {"phi"};
            WriteSingleLevelPlotfile(buf, phi, names, geom, nsteps*dt, nsteps);
            amrex::Print() << "Wrote final " << buf << "\n";
        }
    }
    amrex::Finalize();
    return 0;
}