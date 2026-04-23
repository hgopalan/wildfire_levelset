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
        IntVect dom_lo(0, 0, 0);
        IntVect dom_hi(inputs.n_cell_x-1, inputs.n_cell_y-1, inputs.n_cell_z-1);
        Box domain(dom_lo, dom_hi);

        RealBox rb({inputs.plo_x, inputs.plo_y, inputs.plo_z},
                   {inputs.phi_x, inputs.phi_y, inputs.phi_z});
        

        Array<int,AMREX_SPACEDIM> is_periodic{0, 0, 0};
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
        init_velocity_constant(vel, geom, inputs.ux, inputs.uy, inputs.uz);

        // ---------------- dt from CFL --------------------------
        Real dt = compute_dt(vel, geom, inputs.cfl);
        amrex::Print() << "Computed dt = " << dt << "\n";


        // ---------------- Write initial plotfile ---------------
        {
            Vector<std::string> names = {"phi"};
            WriteSingleLevelPlotfile("plt0000", phi, names, geom, 0.0, 0);
        }

        // ---------------- Time stepping ------------------------
        for (int step = 1; step <= inputs.nsteps; ++step) {
            advect_levelset_weno5z_rk3 (phi, vel, geom, dt);
            dt = compute_dt(vel, geom, inputs.cfl);
            Real philomax = phi.min(0);
            Real phihimax = phi.max(0);

            amrex::Print() << "Step " << step
                        << " : phi_min = " << philomax
                        << " , phi_max = " << phihimax << "\n";
            
            #include "plot_results.H"
        }

        // ---------------- Final write --------------------------
        {
            char buf[64];
            std::snprintf(buf, sizeof(buf), "plt%04d", inputs.nsteps);
            Vector<std::string> names = {"phi"};
            WriteSingleLevelPlotfile(buf, phi, names, geom, inputs.nsteps*dt, inputs.nsteps);
            amrex::Print() << "Wrote final " << buf << "\n";
        }
    }
    amrex::Finalize();
    return 0;
}