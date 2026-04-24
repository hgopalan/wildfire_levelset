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
#include "farsite_ellipse.H"
#include "terrain_slope.H"



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
        
    // FARSITE spread field: stores x,y,z displacement (3 components in 3D, 2 in 2D)
    MultiFab farsite_spread(ba, dm, AMREX_SPACEDIM, 0);


    // ---------------- Initialize ---------------------------
    // When FARSITE is enabled and level set is skipped, use indicator initialization
    // (phi = 1 inside, phi = 0 outside) instead of signed distance
    bool use_indicator = (inputs.farsite.enable == 1 && inputs.skip_levelset == 1);
        
    if (inputs.source_type == "sphere") {
      if (use_indicator) {
	init_phi_sphere_indicator(phi, geom, inputs.cx, inputs.cy, inputs.cz, inputs.radius);
      } else {
	init_phi_sphere(phi, geom, inputs.cx, inputs.cy, inputs.cz, inputs.radius);
      }
    }
    else if(inputs.source_type == "box") {
      if (use_indicator) {
	init_phi_box_indicator(phi, geom, inputs.cx, inputs.cy, inputs.cz, inputs.bx, inputs.by, inputs.bz);
      } else {
	init_phi_box(phi, geom, inputs.cx, inputs.cy, inputs.cz, inputs.bx, inputs.by, inputs.bz);
      }
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

    // Initialize terrain slopes if terrain file is provided
    std::unique_ptr<MultiFab> terrain_slopes;
    if (!inputs.rothermel.terrain_file.empty()) {
      // Create MultiFab for slopes (2 components: slope_x, slope_y)
      terrain_slopes = std::make_unique<MultiFab>(ba, dm, 2, 0);
      compute_slopes_from_terrain(*terrain_slopes, geom, inputs.rothermel.terrain_file);
      amrex::Print() << "Initialized terrain slopes from file: " 
		     << inputs.rothermel.terrain_file << "\n";
    }

    // ---------------- dt from CFL --------------------------
    Real dt=10;
    if (inputs.skip_levelset == 0)
      {
        dt = compute_dt(vel, geom, inputs.cfl);
        amrex::Print() << "Computed dt = " << dt << "\n";
      } else {
      amrex::Print() << "Skipping level set advection; using dt = " << dt << " for FARSITE spread\n";
    }
    Real time = 0.0;


    // ---------------- Write initial plotfile ---------------
    {
      Vector<std::string> names = {"phi", "velx", "vely"
#if (AMREX_SPACEDIM == 3)
				   , "velz", "farsite_dx", "farsite_dy", "farsite_dz"
#else
				   , "farsite_dx", "farsite_dy"
#endif
      };
      MultiFab plotmf(ba, dm, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM, 0);
      MultiFab::Copy(plotmf, phi, 0, 0, 1, 0);
      MultiFab::Copy(plotmf, vel, 0, 1, AMREX_SPACEDIM, 0);
      MultiFab::Copy(plotmf, farsite_spread, 0, 1 + AMREX_SPACEDIM, AMREX_SPACEDIM, 0);
      WriteSingleLevelPlotfile("plt0000", plotmf, names, geom, 0.0, 0);
    }

    // ---------------- Time stepping ------------------------
    for (int step = 1; step <= inputs.nsteps; ++step) {
      fill_boundary_extrap(phi, geom);
      const Real dt_step = dt;
      amrex::Print() << "Time:"<< time << std::endl;
      // Run either level set advection OR FARSITE ellipse spread (mutually exclusive)
      if (inputs.skip_levelset == 0) {
	// Traditional level set advection
	advect_levelset_weno5z_rk3 (phi, vel, geom, dt_step, inputs.rothermel, terrain_slopes.get());
	dt = compute_dt(vel, geom, inputs.cfl);
      } else if (inputs.farsite.enable == 1) {
	// FARSITE ellipse spread (only when skip_levelset == 1 and farsite.enable == 1)
	compute_farsite_spread(phi, vel, farsite_spread, geom, dt_step, inputs.rothermel, inputs.farsite, terrain_slopes.get());
      }
      if (inputs.reinit_int > 0 && (step % inputs.reinit_int == 0) && inputs.skip_levelset == 0) {
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
	  amrex::Print() << "Reinitialized phi with dtau = " << dtau << " and niters = " << niters << "\n";
	}
      }

      time += dt_step;

      if (inputs.plot_int > 0 && (step % inputs.plot_int == 0)) {
	char buf[64];
	std::snprintf(buf, sizeof(buf), "plt%04d", step);
	Vector<std::string> names = {"phi", "velx", "vely"
#if (AMREX_SPACEDIM == 3)
				     , "velz", "farsite_dx", "farsite_dy", "farsite_dz"
#else
				     , "farsite_dx", "farsite_dy"
#endif
	};
	MultiFab plotmf(ba, dm, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM, 0);
	MultiFab::Copy(plotmf, phi, 0, 0, 1, 0);
	MultiFab::Copy(plotmf, vel, 0, 1, AMREX_SPACEDIM, 0);
	MultiFab::Copy(plotmf, farsite_spread, 0, 1 + AMREX_SPACEDIM, AMREX_SPACEDIM, 0);
	WriteSingleLevelPlotfile(buf, plotmf, names, geom, time, step);
	amrex::Print() << "Wrote " << buf << "\n";
      }
    }
      // ---------------- Final write --------------------------
      {
	char buf[64];
	std::snprintf(buf, sizeof(buf), "plt%04d", inputs.nsteps);
	Vector<std::string> names = {"phi", "velx", "vely"
#if (AMREX_SPACEDIM == 3)
				     , "velz", "farsite_dx", "farsite_dy", "farsite_dz"
#else
				     , "farsite_dx", "farsite_dy"
#endif
	};
	MultiFab plotmf(ba, dm, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM, 0);
	MultiFab::Copy(plotmf, phi, 0, 0, 1, 0);
	MultiFab::Copy(plotmf, vel, 0, 1, AMREX_SPACEDIM, 0);
	MultiFab::Copy(plotmf, farsite_spread, 0, 1 + AMREX_SPACEDIM, AMREX_SPACEDIM, 0);
	WriteSingleLevelPlotfile(buf, plotmf, names, geom, time, inputs.nsteps);
	amrex::Print() << "Wrote final " << buf << "\n";
      }
    }
  amrex::Finalize();
  return 0;
}
