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
#include "landscape_file.H"
#include "write_xy_data.H"
#include "compute_rothermel_R.H"
#include "firebrand_spotting.H"


// ======================= Main ================================================
// 
// Summary Flow:
// 1. Setup inputs (landscape, fuel, weather, wind)
// 2. Compute surface ROS via Rothermel/Level Set
// 3. Generate elliptical wavelets per vertex
// 4. Merge to new perimeter
// 5. Apply crown/spotting sub-models
// 6. Simulate post-frontal burnout
// 7. Update states, record outputs, step time
//
int main(int argc, char* argv[])
{
  amrex::Initialize(argc, argv);
  {
    // --- Step 1: Setup inputs (landscape, fuel, weather, wind)
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

        // ---------------- Inputs: terrain ----------------------
        bool use_terrain_effects = false;
        bool use_farsite_model = false;
        Real terrain_slope = 0.0;      // slope in degrees
        Real terrain_aspect = 0.0;     // aspect in degrees (0=East, 90=North, 180=West, 270=South)
        amrex::ParmParse pp;
        
        pp.query("use_terrain_effects", use_terrain_effects);
        pp.query("use_farsite_model", use_farsite_model);
        pp.query("terrain_slope", terrain_slope);
        pp.query("terrain_aspect", terrain_aspect);

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

    // Firebrand spotting data field: 4 components (probability, count, distance, active flag)
    MultiFab spotting_data(ba, dm, 4, 0);

    // Rothermel wind speed R field (1 component, no ghost cells)
    MultiFab R_mf(ba, dm, 1, 0);
    
    // Bulk fuel consumption fraction field (1 component, no ghost cells)
    MultiFab fuel_consumption_mf(ba, dm, 1, 0);
    fuel_consumption_mf.setVal(0.0); // Initialize to zero
    
    // Crown fire fraction field (1 component, no ghost cells)
    MultiFab crown_fire_fraction_mf(ba, dm, 1, 0);
    crown_fire_fraction_mf.setVal(0.0); // Initialize to zero


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
	init_phi_box_indicator(phi, geom, inputs.xmin, inputs.ymin, inputs.zmin, inputs.xmax, inputs.ymax, inputs.zmax);
      } else {
	init_phi_box(phi, geom, inputs.xmin, inputs.ymin, inputs.zmin, inputs.xmax, inputs.ymax, inputs.zmax);
      }
    }
    else if(inputs.source_type == "ellipse") {
      if (use_indicator) {
	init_phi_ellipse_indicator(phi, geom, 
	                          inputs.ellipse_center_x, inputs.ellipse_center_y, inputs.ellipse_center_z,
	                          inputs.ellipse_radius_x, inputs.ellipse_radius_y, inputs.ellipse_radius_z);
      } else {
	init_phi_ellipse(phi, geom, 
	                inputs.ellipse_center_x, inputs.ellipse_center_y, inputs.ellipse_center_z,
	                inputs.ellipse_radius_x, inputs.ellipse_radius_y, inputs.ellipse_radius_z);
      }
    }
    else if(inputs.source_type == "eb") {
      // EB implicit function always uses signed distance (no indicator mode)
      init_phi_from_eb_implicit(phi, geom, inputs.eb_type,
                               inputs.eb_param1, inputs.eb_param2, inputs.eb_param3,
                               inputs.eb_param4, inputs.eb_param5, inputs.eb_param6);
    }
    else {
      amrex::Abort("Invalid source_type: " + inputs.source_type);
    }
    fill_boundary_extrap(phi, geom);
        
    // Initialize velocity field
#if (AMREX_SPACEDIM == 2)
    // Storage for time-dependent wind field data
    std::vector<Real> wind_x_data1, wind_y_data1, wind_u_data1, wind_v_data1;
    std::vector<Real> wind_x_data2, wind_y_data2, wind_u_data2, wind_v_data2;
    int current_wind_field_index = -1;
    int next_wind_field_index = -1;
    
    if (!inputs.velocity_file.empty()) {
      if (inputs.use_time_dependent_wind == 1) {
        // Time-dependent wind field: load initial two snapshots
        amrex::Print() << "Using time-dependent wind fields with spacing = " 
                       << inputs.wind_time_spacing << " seconds\n";
        update_time_dependent_velocity(vel, geom, inputs.velocity_file, 0.0, inputs.wind_time_spacing,
                                        wind_x_data1, wind_y_data1, wind_u_data1, wind_v_data1,
                                        wind_x_data2, wind_y_data2, wind_u_data2, wind_v_data2,
                                        current_wind_field_index, next_wind_field_index);
      } else {
        // Static wind field
        init_velocity_from_file(vel, geom, inputs.velocity_file);
      }
    } else {
      init_velocity_constant(vel, geom, inputs.ux, inputs.uy, inputs.uz);
    }
#else
    init_velocity_constant(vel, geom, inputs.ux, inputs.uy, inputs.uz);
#endif

    // Initialize terrain slopes
    // Priority: landscape_file > terrain_file
    // If landscape file is specified, ignore terrain file slope/elevation
    std::unique_ptr<MultiFab> terrain_slopes;
    if (!inputs.rothermel.landscape_file.empty()) {
      // Create MultiFab for slopes (2 components: slope_x, slope_y)
      terrain_slopes = std::make_unique<MultiFab>(ba, dm, 2, 0);
      compute_slopes_from_landscape(*terrain_slopes, geom, inputs.rothermel.landscape_file);
      amrex::Print() << "Initialized terrain slopes from landscape file: " 
		     << inputs.rothermel.landscape_file << "\n";
      amrex::Print() << "NOTE: Ignoring terrain_file (if specified) because landscape_file takes precedence\n";
    } else if (!inputs.rothermel.terrain_file.empty()) {
      // Create MultiFab for slopes (2 components: slope_x, slope_y)
      terrain_slopes = std::make_unique<MultiFab>(ba, dm, 2, 0);
      compute_slopes_from_terrain(*terrain_slopes, geom, inputs.rothermel.terrain_file);
      amrex::Print() << "Initialized terrain slopes from terrain file: " 
		     << inputs.rothermel.terrain_file << "\n";
    }

    // ---------------- dt from CFL --------------------------
    Real dt=10;
    if (inputs.skip_levelset == 0)
      {
        // Compute Rothermel wind speed R
        compute_rothermel_R(R_mf, vel, geom, inputs.rothermel, terrain_slopes.get());
        dt = compute_dt(R_mf, geom, inputs.cfl);
        amrex::Print() << "Computed dt = " << dt << "\n";
      } else {
      amrex::Print() << "Skipping level set advection; using dt = " << dt << " for FARSITE spread\n";
    }
    Real time = 0.0;


    // ---------------- Write initial plotfile ---------------
    {
      Vector<std::string> names = {"phi", "velx", "vely"
#if (AMREX_SPACEDIM == 3)
				   , "velz", "farsite_dx", "farsite_dy", "farsite_dz", "R",
				   "spot_prob", "spot_count", "spot_dist", "spot_active", "fuel_consumption", "crown_fraction"
#else
				   , "farsite_dx", "farsite_dy", "R",
				   "spot_prob", "spot_count", "spot_dist", "spot_active", "fuel_consumption", "crown_fraction"
#endif
      };
      MultiFab plotmf(ba, dm, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1, 0);
      MultiFab::Copy(plotmf, phi, 0, 0, 1, 0);
      MultiFab::Copy(plotmf, vel, 0, 1, AMREX_SPACEDIM, 0);
      MultiFab::Copy(plotmf, farsite_spread, 0, 1 + AMREX_SPACEDIM, AMREX_SPACEDIM, 0);
      MultiFab::Copy(plotmf, R_mf, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM, 1, 0);
      MultiFab::Copy(plotmf, spotting_data, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1, 4, 0);
      MultiFab::Copy(plotmf, fuel_consumption_mf, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4, 1, 0);
      MultiFab::Copy(plotmf, crown_fire_fraction_mf, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1, 1, 0);
      WriteSingleLevelPlotfile("plt0000", plotmf, names, geom, 0.0, 0);
      
      // Write negative phi x-y data files
      write_negative_phi_xy(phi, geom, "phi_negative_0000.dat");
      write_negative_phi_convex_hull(phi, geom, "phi_envelope_0000.dat");
    }

    // ---------------- Time stepping ------------------------
    for (int step = 1; step <= inputs.nsteps; ++step) {
      fill_boundary_extrap(phi, geom);
      const Real dt_step = dt;
      amrex::Print() << "Time:"<< time << " with timestep:" << dt_step <<std::endl;
      
      // Update time-dependent wind field if enabled
#if (AMREX_SPACEDIM == 2)
      if (!inputs.velocity_file.empty() && inputs.use_time_dependent_wind == 1) {
        update_time_dependent_velocity(vel, geom, inputs.velocity_file, time, inputs.wind_time_spacing,
                                        wind_x_data1, wind_y_data1, wind_u_data1, wind_v_data1,
                                        wind_x_data2, wind_y_data2, wind_u_data2, wind_v_data2,
                                        current_wind_field_index, next_wind_field_index);
      }
#endif
      
      // --- Step 2: Compute surface ROS via Rothermel/Level Set
      // Update Rothermel wind speed R and dt
  	compute_rothermel_R(R_mf, vel, geom, inputs.rothermel, terrain_slopes.get());
      if (inputs.skip_levelset == 0) {
	// Traditional level set advection
	advect_levelset_weno5z_rk3 (phi, vel, geom, dt_step, inputs.rothermel, terrain_slopes.get());
	dt = compute_dt(R_mf, geom, inputs.cfl);
      } else if (inputs.farsite.enable == 1) {
	// --- Step 3: Generate elliptical wavelets per vertex
	// --- Step 4: Merge to new perimeter
	// FARSITE ellipse spread (only when skip_levelset == 1 and farsite.enable == 1)
	compute_farsite_spread(phi, vel, farsite_spread, geom, dt_step, inputs.rothermel, inputs.farsite, inputs.crown, terrain_slopes.get(), &fuel_consumption_mf, &crown_fire_fraction_mf);
	
	// --- Step 5: Apply crown/spotting sub-models
	// Add firebrand spotting model
	if (inputs.spotting.enable == 1 && (step % inputs.spotting.check_interval == 0)) {
	  compute_spotting_probability(spotting_data, phi, vel, geom, inputs.rothermel, inputs.spotting, terrain_slopes.get());
	  generate_firebrand_spots(phi, spotting_data, vel, geom, inputs.spotting, step);
	}
	
	// --- Step 6: Simulate post-frontal burnout
	// (Bulk fuel consumption is computed within compute_farsite_spread)
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

      // --- Step 7: Update states, record outputs, step time
      if (inputs.plot_int > 0 && (step % inputs.plot_int == 0)) {
	char buf[64];
	std::snprintf(buf, sizeof(buf), "plt%04d", step);
	Vector<std::string> names = {"phi", "velx", "vely"
#if (AMREX_SPACEDIM == 3)
				     , "velz", "farsite_dx", "farsite_dy", "farsite_dz", "R",
				     "spot_prob", "spot_count", "spot_dist", "spot_active", "fuel_consumption", "crown_fraction"
#else
				     , "farsite_dx", "farsite_dy", "R",
				     "spot_prob", "spot_count", "spot_dist", "spot_active", "fuel_consumption", "crown_fraction"
#endif
	};
	MultiFab plotmf(ba, dm, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1, 0);
	MultiFab::Copy(plotmf, phi, 0, 0, 1, 0);
	MultiFab::Copy(plotmf, vel, 0, 1, AMREX_SPACEDIM, 0);
	MultiFab::Copy(plotmf, farsite_spread, 0, 1 + AMREX_SPACEDIM, AMREX_SPACEDIM, 0);
	MultiFab::Copy(plotmf, R_mf, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM, 1, 0);
	MultiFab::Copy(plotmf, spotting_data, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1, 4, 0);
	MultiFab::Copy(plotmf, fuel_consumption_mf, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4, 1, 0);
	MultiFab::Copy(plotmf, crown_fire_fraction_mf, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1, 1, 0);
	WriteSingleLevelPlotfile(buf, plotmf, names, geom, time, step);
	amrex::Print() << "Wrote " << buf << "\n";
	
	// Write negative phi x-y data files
	char xy_buf[64];
	std::snprintf(xy_buf, sizeof(xy_buf), "phi_negative_%04d.dat", step);
	write_negative_phi_xy(phi, geom, xy_buf);
	
	std::snprintf(xy_buf, sizeof(xy_buf), "phi_envelope_%04d.dat", step);
	write_negative_phi_convex_hull(phi, geom, xy_buf);
      }
    }
      // ---------------- Final write --------------------------
      // Only write final if it wasn't already written at plot_int
      bool should_write_final = (inputs.plot_int <= 0);
      if (inputs.plot_int > 0) {
          should_write_final = (inputs.nsteps % inputs.plot_int != 0);
      }
      if (should_write_final)
      {
	char buf[64];
	std::snprintf(buf, sizeof(buf), "plt%04d", inputs.nsteps);
	Vector<std::string> names = {"phi", "velx", "vely"
#if (AMREX_SPACEDIM == 3)
				     , "velz", "farsite_dx", "farsite_dy", "farsite_dz", "R",
				     "spot_prob", "spot_count", "spot_dist", "spot_active", "fuel_consumption", "crown_fraction"
#else
				     , "farsite_dx", "farsite_dy", "R",
				     "spot_prob", "spot_count", "spot_dist", "spot_active", "fuel_consumption", "crown_fraction"
#endif
	};
	MultiFab plotmf(ba, dm, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1, 0);
	MultiFab::Copy(plotmf, phi, 0, 0, 1, 0);
	MultiFab::Copy(plotmf, vel, 0, 1, AMREX_SPACEDIM, 0);
	MultiFab::Copy(plotmf, farsite_spread, 0, 1 + AMREX_SPACEDIM, AMREX_SPACEDIM, 0);
	MultiFab::Copy(plotmf, R_mf, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM, 1, 0);
	MultiFab::Copy(plotmf, spotting_data, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1, 4, 0);
	MultiFab::Copy(plotmf, fuel_consumption_mf, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4, 1, 0);
	MultiFab::Copy(plotmf, crown_fire_fraction_mf, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1, 1, 0);
	WriteSingleLevelPlotfile(buf, plotmf, names, geom, time, inputs.nsteps);
	amrex::Print() << "Wrote final " << buf << "\n";
	
	// Write negative phi x-y data files for final step
	char xy_buf[64];
	std::snprintf(xy_buf, sizeof(xy_buf), "phi_negative_%04d.dat", inputs.nsteps);
	write_negative_phi_xy(phi, geom, xy_buf);
	
	std::snprintf(xy_buf, sizeof(xy_buf), "phi_envelope_%04d.dat", inputs.nsteps);
	write_negative_phi_convex_hull(phi, geom, xy_buf);
      }
    }
  amrex::Finalize();
  return 0;
}
