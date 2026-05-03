#include <AMReX.H>
#include <AMReX_Array4.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_BoxArray.H>
#include <AMReX_VisMF.H>
#include <cmath>
#include <vector>
#include <string>
#include <memory>
#include <fstream>
#include <sstream>

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
#include "albini_spotting.H"
#include "compute_fire_behavior.H"


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

    // Albini spotting diagnostic field: 4 components
    //   0 – lofting height H_z [m] at fire-front source cells
    //   1 – number of firebrands launched from each source cell
    //   2 – maximum landing distance [m] from each source cell
    //   3 – active flag at cells that received a spot ignition
    MultiFab albini_data(ba, dm, 4, 0);
    albini_data.setVal(0.0);

    // Fire behavior diagnostic fields
    MultiFab fireline_intensity_mf(ba, dm, 1, 0);  // Byram fireline intensity [kW/m]
    MultiFab flame_length_mf(ba, dm, 1, 0);         // Byram flame length [m]
    fireline_intensity_mf.setVal(0.0);
    flame_length_mf.setVal(0.0);


    // ---------------- Initialize ---------------------------
    // When FARSITE is enabled and level set is skipped, use indicator initialization
    // (phi = 1 inside, phi = 0 outside) instead of signed distance
    bool use_indicator = (inputs.farsite.enable == 1 && inputs.skip_levelset == 1);

    // Restart state (overridden by checkpoint when restart_chkfile is non-empty)
    Real time = 0.0;
    int restart_step = 0;

    if (!inputs.restart_chkfile.empty()) {
      // Restart from checkpoint: skip normal initialization and load phi from file
      read_checkpoint(inputs.restart_chkfile, phi, restart_step, time);
      fill_boundary_extrap(phi, geom);
    } else if (!inputs.fire_points_file.empty()) {
      // CSV fire-points takes precedence over all other source_type options
      init_phi_from_fire_points_csv(phi, geom, inputs.fire_points_file, inputs.fire_gaussian_sigma);
      fill_boundary_extrap(phi, geom);
    }
    else if (inputs.source_type == "sphere") {
      if (use_indicator) {
	init_phi_sphere_indicator(phi, geom, inputs.cx, inputs.cy, inputs.cz, inputs.radius);
      } else {
	init_phi_sphere(phi, geom, inputs.cx, inputs.cy, inputs.cz, inputs.radius);
      }
      fill_boundary_extrap(phi, geom);
    }
    else if(inputs.source_type == "box") {
      if (use_indicator) {
	init_phi_box_indicator(phi, geom, inputs.xmin, inputs.ymin, inputs.zmin, inputs.xmax, inputs.ymax, inputs.zmax);
      } else {
	init_phi_box(phi, geom, inputs.xmin, inputs.ymin, inputs.zmin, inputs.xmax, inputs.ymax, inputs.zmax);
      }
      fill_boundary_extrap(phi, geom);
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
      fill_boundary_extrap(phi, geom);
    }
    else if(inputs.source_type == "eb") {
      // EB implicit function: use indicator (-1/0) for FARSITE, SDF for level set
      if (use_indicator) {
        init_phi_from_eb_implicit_indicator(phi, geom, inputs.eb_type,
                                            inputs.eb_param1, inputs.eb_param2, inputs.eb_param3,
                                            inputs.eb_param4, inputs.eb_param5, inputs.eb_param6);
      } else {
        init_phi_from_eb_implicit(phi, geom, inputs.eb_type,
                                  inputs.eb_param1, inputs.eb_param2, inputs.eb_param3,
                                  inputs.eb_param4, inputs.eb_param5, inputs.eb_param6);
      }
      fill_boundary_extrap(phi, geom);
    }
    else {
      amrex::Abort("Invalid source_type: " + inputs.source_type);
    }
        
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
    // Priority for elevation/slope/aspect: terrain_file > landscape_file
    // Fuel model always comes from landscape_file when present (terrain file has no fuel data)
    std::unique_ptr<MultiFab> terrain_slopes;
    if (!inputs.rothermel.terrain_file.empty()) {
      // Create MultiFab for slopes (2 components: slope_x, slope_y)
      terrain_slopes = std::make_unique<MultiFab>(ba, dm, 2, 0);
      compute_slopes_from_terrain(*terrain_slopes, geom, inputs.rothermel.terrain_file);
      amrex::Print() << "Initialized terrain slopes from terrain file: " 
		     << inputs.rothermel.terrain_file << "\n";
      if (!inputs.rothermel.landscape_file.empty()) {
        amrex::Print() << "NOTE: terrain_file takes precedence for elevation/slope/aspect; "
                          "landscape_file used for fuel model only\n";
      }
    } else if (!inputs.rothermel.landscape_file.empty()) {
      // Create MultiFab for slopes (2 components: slope_x, slope_y)
      terrain_slopes = std::make_unique<MultiFab>(ba, dm, 2, 0);
      compute_slopes_from_landscape(*terrain_slopes, geom, inputs.rothermel.landscape_file);
      amrex::Print() << "Initialized terrain slopes from landscape file: " 
		     << inputs.rothermel.landscape_file << "\n";
    }

    // Create elevation MultiFab (1 component, no ghost cells).
    // Populated from terrain or landscape file when available; zero otherwise.
    // terrain_file takes precedence over landscape_file for elevation.
    MultiFab elevation_mf(ba, dm, 1, 0);
    elevation_mf.setVal(0.0);
    if (!inputs.rothermel.terrain_file.empty()) {
      compute_elevation_from_terrain(elevation_mf, geom, inputs.rothermel.terrain_file);
    } else if (!inputs.rothermel.landscape_file.empty()) {
      compute_elevation_from_landscape(elevation_mf, geom, inputs.rothermel.landscape_file);
    }

    // Create slope (degrees), aspect (degrees), and fuel model MultiFabs.
    // Populated from terrain or landscape file when available; zero otherwise.
    // terrain_file takes precedence for slope/aspect; fuel model always comes
    // from landscape_file when present (terrain file carries no fuel data).
    MultiFab slope_mf(ba, dm, 1, 0);
    MultiFab aspect_mf(ba, dm, 1, 0);
    MultiFab fuel_model_mf(ba, dm, 1, 0);
    slope_mf.setVal(0.0);
    aspect_mf.setVal(0.0);
    fuel_model_mf.setVal(0.0);
    if (!inputs.rothermel.terrain_file.empty() && terrain_slopes) {
      compute_slope_aspect_from_slopes(*terrain_slopes, slope_mf, aspect_mf);
      // fuel_model is not available from terrain file; use landscape file if present
      if (!inputs.rothermel.landscape_file.empty()) {
        compute_fuel_model_from_landscape(fuel_model_mf, geom, inputs.rothermel.landscape_file);
      }
    } else if (!inputs.rothermel.landscape_file.empty()) {
      compute_landscape_slope_aspect_fuel(slope_mf, aspect_mf, fuel_model_mf,
                                          geom, inputs.rothermel.landscape_file);
    }

    // ---------------- Fuel model lookup table for per-cell spread --------
    // When a landscape file is present with fuel model data, precompute a
    // RothermelComputed lookup table indexed by fuel code so the GPU kernel
    // can use per-cell fuel parameters rather than the single global model.
    Gpu::DeviceVector<RothermelComputed> d_fuel_table;
    const RothermelComputed* d_fuel_table_ptr = nullptr;
    int fuel_table_size = 0;

    if (!inputs.rothermel.landscape_file.empty()) {
        std::vector<RothermelComputed> h_fuel_table =
            build_fuel_rothermel_table(inputs.rothermel,
                                       inputs.rothermel.landscape_fuel_type);
        fuel_table_size = static_cast<int>(h_fuel_table.size());
        d_fuel_table.resize(fuel_table_size);
        Gpu::copy(Gpu::hostToDevice,
                  h_fuel_table.begin(), h_fuel_table.end(),
                  d_fuel_table.begin());
        d_fuel_table_ptr = d_fuel_table.data();
        amrex::Print() << "Built per-cell Rothermel lookup table: "
                       << fuel_table_size << " entries ("
                       << (inputs.rothermel.landscape_fuel_type == "40"
                           ? "FBFM40" : "FBFM13")
                       << ")\n";
    }

    // ---------------- dt from CFL --------------------------
    Real dt=10;
    if (inputs.skip_levelset == 0)
      {
        // Compute Rothermel wind speed R
        compute_rothermel_R(R_mf, vel, geom, inputs.rothermel,
                             terrain_slopes.get(),
                             !inputs.rothermel.landscape_file.empty() ? &fuel_model_mf : nullptr,
                             d_fuel_table_ptr, fuel_table_size);
        dt = compute_dt(R_mf, geom, inputs.cfl);
        amrex::Print() << "Computed dt = " << dt << "\n";
      } else {
      amrex::Print() << "Skipping level set advection; using dt = " << dt << " for FARSITE spread\n";
    }
    compute_fire_behavior(fireline_intensity_mf, flame_length_mf, R_mf, inputs.rothermel);


    // ---------------- Write initial plotfile ---------------
    {
      Vector<std::string> names = {"phi", "velx", "vely"
#if (AMREX_SPACEDIM == 3)
				   , "velz", "farsite_dx", "farsite_dy", "farsite_dz", "R",
				   "spot_prob", "spot_count", "spot_dist", "spot_active", "fuel_consumption", "crown_fraction",
				   "albini_Hz", "albini_count", "albini_dist", "albini_active",
				   "elevation", "slope", "aspect", "fuel_model",
				   "fireline_intensity", "flame_length"
#else
				   , "farsite_dx", "farsite_dy", "R",
				   "spot_prob", "spot_count", "spot_dist", "spot_active", "fuel_consumption", "crown_fraction",
				   "albini_Hz", "albini_count", "albini_dist", "albini_active",
				   "elevation", "slope", "aspect", "fuel_model",
				   "fireline_intensity", "flame_length"
#endif
      };
      MultiFab plotmf(ba, dm, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 1 + 5, 0);
      MultiFab::Copy(plotmf, phi, 0, 0, 1, 0);
      MultiFab::Copy(plotmf, vel, 0, 1, AMREX_SPACEDIM, 0);
      MultiFab::Copy(plotmf, farsite_spread, 0, 1 + AMREX_SPACEDIM, AMREX_SPACEDIM, 0);
      MultiFab::Copy(plotmf, R_mf, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM, 1, 0);
      MultiFab::Copy(plotmf, spotting_data, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1, 4, 0);
      MultiFab::Copy(plotmf, fuel_consumption_mf, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4, 1, 0);
      MultiFab::Copy(plotmf, crown_fire_fraction_mf, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1, 1, 0);
      MultiFab::Copy(plotmf, albini_data, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1, 4, 0);
      MultiFab::Copy(plotmf, elevation_mf, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4, 1, 0);
      MultiFab::Copy(plotmf, slope_mf, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 1, 1, 0);
      MultiFab::Copy(plotmf, aspect_mf, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 2, 1, 0);
      MultiFab::Copy(plotmf, fuel_model_mf, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 3, 1, 0);
      MultiFab::Copy(plotmf, fireline_intensity_mf, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 4, 1, 0);
      MultiFab::Copy(plotmf, flame_length_mf, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 5, 1, 0);
      {
        char buf[64];
        std::snprintf(buf, sizeof(buf), "plt%04d", restart_step);
        WriteSingleLevelPlotfile(buf, plotmf, names, geom, time, restart_step);
      }
      
      // Write negative phi x-y data files
      {
        char xy_buf[64];
        std::snprintf(xy_buf, sizeof(xy_buf), "phi_negative_%04d.dat", restart_step);
        write_negative_phi_xy(phi, geom, xy_buf);
        std::snprintf(xy_buf, sizeof(xy_buf), "phi_envelope_%04d.dat", restart_step);
        write_negative_phi_convex_hull(phi, geom, xy_buf);
      }
    }

    // ---------------- Time stepping ------------------------
    // Run until final_time (if > 0) or nsteps steps (backward-compatible fallback)
    const bool use_final_time = (inputs.final_time > 0.0);
    int step = restart_step;
    while ((use_final_time && time < inputs.final_time) ||
           (!use_final_time && step < restart_step + inputs.nsteps)) {
      ++step;
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
  	compute_rothermel_R(R_mf, vel, geom, inputs.rothermel,
                          terrain_slopes.get(),
                          !inputs.rothermel.landscape_file.empty() ? &fuel_model_mf : nullptr,
                          d_fuel_table_ptr, fuel_table_size);
      compute_fire_behavior(fireline_intensity_mf, flame_length_mf, R_mf, inputs.rothermel);
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
	// Albini (1983) firebrand spotting with 2-D trajectory integration
	if (inputs.albini_spotting.enable == 1 && (step % inputs.albini_spotting.check_interval == 0)) {
	  compute_albini_spotting(phi, albini_data, vel, R_mf, geom,
	                          inputs.rothermel, inputs.albini_spotting, step);
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

      // --- Dynamic fire points: check if a new ignition file has appeared
      if (!inputs.dynamic_fire_points_file.empty()) {
        apply_dynamic_fire_points(phi, geom,
                                  inputs.dynamic_fire_points_file,
                                  inputs.fire_gaussian_sigma);
      }

      // --- Write checkpoint if requested
      if (inputs.chk_int > 0 && (step % inputs.chk_int == 0)) {
        char chk_buf[64];
        std::snprintf(chk_buf, sizeof(chk_buf), "chk%04d", step);
        write_checkpoint(chk_buf, phi, geom, step, time);
      }

      // --- Step 7: Update states, record outputs, step time
      if (inputs.plot_int > 0 && (step % inputs.plot_int == 0)) {
	char buf[64];
	std::snprintf(buf, sizeof(buf), "plt%04d", step);
	Vector<std::string> names = {"phi", "velx", "vely"
#if (AMREX_SPACEDIM == 3)
				     , "velz", "farsite_dx", "farsite_dy", "farsite_dz", "R",
				     "spot_prob", "spot_count", "spot_dist", "spot_active", "fuel_consumption", "crown_fraction",
				     "albini_Hz", "albini_count", "albini_dist", "albini_active",
				     "elevation", "slope", "aspect", "fuel_model",
				     "fireline_intensity", "flame_length"
#else
				     , "farsite_dx", "farsite_dy", "R",
				     "spot_prob", "spot_count", "spot_dist", "spot_active", "fuel_consumption", "crown_fraction",
				     "albini_Hz", "albini_count", "albini_dist", "albini_active",
				     "elevation", "slope", "aspect", "fuel_model",
				     "fireline_intensity", "flame_length"
#endif
	};
	MultiFab plotmf(ba, dm, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 1 + 5, 0);
	MultiFab::Copy(plotmf, phi, 0, 0, 1, 0);
	MultiFab::Copy(plotmf, vel, 0, 1, AMREX_SPACEDIM, 0);
	MultiFab::Copy(plotmf, farsite_spread, 0, 1 + AMREX_SPACEDIM, AMREX_SPACEDIM, 0);
	MultiFab::Copy(plotmf, R_mf, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM, 1, 0);
	MultiFab::Copy(plotmf, spotting_data, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1, 4, 0);
	MultiFab::Copy(plotmf, fuel_consumption_mf, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4, 1, 0);
	MultiFab::Copy(plotmf, crown_fire_fraction_mf, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1, 1, 0);
	MultiFab::Copy(plotmf, albini_data, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1, 4, 0);
	MultiFab::Copy(plotmf, elevation_mf, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4, 1, 0);
	MultiFab::Copy(plotmf, slope_mf, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 1, 1, 0);
	MultiFab::Copy(plotmf, aspect_mf, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 2, 1, 0);
	MultiFab::Copy(plotmf, fuel_model_mf, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 3, 1, 0);
	MultiFab::Copy(plotmf, fireline_intensity_mf, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 4, 1, 0);
	MultiFab::Copy(plotmf, flame_length_mf, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 5, 1, 0);
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
      const int final_step = step;
      bool should_write_final = (inputs.plot_int <= 0);
      if (inputs.plot_int > 0) {
          should_write_final = (final_step % inputs.plot_int != 0);
      }
      if (should_write_final)
      {
	char buf[64];
	std::snprintf(buf, sizeof(buf), "plt%04d", final_step);
	Vector<std::string> names = {"phi", "velx", "vely"
#if (AMREX_SPACEDIM == 3)
				     , "velz", "farsite_dx", "farsite_dy", "farsite_dz", "R",
				     "spot_prob", "spot_count", "spot_dist", "spot_active", "fuel_consumption", "crown_fraction",
				     "albini_Hz", "albini_count", "albini_dist", "albini_active",
				     "elevation", "slope", "aspect", "fuel_model",
				     "fireline_intensity", "flame_length"
#else
				     , "farsite_dx", "farsite_dy", "R",
				     "spot_prob", "spot_count", "spot_dist", "spot_active", "fuel_consumption", "crown_fraction",
				     "albini_Hz", "albini_count", "albini_dist", "albini_active",
				     "elevation", "slope", "aspect", "fuel_model",
				     "fireline_intensity", "flame_length"
#endif
	};
	MultiFab plotmf(ba, dm, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 1 + 5, 0);
	MultiFab::Copy(plotmf, phi, 0, 0, 1, 0);
	MultiFab::Copy(plotmf, vel, 0, 1, AMREX_SPACEDIM, 0);
	MultiFab::Copy(plotmf, farsite_spread, 0, 1 + AMREX_SPACEDIM, AMREX_SPACEDIM, 0);
	MultiFab::Copy(plotmf, R_mf, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM, 1, 0);
	MultiFab::Copy(plotmf, spotting_data, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1, 4, 0);
	MultiFab::Copy(plotmf, fuel_consumption_mf, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4, 1, 0);
	MultiFab::Copy(plotmf, crown_fire_fraction_mf, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1, 1, 0);
	MultiFab::Copy(plotmf, albini_data, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1, 4, 0);
	MultiFab::Copy(plotmf, elevation_mf, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4, 1, 0);
	MultiFab::Copy(plotmf, slope_mf, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 1, 1, 0);
	MultiFab::Copy(plotmf, aspect_mf, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 2, 1, 0);
	MultiFab::Copy(plotmf, fuel_model_mf, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 3, 1, 0);
	MultiFab::Copy(plotmf, fireline_intensity_mf, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 4, 1, 0);
	MultiFab::Copy(plotmf, flame_length_mf, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 5, 1, 0);
	WriteSingleLevelPlotfile(buf, plotmf, names, geom, time, final_step);
	amrex::Print() << "Wrote final " << buf << "\n";
	
	// Write negative phi x-y data files for final step
	char xy_buf[64];
	std::snprintf(xy_buf, sizeof(xy_buf), "phi_negative_%04d.dat", final_step);
	write_negative_phi_xy(phi, geom, xy_buf);
	
	std::snprintf(xy_buf, sizeof(xy_buf), "phi_envelope_%04d.dat", final_step);
	write_negative_phi_convex_hull(phi, geom, xy_buf);
      }
    }
  amrex::Finalize();
  return 0;
}
