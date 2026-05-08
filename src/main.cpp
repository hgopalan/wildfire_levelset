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
#include <algorithm>
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
#include "balbi_model.H"
#include "compute_balbi_R.H"
#include "andrews_model.H"
#include "cheney_gould_model.H"
#include "compute_cheney_gould_R.H"
#include "weise_biging_whirl.H"
#include "viegas_model.H"
#include "wind_terrain_models.H"
#include "heat_flux_model.H"
#include "cruz_crown_model.H"
#include "compute_cruz_crown_R.H"
#include "fire_ecology_model.H"
#include "fire_emissions_model.H"
#include "heat_per_unit_area.H"
#include "burnout_time.H"
#include "vorticity_model.H"
#include "turb_wind.H"
#include "fuel_adjustment.H"
#include "fuel_moisture_scheduler.H"
#include "scott_spotting_table.H"
#include "mtt_propagation.H"
#include "barrier_polygons.H"
#include "compute_reaction_intensity.H"
#include "albini_torching_spotting.H"
#include "fire_acceleration.H"
#include "fmc_schedule.H"
#include "precipitation_moisture.H"
#include "polygon_ignition.H"
#include "fms_moisture.H"
#include "wind_dir_schedule.H"
#include "fbp_model.H"
#include "compute_fbp_R.H"
#include "lautenberger_model.H"
#include "compute_lautenberger_R.H"
#include "solar_radiation.H"
#include "ignition_schedule.H"
#include "wtr_weather.H"
#include "retardant_drop.H"
#include "herb_moisture_schedule.H"
#include "multi_wtr_weather.H"


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

        // ---------------- Grids & distribution -----------------
    BoxArray ba(domain);
    ba.maxSize(inputs.max_grid);
    DistributionMapping dm(ba);

    // ---------------- Fields: phi (1 comp), vel (3 comps) --
    const int ng_phi = 3; // 3 ghost cells for stencil operations (WENO5-Z flux divergence uses up to ±3 cells)
    MultiFab phi(ba, dm, 1, ng_phi);
    MultiFab vel(ba, dm, 3, 1);

    // Effective wind field for wind-terrain feedback models (Options 3–6).
    // For "none" and "viegas_ros" this is unused; Rothermel sees vel directly.
    MultiFab vel_effective(ba, dm, 3, 1);
        
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

    // Fire ecology diagnostics (ECOLOGY_NCOMP = 7 components):
    //   0 – scorch_height            [m]  (Van Wagner 1973)
    //   1 – prob_ignition            [-]  (Anderson 1970 / Rothermel 1983)
    //   2 – tree_mortality           [-]  (Ryan-Reinhardt 1988 style logistic)
    //   3 – crown_activity           [-]  (0=surface, 1=passive, 2=active crown fire)
    //   4 – torching_ratio           [-]  (Scott & Reinhardt 2001 TI proxy)
    //   5 – crowning_ratio           [-]  (Scott & Reinhardt 2001 CI proxy)
    //   6 – energy_release_component [-]  (Deeming et al. 1977 NFDRS ERC)
    MultiFab ecology_mf(ba, dm, ECOLOGY_NCOMP, 0);
    ecology_mf.setVal(0.0);

    // Fire emissions (EMISSIONS_NCOMP = 3 components):
    //   0 – co2_emissions   [kg CO₂/m²]
    //   1 – co_emissions    [kg CO/m²]
    //   2 – pm25_emissions  [kg PM₂.₅/m²]
    MultiFab emissions_mf(ba, dm, EMISSIONS_NCOMP, 0);
    emissions_mf.setVal(0.0);

    // Arrival time field [s]: simulation time when each cell first ignited.
    // Initialized to -1.0 (sentinel for "not yet burned").
    // Set to the current simulation time on the first timestep phi < 0.
    MultiFab arrival_time_mf(ba, dm, 1, 0);
    arrival_time_mf.setVal(Real(-1.0));

    // Burnout time field [s]: arrival_time + tau_residence per burned cell.
    // -1.0 sentinel for cells not yet burned.
    MultiFab burnout_time_mf(ba, dm, 1, 0);
    burnout_time_mf.setVal(Real(-1.0));

    // Feature 9: Residual fuel fraction [-] per cell (1.0 = fully loaded, 0.0 = exhausted).
    // Decays exponentially after fire passage at rate exp(-dt/tau_burnout).
    // Written to every plotfile as "residual_fuel".
    MultiFab residual_fuel_mf(ba, dm, 1, 0);
    residual_fuel_mf.setVal(Real(1.0));

    // Heat per unit area [BTU/ft²]: I_R × residence_time [min] for burned cells.
    // Zero for unburned cells. Recomputed before each plotfile write.
    MultiFab heat_per_unit_area_mf(ba, dm, 1, 0);
    heat_per_unit_area_mf.setVal(0.0);

    // Vertical vorticity ω_z = ∂v/∂x − ∂u/∂y [s⁻¹].
    // Recomputed from the current velocity field before each plotfile write.
    MultiFab vorticity_mf(ba, dm, 1, 0);
    vorticity_mf.setVal(0.0);

    // Scott & Reinhardt (2001) full TI and CI [km/h]
    // Computed once at initialization (fuel-property-based, not fire-state-based).
    MultiFab ti_full_mf(ba, dm, 1, 0);
    MultiFab ci_full_mf(ba, dm, 1, 0);
    ti_full_mf.setVal(Real(-1.0));  // -1 = not computed / disabled
    ci_full_mf.setVal(Real(-1.0));

    // Reaction intensity I_R [kW/m²] (Rothermel 1972).
    // Per-cell when a landscape fuel table is active; uniform otherwise.
    // Recomputed before each plotfile write.
    MultiFab reaction_intensity_mf(ba, dm, 1, 0);
    reaction_intensity_mf.setVal(0.0);

    // Weise & Biging (1996) fire whirl diagnostic fields (WEISE_NCOMP components)
    //   0 – weise_flame_height     [m]
    //   1 – weise_flame_tilt       [rad]
    //   2 – weise_whirl_height     [m]
    //   3 – weise_whirl_radius     [m]
    //   4 – weise_angular_velocity [rad/s]
    //   5 – weise_max_tang_vel     [m/s]
    MultiFab weise_data(ba, dm, WEISE_NCOMP, 0);
    weise_data.setVal(0.0);

    // Viegas (2004) eruptive fire diagnostics and fire-terrain feedback
    //   0 – viegas_ROS          [m/s]
    //   1 – viegas_eruptive_flag [-]   (1 when eruptive conditions met)
    //   2 – viegas_ROS_excess   [-]   (R_V - R_primary)/R_primary
    //   3 – viegas_flame_tilt   [rad]  (Viegas flame-tilt angle from vertical)
    //   4 – viegas_slope_factor [-]   (Viegas slope enhancement factor Phi_s_V)
    // Number of output components in the viegas_data MultiFab.
    static constexpr int VIEGAS_NCOMP = 5;
    MultiFab viegas_data(ba, dm, VIEGAS_NCOMP, 0);
    viegas_data.setVal(0.0);

    // Heat flux field [W/m²] (1 component, no ghost cells).
    // Initialized from heat_flux.value (uniform) or heat_flux.file (spatially varying).
    MultiFab heat_flux_mf(ba, dm, 1, 0);
    heat_flux_mf.setVal(0.0);

    // Per-cell shade fraction [-] (0 = fully insolated, 1 = fully shaded).
    // Updated each timestep when solar_radiation.enable = 1.
    // 0.0 default means no shading (no solar adjustment unless enabled).
    MultiFab shade_fraction_mf(ba, dm, 1, 0);
    shade_fraction_mf.setVal(0.0);


    // ---------------- Initialize ---------------------------
    // When FARSITE propagation is used, initialize phi as indicator (1 inside, 0 outside)
    bool use_indicator = (inputs.propagation_method == "farsite");

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
    else if (inputs.source_type == "polygon") {
      // Closed polygon ignition rasterizer
      std::vector<Real> poly_xs, poly_ys;
      read_polygon_vertices(inputs.fire_polygon_file, poly_xs, poly_ys);
      init_phi_from_polygon(phi, geom, poly_xs, poly_ys, inputs.fire_polygon_z_level);
      fill_boundary_extrap(phi, geom);
    }
    else if (inputs.source_type == "polyline") {
      // Polyline (line-fire) ignition rasterizer
      std::vector<Real> poly_xs, poly_ys;
      read_polygon_vertices(inputs.fire_polygon_file, poly_xs, poly_ys);
      init_phi_from_polyline(phi, geom, poly_xs, poly_ys,
                             inputs.polyline_width, inputs.fire_polygon_z_level);
      fill_boundary_extrap(phi, geom);
    }
    else {
      amrex::Abort("Invalid source_type: " + inputs.source_type);
    }

    // Mark cells that are already burned in the initial state (arrival_time = 0.0)
    {
      const Real t0 = Real(0.0);
      for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const p  = phi.const_array(mfi);
        auto       at = arrival_time_mf.array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
          if (p(i, j, k) < Real(0.0) && at(i, j, k) < Real(0.0)) {
            at(i, j, k) = t0;
          }
        });
      }
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
    // File-based velocity initialization is only supported in 2D.
    // In 3D, use constant velocity specified via inputs.ux/uy/uz.
    init_velocity_constant(vel, geom, inputs.ux, inputs.uy, inputs.uz);
#endif

    // ---------------- Turbulent wind perturbation setup ----------------
    // vel_base stores the unperturbed background wind (constant or time-dependent).
    // ou_state_mf stores the per-cell OU state (u', v') only when model=ou_process
    // and L_c > 0.  For domain-uniform OU, spectral_noise, and direction_walk the
    // scalar / spectral state in turb_state is sufficient; ou_state_mf is null.
    const bool turb_wind_active = (inputs.turb_wind.model != "none");
    std::unique_ptr<MultiFab> vel_base;
    std::unique_ptr<MultiFab> ou_state_mf;
    TurbWindState turb_state;
    if (turb_wind_active) {
        vel_base = std::make_unique<MultiFab>(ba, dm, 3, 1);
        MultiFab::Copy(*vel_base, vel, 0, 0, 3, 1);
        init_turb_wind_state(turb_state, inputs.turb_wind);
        if (inputs.turb_wind.model == "ou_process" && inputs.turb_wind.L_c > amrex::Real(0.0)) {
            ou_state_mf = std::make_unique<MultiFab>(ba, dm, 2, 0);
            ou_state_mf->setVal(amrex::Real(0.0));
            amrex::Print() << "Turbulent wind: per-cell OU with L_c="
                           << inputs.turb_wind.L_c << " m"
                           << "  sigma_k=" << inputs.turb_wind.L_c / geom.CellSize(0)
                           << " cells\n";
        }
    }

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

    // ---------------- Crown spatial layers from binary LCP ------------------
    // Per-cell CBH, CBD, and canopy cover MultiFabs are populated when:
    //   (a) a binary LCP with crown fuel layers is specified, AND
    //   (b) use_spatial_crown = 1 (default).
    // The global crown.CBH / crown.CBD scalars remain as fallbacks.
    MultiFab cbh_mf(ba, dm, 1, 0);   // crown base height [m]
    MultiFab cbd_mf(ba, dm, 1, 0);   // crown bulk density [kg/m³]
    MultiFab cc_mf(ba, dm, 1, 0);    // canopy cover [%]
    MultiFab canopy_height_mf(ba, dm, 1, 0);  // stand (canopy) height [m]
    cbh_mf.setVal(inputs.crown.CBH);
    cbd_mf.setVal(inputs.crown.CBD);
    cc_mf.setVal(0.0);
    canopy_height_mf.setVal(0.0);
    bool has_spatial_crown = false;
    if (inputs.use_spatial_crown == 1 &&
        !inputs.rothermel.landscape_file.empty()) {
      has_spatial_crown = compute_crown_layers_from_lcp(
          &cbh_mf, &cbd_mf, &cc_mf,
          geom, inputs.rothermel.landscape_file,
          &canopy_height_mf);
    }

    // ---------------- Per-cell live canopy moisture from FMS file -----------
    // Spatially-varying dead and live moisture per fuel model code.
    // Populated from a FARSITE .fms scenario file when fms_file is non-empty.
    // When absent, the global RothermelParams moisture values are used.
    const int FMS_NCOMP = FMS_MOISTURE_NCOMP;  // 5 components
    MultiFab spatial_moisture_mf(ba, dm, FMS_NCOMP, 0);
    bool has_spatial_moisture = false;
    // Initialise to global Rothermel moisture values as fallback
    spatial_moisture_mf.setVal(static_cast<Real>(inputs.rothermel.M_d1),   0, 1, 0);
    spatial_moisture_mf.setVal(static_cast<Real>(inputs.rothermel.M_d10),  1, 1, 0);
    spatial_moisture_mf.setVal(static_cast<Real>(inputs.rothermel.M_d100), 2, 1, 0);
    spatial_moisture_mf.setVal(static_cast<Real>(inputs.rothermel.M_lh),   3, 1, 0);
    spatial_moisture_mf.setVal(static_cast<Real>(inputs.rothermel.M_lw),   4, 1, 0);
    if (!inputs.fms_file.empty() && !inputs.rothermel.landscape_file.empty()) {
        load_fms_spatial_moisture(spatial_moisture_mf, fuel_model_mf,
                                  geom, inputs.fms_file, inputs.rothermel);
        has_spatial_moisture = true;
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
        // Apply fuel adjustment to the per-fuel lookup table
        if (!inputs.fuel_adj_file.empty()) {
            auto adjs = parse_fuel_adjustment_file(inputs.fuel_adj_file);
            apply_fuel_adjustment_to_table(h_fuel_table, adjs);
        }
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
    } else if (!inputs.fuel_adj_file.empty()) {
        // No landscape file: apply adjustment to global RothermelParams
        auto adjs = parse_fuel_adjustment_file(inputs.fuel_adj_file);
        apply_fuel_adjustment_to_params(inputs.rothermel, adjs, inputs.fuel_adj_model);
    }

    // ---------------- Per-fuel burnout (residence) time table ---------------
    // Build a device-side table of per-fuel-model residence times used by
    // compute_burnout_time.  When no landscape file is present the global
    // tau_residence scalar is used instead (d_tau_table_ptr stays null).
    Gpu::DeviceVector<FuelResidenceTime> d_tau_table;
    const FuelResidenceTime* d_tau_table_ptr = nullptr;
    int tau_table_size = 0;

    if (d_fuel_table_ptr != nullptr && fuel_table_size > 0) {
        std::vector<FuelResidenceTime> h_tau_table =
            build_fuel_tau_table(fuel_table_size);
        tau_table_size = static_cast<int>(h_tau_table.size());
        d_tau_table.resize(tau_table_size);
        Gpu::copy(Gpu::hostToDevice,
                  h_tau_table.begin(), h_tau_table.end(),
                  d_tau_table.begin());
        d_tau_table_ptr = d_tau_table.data();
        amrex::Print() << "Built per-fuel burnout time table: "
                       << tau_table_size << " entries\n";
    }

    // ---------------- Time-varying fuel moisture schedule ------------------
    FuelMoistureSchedule fmd_sched;
    if (!inputs.fmd_file.empty()) {
        load_fuel_moisture_schedule(inputs.fmd_file, fmd_sched,
                                    inputs.fmd_start_year,
                                    inputs.fmd_start_month,
                                    inputs.fmd_start_day,
                                    inputs.fmd_start_hour);
    }

    // Helper lambda: update RothermelParams from the FMD schedule at time t_s.
    // When no FMD schedule is provided but the diurnal moisture model is enabled,
    // moisture is computed from the Nelson (2000) EMC diurnal cycle instead.
    // Rebuilds fuel and Balbi tables to reflect the new moisture values.
    auto apply_fmd_moisture = [&](Real t_s) {
        RothermelMoistures m;
        bool moisture_updated = false;

        if (!fmd_sched.empty()) {
            // FMD schedule takes priority
            m = get_moisture_at_time(
                fmd_sched, static_cast<double>(t_s), inputs.fmd_fuel_model,
                {static_cast<float>(inputs.rothermel.M_d1),
                 static_cast<float>(inputs.rothermel.M_d10),
                 static_cast<float>(inputs.rothermel.M_d100),
                 static_cast<float>(inputs.rothermel.M_d1000),
                 static_cast<float>(inputs.rothermel.M_lh),
                 static_cast<float>(inputs.rothermel.M_lw)});
            moisture_updated = true;
        } else if (inputs.diurnal_moisture.enable == 1) {
            // Diurnal EMC model (Nelson 2000)
            m = compute_diurnal_emc(
                inputs.diurnal_moisture,
                static_cast<double>(t_s),
                {static_cast<float>(inputs.rothermel.M_d1),
                 static_cast<float>(inputs.rothermel.M_d10),
                 static_cast<float>(inputs.rothermel.M_d100),
                 static_cast<float>(inputs.rothermel.M_d1000),
                 static_cast<float>(inputs.rothermel.M_lh),
                 static_cast<float>(inputs.rothermel.M_lw)});
            moisture_updated = true;
        }

        if (!moisture_updated) return;

        inputs.rothermel.M_d1    = static_cast<amrex::Real>(m.M_d1);
        inputs.rothermel.M_d10   = static_cast<amrex::Real>(m.M_d10);
        inputs.rothermel.M_d100  = static_cast<amrex::Real>(m.M_d100);
        inputs.rothermel.M_d1000 = static_cast<amrex::Real>(m.M_d1000);
        inputs.rothermel.M_lh    = static_cast<amrex::Real>(m.M_lh);
        inputs.rothermel.M_lw    = static_cast<amrex::Real>(m.M_lw);
        // Also keep the single-class M_f in sync with the 1-hr dead value
        inputs.rothermel.M_f    = static_cast<amrex::Real>(m.M_d1);

        // Rebuild per-cell Rothermel table with updated moisture
        if (!inputs.rothermel.landscape_file.empty() && fuel_table_size > 0) {
            std::vector<RothermelComputed> h_table =
                build_fuel_rothermel_table(inputs.rothermel,
                                           inputs.rothermel.landscape_fuel_type);
            if (!inputs.fuel_adj_file.empty()) {
                auto adjs = parse_fuel_adjustment_file(inputs.fuel_adj_file);
                apply_fuel_adjustment_to_table(h_table, adjs);
            }
            Gpu::copy(Gpu::hostToDevice,
                      h_table.begin(), h_table.end(),
                      d_fuel_table.begin());
        }
    };

    // Apply FMD at t=0 so initial plotfile uses correct moisture
    apply_fmd_moisture(Real(0.0));

    // ---- Solar radiation shading at t=0 (initial state) ----
    // Apply shade-adjusted EMC to spatial_moisture_mf for the initial plotfile.
    if (inputs.solar_radiation.enable == 1) {
        int sol_year0  = inputs.solar_radiation.year;
        int sol_month0 = inputs.solar_radiation.month;
        int sol_day0   = inputs.solar_radiation.day;
        // At t=0 no days have elapsed; hour_local = sim_start_hour
        amrex::Real hour_local0 = inputs.solar_radiation.sim_start_hour;
        SolarPosition sun0 = compute_solar_position(
            inputs.solar_radiation.latitude,
            inputs.solar_radiation.longitude,
            sol_year0, sol_month0, sol_day0,
            hour_local0,
            inputs.solar_radiation.timezone_offset);
        amrex::Print() << "Solar position at t=0: elevation="
                       << sun0.elevation_rad * (180.0 / WildfireConst::PI) << " deg"
                       << "  azimuth="
                       << sun0.azimuth_rad   * (180.0 / WildfireConst::PI) << " deg\n";
        const MultiFab* cc_for_shade0 = (inputs.solar_radiation.use_canopy_shading == 1
                                         && has_spatial_crown) ? &cc_mf : nullptr;
        compute_shade_fraction_mf(shade_fraction_mf,
                                  slope_mf, aspect_mf,
                                  cc_for_shade0,
                                  sun0,
                                  inputs.solar_radiation.use_canopy_shading,
                                  static_cast<amrex::Real>(inputs.solar_radiation.cloud_cover));
        if (inputs.diurnal_moisture.enable == 1) {
            const double phase0 = WildfireConst::TWO_PI *
                (inputs.diurnal_moisture.t_start_s - inputs.diurnal_moisture.t_T_peak_s)
                / WildfireConst::DAY_S;
            double T_mean0  = 0.5 * (inputs.diurnal_moisture.T_max + inputs.diurnal_moisture.T_min);
            double A_T0     = 0.5 * (inputs.diurnal_moisture.T_max - inputs.diurnal_moisture.T_min);
            double RH_mean0 = 0.5 * (inputs.diurnal_moisture.RH_max + inputs.diurnal_moisture.RH_min);
            double A_RH0    = 0.5 * (inputs.diurnal_moisture.RH_max - inputs.diurnal_moisture.RH_min);
            double T_C0  = T_mean0  + A_T0  * std::sin(phase0);
            double RH_p0 = RH_mean0 - A_RH0 * std::sin(phase0);
            T_C0  = std::max(-40.0, std::min(60.0,  T_C0));
            RH_p0 = std::max(1.0,   std::min(100.0, RH_p0));
            apply_solar_emc_to_spatial_moisture(
                spatial_moisture_mf,
                shade_fraction_mf,
                static_cast<amrex::Real>(T_C0),
                static_cast<amrex::Real>(RH_p0),
                inputs.solar_radiation.solar_heating_C);
        }
    }

    // ---------------- FMC seasonal schedule --------------------------------
    FMCSchedule fmc_sched;
    if (inputs.fmc_schedule.enable == 1) {
        if (!inputs.fmc_schedule.file.empty()) {
            load_fmc_schedule(inputs.fmc_schedule.file, fmc_sched,
                              inputs.fmc_schedule.start_doy);
        } else if (inputs.fmc_schedule.use_farsite_curve == 1) {
            build_farsite_fmc_curve(fmc_sched,
                                    inputs.fmc_schedule.start_doy,
                                    inputs.fmc_schedule.spring_start,
                                    inputs.fmc_schedule.summer_peak,
                                    inputs.fmc_schedule.fall_start,
                                    inputs.fmc_schedule.fall_end,
                                    static_cast<double>(inputs.fmc_schedule.fmc_min),
                                    static_cast<double>(inputs.fmc_schedule.fmc_max));
        }
        // Apply at t=0 to set initial FMC in crown params
        if (!fmc_sched.empty()) {
            inputs.crown.FMC = static_cast<amrex::Real>(
                get_fmc_at_time(fmc_sched, 0.0));
            amrex::Print() << "FMC at t=0: " << inputs.crown.FMC << " %\n";
        }
    }

    // ---------------- Live herbaceous moisture schedule ---------------------
    HerbMoistureSchedule herb_sched;
    if (inputs.herb_moisture_schedule.enable == 1) {
        if (!inputs.herb_moisture_schedule.file.empty()) {
            load_herb_moisture_schedule(inputs.herb_moisture_schedule.file,
                                        herb_sched,
                                        inputs.herb_moisture_schedule.start_doy);
        } else if (inputs.herb_moisture_schedule.use_curing_curve == 1) {
            build_herb_curing_curve(herb_sched,
                                    inputs.herb_moisture_schedule.start_doy,
                                    inputs.herb_moisture_schedule.spring_start,
                                    inputs.herb_moisture_schedule.summer_peak,
                                    inputs.herb_moisture_schedule.fall_start,
                                    inputs.herb_moisture_schedule.fall_end,
                                    static_cast<double>(inputs.herb_moisture_schedule.m_lh_min),
                                    static_cast<double>(inputs.herb_moisture_schedule.m_lh_max));
        }
        // Apply at t=0 to set initial M_lh
        if (!herb_sched.empty()) {
            double m_lh_pct0 = get_herb_moisture_at_time(herb_sched, 0.0);
            inputs.rothermel.M_lh = static_cast<amrex::Real>(m_lh_pct0 / 100.0);
            amrex::Print() << "Herb moisture at t=0: " << m_lh_pct0 << " % (M_lh="
                           << inputs.rothermel.M_lh << " fraction)\n";
        }
    }

    // ---------------- Compact wind direction schedule -----------------------
    WindDirSchedule wind_dir_sched;
    if (!inputs.wind_dir_schedule_file.empty()) {
        load_wind_dir_schedule(inputs.wind_dir_schedule_file, wind_dir_sched);
        // Apply at t=0
        auto [ux0, uy0] = get_wind_at_time(wind_dir_sched, 0.0);
        inputs.ux = static_cast<amrex::Real>(ux0);
        inputs.uy = static_cast<amrex::Real>(uy0);
        init_velocity_constant(vel, geom, inputs.ux, inputs.uy, inputs.uz);
    }

    // ---------------- FARSITE .wtr single-file hourly weather ---------------
    // When wtr_file is provided, parse it once and populate the wind schedule
    // and a WtrWeatherData structure.  At each timestep, T/RH/precip from the
    // .wtr data override the diurnal_moisture parameters.
    WtrWeatherData wtr_data;
    const bool wtr_active = !inputs.wtr_file.empty();
    if (wtr_active) {
        load_wtr_weather(inputs.wtr_file, wtr_data,
                         inputs.wtr_start_year, inputs.wtr_start_month,
                         inputs.wtr_start_day,  inputs.wtr_start_hour);
        // Use the .wtr wind schedule as the wind direction schedule (overrides
        // wind_dir_schedule_file if both are set, with a warning).
        if (!wind_dir_sched.empty()) {
            amrex::Print() << "WARNING: both wind_dir_schedule_file and wtr_file set; "
                              "wtr_file wind schedule takes precedence.\n";
        }
        wind_dir_sched = wtr_data.wind_sched;
        auto [ux_wtr0, uy_wtr0] = get_wind_at_time(wind_dir_sched, 0.0);
        inputs.ux = static_cast<amrex::Real>(ux_wtr0);
        inputs.uy = static_cast<amrex::Real>(uy_wtr0);
        init_velocity_constant(vel, geom, inputs.ux, inputs.uy, inputs.uz);
    }

    // ---------------- Aerial retardant drop zones ---------------------------
    RetardantDropList retardant_drops;
    if (!inputs.retardant_file.empty()) {
        load_retardant_drops(inputs.retardant_file, retardant_drops);
    }

    // ---------------- Multiple weather stations (IDW spatial interpolation) -
    // When multi_wtr_file is set, load all station .wtr files and apply IDW
    // wind interpolation at each timestep.  This supersedes the single wtr_file
    // wind schedule; T/RH/precip are taken as the domain-mean IDW centroid.
    MultiWtrWeather multi_wtr;
    const bool multi_wtr_active = !inputs.multi_wtr_file.empty();
    if (multi_wtr_active) {
        load_multi_wtr_stations(inputs.multi_wtr_file, multi_wtr,
                                inputs.wtr_start_year, inputs.wtr_start_month,
                                inputs.wtr_start_day,  inputs.wtr_start_hour,
                                static_cast<double>(inputs.multi_wtr_idw_power));
        // Apply IDW wind at t=0
        apply_multi_wtr_to_vel(vel, geom, multi_wtr, 0.0);
    }

    // ---------------- Multiple scheduled ignitions -------------------------
    IgnitionSchedule ignition_sched;
    if (!inputs.ignition_schedule_file.empty()) {
        load_ignition_schedule(inputs.ignition_schedule_file, ignition_sched);
        // Apply any t=0 events (time_s = 0)
        apply_scheduled_ignitions(phi, geom, ignition_sched,
                                  Real(0.0), Real(0.0), use_indicator);
        fill_boundary_extrap(phi, geom);
    }

    // ---------------- Precipitation wetting state --------------------------
    // Must be declared before the conditioning block which uses precip_state.
    PrecipState precip_state;
    precip_state.M_d1   = static_cast<float>(inputs.rothermel.M_d1);
    precip_state.M_d10  = static_cast<float>(inputs.rothermel.M_d10);
    precip_state.M_d100 = static_cast<float>(inputs.rothermel.M_d100);
    precip_state.M_d1000 = static_cast<float>(inputs.rothermel.M_d1000);
    precip_state.initialized = true;

    // ---------------- Fuel moisture conditioning period --------------------
    // Pre-run the Nelson (2000) diurnal EMC model for conditioning.n_days of
    // synthetic hourly weather to arrive at realistic starting moisture values.
    if (inputs.conditioning.n_days > 0 && inputs.diurnal_moisture.enable == 1) {
        amrex::Print() << "Fuel moisture conditioning: spinning up "
                       << inputs.conditioning.n_days << " day(s)...\n";
        const double cond_dt = 3600.0;  // 1-hour steps
        const int    n_steps_cond = inputs.conditioning.n_days * 24;

        // Use conditioning.wtr_file if provided, otherwise inputs.wtr_file,
        // otherwise fall back to diurnal_moisture parameters already loaded.
        WtrWeatherData cond_wtr;
        const bool cond_wtr_active =
            (!inputs.conditioning.wtr_file.empty()) ||
            (!inputs.wtr_file.empty() && wtr_active);
        if (!inputs.conditioning.wtr_file.empty()) {
            load_wtr_weather(inputs.conditioning.wtr_file, cond_wtr,
                             inputs.wtr_start_year, inputs.wtr_start_month,
                             inputs.wtr_start_day,  inputs.wtr_start_hour);
        } else if (wtr_active) {
            cond_wtr = wtr_data;
        }

        for (int ci = 0; ci < n_steps_cond; ++ci) {
            const double ct = static_cast<double>(ci) * cond_dt;
            float rain_rate = static_cast<float>(inputs.precip_rain_rate_mm_hr);
            double T_cond = 0.5 * (inputs.diurnal_moisture.T_max + inputs.diurnal_moisture.T_min);
            double RH_cond = 0.5 * (inputs.diurnal_moisture.RH_max + inputs.diurnal_moisture.RH_min);

            if (cond_wtr_active && !cond_wtr.empty()) {
                auto [T_wtr, RH_wtr] = cond_wtr.get_TRH_at_time(ct);
                T_cond  = T_wtr;
                RH_cond = RH_wtr;
                rain_rate = static_cast<float>(cond_wtr.get_precip_at_time(ct));
            }
            // Use Nelson EMC formula inline (Simard 1968)
            double rh_f = std::max(1.0, std::min(100.0, RH_cond)) / 100.0;
            double T_C  = std::max(-40.0, std::min(60.0, T_cond));
            double emc_pct;
            if (RH_cond < 10.0)
                emc_pct = 0.03229 + 0.281073*rh_f - 0.000578*T_C*rh_f;
            else if (RH_cond < 50.0)
                emc_pct = 2.22749 + 0.160107*rh_f - 0.014784*T_C;
            else
                emc_pct = 21.0606 + 0.005565*rh_f*rh_f - 0.00035*T_C*rh_f - 0.483199*rh_f;
            emc_pct = std::max(0.5, emc_pct);
            RothermelMoistures emc_cond;
            emc_cond.M_d1   = static_cast<float>(emc_pct / 100.0);
            emc_cond.M_d10  = static_cast<float>(emc_pct * 1.10 / 100.0);
            emc_cond.M_d100 = static_cast<float>(emc_pct * 1.30 / 100.0);
            // Live fuel moisture during conditioning: linearly ramp from initial
            // M_lh / M_lw toward the current dead-fuel equilibrium over the
            // conditioning period (FARSITE-style linear schedule, Finney 1998).
            // frac = ci / max(n_steps_cond - 1, 1)  → 0 at start, 1 at end.
            const double cond_frac = (n_steps_cond > 1)
                ? static_cast<double>(ci) / static_cast<double>(n_steps_cond - 1)
                : 1.0;
            const float M_lh_init = static_cast<float>(inputs.rothermel.M_lh);
            const float M_lw_init = static_cast<float>(inputs.rothermel.M_lw);
            const float M_lh_target = emc_cond.M_d100 * static_cast<float>(WildfireConst::COND_LH_EMC_MULT);
            const float M_lw_target = emc_cond.M_d100 * static_cast<float>(WildfireConst::COND_LW_EMC_MULT);
            emc_cond.M_lh = M_lh_init + static_cast<float>(cond_frac) * (M_lh_target - M_lh_init);
            emc_cond.M_lw = M_lw_init + static_cast<float>(cond_frac) * (M_lw_target - M_lw_init);
            apply_precipitation_moisture(precip_state, emc_cond, rain_rate,
                                         static_cast<float>(cond_dt),
                                         static_cast<float>(inputs.precip_threshold_mm_hr),
                                         static_cast<float>(inputs.M_sat));
        }
        // Apply conditioned moisture to Rothermel params
        inputs.rothermel.M_d1   = static_cast<amrex::Real>(precip_state.M_d1);
        inputs.rothermel.M_d10  = static_cast<amrex::Real>(precip_state.M_d10);
        inputs.rothermel.M_d100 = static_cast<amrex::Real>(precip_state.M_d100);
        inputs.rothermel.M_d1000 = static_cast<amrex::Real>(precip_state.M_d1000);
        inputs.rothermel.M_f    = static_cast<amrex::Real>(precip_state.M_d1);
        // Live fuel moisture from conditioning (linear ramp result)
        inputs.rothermel.M_lh   = static_cast<amrex::Real>(precip_state.M_d100 * static_cast<float>(WildfireConst::COND_LH_EMC_MULT));
        inputs.rothermel.M_lw   = static_cast<amrex::Real>(precip_state.M_d100 * static_cast<float>(WildfireConst::COND_LW_EMC_MULT));
        amrex::Print() << "Conditioning complete: M_d1=" << precip_state.M_d1
                       << " M_d10=" << precip_state.M_d10
                       << " M_d100=" << precip_state.M_d100
                       << " M_d1000=" << precip_state.M_d1000
                       << " M_lh=" << inputs.rothermel.M_lh
                       << " M_lw=" << inputs.rothermel.M_lw << "\n";
    }

    // Load precipitation time series if provided
    std::vector<std::pair<double,double>> precip_schedule;  // (time_s, rain_mm_hr)
    if (!inputs.precip_schedule_file.empty()) {
        std::ifstream pf(inputs.precip_schedule_file);
        if (!pf.is_open()) {
            amrex::Abort("Cannot open precip_schedule_file: " + inputs.precip_schedule_file);
        }
        std::string pline;
        while (std::getline(pf, pline)) {
            if (pline.empty() || pline[0] == '#' || pline[0] == '!') continue;
            for (char& c : pline) if (c == ',') c = ' ';
            std::istringstream piss(pline);
            double ts, rate;
            if (piss >> ts >> rate) precip_schedule.push_back({ts, rate});
        }
        amrex::Print() << "Precipitation schedule: loaded " << precip_schedule.size()
                       << " records from " << inputs.precip_schedule_file << "\n";
    }

    // ---------------- Cruz crown pre-computed coefficients ------------------
    // Pre-compute Cruz, Alexander & Wakimoto (2005) crown fire ROS coefficients
    // once at startup.  These are used both in the FARSITE ellipse kernel and in
    // the Cap-8 crown ROS override when crown.use_cruz_crown = 1.
    CruzCrownComputed ccc_global = compute_cruz_crown_params(inputs.cruz_crown);
    const CruzCrownComputed* ccc_ptr =
        (inputs.crown.use_cruz_crown == 1) ? &ccc_global : nullptr;

    // ---------------- Balbi (2009) lookup table (if enabled) -------------
    Gpu::DeviceVector<BalbiComputed> d_balbi_table;
    const BalbiComputed* d_balbi_table_ptr = nullptr;
    int balbi_table_size = 0;
    // Store global default Balbi coefficients for Viegas+Balbi coupling
    BalbiComputed bc_global_default;
    bc_global_default.A_coeff = 0.0;
    bc_global_default.v_b     = 1.0;

    if (inputs.fire_spread_model == "balbi") {
        // Print Balbi fuel parameter table
        print_balbi_fuel_table(inputs.rothermel, inputs.balbi,
                               inputs.rothermel.landscape_fuel_type);
        bc_global_default = compute_balbi_params(inputs.rothermel, inputs.balbi);

        if (!inputs.rothermel.landscape_file.empty()) {
            std::vector<BalbiComputed> h_balbi_table =
                build_fuel_balbi_table(inputs.rothermel, inputs.balbi,
                                       inputs.rothermel.landscape_fuel_type);
            balbi_table_size = static_cast<int>(h_balbi_table.size());
            d_balbi_table.resize(balbi_table_size);
            Gpu::copy(Gpu::hostToDevice,
                      h_balbi_table.begin(), h_balbi_table.end(),
                      d_balbi_table.begin());
            d_balbi_table_ptr = d_balbi_table.data();
            amrex::Print() << "Built per-cell Balbi lookup table: "
                           << balbi_table_size << " entries ("
                           << (inputs.rothermel.landscape_fuel_type == "40"
                               ? "FBFM40" : "FBFM13")
                           << ")\n";
        }
    }

    // ---------------- Heat flux MultiFab initialization ------------------
    // Initialize from file (2D only) or uniform value
    {
        const bool hf_active = (inputs.heat_flux.enable_upward == 1 ||
                                 inputs.heat_flux.enable_induced == 1);
        if (hf_active) {
#if (AMREX_SPACEDIM == 2)
            if (!inputs.heat_flux.heat_flux_file.empty()) {
                init_heat_flux_from_file(heat_flux_mf, geom,
                                         inputs.heat_flux.heat_flux_file);
                amrex::Print() << "Initialized heat flux from file: "
                               << inputs.heat_flux.heat_flux_file << "\n";
            } else {
                init_heat_flux_from_value(heat_flux_mf,
                                          inputs.heat_flux.heat_flux_value);
                amrex::Print() << "Initialized uniform heat flux: "
                               << inputs.heat_flux.heat_flux_value << " W/m2\n";
            }
#else
            // In 3D, only uniform value is supported for now
            init_heat_flux_from_value(heat_flux_mf,
                                      inputs.heat_flux.heat_flux_value);
            if (!inputs.heat_flux.heat_flux_file.empty()) {
                amrex::Print() << "WARNING: heat_flux.file is only supported in 2D builds; "
                                  "using heat_flux.value instead.\n";
            } else {
                amrex::Print() << "Initialized uniform heat flux: "
                               << inputs.heat_flux.heat_flux_value << " W/m2\n";
            }
#endif
        }
    }
    const bool heat_flux_active = (inputs.heat_flux.enable_upward == 1 ||
                                   inputs.heat_flux.enable_induced == 1);

    // ---------------- Wind-terrain model setup ------------------
    // wind_terrain_modifies_vel: true for Options 3-7 which produce vel_effective.
    // For "none" (Option 1) and "viegas_ros" (Option 2) the original vel is used.
    const bool wind_terrain_modifies_vel =
        (inputs.wind_terrain.model == "viegas_wind" ||
         inputs.wind_terrain.model == "canyon_wind"  ||
         inputs.wind_terrain.model == "viegas_neto"  ||
         inputs.wind_terrain.model == "pimont"        ||
         inputs.wind_terrain.model == "windninja_ridge_canyon");

    // use_precomp_R_for_advection: when a wind-terrain model or a non-Rothermel
    // spread model is active, pass R_mf as the pre-computed ROS to advection so
    // that the advection kernel does not internally recompute Rothermel with the
    // unmodified vel.  This ensures the advected ROS matches what was computed
    // above (including any terrain-corrected velocity or Viegas ROS override).
    const bool use_precomp_R_for_advection =
        (wind_terrain_modifies_vel ||
         heat_flux_active                              ||
         inputs.wind_terrain.model == "viegas_ros"    ||
         inputs.fire_spread_model  == "balbi"         ||
         inputs.fire_spread_model  == "cheney_gould"  ||
         inputs.fire_spread_model  == "cruz_crown"    ||
         inputs.fire_spread_model  == "fbp_o1a"       ||
         inputs.fire_spread_model  == "fbp_o1b"       ||
         inputs.fire_spread_model  == "fbp_s1"        ||
         inputs.fire_spread_model  == "fbp_s2"        ||
         inputs.fire_spread_model  == "fbp_s3"        ||
         inputs.fire_spread_model  == "lautenberger");

    // Helper flag for Viegas+Balbi coupling
    const bool use_balbi_for_viegas = (inputs.fire_spread_model == "balbi");

    // ---------------- dt from CFL --------------------------
    Real dt=10;
    const bool use_levelset = (inputs.propagation_method == "levelset");
    const bool use_mtt      = (inputs.propagation_method == "mtt");
    if (use_levelset)
      {
        // Apply wind-terrain velocity modification (Options 3-7) before ROS computation
        if (wind_terrain_modifies_vel) {
            apply_wind_terrain_velocity(vel_effective, vel, terrain_slopes.get(), inputs);
        } else {
            // Copy ambient wind into vel_effective (needed for heat flux application below)
            MultiFab::Copy(vel_effective, vel, 0, 0, 3, 0);
        }

        // Apply heat flux wind corrections (upward velocity + induced inflow)
        if (heat_flux_active) {
            apply_heatflux_wind(vel_effective, vel, heat_flux_mf, &phi,
                                inputs.heat_flux);
        }

        const MultiFab& vel_for_model = (wind_terrain_modifies_vel || heat_flux_active)
                                        ? vel_effective : vel;

        if (inputs.fire_spread_model == "balbi") {
            compute_balbi_R(R_mf, vel_for_model, geom, inputs.rothermel, inputs.balbi,
                             terrain_slopes.get(),
                             !inputs.rothermel.landscape_file.empty() ? &fuel_model_mf : nullptr,
                             d_balbi_table_ptr, balbi_table_size,
                             heat_flux_active ? &heat_flux_mf : nullptr,
                             heat_flux_active ? &inputs.heat_flux : nullptr);
        } else if (inputs.fire_spread_model == "cheney_gould") {
            compute_cheney_gould_R(R_mf, vel_for_model, inputs.cheney_gould);
        } else if (inputs.fire_spread_model == "cruz_crown") {
            compute_cruz_crown_R(R_mf, vel_for_model, inputs.cruz_crown);
        } else if (inputs.fire_spread_model == "fbp_o1a" ||
                   inputs.fire_spread_model == "fbp_o1b" ||
                   inputs.fire_spread_model == "fbp_s1"  ||
                   inputs.fire_spread_model == "fbp_s2"  ||
                   inputs.fire_spread_model == "fbp_s3") {
            compute_fbp_R(R_mf, vel_for_model, inputs.fbp);
        } else if (inputs.fire_spread_model == "lautenberger") {
            compute_lautenberger_R(R_mf, vel_for_model, inputs.rothermel, inputs.lautenberger,
                                    terrain_slopes.get());
        } else {
            // Compute Rothermel wind speed R
            compute_rothermel_R(R_mf, vel_for_model, geom, inputs.rothermel,
                                 terrain_slopes.get(),
                                 !inputs.rothermel.landscape_file.empty() ? &fuel_model_mf : nullptr,
                                 d_fuel_table_ptr, fuel_table_size,
                                 has_spatial_crown ? &cc_mf : nullptr,
                                 has_spatial_crown ? &canopy_height_mf : nullptr);
        }
        dt = compute_dt(R_mf, geom, inputs.cfl);
        amrex::Print() << "Computed dt = " << dt << "\n";
      } else if (use_mtt) {
      // For MTT: compute the ROS field now (same as levelset path),
      // then run the Dijkstra fast-march to fill arrival_time_mf,
      // and set a sensible dt based on the ROS field.
      {
        if (wind_terrain_modifies_vel) {
            apply_wind_terrain_velocity(vel_effective, vel, terrain_slopes.get(), inputs);
        } else {
            MultiFab::Copy(vel_effective, vel, 0, 0, 3, 0);
        }
        if (heat_flux_active) {
            apply_heatflux_wind(vel_effective, vel, heat_flux_mf, &phi, inputs.heat_flux);
        }
        const MultiFab& vel_for_mtt = (wind_terrain_modifies_vel || heat_flux_active)
                                      ? vel_effective : vel;
        if (inputs.fire_spread_model == "balbi") {
            compute_balbi_R(R_mf, vel_for_mtt, geom, inputs.rothermel, inputs.balbi,
                             terrain_slopes.get(),
                             !inputs.rothermel.landscape_file.empty() ? &fuel_model_mf : nullptr,
                             d_balbi_table_ptr, balbi_table_size,
                             heat_flux_active ? &heat_flux_mf : nullptr,
                             heat_flux_active ? &inputs.heat_flux : nullptr);
        } else if (inputs.fire_spread_model == "cheney_gould") {
            compute_cheney_gould_R(R_mf, vel_for_mtt, inputs.cheney_gould);
        } else if (inputs.fire_spread_model == "cruz_crown") {
            compute_cruz_crown_R(R_mf, vel_for_mtt, inputs.cruz_crown);
        } else if (inputs.fire_spread_model == "fbp_o1a" ||
                   inputs.fire_spread_model == "fbp_o1b" ||
                   inputs.fire_spread_model == "fbp_s1"  ||
                   inputs.fire_spread_model == "fbp_s2"  ||
                   inputs.fire_spread_model == "fbp_s3") {
            compute_fbp_R(R_mf, vel_for_mtt, inputs.fbp);
        } else if (inputs.fire_spread_model == "lautenberger") {
            compute_lautenberger_R(R_mf, vel_for_mtt, inputs.rothermel, inputs.lautenberger,
                                    terrain_slopes.get());
        } else {
            compute_rothermel_R(R_mf, vel_for_mtt, geom, inputs.rothermel,
                                 terrain_slopes.get(),
                                 !inputs.rothermel.landscape_file.empty() ? &fuel_model_mf : nullptr,
                                 d_fuel_table_ptr, fuel_table_size,
                                 has_spatial_crown ? &cc_mf : nullptr,
                                 has_spatial_crown ? &canopy_height_mf : nullptr);
        }
        dt = compute_dt(R_mf, geom, inputs.cfl);
        amrex::Print() << "MTT: pre-computing arrival times (dt = " << dt << ") ...\n";
        compute_mtt_arrival_times(arrival_time_mf, phi, R_mf, geom);
        // Set phi from arrival times at t=0
        apply_mtt_phi_update(phi, arrival_time_mf, Real(0.0));
        fill_boundary_extrap(phi, geom);
      }
      } else {
      amrex::Print() << "Using FARSITE propagation; dt = " << dt << "\n";
    }
    compute_fire_behavior(fireline_intensity_mf, flame_length_mf, R_mf, inputs.rothermel);
    if (inputs.weise_biging.enable == 1) {
        compute_weise_biging_whirl(weise_data, fireline_intensity_mf, flame_length_mf,
                                   vel, terrain_slopes.get(), inputs.weise_biging);
    }
    if (inputs.viegas.enable == 1) {
        // For Balbi+Viegas: pass Balbi table so diagnostic uses Balbi amplitude
        compute_viegas_diagnostics(viegas_data, R_mf, vel, inputs.rothermel, inputs.viegas,
                                   terrain_slopes.get(),
                                   !inputs.rothermel.landscape_file.empty() ? &fuel_model_mf : nullptr,
                                   d_fuel_table_ptr, fuel_table_size,
                                   use_balbi_for_viegas,
                                   use_balbi_for_viegas ? &bc_global_default : nullptr,
                                   use_balbi_for_viegas ? d_balbi_table_ptr : nullptr,
                                   use_balbi_for_viegas ? balbi_table_size : 0);
        // Option 2: override R_mf with max(primary, R_viegas) in eruptive cells
        if (inputs.wind_terrain.model == "viegas_ros") {
            apply_viegas_ros_override(R_mf, viegas_data);
            if (use_levelset) dt = compute_dt(R_mf, geom, inputs.cfl);
        }
    }
    // Fire ecology diagnostics (scorch height, prob. ignition, tree mortality,
    // crown activity) – always computed, always written to plotfile
    compute_fire_ecology(ecology_mf, fireline_intensity_mf, R_mf,
                         inputs.rothermel, inputs.fire_ecology, inputs.crown);

    // Ecology → propagation coupling (initial step):
    //   (a) P_ignition scales R_mf when fire_ecology.couple_to_ros = 1
    //   (b) Active crown fire (crown_activity == 2) overrides R_mf with the
    //       Cruz et al. (2005) crown ROS (when crown.use_cruz_crown = 1) or
    //       the Van Wagner 3/CBD proxy (when crown.use_cruz_crown = 0).
    apply_ecology_p_ignition_feedback(R_mf, ecology_mf, phi, inputs.fire_ecology);

    // Cap 8: Active crown fire ROS feedback (initial dt computation).
    // Apply the same crown-ROS override used in the time loop so the initial
    // CFL dt correctly accounts for any active crown fire cells.
    if (inputs.crown.enable == 1 && use_levelset) {
        const Real CBD_global_i = Real(inputs.crown.CBD);
        const Real FMC_global_i = Real(inputs.crown.FMC);
        const Real mf_i = amrex::max(Real(0.3),
                              amrex::min(Real(1.0),
                              Real(1.0) - (FMC_global_i - Real(100.0)) / Real(200.0)));
        // Global crown ROS (Van Wagner proxy) for dt guard
        const Real R_crown_g_ms_init = (Real(3.0) / amrex::max(CBD_global_i, Real(0.01)))
                                       * mf_i / Real(60.0);
        const bool use_cruz_init          = (inputs.crown.use_cruz_crown == 1);
        const bool use_roth1991_init      = (inputs.crown.use_rothermel1991_crown == 1);
        const bool use_passive_blend_init = (inputs.crown.use_passive_blend == 1);
        const Real CBH_init = Real(inputs.crown.CBH);
        const Real FMC_init = FMC_global_i;
        const Real I_o_init = Real(0.010) * CBH_init * (Real(460.0) + Real(25.9) * FMC_init);
        CruzCrownComputed ccc_init;
        if (use_cruz_init) { ccc_init = ccc_global; }

        for (MFIter mfi(R_mf); mfi.isValid(); ++mfi) {
            const Box& bx  = mfi.validbox();
            auto       R   = R_mf.array(mfi);
            auto const eco = ecology_mf.const_array(mfi);
            auto const v   = vel.const_array(mfi);
            auto const fi  = fireline_intensity_mf.const_array(mfi);
            const bool use_sp_i = has_spatial_crown;
            Array4<const Real> cbd_arr_i;
            if (use_sp_i) {
                cbd_arr_i = cbd_mf.const_array(mfi);
            }
            const Real CBD_g_i  = CBD_global_i;
            const Real mf_val_i = mf_i;
            const Real MC10_i   = Real(inputs.cruz_crown.MC10);
            const bool use_cruz_i    = use_cruz_init;
            const bool use_roth91_i  = use_roth1991_init;
            const bool use_passive_i = use_passive_blend_init;
            const Real I_o_i = I_o_init;
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                if (eco(i, j, k, 3) < Real(1.5)) return; // surface or passive crown
                Real R_surface = R(i, j, k);
                Real R_crown_ms;
                if (use_roth91_i) {
                    // Rothermel (1991): R_crown = 3.34 x R_surface
                    R_crown_ms = compute_rothermel_1991_crown_ros(R_surface);
                } else if (use_cruz_i) {
                    // Cruz, Alexander & Wakimoto (2005) active crown ROS
                    Real ux = v(i,j,k,0);
                    Real uy = v(i,j,k,1);
                    Real wind_mag = std::sqrt(ux*ux + uy*uy);
                    Real CBD_c = use_sp_i ? cbd_arr_i(i,j,k) : CBD_g_i;
                    R_crown_ms = compute_crown_fire_spread_rate_cruz(wind_mag, CBD_c, MC10_i);
                } else {
                    // Van Wagner (1977) simplified proxy
                    Real CBD_c = use_sp_i ? cbd_arr_i(i,j,k) : CBD_g_i;
                    CBD_c = amrex::max(CBD_c, Real(0.01));
                    R_crown_ms = (Real(3.0) / CBD_c) * mf_val_i / Real(60.0);
                }
                if (use_passive_i) {
                    const Real I_B_kwm = fi(i, j, k);
                    R(i, j, k) = compute_van_wagner_passive_blend(
                        R_surface, R_crown_ms, I_B_kwm, I_o_i);
                } else {
                    R(i, j, k) = amrex::max(R_surface, R_crown_ms);
                }
            });
        }
        // Recompute dt when crown ROS is positive and could tighten the CFL.
        if (R_crown_g_ms_init > Real(0.0)) {
            dt = compute_dt(R_mf, geom, inputs.cfl);
        }
    }

    // Fire emissions (CO₂, CO, PM₂.₅) from fuel load × consumption fraction
    compute_fire_emissions(emissions_mf, phi, fuel_consumption_mf,
                           inputs.rothermel, inputs.emissions);

    // Scott & Reinhardt (2001) full bisection-based TI/CI (optional, host-only)
    if (inputs.scott_reinhardt_full.enable == 1) {
        compute_full_scott_reinhardt(
            ti_full_mf, ci_full_mf,
            inputs.rothermel, inputs.crown,
            !inputs.rothermel.landscape_file.empty() ? &fuel_model_mf : nullptr,
            nullptr, 0,
            amrex::Real(inputs.scott_reinhardt_full.U_max_kmh));
    }


    // ---------------- Barrier polygon cells (firebreaks) ------------------
    // Load once at startup; applied every time step inside the time loop.
    std::vector<IntVect> barrier_cells =
        load_barrier_cells(inputs.barrier_files, geom);

    // Apply barriers to the initial phi (they are present from t=0)
    if (!barrier_cells.empty())
        apply_barrier_polygons(phi, barrier_cells, geom);

    // ---------------- Initialize fire statistics CSV -------------------
    if (!inputs.fire_stats_file.empty())
        write_fire_stats_header(inputs.fire_stats_file);

    // ---------------- Write initial plotfile ---------------
    {
      // Compute per-plotfile diagnostics
      {
        const RothermelComputed rc_plt = compute_rothermel_params(inputs.rothermel);
        compute_heat_per_unit_area(heat_per_unit_area_mf, phi, rc_plt.I_R,
                                   inputs.farsite.tau_residence);
      }
      compute_burnout_time(burnout_time_mf, arrival_time_mf, inputs.farsite.tau_residence,
                       !inputs.rothermel.landscape_file.empty() ? &fuel_model_mf : nullptr,
                       d_tau_table_ptr, tau_table_size);
      compute_vorticity(vorticity_mf, vel, geom);
      compute_reaction_intensity(reaction_intensity_mf,
                                 !inputs.rothermel.landscape_file.empty() ? &fuel_model_mf : nullptr,
                                 inputs.rothermel, d_fuel_table_ptr, fuel_table_size);
      Vector<std::string> names = {"phi", "velx", "vely"
#if (AMREX_SPACEDIM == 3)
				   , "velz", "farsite_dx", "farsite_dy", "farsite_dz", "R",
				   "spot_prob", "spot_count", "spot_dist", "spot_active", "fuel_consumption", "crown_fraction",
				   "albini_Hz", "albini_count", "albini_dist", "albini_active",
				   "elevation", "slope", "aspect", "fuel_model",
				   "fireline_intensity", "flame_length",
				   "weise_flame_height", "weise_flame_tilt",
				   "weise_whirl_height", "weise_whirl_radius",
				   "weise_angular_velocity", "weise_max_tang_vel",
				   "viegas_ROS", "viegas_eruptive_flag",
				   "viegas_ROS_excess", "viegas_flame_tilt", "viegas_slope_factor",
				   "scorch_height", "prob_ignition", "tree_mortality", "crown_activity",
					   "torching_ratio", "crowning_ratio", "energy_release_component",
				   "co2_emissions", "co_emissions", "pm25_emissions",
				   "arrival_time", "heat_per_unit_area", "vorticity_z",
				   "cbh", "cbd", "canopy_cover", "canopy_height", "burnout_time",
   "reaction_intensity", "residual_fuel", "shade_fraction",
   "torching_index_kmh", "crowning_index_kmh",
   "moisture_d1", "moisture_d10", "moisture_d100", "moisture_lh", "moisture_lw"
#else
				   , "farsite_dx", "farsite_dy", "R",
				   "spot_prob", "spot_count", "spot_dist", "spot_active", "fuel_consumption", "crown_fraction",
				   "albini_Hz", "albini_count", "albini_dist", "albini_active",
				   "elevation", "slope", "aspect", "fuel_model",
				   "fireline_intensity", "flame_length",
				   "weise_flame_height", "weise_flame_tilt",
				   "weise_whirl_height", "weise_whirl_radius",
				   "weise_angular_velocity", "weise_max_tang_vel",
				   "viegas_ROS", "viegas_eruptive_flag",
				   "viegas_ROS_excess", "viegas_flame_tilt", "viegas_slope_factor",
				   "scorch_height", "prob_ignition", "tree_mortality", "crown_activity",
					   "torching_ratio", "crowning_ratio", "energy_release_component",
				   "co2_emissions", "co_emissions", "pm25_emissions",
				   "arrival_time", "heat_per_unit_area", "vorticity_z",
				   "cbh", "cbd", "canopy_cover", "canopy_height", "burnout_time",
   "reaction_intensity", "residual_fuel", "shade_fraction",
   "torching_index_kmh", "crowning_index_kmh",
   "moisture_d1", "moisture_d10", "moisture_d100", "moisture_lh", "moisture_lw"
#endif
      };
      MultiFab plotmf(ba, dm, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 1 + 5 + WEISE_NCOMP + VIEGAS_NCOMP + ECOLOGY_NCOMP + EMISSIONS_NCOMP + 13 + FMS_NCOMP, 0);
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
      MultiFab::Copy(plotmf, weise_data, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 1 + 5, WEISE_NCOMP, 0);
      MultiFab::Copy(plotmf, viegas_data, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 1 + 5 + WEISE_NCOMP, VIEGAS_NCOMP, 0);
      MultiFab::Copy(plotmf, ecology_mf,   0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 1 + 5 + WEISE_NCOMP + VIEGAS_NCOMP, ECOLOGY_NCOMP, 0);
      MultiFab::Copy(plotmf, emissions_mf, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 1 + 5 + WEISE_NCOMP + VIEGAS_NCOMP + ECOLOGY_NCOMP, EMISSIONS_NCOMP, 0);
      MultiFab::Copy(plotmf, arrival_time_mf,       0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 1 + 5 + WEISE_NCOMP + VIEGAS_NCOMP + ECOLOGY_NCOMP + EMISSIONS_NCOMP, 1, 0);
      MultiFab::Copy(plotmf, heat_per_unit_area_mf, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 1 + 5 + WEISE_NCOMP + VIEGAS_NCOMP + ECOLOGY_NCOMP + EMISSIONS_NCOMP + 1, 1, 0);
      MultiFab::Copy(plotmf, vorticity_mf,          0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 1 + 5 + WEISE_NCOMP + VIEGAS_NCOMP + ECOLOGY_NCOMP + EMISSIONS_NCOMP + 2, 1, 0);
      MultiFab::Copy(plotmf, cbh_mf, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 1 + 5 + WEISE_NCOMP + VIEGAS_NCOMP + ECOLOGY_NCOMP + EMISSIONS_NCOMP + 3, 1, 0);
      MultiFab::Copy(plotmf, cbd_mf, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 1 + 5 + WEISE_NCOMP + VIEGAS_NCOMP + ECOLOGY_NCOMP + EMISSIONS_NCOMP + 4, 1, 0);
      MultiFab::Copy(plotmf, cc_mf,            0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 1 + 5 + WEISE_NCOMP + VIEGAS_NCOMP + ECOLOGY_NCOMP + EMISSIONS_NCOMP + 5, 1, 0);
      MultiFab::Copy(plotmf, canopy_height_mf, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 1 + 5 + WEISE_NCOMP + VIEGAS_NCOMP + ECOLOGY_NCOMP + EMISSIONS_NCOMP + 6, 1, 0);
      MultiFab::Copy(plotmf, burnout_time_mf,  0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 1 + 5 + WEISE_NCOMP + VIEGAS_NCOMP + ECOLOGY_NCOMP + EMISSIONS_NCOMP + 7, 1, 0);
      MultiFab::Copy(plotmf, reaction_intensity_mf, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 1 + 5 + WEISE_NCOMP + VIEGAS_NCOMP + ECOLOGY_NCOMP + EMISSIONS_NCOMP + 8, 1, 0);
	MultiFab::Copy(plotmf, residual_fuel_mf, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 1 + 5 + WEISE_NCOMP + VIEGAS_NCOMP + ECOLOGY_NCOMP + EMISSIONS_NCOMP + 9, 1, 0);
	MultiFab::Copy(plotmf, shade_fraction_mf, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 1 + 5 + WEISE_NCOMP + VIEGAS_NCOMP + ECOLOGY_NCOMP + EMISSIONS_NCOMP + 10, 1, 0);
	MultiFab::Copy(plotmf, ti_full_mf, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 1 + 5 + WEISE_NCOMP + VIEGAS_NCOMP + ECOLOGY_NCOMP + EMISSIONS_NCOMP + 11, 1, 0);
	MultiFab::Copy(plotmf, ci_full_mf, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 1 + 5 + WEISE_NCOMP + VIEGAS_NCOMP + ECOLOGY_NCOMP + EMISSIONS_NCOMP + 12, 1, 0);
      MultiFab::Copy(plotmf, spatial_moisture_mf, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 1 + 5 + WEISE_NCOMP + VIEGAS_NCOMP + ECOLOGY_NCOMP + EMISSIONS_NCOMP + 13, FMS_NCOMP, 0);
      {
        char buf[64];
        std::snprintf(buf, sizeof(buf), "plt%04d", restart_step);
        WriteSingleLevelPlotfile(buf, plotmf, names, geom, time, restart_step);
      }
      // Print burned area and perimeter statistics
      {
        amrex::Long n_burned = 0, n_perim = 0;
        const auto dxx = geom.CellSize();
        for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
          const Box& bx = mfi.validbox();
          auto const p = phi.const_array(mfi);
          // Bounds-guarded neighbor access (valid box only; ghost cells
          // at MPI subdomain edges are filled but domain-edge ghosts may
          // not be meaningful for the perimeter check, so we stay in-box).
          const IntVect bxlo = bx.smallEnd();
          const IntVect bxhi = bx.bigEnd();
          amrex::LoopOnCpu(bx, [&](int i, int j, int k) {
            if (p(i,j,k,0) < Real(0.0)) {
              ++n_burned;
              bool on_perim = false;
              if (i > bxlo[0]) on_perim |= (p(i-1,j,k,0) >= Real(0.0));
              if (i < bxhi[0]) on_perim |= (p(i+1,j,k,0) >= Real(0.0));
              if (j > bxlo[1]) on_perim |= (p(i,j-1,k,0) >= Real(0.0));
              if (j < bxhi[1]) on_perim |= (p(i,j+1,k,0) >= Real(0.0));
              if (on_perim) ++n_perim;
            }
          });
        }
        ParallelDescriptor::ReduceLongSum(n_burned);
        ParallelDescriptor::ReduceLongSum(n_perim);
        Real cell_area = dxx[0] * dxx[1];
        amrex::Print() << "  Burned area: " << Real(n_burned)*cell_area/1.0e4 << " ha"
                       << "  Perimeter: "   << Real(n_perim)*std::sqrt(cell_area)/1000.0 << " km\n";
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

    // ---- Timed isochrone state ----
    // Track the last isochrone index written so we fire exactly once per interval.
    int last_isochrone_idx = -1;
    const bool isochrone_active = (inputs.isochrone_interval_s > amrex::Real(0.0));
    while ((use_final_time && time < inputs.final_time) ||
           (!use_final_time && step < restart_step + inputs.nsteps)) {
      ++step;
      fill_boundary_extrap(phi, geom);
      const Real dt_step = dt;
      amrex::Print() << "Time:"<< time << " with timestep:" << dt_step <<std::endl;
      
      // Update time-dependent wind field if enabled
#if (AMREX_SPACEDIM == 2)
      if (!inputs.velocity_file.empty() && inputs.use_time_dependent_wind == 1) {
        // When turbulence is active, reload into vel_base; apply_turb_wind will
        // compute vel = vel_base + perturbation immediately below.
        MultiFab& wind_target = turb_wind_active ? *vel_base : vel;
        update_time_dependent_velocity(wind_target, geom, inputs.velocity_file, time, inputs.wind_time_spacing,
                                        wind_x_data1, wind_y_data1, wind_u_data1, wind_v_data1,
                                        wind_x_data2, wind_y_data2, wind_u_data2, wind_v_data2,
                                        current_wind_field_index, next_wind_field_index);
      }
#endif

      // Apply turbulent wind perturbation (vel = vel_base + stochastic perturbation)
      if (turb_wind_active) {
          apply_turb_wind(vel, *vel_base, ou_state_mf.get(), turb_state,
                          dt_step, inputs.turb_wind, geom);
      }

      // Update time-varying fuel moisture from FMD schedule
      apply_fmd_moisture(time);

      // ---- Solar radiation shading → per-cell shade-adjusted EMC ----
      // Computes the sun's position for the current simulation time, derives
      // per-cell terrain and canopy shade fractions, then writes shade-adjusted
      // dead fuel moistures into spatial_moisture_mf (components 0-2).
      // Requires diurnal_moisture.enable = 1 to supply T_air and RH.
      if (inputs.solar_radiation.enable == 1) {
          // Advance calendar date by elapsed simulation time
          int sol_year  = inputs.solar_radiation.year;
          int sol_month = inputs.solar_radiation.month;
          int sol_day   = inputs.solar_radiation.day;
          amrex::Real hour_local = advance_solar_date(
              sol_year, sol_month, sol_day,
              inputs.solar_radiation.sim_start_hour,
              static_cast<amrex::Real>(time));

          // Compute solar position
          SolarPosition sun = compute_solar_position(
              inputs.solar_radiation.latitude,
              inputs.solar_radiation.longitude,
              sol_year, sol_month, sol_day,
              hour_local,
              inputs.solar_radiation.timezone_offset);

          // Per-cell terrain + canopy shade fraction
          const MultiFab* cc_for_shade = (inputs.solar_radiation.use_canopy_shading == 1
                                          && has_spatial_crown) ? &cc_mf : nullptr;
          compute_shade_fraction_mf(shade_fraction_mf,
                                    slope_mf, aspect_mf,
                                    cc_for_shade,
                                    sun,
                                    inputs.solar_radiation.use_canopy_shading,
                                    static_cast<amrex::Real>(inputs.solar_radiation.cloud_cover));

          // Per-cell shade-adjusted EMC (requires diurnal T_air / RH)
          if (inputs.diurnal_moisture.enable == 1) {
              // Compute diurnal T_air and RH at current time (same formulas as
              // compute_diurnal_emc in fuel_moisture_scheduler.H)
              const double phase = WildfireConst::TWO_PI *
                  (inputs.diurnal_moisture.t_start_s +
                   static_cast<double>(time) -
                   inputs.diurnal_moisture.t_T_peak_s)
                  / WildfireConst::DAY_S;
              const double T_mean  = 0.5 * (inputs.diurnal_moisture.T_max +
                                            inputs.diurnal_moisture.T_min);
              const double A_T     = 0.5 * (inputs.diurnal_moisture.T_max -
                                            inputs.diurnal_moisture.T_min);
              const double RH_mean = 0.5 * (inputs.diurnal_moisture.RH_max +
                                            inputs.diurnal_moisture.RH_min);
              const double A_RH    = 0.5 * (inputs.diurnal_moisture.RH_max -
                                            inputs.diurnal_moisture.RH_min);
              double T_C  = T_mean  + A_T  * std::sin(phase);
              double RH_p = RH_mean - A_RH * std::sin(phase);
              T_C  = std::max(-40.0, std::min(60.0,  T_C));
              RH_p = std::max(1.0,   std::min(100.0, RH_p));

              apply_solar_emc_to_spatial_moisture(
                  spatial_moisture_mf,
                  shade_fraction_mf,
                  static_cast<amrex::Real>(T_C),
                  static_cast<amrex::Real>(RH_p),
                  inputs.solar_radiation.solar_heating_C);
          }
      }

      // Update FMC seasonal schedule (updates crown.FMC used by Van Wagner model)
      if (inputs.fmc_schedule.enable == 1 && !fmc_sched.empty()) {
          inputs.crown.FMC = static_cast<amrex::Real>(
              get_fmc_at_time(fmc_sched, static_cast<double>(time)));
      }

      // Update live herbaceous moisture from curing schedule
      if (inputs.herb_moisture_schedule.enable == 1 && !herb_sched.empty()) {
          double m_lh_pct = get_herb_moisture_at_time(herb_sched, static_cast<double>(time));
          inputs.rothermel.M_lh = static_cast<amrex::Real>(m_lh_pct / 100.0);
          // Propagate into the spatial moisture MultiFab component 3 (M_lh)
          spatial_moisture_mf.setVal(inputs.rothermel.M_lh, 3, 1, 0);
          // Rebuild fuel lookup table if per-cell landscape is active
          if (!inputs.rothermel.landscape_file.empty() && fuel_table_size > 0) {
              std::vector<RothermelComputed> h_herb_table =
                  build_fuel_rothermel_table(inputs.rothermel,
                                             inputs.rothermel.landscape_fuel_type);
              if (!inputs.fuel_adj_file.empty()) {
                  auto adjs = parse_fuel_adjustment_file(inputs.fuel_adj_file);
                  apply_fuel_adjustment_to_table(h_herb_table, adjs);
              }
              Gpu::copy(Gpu::hostToDevice,
                        h_herb_table.begin(), h_herb_table.end(),
                        d_fuel_table.begin());
          }
      }

      // Update wind from compact direction schedule (overrides constant wind)
      if (!wind_dir_sched.empty()) {
          MultiFab& wind_target = turb_wind_active ? *vel_base : vel;
          auto [ux_sched, uy_sched] = get_wind_at_time(wind_dir_sched,
                                                        static_cast<double>(time));
          inputs.ux = static_cast<amrex::Real>(ux_sched);
          inputs.uy = static_cast<amrex::Real>(uy_sched);
          init_velocity_constant(wind_target, geom, inputs.ux, inputs.uy, inputs.uz);
      }

      // Multi-station IDW wind interpolation (overrides single station / schedule)
      if (multi_wtr_active) {
          MultiFab& wind_target = turb_wind_active ? *vel_base : vel;
          apply_multi_wtr_to_vel(wind_target, geom, multi_wtr,
                                 static_cast<double>(time));
          // Domain-mean T/RH for global moisture model
          if (inputs.diurnal_moisture.enable == 1) {
              auto [T_mww, RH_mww] = multi_wtr.get_domain_TRH_at_time(
                                       static_cast<double>(time));
              inputs.diurnal_moisture.T_min  = static_cast<amrex::Real>(T_mww);
              inputs.diurnal_moisture.T_max  = static_cast<amrex::Real>(T_mww);
              inputs.diurnal_moisture.RH_min = static_cast<amrex::Real>(RH_mww);
              inputs.diurnal_moisture.RH_max = static_cast<amrex::Real>(RH_mww);
          }
      }

      // Apply precipitation wetting to dead fuel moisture (extends diurnal model)
      if (inputs.diurnal_moisture.enable == 1 && precip_state.initialized) {
          // Determine current rain rate (wtr_file overrides precip_schedule)
          float rain_rate = static_cast<float>(inputs.precip_rain_rate_mm_hr);
          if (wtr_active && !wtr_data.empty()) {
              rain_rate = static_cast<float>(wtr_data.get_precip_at_time(
                  static_cast<double>(time)));
          } else if (!precip_schedule.empty()) {
              // Linearly interpolate from schedule using binary search (O(log n))
              double t_d = static_cast<double>(time);
              if (t_d <= precip_schedule.front().first) {
                  rain_rate = static_cast<float>(precip_schedule.front().second);
              } else if (t_d >= precip_schedule.back().first) {
                  rain_rate = static_cast<float>(precip_schedule.back().second);
              } else {
                  auto it = std::upper_bound(precip_schedule.begin(), precip_schedule.end(),
                                             t_d,
                                             [](double t, const std::pair<double,double>& r){
                                                 return t < r.first; });
                  auto prev = std::prev(it);
                  const double alpha = (t_d - prev->first) / (it->first - prev->first);
                  rain_rate = static_cast<float>(prev->second +
                                                  alpha * (it->second - prev->second));
              }
          }
          // If wtr_file is active, override diurnal T/RH from .wtr data
          if (wtr_active && !wtr_data.empty()) {
              auto [T_wtr, RH_wtr] = wtr_data.get_TRH_at_time(static_cast<double>(time));
              // Update diurnal params in-place for this timestep (T_min=T_max and
              // RH_min=RH_max collapses the sinusoid to a constant – i.e. the .wtr
              // hourly value is used directly without a diurnal cycle on top).
              inputs.diurnal_moisture.T_min  = static_cast<amrex::Real>(T_wtr);
              inputs.diurnal_moisture.T_max  = static_cast<amrex::Real>(T_wtr);
              inputs.diurnal_moisture.RH_min = static_cast<amrex::Real>(RH_wtr);
              inputs.diurnal_moisture.RH_max = static_cast<amrex::Real>(RH_wtr);
          }
          // Build EMC from diurnal model as drying target
          RothermelMoistures emc_target = compute_diurnal_emc(
              inputs.diurnal_moisture,
              static_cast<double>(time),
              {static_cast<float>(inputs.rothermel.M_d1),
               static_cast<float>(inputs.rothermel.M_d10),
               static_cast<float>(inputs.rothermel.M_d100),
               static_cast<float>(inputs.rothermel.M_d1000),
               static_cast<float>(inputs.rothermel.M_lh),
               static_cast<float>(inputs.rothermel.M_lw)});
          apply_precipitation_moisture(precip_state, emc_target, rain_rate,
                                       static_cast<float>(dt_step),
                                       static_cast<float>(inputs.precip_threshold_mm_hr),
                                       static_cast<float>(inputs.M_sat));
          // Apply wetting result to Rothermel params (dead fuels only)
          inputs.rothermel.M_d1    = static_cast<amrex::Real>(precip_state.M_d1);
          inputs.rothermel.M_d10   = static_cast<amrex::Real>(precip_state.M_d10);
          inputs.rothermel.M_d100  = static_cast<amrex::Real>(precip_state.M_d100);
          inputs.rothermel.M_d1000 = static_cast<amrex::Real>(precip_state.M_d1000);
          inputs.rothermel.M_f     = static_cast<amrex::Real>(precip_state.M_d1);
          // Rebuild fuel table with updated moisture
          if (!inputs.rothermel.landscape_file.empty() && fuel_table_size > 0) {
              std::vector<RothermelComputed> h_table =
                  build_fuel_rothermel_table(inputs.rothermel,
                                             inputs.rothermel.landscape_fuel_type);
              if (!inputs.fuel_adj_file.empty()) {
                  auto adjs = parse_fuel_adjustment_file(inputs.fuel_adj_file);
                  apply_fuel_adjustment_to_table(h_table, adjs);
              }
              Gpu::copy(Gpu::hostToDevice,
                        h_table.begin(), h_table.end(),
                        d_fuel_table.begin());
          }
      }

      // ---- Live fuel moisture FMC seasonal link ----
      // When live_fuel_seasonal.enable = 1 and fmc_schedule is active, scale
      // M_lh and M_lw between their winter and summer values in proportion to
      // the current FMC fraction (0 = dormant FMC_min, 1 = peak FMC_max).
      if (inputs.live_fuel_seasonal.enable == 1 &&
          inputs.fmc_schedule.enable == 1 && !fmc_sched.empty()) {
          const amrex::Real fmc_now = inputs.crown.FMC;  // already updated above
          const amrex::Real fmc_min = static_cast<amrex::Real>(inputs.fmc_schedule.fmc_min);
          const amrex::Real fmc_max = static_cast<amrex::Real>(inputs.fmc_schedule.fmc_max);
          const amrex::Real frac = (fmc_max > fmc_min)
              ? std::max(amrex::Real(0.0), std::min(amrex::Real(1.0),
                          (fmc_now - fmc_min) / (fmc_max - fmc_min)))
              : amrex::Real(0.5);
          inputs.rothermel.M_lh = static_cast<amrex::Real>(
              inputs.live_fuel_seasonal.M_lh_winter +
              frac * (inputs.live_fuel_seasonal.M_lh_summer -
                      inputs.live_fuel_seasonal.M_lh_winter));
          inputs.rothermel.M_lw = static_cast<amrex::Real>(
              inputs.live_fuel_seasonal.M_lw_winter +
              frac * (inputs.live_fuel_seasonal.M_lw_summer -
                      inputs.live_fuel_seasonal.M_lw_winter));
      }

      // ---- Elevation lapse-rate T/RH correction for per-cell solar EMC ----
      // Apply lapse-rate T/RH adjustment to spatial_moisture_mf before
      // the solar EMC pass.  This pass is inline: for each cell read elevation
      // and adjust T/RH before feeding into apply_solar_emc_to_spatial_moisture.
      // When use_elevation_lapse = 0 (default) this block is skipped.
      if (inputs.use_elevation_lapse == 1 && inputs.diurnal_moisture.enable == 1) {
          // Magnus formula empirical constants (Murray 1967, J. Applied Meteorology 6:203):
          //   a = 17.67  (dimensionless)
          //   b = 243.5  [°C]
          // Used in Clausius–Clapeyron RH approximation:
          //   RH_cell ≈ RH_ref * exp(a * (T_cell - T_ref) / (T_ref + b))
          constexpr amrex::Real MAGNUS_A = amrex::Real(17.67);
          constexpr amrex::Real MAGNUS_B = amrex::Real(243.5);  // °C

          // Compute domain-mean T and RH at this timestep (same formulas as diurnal)
          const double phase_lr = WildfireConst::TWO_PI *
              (inputs.diurnal_moisture.t_start_s +
               static_cast<double>(time) -
               inputs.diurnal_moisture.t_T_peak_s)
              / WildfireConst::DAY_S;
          const double T_mean_lr  = 0.5 * (inputs.diurnal_moisture.T_max + inputs.diurnal_moisture.T_min);
          const double A_T_lr     = 0.5 * (inputs.diurnal_moisture.T_max - inputs.diurnal_moisture.T_min);
          const double RH_mean_lr = 0.5 * (inputs.diurnal_moisture.RH_max + inputs.diurnal_moisture.RH_min);
          const double A_RH_lr    = 0.5 * (inputs.diurnal_moisture.RH_max - inputs.diurnal_moisture.RH_min);
          double T_ref_lr  = T_mean_lr  + A_T_lr  * std::sin(phase_lr);
          double RH_ref_lr = RH_mean_lr - A_RH_lr * std::sin(phase_lr);
          T_ref_lr  = std::max(-40.0, std::min(60.0,  T_ref_lr));
          RH_ref_lr = std::max(1.0,   std::min(100.0, RH_ref_lr));

          const amrex::Real lapse   = static_cast<amrex::Real>(inputs.lapse_rate_C_per_m);
          const amrex::Real elev0   = static_cast<amrex::Real>(inputs.lapse_ref_elevation_m);
          const amrex::Real T_ref_r = static_cast<amrex::Real>(T_ref_lr);
          const amrex::Real RH_ref_r= static_cast<amrex::Real>(RH_ref_lr);
          const amrex::Real sol_C   = static_cast<amrex::Real>(inputs.solar_radiation.solar_heating_C);

          // Per-cell: apply lapse-rate correction to T, then Clausius-Clapeyron to RH,
          // then Nelson EMC, then update spatial_moisture_mf components 0-2.
          for (MFIter mfi(spatial_moisture_mf); mfi.isValid(); ++mfi) {
              const Box& bx = mfi.validbox();
              auto sm = spatial_moisture_mf.array(mfi);
              auto const elev = elevation_mf.const_array(mfi);
              auto const shade = shade_fraction_mf.const_array(mfi);
              ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                  const amrex::Real dT_lapse = lapse * (elev(i,j,k) - elev0);
                  amrex::Real T_cell = T_ref_r - dT_lapse;
                  // Clausius–Clapeyron RH correction using the Magnus formula approximation:
                  //   RH_cell ≈ RH_ref * exp(a * dT / (T_ref + b))
                  // where dT = T_cell - T_ref, a = MAGNUS_A (17.67, dimensionless),
                  // b = MAGNUS_B (243.5 °C).  Reference: Murray (1967).
                  // At higher elevations (dT < 0, T_cell < T_ref) RH increases.
                  const amrex::Real dT_cc  = -(dT_lapse);  // T_cell - T_ref
                  // Magnus formula constants (Murray 1967): a=17.67, b=243.5 °C
                  const amrex::Real RH_adj = RH_ref_r * std::exp(
                      MAGNUS_A * dT_cc / (T_ref_r + MAGNUS_B));
                  amrex::Real RH_cell = amrex::max(amrex::Real(1.0),
                                                   amrex::min(amrex::Real(100.0), RH_adj));
                  // Apply solar heating to unshaded T
                  amrex::Real T_fuel = T_cell + sol_C * (amrex::Real(1.0) - shade(i,j,k));
                  T_fuel = amrex::max(amrex::Real(-40.0), amrex::min(amrex::Real(60.0), T_fuel));
                  // Nelson / Simard EMC
                  amrex::Real rh_f = RH_cell / amrex::Real(100.0);
                  amrex::Real emc_pct;
                  if (RH_cell < amrex::Real(10.0)) {
                      emc_pct = amrex::Real(0.03229) + amrex::Real(0.281073)*rh_f
                                - amrex::Real(0.000578)*T_fuel*rh_f;
                  } else if (RH_cell < amrex::Real(50.0)) {
                      emc_pct = amrex::Real(2.22749) + amrex::Real(0.160107)*rh_f
                                - amrex::Real(0.014784)*T_fuel;
                  } else {
                      emc_pct = amrex::Real(21.0606) + amrex::Real(0.005565)*rh_f*rh_f
                                - amrex::Real(0.00035)*T_fuel*rh_f - amrex::Real(0.483199)*rh_f;
                  }
                  emc_pct = amrex::max(amrex::Real(0.005), emc_pct);
                  sm(i,j,k,0) = emc_pct / amrex::Real(100.0);            // M_d1
                  sm(i,j,k,1) = emc_pct * amrex::Real(1.10) / amrex::Real(100.0);  // M_d10
                  sm(i,j,k,2) = emc_pct * amrex::Real(1.30) / amrex::Real(100.0);  // M_d100
              });
          }
      }

      // --- Step 2: Compute surface ROS via selected fire spread model
      // Apply wind-terrain velocity modification (Options 3-7) before ROS computation
      if (wind_terrain_modifies_vel) {
          apply_wind_terrain_velocity(vel_effective, vel, terrain_slopes.get(), inputs);
      } else {
          MultiFab::Copy(vel_effective, vel, 0, 0, 3, 0);
      }

      // Apply heat flux wind corrections (upward velocity + induced inflow)
      if (heat_flux_active) {
          apply_heatflux_wind(vel_effective, vel, heat_flux_mf, &phi,
                              inputs.heat_flux);
      }

      const MultiFab& vel_for_model = (wind_terrain_modifies_vel || heat_flux_active)
                                      ? vel_effective : vel;

      if (inputs.fire_spread_model == "balbi") {
          compute_balbi_R(R_mf, vel_for_model, geom, inputs.rothermel, inputs.balbi,
                           terrain_slopes.get(),
                           !inputs.rothermel.landscape_file.empty() ? &fuel_model_mf : nullptr,
                           d_balbi_table_ptr, balbi_table_size,
                           heat_flux_active ? &heat_flux_mf : nullptr,
                           heat_flux_active ? &inputs.heat_flux : nullptr);
      } else if (inputs.fire_spread_model == "cheney_gould") {
          compute_cheney_gould_R(R_mf, vel_for_model, inputs.cheney_gould);
      } else if (inputs.fire_spread_model == "cruz_crown") {
          compute_cruz_crown_R(R_mf, vel_for_model, inputs.cruz_crown);
      } else if (inputs.fire_spread_model == "fbp_o1a" ||
                 inputs.fire_spread_model == "fbp_o1b" ||
                 inputs.fire_spread_model == "fbp_s1"  ||
                 inputs.fire_spread_model == "fbp_s2"  ||
                 inputs.fire_spread_model == "fbp_s3") {
          compute_fbp_R(R_mf, vel_for_model, inputs.fbp);
      } else if (inputs.fire_spread_model == "lautenberger") {
          compute_lautenberger_R(R_mf, vel_for_model, inputs.rothermel, inputs.lautenberger,
                                  terrain_slopes.get());
      } else {
          compute_rothermel_R(R_mf, vel_for_model, geom, inputs.rothermel,
                               terrain_slopes.get(),
                               !inputs.rothermel.landscape_file.empty() ? &fuel_model_mf : nullptr,
                               d_fuel_table_ptr, fuel_table_size,
                               has_spatial_crown ? &cc_mf : nullptr,
                               has_spatial_crown ? &canopy_height_mf : nullptr);
      }
      compute_fire_behavior(fireline_intensity_mf, flame_length_mf, R_mf, inputs.rothermel);
      if (inputs.weise_biging.enable == 1) {
          compute_weise_biging_whirl(weise_data, fireline_intensity_mf, flame_length_mf,
                                     vel, terrain_slopes.get(), inputs.weise_biging);
      }
      if (inputs.viegas.enable == 1) {
          // For Balbi+Viegas: pass Balbi table so diagnostic uses Balbi amplitude
          compute_viegas_diagnostics(viegas_data, R_mf, vel, inputs.rothermel, inputs.viegas,
                                     terrain_slopes.get(),
                                     !inputs.rothermel.landscape_file.empty() ? &fuel_model_mf : nullptr,
                                     d_fuel_table_ptr, fuel_table_size,
                                     use_balbi_for_viegas,
                                     use_balbi_for_viegas ? &bc_global_default : nullptr,
                                     use_balbi_for_viegas ? d_balbi_table_ptr : nullptr,
                                     use_balbi_for_viegas ? balbi_table_size : 0);
          // Option 2: override R_mf with max(primary, R_viegas) in eruptive cells
          if (inputs.wind_terrain.model == "viegas_ros") {
              apply_viegas_ros_override(R_mf, viegas_data);
          }
      }
      // Fire ecology diagnostics (always computed)
      compute_fire_ecology(ecology_mf, fireline_intensity_mf, R_mf,
                           inputs.rothermel, inputs.fire_ecology, inputs.crown);

      // Scott & Reinhardt (2001) full bisection-based TI/CI (optional, host-only)
      if (inputs.scott_reinhardt_full.enable == 1) {
          compute_full_scott_reinhardt(
              ti_full_mf, ci_full_mf,
              inputs.rothermel, inputs.crown,
              !inputs.rothermel.landscape_file.empty() ? &fuel_model_mf : nullptr,
              nullptr, 0,
              amrex::Real(inputs.scott_reinhardt_full.U_max_kmh));
      }

      // Ecology → propagation coupling:
      //   (a) P_ignition scales R_mf in unburned cells when couple_to_ros = 1
      apply_ecology_p_ignition_feedback(R_mf, ecology_mf, phi, inputs.fire_ecology);

      // Cap 8: Active crown fire ROS feedback into the level-set ROS field.
      // When crown.enable = 1 and a cell is classified as active crown fire
      // (crown_activity == 2), override R_mf with the active crown fire ROS so
      // the level-set front propagates at the faster crown speed.
      // Route: Cruz, Alexander & Wakimoto (2005) when use_cruz_crown = 1;
      //        Van Wagner (1977) 3/CBD proxy otherwise.
      if (inputs.crown.enable == 1 && use_levelset) {
          const Real CBD_global   = Real(inputs.crown.CBD);
          const Real FMC_global   = Real(inputs.crown.FMC);
          const Real m_factor_g   = amrex::max(Real(0.3),
                                        amrex::min(Real(1.0),
                                        Real(1.0) - (FMC_global - Real(100.0)) / Real(200.0)));
          const bool use_cruz_tl          = (inputs.crown.use_cruz_crown == 1);
          const bool use_roth1991_tl      = (inputs.crown.use_rothermel1991_crown == 1);
          const bool use_passive_blend_tl = (inputs.crown.use_passive_blend == 1);
          const Real CBH_tl = Real(inputs.crown.CBH);
          const Real I_o_tl = Real(0.010) * CBH_tl * (Real(460.0) + Real(25.9) * FMC_global);
          CruzCrownComputed ccc_tl;
          if (use_cruz_tl) { ccc_tl = ccc_global; }
          const Real MC10_tl = Real(inputs.cruz_crown.MC10);

          for (MFIter mfi(R_mf); mfi.isValid(); ++mfi) {
              const Box& bx   = mfi.validbox();
              auto       R    = R_mf.array(mfi);
              auto const eco  = ecology_mf.const_array(mfi);
              auto const v    = vel.const_array(mfi);
              auto const fi   = fireline_intensity_mf.const_array(mfi);
              const bool use_sp = has_spatial_crown;
              Array4<const Real> cbd_arr;
              if (use_sp) {
                  cbd_arr = cbd_mf.const_array(mfi);
              }
              const Real CBD_g   = CBD_global;
              const Real mf_val  = m_factor_g;
              const bool use_passive_tl = use_passive_blend_tl;
              const Real I_o_tl_k = I_o_tl;

              ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                  if (eco(i, j, k, 3) < Real(1.5)) return; // surface or passive
                  Real R_surface = R(i, j, k);
                  Real R_crown_ms;
                  if (use_roth1991_tl) {
                      R_crown_ms = compute_rothermel_1991_crown_ros(R_surface);
                  } else if (use_cruz_tl) {
                      Real ux = v(i,j,k,0);
                      Real uy = v(i,j,k,1);
                      Real wind_mag = std::sqrt(ux*ux + uy*uy);
                      Real CBD_c = use_sp ? cbd_arr(i,j,k) : CBD_g;
                      R_crown_ms = compute_crown_fire_spread_rate_cruz(wind_mag, CBD_c, MC10_tl);
                  } else {
                      Real CBD_c = use_sp ? cbd_arr(i,j,k) : CBD_g;
                      CBD_c = amrex::max(CBD_c, Real(0.01));
                      R_crown_ms = (Real(3.0) / CBD_c) * mf_val / Real(60.0);
                  }
                  if (use_passive_tl) {
                      const Real I_B_kwm = fi(i, j, k);
                      R(i, j, k) = compute_van_wagner_passive_blend(
                          R_surface, R_crown_ms, I_B_kwm, I_o_tl_k);
                  } else {
                      R(i, j, k) = amrex::max(R_surface, R_crown_ms);
                  }
              });
          }
          amrex::Print() << "  Crown ROS override applied (crown.enable=1"
                         << (use_roth1991_tl ? ", Rothermel1991"
                             : use_cruz_tl   ? ", Cruz 2005"
                                             : ", Van Wagner 1977") << ").\n";
      }

      // Fire emissions (CO₂, CO, PM₂.₅)
      compute_fire_emissions(emissions_mf, phi, fuel_consumption_mf,
                             inputs.rothermel, inputs.emissions);

      // Feature 10: Fire acceleration model (Rothermel 1983 / Catchpole et al. 1992)
      // Scales R_mf by alpha = 1 - exp(-r_fire / L_acc) to account for the
      // slower spread of small fires before quasi-steady-state is reached.
      // When disabled (default) or fire is large, this is a no-op.
      apply_fire_acceleration(R_mf, phi, geom, inputs.acceleration);

      // ---- Aerial retardant suppression: scale R_mf in active drop zones ----
      if (!retardant_drops.empty()) {
          apply_retardant_to_ros(R_mf, retardant_drops, geom, time);
          // Suppress spotting probability inside active retardant zones
          apply_retardant_to_spotting_probability(spotting_data, retardant_drops, geom, time);
      }

      // ---- Burn-period gate: zero R_mf outside the active daily spread window ----
      // When burn_period.enable = 1, fire spread is suppressed to zero outside
      // the [start_hour, end_hour) local clock window.  Moisture evolution and all
      // other diagnostics continue normally.  The FARSITE propagation and level-set
      // advection each read R_mf, so zeroing it here silences all spread paths.
      bool burn_period_active = true;
      if (inputs.burn_period.enable == 1) {
          // Compute current local clock hour (wraps at 24)
          const double clock_hour = std::fmod(
              inputs.burn_period.sim_start_hour +
              static_cast<double>(time) / 3600.0,
              24.0);
          const double sh = static_cast<double>(inputs.burn_period.start_hour);
          const double eh = static_cast<double>(inputs.burn_period.end_hour);

          // Normal window (no midnight crossing): active when sh <= clock < eh
          // Overnight window (midnight crossing):  active when clock >= sh OR clock < eh
          if (sh < eh) {
              burn_period_active = (clock_hour >= sh && clock_hour < eh);
          } else {
              // Window crosses midnight (e.g. 22:00–06:00)
              burn_period_active = (clock_hour >= sh || clock_hour < eh);
          }

          if (!burn_period_active) {
              R_mf.setVal(amrex::Real(0.0));
              amrex::Print() << "  Burn period inactive at hour " << clock_hour
                             << " (window " << sh << ":00-" << eh
                             << ":00) – spread suppressed.\n";
          }
      }

      if (use_levelset) {
	// Traditional level set advection.
	// Pass pre-computed R_mf when a wind-terrain model or non-Rothermel spread
	// model is active (see use_precomp_R_for_advection defined above).
	advect_levelset_weno5z_rk3(phi, vel, geom, dt_step, inputs.rothermel,
                                   terrain_slopes.get(),
                                   use_precomp_R_for_advection ? &R_mf : nullptr);
	dt = compute_dt(R_mf, geom, inputs.cfl);
      } else if (use_mtt) {
	// --- MTT: phi is updated analytically from pre-computed arrival times.
	// No re-advection needed; just set phi = arrival_time - current_time
	// where current_time = time + dt_step (already advanced below).
	apply_mtt_phi_update(phi, arrival_time_mf, time + dt_step);
	fill_boundary_extrap(phi, geom);
	// dt stays constant for MTT (fixed by initial ROS computation)
      } else {
	// --- Step 3: FARSITE elliptical wavelet propagation (Richards 1990)
	// --- Step 4: Merge to new perimeter
	// For wind-terrain models, pass vel_for_model so FARSITE ellipse
	// orientation and ROS also reflect the terrain-corrected wind.
	compute_farsite_spread(phi, vel_for_model, farsite_spread, geom, dt_step, inputs.rothermel, inputs.farsite, inputs.crown, terrain_slopes.get(), &fuel_consumption_mf, &crown_fire_fraction_mf, has_spatial_crown ? &cc_mf : nullptr, has_spatial_crown ? &canopy_height_mf : nullptr, ccc_ptr);
	
	// --- Step 5: Apply crown/spotting sub-models
	if (inputs.spotting.enable == 1 && (step % inputs.spotting.check_interval == 0)) {
	  compute_spotting_probability(spotting_data, phi, vel, geom, inputs.rothermel, inputs.spotting, terrain_slopes.get());
	  generate_firebrand_spots(phi, spotting_data, vel, geom, inputs.spotting, step,
	                           !inputs.rothermel.landscape_file.empty() ? &fuel_model_mf : nullptr,
	                           !inputs.rothermel.landscape_file.empty() ? &inputs.rothermel.landscape_fuel_type : nullptr);
	}
	// Albini (1983) firebrand spotting with 2-D trajectory integration
	if (inputs.albini_spotting.enable == 1 && (step % inputs.albini_spotting.check_interval == 0)) {
	  compute_albini_spotting(phi, albini_data, vel, R_mf, geom,
	                          inputs.rothermel, inputs.albini_spotting, step,
	                          !inputs.rothermel.landscape_file.empty() ? &fuel_model_mf : nullptr,
	                          !inputs.rothermel.landscape_file.empty() ? &inputs.rothermel.landscape_fuel_type : nullptr);
	  // Scott/Albini (1979) maximum spotting distance table diagnostic:
	  // Print the table maximum for the dominant global fuel model and
	  // the current mean wind speed to help the user assess whether
	  // trajectory-integrated distances are physically reasonable.
	  if (ParallelDescriptor::IOProcessor()) {
	    Real wind_speed_ms = std::sqrt(inputs.ux * inputs.ux + inputs.uy * inputs.uy);
	    float d_max_base   = get_max_spot_dist_m(inputs.fuel_adj_model > 0
	                             ? inputs.fuel_adj_model : 4,  // fallback to FM4 (chaparral)
	                             inputs.rothermel.landscape_fuel_type);
	    float d_max_scaled = get_max_spot_dist_scaled_m(
	                             inputs.fuel_adj_model > 0 ? inputs.fuel_adj_model : 4,
	                             static_cast<float>(wind_speed_ms),
	                             inputs.rothermel.landscape_fuel_type);
	    amrex::Print() << "  Scott/Albini max spotting: base=" << d_max_base
	                   << " m  scaled(U=" << wind_speed_ms << " m/s)="
	                   << d_max_scaled << " m\n";
	  }
	}
	
	// --- Step 6: Simulate post-frontal burnout
	// (Bulk fuel consumption is computed within compute_farsite_spread)

	// Feature 8: Albini (1979) torching-tree spotting from crown-fire cells
	if (inputs.torching_spotting.enable == 1 &&
	    (step % inputs.torching_spotting.check_interval == 0)) {
	  compute_albini_torching_spots(phi, ecology_mf, fireline_intensity_mf,
	                                flame_length_mf, canopy_height_mf,
	                                vel, geom, inputs.torching_spotting, step);
	}
      }

      // --- Apply barrier polygons (firebreaks): extinguish burning cells
      //     that coincide with barrier locations, regardless of propagation method.
      if (!barrier_cells.empty()) {
	apply_barrier_polygons(phi, barrier_cells, geom);
      }
      if (inputs.reinit_int > 0 && (step % inputs.reinit_int == 0) && use_levelset) {
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

      // --- Update arrival time: mark cells that first became burned this step
      {
        const Real cur_time = time;
        for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
          const Box& bx = mfi.validbox();
          auto const p  = phi.const_array(mfi);
          auto       at = arrival_time_mf.array(mfi);
          ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            if (p(i, j, k) < Real(0.0) && at(i, j, k) < Real(0.0)) {
              at(i, j, k) = cur_time;
            }
          });
        }
      }

      // Feature 9: Update per-cell residual fuel load.
      // For each cell that has burned (arrival_time >= 0), compute the fraction of
      // fuel remaining: f_r = exp(-(time - arrival_time) / tau_burnout).
      // Optionally scale R_mf in re-entry cells by f_r (fuel_depletion.couple_to_ros).
      if (inputs.fuel_depletion.enable == 1) {
        const Real tau_b    = Real(inputs.fuel_depletion.tau_burnout);
        const Real cur_time = time;
        for (MFIter mfi(residual_fuel_mf); mfi.isValid(); ++mfi) {
          const Box& bx = mfi.validbox();
          auto const at = arrival_time_mf.const_array(mfi);
          auto       rf = residual_fuel_mf.array(mfi);
          ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const Real t_arrive = at(i, j, k);
            if (t_arrive >= Real(0.0)) {
              const Real elapsed = cur_time - t_arrive;
              rf(i, j, k) = std::exp(-elapsed / tau_b);
            }
          });
        }
        // Optional: scale ROS by residual fuel in re-entry cells (fuel previously burned)
        if (inputs.fuel_depletion.couple_to_ros == 1) {
          for (MFIter mfi(R_mf); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            auto const phi_arr = phi.const_array(mfi);
            auto const at      = arrival_time_mf.const_array(mfi);
            auto const rf      = residual_fuel_mf.const_array(mfi);
            auto       R       = R_mf.array(mfi);
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
              // Scale only unburned cells where fuel was previously depleted (re-entry)
              if (at(i, j, k) >= Real(0.0) && phi_arr(i, j, k) >= Real(0.0)) {
                R(i, j, k) *= rf(i, j, k);
              }
            });
          }
        }
      }

      // --- Dynamic fire points: check if a new ignition file has appeared
      if (!inputs.dynamic_fire_points_file.empty()) {
        apply_dynamic_fire_points(phi, geom,
                                  inputs.dynamic_fire_points_file,
                                  inputs.fire_gaussian_sigma);
      }

      // ---- Multiple scheduled ignitions: fire any events due this timestep ----
      if (!ignition_sched.empty()) {
          apply_scheduled_ignitions(phi, geom, ignition_sched,
                                    time, time - dt_step, use_indicator);
          fill_boundary_extrap(phi, geom);
      }

      // --- Write checkpoint if requested
      if (inputs.chk_int > 0 && ((step - restart_step) % inputs.chk_int == 0)) {
        char chk_buf[64];
        std::snprintf(chk_buf, sizeof(chk_buf), "chk%04d", step);
        write_checkpoint(chk_buf, phi, geom, step, time);
      }

      // ---- Timed isochrone output ----
      // Fire at each clock time that is a new multiple of isochrone_interval_s.
      if (isochrone_active) {
        int cur_iso_idx = static_cast<int>(std::floor(
            static_cast<double>(time) /
            static_cast<double>(inputs.isochrone_interval_s)));
        if (cur_iso_idx > last_isochrone_idx) {
          // Write all intervals that elapsed since the last check
          // (handles large dt that skips multiple intervals, though unusual)
          for (int iso_i = last_isochrone_idx + 1; iso_i <= cur_iso_idx; ++iso_i) {
            char iso_csv_buf[128];
            std::snprintf(iso_csv_buf, sizeof(iso_csv_buf),
                          "isochrone_%06d.csv", iso_i);
            write_fire_perimeter_csv(phi, geom, iso_csv_buf);
            amrex::Print() << "Isochrone " << iso_i
                           << " written at t=" << time << " s -> "
                           << iso_csv_buf << "\n";
            if (inputs.write_perimeter_geojson == 1) {
              char iso_gj_buf[128];
              std::snprintf(iso_gj_buf, sizeof(iso_gj_buf),
                            "isochrone_%06d.geojson", iso_i);
              write_fire_perimeter_geojson(phi, geom, iso_gj_buf, iso_i, time);
            }
          }
          last_isochrone_idx = cur_iso_idx;
        }
      }

      // --- Step 7: Update states, record outputs, step time
      if (inputs.plot_int > 0 && (step % inputs.plot_int == 0)) {
	// Compute per-plotfile diagnostics
	{
	  const RothermelComputed rc_plt = compute_rothermel_params(inputs.rothermel);
	  compute_heat_per_unit_area(heat_per_unit_area_mf, phi, rc_plt.I_R,
	                             inputs.farsite.tau_residence);
	}
	compute_burnout_time(burnout_time_mf, arrival_time_mf, inputs.farsite.tau_residence,
                       !inputs.rothermel.landscape_file.empty() ? &fuel_model_mf : nullptr,
                       d_tau_table_ptr, tau_table_size);
	compute_vorticity(vorticity_mf, vel, geom);
	compute_reaction_intensity(reaction_intensity_mf,
	                           !inputs.rothermel.landscape_file.empty() ? &fuel_model_mf : nullptr,
	                           inputs.rothermel, d_fuel_table_ptr, fuel_table_size);
	char buf[64];
	std::snprintf(buf, sizeof(buf), "plt%04d", step);
	Vector<std::string> names = {"phi", "velx", "vely"
#if (AMREX_SPACEDIM == 3)
				     , "velz", "farsite_dx", "farsite_dy", "farsite_dz", "R",
				     "spot_prob", "spot_count", "spot_dist", "spot_active", "fuel_consumption", "crown_fraction",
				     "albini_Hz", "albini_count", "albini_dist", "albini_active",
				     "elevation", "slope", "aspect", "fuel_model",
				     "fireline_intensity", "flame_length",
				     "weise_flame_height", "weise_flame_tilt",
				     "weise_whirl_height", "weise_whirl_radius",
				     "weise_angular_velocity", "weise_max_tang_vel",
				     "viegas_ROS", "viegas_eruptive_flag",
				     "viegas_ROS_excess", "viegas_flame_tilt", "viegas_slope_factor",
				     "scorch_height", "prob_ignition", "tree_mortality", "crown_activity",
					   "torching_ratio", "crowning_ratio", "energy_release_component",
				     "co2_emissions", "co_emissions", "pm25_emissions",
				   "arrival_time", "heat_per_unit_area", "vorticity_z",
				   "cbh", "cbd", "canopy_cover", "canopy_height", "burnout_time",
   "reaction_intensity", "residual_fuel", "shade_fraction",
   "torching_index_kmh", "crowning_index_kmh",
   "moisture_d1", "moisture_d10", "moisture_d100", "moisture_lh", "moisture_lw"
#else
				     , "farsite_dx", "farsite_dy", "R",
				     "spot_prob", "spot_count", "spot_dist", "spot_active", "fuel_consumption", "crown_fraction",
				     "albini_Hz", "albini_count", "albini_dist", "albini_active",
				     "elevation", "slope", "aspect", "fuel_model",
				     "fireline_intensity", "flame_length",
				     "weise_flame_height", "weise_flame_tilt",
				     "weise_whirl_height", "weise_whirl_radius",
				     "weise_angular_velocity", "weise_max_tang_vel",
				     "viegas_ROS", "viegas_eruptive_flag",
				     "viegas_ROS_excess", "viegas_flame_tilt", "viegas_slope_factor",
				     "scorch_height", "prob_ignition", "tree_mortality", "crown_activity",
					   "torching_ratio", "crowning_ratio", "energy_release_component",
				     "co2_emissions", "co_emissions", "pm25_emissions",
				   "arrival_time", "heat_per_unit_area", "vorticity_z",
				   "cbh", "cbd", "canopy_cover", "canopy_height", "burnout_time",
   "reaction_intensity", "residual_fuel", "shade_fraction",
   "torching_index_kmh", "crowning_index_kmh",
   "moisture_d1", "moisture_d10", "moisture_d100", "moisture_lh", "moisture_lw"
#endif
	};
	MultiFab plotmf(ba, dm, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 1 + 5 + WEISE_NCOMP + VIEGAS_NCOMP + ECOLOGY_NCOMP + EMISSIONS_NCOMP + 13 + FMS_NCOMP, 0);
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
	MultiFab::Copy(plotmf, weise_data, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 1 + 5, WEISE_NCOMP, 0);
	MultiFab::Copy(plotmf, viegas_data, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 1 + 5 + WEISE_NCOMP, VIEGAS_NCOMP, 0);
	MultiFab::Copy(plotmf, ecology_mf,   0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 1 + 5 + WEISE_NCOMP + VIEGAS_NCOMP, ECOLOGY_NCOMP, 0);
	MultiFab::Copy(plotmf, emissions_mf, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 1 + 5 + WEISE_NCOMP + VIEGAS_NCOMP + ECOLOGY_NCOMP, EMISSIONS_NCOMP, 0);
	MultiFab::Copy(plotmf, arrival_time_mf,       0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 1 + 5 + WEISE_NCOMP + VIEGAS_NCOMP + ECOLOGY_NCOMP + EMISSIONS_NCOMP, 1, 0);
	MultiFab::Copy(plotmf, heat_per_unit_area_mf, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 1 + 5 + WEISE_NCOMP + VIEGAS_NCOMP + ECOLOGY_NCOMP + EMISSIONS_NCOMP + 1, 1, 0);
	MultiFab::Copy(plotmf, vorticity_mf,          0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 1 + 5 + WEISE_NCOMP + VIEGAS_NCOMP + ECOLOGY_NCOMP + EMISSIONS_NCOMP + 2, 1, 0);
	MultiFab::Copy(plotmf, cbh_mf, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 1 + 5 + WEISE_NCOMP + VIEGAS_NCOMP + ECOLOGY_NCOMP + EMISSIONS_NCOMP + 3, 1, 0);
	MultiFab::Copy(plotmf, cbd_mf, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 1 + 5 + WEISE_NCOMP + VIEGAS_NCOMP + ECOLOGY_NCOMP + EMISSIONS_NCOMP + 4, 1, 0);
	MultiFab::Copy(plotmf, cc_mf,            0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 1 + 5 + WEISE_NCOMP + VIEGAS_NCOMP + ECOLOGY_NCOMP + EMISSIONS_NCOMP + 5, 1, 0);
	MultiFab::Copy(plotmf, canopy_height_mf, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 1 + 5 + WEISE_NCOMP + VIEGAS_NCOMP + ECOLOGY_NCOMP + EMISSIONS_NCOMP + 6, 1, 0);
	MultiFab::Copy(plotmf, burnout_time_mf,  0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 1 + 5 + WEISE_NCOMP + VIEGAS_NCOMP + ECOLOGY_NCOMP + EMISSIONS_NCOMP + 7, 1, 0);
	MultiFab::Copy(plotmf, reaction_intensity_mf, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 1 + 5 + WEISE_NCOMP + VIEGAS_NCOMP + ECOLOGY_NCOMP + EMISSIONS_NCOMP + 8, 1, 0);
	MultiFab::Copy(plotmf, residual_fuel_mf, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 1 + 5 + WEISE_NCOMP + VIEGAS_NCOMP + ECOLOGY_NCOMP + EMISSIONS_NCOMP + 9, 1, 0);
	MultiFab::Copy(plotmf, shade_fraction_mf, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 1 + 5 + WEISE_NCOMP + VIEGAS_NCOMP + ECOLOGY_NCOMP + EMISSIONS_NCOMP + 10, 1, 0);
	MultiFab::Copy(plotmf, ti_full_mf, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 1 + 5 + WEISE_NCOMP + VIEGAS_NCOMP + ECOLOGY_NCOMP + EMISSIONS_NCOMP + 11, 1, 0);
	MultiFab::Copy(plotmf, ci_full_mf, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 1 + 5 + WEISE_NCOMP + VIEGAS_NCOMP + ECOLOGY_NCOMP + EMISSIONS_NCOMP + 12, 1, 0);
	MultiFab::Copy(plotmf, spatial_moisture_mf, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 1 + 5 + WEISE_NCOMP + VIEGAS_NCOMP + ECOLOGY_NCOMP + EMISSIONS_NCOMP + 13, FMS_NCOMP, 0);
	WriteSingleLevelPlotfile(buf, plotmf, names, geom, time, step);
	amrex::Print() << "Wrote " << buf << "\n";
	{
	  amrex::Long n_burned = 0, n_perim = 0;
	  const auto dxx = geom.CellSize();
	  for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
	    const Box& bx = mfi.validbox();
	    auto const p = phi.const_array(mfi);
	    const IntVect bxlo = bx.smallEnd();
	    const IntVect bxhi = bx.bigEnd();
	    amrex::LoopOnCpu(bx, [&](int i, int j, int k) {
	      if (p(i,j,k,0) < Real(0.0)) {
	        ++n_burned;
	        bool on_perim = false;
	        if (i > bxlo[0]) on_perim |= (p(i-1,j,k,0) >= Real(0.0));
	        if (i < bxhi[0]) on_perim |= (p(i+1,j,k,0) >= Real(0.0));
	        if (j > bxlo[1]) on_perim |= (p(i,j-1,k,0) >= Real(0.0));
	        if (j < bxhi[1]) on_perim |= (p(i,j+1,k,0) >= Real(0.0));
	        if (on_perim) ++n_perim;
	      }
	    });
	  }
	  ParallelDescriptor::ReduceLongSum(n_burned);
	  ParallelDescriptor::ReduceLongSum(n_perim);
	  Real cell_area = dxx[0] * dxx[1];
	  amrex::Print() << "  Burned area: " << Real(n_burned)*cell_area/1.0e4 << " ha"
	                 << "  Perimeter: " << Real(n_perim)*std::sqrt(cell_area)/1000.0 << " km\n";
	}
	// Write negative phi x-y data files
	char xy_buf[64];
	std::snprintf(xy_buf, sizeof(xy_buf), "phi_negative_%04d.dat", step);
	write_negative_phi_xy(phi, geom, xy_buf);
	
	std::snprintf(xy_buf, sizeof(xy_buf), "phi_envelope_%04d.dat", step);
	write_negative_phi_convex_hull(phi, geom, xy_buf);
	if (!inputs.fire_stats_file.empty()) append_fire_stats(phi, geom, &emissions_mf, step, time, inputs.fire_stats_file, &R_mf);
	if (inputs.write_perimeter_csv == 1) { char csv_buf[64]; std::snprintf(csv_buf, sizeof(csv_buf), "perimeter_%04d.csv", step); write_fire_perimeter_csv(phi, geom, csv_buf); }
	if (inputs.write_perimeter_geojson == 1) { char gjson_buf[64]; std::snprintf(gjson_buf, sizeof(gjson_buf), "perimeter_%04d.geojson", step); write_fire_perimeter_geojson(phi, geom, gjson_buf, step, time); }
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
	// Compute per-plotfile diagnostics
	{
	  const RothermelComputed rc_plt = compute_rothermel_params(inputs.rothermel);
	  compute_heat_per_unit_area(heat_per_unit_area_mf, phi, rc_plt.I_R,
	                             inputs.farsite.tau_residence);
	}
	compute_burnout_time(burnout_time_mf, arrival_time_mf, inputs.farsite.tau_residence,
                       !inputs.rothermel.landscape_file.empty() ? &fuel_model_mf : nullptr,
                       d_tau_table_ptr, tau_table_size);
	compute_vorticity(vorticity_mf, vel, geom);
	compute_reaction_intensity(reaction_intensity_mf,
	                           !inputs.rothermel.landscape_file.empty() ? &fuel_model_mf : nullptr,
	                           inputs.rothermel, d_fuel_table_ptr, fuel_table_size);
	char buf[64];
	std::snprintf(buf, sizeof(buf), "plt%04d", final_step);
	Vector<std::string> names = {"phi", "velx", "vely"
#if (AMREX_SPACEDIM == 3)
				     , "velz", "farsite_dx", "farsite_dy", "farsite_dz", "R",
				     "spot_prob", "spot_count", "spot_dist", "spot_active", "fuel_consumption", "crown_fraction",
				     "albini_Hz", "albini_count", "albini_dist", "albini_active",
				     "elevation", "slope", "aspect", "fuel_model",
				     "fireline_intensity", "flame_length",
				     "weise_flame_height", "weise_flame_tilt",
				     "weise_whirl_height", "weise_whirl_radius",
				     "weise_angular_velocity", "weise_max_tang_vel",
				     "viegas_ROS", "viegas_eruptive_flag",
				     "viegas_ROS_excess", "viegas_flame_tilt", "viegas_slope_factor",
				     "scorch_height", "prob_ignition", "tree_mortality", "crown_activity",
					   "torching_ratio", "crowning_ratio", "energy_release_component",
				     "co2_emissions", "co_emissions", "pm25_emissions",
				   "arrival_time", "heat_per_unit_area", "vorticity_z",
				   "cbh", "cbd", "canopy_cover", "canopy_height", "burnout_time",
   "reaction_intensity", "residual_fuel", "shade_fraction",
   "torching_index_kmh", "crowning_index_kmh",
   "moisture_d1", "moisture_d10", "moisture_d100", "moisture_lh", "moisture_lw"
#else
				     , "farsite_dx", "farsite_dy", "R",
				     "spot_prob", "spot_count", "spot_dist", "spot_active", "fuel_consumption", "crown_fraction",
				     "albini_Hz", "albini_count", "albini_dist", "albini_active",
				     "elevation", "slope", "aspect", "fuel_model",
				     "fireline_intensity", "flame_length",
				     "weise_flame_height", "weise_flame_tilt",
				     "weise_whirl_height", "weise_whirl_radius",
				     "weise_angular_velocity", "weise_max_tang_vel",
				     "viegas_ROS", "viegas_eruptive_flag",
				     "viegas_ROS_excess", "viegas_flame_tilt", "viegas_slope_factor",
				     "scorch_height", "prob_ignition", "tree_mortality", "crown_activity",
					   "torching_ratio", "crowning_ratio", "energy_release_component",
				     "co2_emissions", "co_emissions", "pm25_emissions",
				   "arrival_time", "heat_per_unit_area", "vorticity_z",
				   "cbh", "cbd", "canopy_cover", "canopy_height", "burnout_time",
   "reaction_intensity", "residual_fuel", "shade_fraction",
   "torching_index_kmh", "crowning_index_kmh",
   "moisture_d1", "moisture_d10", "moisture_d100", "moisture_lh", "moisture_lw"
#endif
	};
	MultiFab plotmf(ba, dm, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 1 + 5 + WEISE_NCOMP + VIEGAS_NCOMP + ECOLOGY_NCOMP + EMISSIONS_NCOMP + 13 + FMS_NCOMP, 0);
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
	MultiFab::Copy(plotmf, weise_data, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 1 + 5, WEISE_NCOMP, 0);
	MultiFab::Copy(plotmf, viegas_data, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 1 + 5 + WEISE_NCOMP, VIEGAS_NCOMP, 0);
	MultiFab::Copy(plotmf, ecology_mf,   0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 1 + 5 + WEISE_NCOMP + VIEGAS_NCOMP, ECOLOGY_NCOMP, 0);
	MultiFab::Copy(plotmf, emissions_mf, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 1 + 5 + WEISE_NCOMP + VIEGAS_NCOMP + ECOLOGY_NCOMP, EMISSIONS_NCOMP, 0);
	MultiFab::Copy(plotmf, arrival_time_mf,       0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 1 + 5 + WEISE_NCOMP + VIEGAS_NCOMP + ECOLOGY_NCOMP + EMISSIONS_NCOMP, 1, 0);
	MultiFab::Copy(plotmf, heat_per_unit_area_mf, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 1 + 5 + WEISE_NCOMP + VIEGAS_NCOMP + ECOLOGY_NCOMP + EMISSIONS_NCOMP + 1, 1, 0);
	MultiFab::Copy(plotmf, vorticity_mf,          0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 1 + 5 + WEISE_NCOMP + VIEGAS_NCOMP + ECOLOGY_NCOMP + EMISSIONS_NCOMP + 2, 1, 0);
	MultiFab::Copy(plotmf, cbh_mf, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 1 + 5 + WEISE_NCOMP + VIEGAS_NCOMP + ECOLOGY_NCOMP + EMISSIONS_NCOMP + 3, 1, 0);
	MultiFab::Copy(plotmf, cbd_mf, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 1 + 5 + WEISE_NCOMP + VIEGAS_NCOMP + ECOLOGY_NCOMP + EMISSIONS_NCOMP + 4, 1, 0);
	MultiFab::Copy(plotmf, cc_mf,            0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 1 + 5 + WEISE_NCOMP + VIEGAS_NCOMP + ECOLOGY_NCOMP + EMISSIONS_NCOMP + 5, 1, 0);
	MultiFab::Copy(plotmf, canopy_height_mf, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 1 + 5 + WEISE_NCOMP + VIEGAS_NCOMP + ECOLOGY_NCOMP + EMISSIONS_NCOMP + 6, 1, 0);
	MultiFab::Copy(plotmf, burnout_time_mf,  0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 1 + 5 + WEISE_NCOMP + VIEGAS_NCOMP + ECOLOGY_NCOMP + EMISSIONS_NCOMP + 7, 1, 0);
	MultiFab::Copy(plotmf, reaction_intensity_mf, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 1 + 5 + WEISE_NCOMP + VIEGAS_NCOMP + ECOLOGY_NCOMP + EMISSIONS_NCOMP + 8, 1, 0);
	MultiFab::Copy(plotmf, residual_fuel_mf, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 1 + 5 + WEISE_NCOMP + VIEGAS_NCOMP + ECOLOGY_NCOMP + EMISSIONS_NCOMP + 9, 1, 0);
	MultiFab::Copy(plotmf, shade_fraction_mf, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 1 + 5 + WEISE_NCOMP + VIEGAS_NCOMP + ECOLOGY_NCOMP + EMISSIONS_NCOMP + 10, 1, 0);
	MultiFab::Copy(plotmf, ti_full_mf, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 1 + 5 + WEISE_NCOMP + VIEGAS_NCOMP + ECOLOGY_NCOMP + EMISSIONS_NCOMP + 11, 1, 0);
	MultiFab::Copy(plotmf, ci_full_mf, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 1 + 5 + WEISE_NCOMP + VIEGAS_NCOMP + ECOLOGY_NCOMP + EMISSIONS_NCOMP + 12, 1, 0);
	MultiFab::Copy(plotmf, spatial_moisture_mf, 0, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 1 + 5 + WEISE_NCOMP + VIEGAS_NCOMP + ECOLOGY_NCOMP + EMISSIONS_NCOMP + 13, FMS_NCOMP, 0);
	WriteSingleLevelPlotfile(buf, plotmf, names, geom, time, final_step);
	amrex::Print() << "Wrote final " << buf << "\n";
	{
	  amrex::Long n_burned = 0, n_perim = 0;
	  const auto dxx = geom.CellSize();
	  for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
	    const Box& bx = mfi.validbox();
	    auto const p = phi.const_array(mfi);
	    const IntVect bxlo = bx.smallEnd();
	    const IntVect bxhi = bx.bigEnd();
	    amrex::LoopOnCpu(bx, [&](int i, int j, int k) {
	      if (p(i,j,k,0) < Real(0.0)) {
	        ++n_burned;
	        bool on_perim = false;
	        if (i > bxlo[0]) on_perim |= (p(i-1,j,k,0) >= Real(0.0));
	        if (i < bxhi[0]) on_perim |= (p(i+1,j,k,0) >= Real(0.0));
	        if (j > bxlo[1]) on_perim |= (p(i,j-1,k,0) >= Real(0.0));
	        if (j < bxhi[1]) on_perim |= (p(i,j+1,k,0) >= Real(0.0));
	        if (on_perim) ++n_perim;
	      }
	    });
	  }
	  ParallelDescriptor::ReduceLongSum(n_burned);
	  ParallelDescriptor::ReduceLongSum(n_perim);
	  Real cell_area = dxx[0] * dxx[1];
	  amrex::Print() << "  Burned area: " << Real(n_burned)*cell_area/1.0e4 << " ha"
	                 << "  Perimeter: " << Real(n_perim)*std::sqrt(cell_area)/1000.0 << " km\n";
	}
	// Write negative phi x-y data files for final step
	char xy_buf[64];
	std::snprintf(xy_buf, sizeof(xy_buf), "phi_negative_%04d.dat", final_step);
	write_negative_phi_xy(phi, geom, xy_buf);
	
	std::snprintf(xy_buf, sizeof(xy_buf), "phi_envelope_%04d.dat", final_step);
	write_negative_phi_convex_hull(phi, geom, xy_buf);
	if (!inputs.fire_stats_file.empty()) append_fire_stats(phi, geom, &emissions_mf, final_step, time, inputs.fire_stats_file, &R_mf);
	if (inputs.write_perimeter_csv == 1) { char csv_buf[64]; std::snprintf(csv_buf, sizeof(csv_buf), "perimeter_%04d.csv", final_step); write_fire_perimeter_csv(phi, geom, csv_buf); }
	if (inputs.write_perimeter_geojson == 1) { char gjson_buf[64]; std::snprintf(gjson_buf, sizeof(gjson_buf), "perimeter_%04d.geojson", final_step); write_fire_perimeter_geojson(phi, geom, gjson_buf, final_step, time); }
      }
    }
  amrex::Finalize();
  return 0;
}
