#include "wildfire_includes.H"

static_assert(AMREX_SPACEDIM == 2,
              "wildfire_levelset requires a 2D AMReX build. "
              "Configure with -DLEVELSET_DIM_2D=ON (the default).");


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

    // ---------------- Fields -----------------------------------------------
    // All MultiFabs allocated and zero/sentinel-initialised in one call.
    WildfireFields f(ba, dm, inputs);
    // Convenience references so the rest of main.cpp is unchanged.
    MultiFab& phi                   = f.phi;
    MultiFab& vel                   = f.vel;
    MultiFab& vel_effective         = f.vel_effective;
    MultiFab& farsite_spread        = f.farsite_spread;
    MultiFab& spotting_data         = f.spotting_data;
    MultiFab& albini_data           = f.albini_data;
    MultiFab& ember_cascade_mf      = f.ember_cascade_mf;
    MultiFab& R_mf                  = f.R_mf;
    MultiFab& fuel_consumption_mf   = f.fuel_consumption_mf;
    MultiFab& crown_fire_fraction_mf= f.crown_fire_fraction_mf;
    MultiFab& fireline_intensity_mf = f.fireline_intensity_mf;
    MultiFab& flame_length_mf       = f.flame_length_mf;
    MultiFab& ecology_mf            = f.ecology_mf;
    MultiFab& emissions_mf          = f.emissions_mf;
    MultiFab& arrival_time_mf       = f.arrival_time_mf;
    MultiFab& residual_fuel_mf      = f.residual_fuel_mf;
    MultiFab& heat_flux_mf          = f.heat_flux_mf;
    MultiFab& ti_full_mf            = f.ti_full_mf;
    MultiFab& ci_full_mf            = f.ci_full_mf;
    MultiFab& weise_data            = f.weise_data;
    MultiFab& viegas_data           = f.viegas_data;
    MultiFab& shade_fraction_mf     = f.shade_fraction_mf;
    MultiFab& elevation_mf          = f.elevation_mf;
    MultiFab& slope_mf              = f.slope_mf;
    MultiFab& aspect_mf             = f.aspect_mf;
    MultiFab& fuel_model_mf         = f.fuel_model_mf;
    MultiFab& cbh_mf                = f.cbh_mf;
    MultiFab& cbd_mf                = f.cbd_mf;
    MultiFab& cc_mf                 = f.cc_mf;
    MultiFab& canopy_height_mf      = f.canopy_height_mf;
    MultiFab& spatial_moisture_mf   = f.spatial_moisture_mf;
    MultiFab& temperature_mf        = f.temperature_mf;
    MultiFab& humidity_mf           = f.humidity_mf;
    MultiFab& plume_rise_mf         = f.plume_rise_mf;
    std::unique_ptr<MultiFab>& terrain_slopes = f.terrain_slopes;

    // New diagnostic fields
    MultiFab& spread_dir_mf         = f.spread_dir_mf;
    MultiFab& spot_catch_prob_mf    = f.spot_catch_prob_mf;
    MultiFab& spotting_lineage_mf   = f.spotting_lineage_mf;
    MultiFab& ros_at_arrival_mf     = f.ros_at_arrival_mf;
    MultiFab& fl_exceedance_mf      = f.fl_exceedance_mf;

    // ---------------- Initialise phi and mark initially-burned cells --------
    bool use_indicator = (inputs.propagation_method == "farsite");
    Real time = 0.0;
    int  restart_step = 0;
    init_phi_source(phi, arrival_time_mf, geom, inputs,
                    use_indicator, restart_step, time);

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

    // ---- 3-D wind from massconsistent_amr plt file ----
    // When albini_spotting.use_3d_wind = 1, read the plt file once before the
    // time loop.  The PltWindData struct stores flat 1-D GPU arrays of
    // (x, y, z, u, v, w) and a precomputed height-averaged 2-D wind field.
    PltWindData albini_plt_wind;
    if (inputs.albini_spotting.enable == 1 &&
        inputs.albini_spotting.use_3d_wind == 1 &&
        !inputs.albini_spotting.plt_wind_file.empty()) {
        if (!read_plt_wind_file(inputs.albini_spotting.plt_wind_file, albini_plt_wind)) {
            amrex::Abort("Failed to read 3-D wind plt file for Albini spotting: "
                         + inputs.albini_spotting.plt_wind_file);
        }
    }

    // ---- Fire acceleration state initialization ----
    // For FARSITE temporal model, allocate per-cell state tracking
    if (inputs.acceleration.enable == 1 && inputs.acceleration.use_temporal == 1) {
        f.accel_state_mf = std::make_unique<MultiFab>(ba, dm, 3, 0);
        f.accel_state_mf->setVal(Real(0.0));
        amrex::Print() << "Fire acceleration: allocated temporal state tracking (3 components per cell)\n";
    }

    // ---- 3-D wind for flux-based ember cascade model ----
    // When ember_cascade.use_3d_wind = 1 or ember_cascade.require_3d_wind = 1,
    // read the plt file once before the time loop.  A missing or unreadable
    // file is a fatal error when require_3d_wind = 1.
    PltWindData ember_cascade_plt_wind;
    if (inputs.ember_cascade.enable == 1 &&
        (inputs.ember_cascade.use_3d_wind == 1 || inputs.ember_cascade.require_3d_wind == 1) &&
        !inputs.ember_cascade.plt_wind_file.empty()) {
        if (!read_plt_wind_file(inputs.ember_cascade.plt_wind_file, ember_cascade_plt_wind)) {
            amrex::Abort("Failed to read 3-D wind plt file for ember cascade: "
                         + inputs.ember_cascade.plt_wind_file);
        }
    }

    // ---- Terrain, landscape, crown layers, and spatial moisture ------------
    bool has_spatial_crown    = false;
    bool has_spatial_moisture = false;
    setup_terrain(terrain_slopes, elevation_mf, slope_mf, aspect_mf,
                  fuel_model_mf, cbh_mf, cbd_mf, cc_mf, canopy_height_mf,
                  has_spatial_crown, spatial_moisture_mf, has_spatial_moisture,
                  ba, dm, geom, inputs);

    // ---- Topographic horizon angles (FARSITE 8-direction scan) -------------
    // Precomputed once here because elevation never changes during a run.
    // Guarded by use_topographic_horizon so users can skip the expensive
    // MPI global-gather + O(N^2) CPU sweep when it is not needed.
    std::unique_ptr<MultiFab> horizon_mf;
    if (inputs.solar_radiation.enable == 1 &&
        inputs.solar_radiation.use_topographic_horizon == 1) {
        horizon_mf = std::make_unique<MultiFab>(ba, dm, 8, 0);
        compute_topographic_horizon_angles(
            *horizon_mf, elevation_mf, geom,
            inputs.solar_radiation.horizon_scan_max_dist_m);
    }

    // ---- GPU lookup tables and crown/Balbi coefficients -------------------
    FuelTableData ftd = setup_fuel_tables(inputs, fuel_model_mf);
    // Unpack frequently used pointers for backward compatibility with time loop.
    Gpu::DeviceVector<RothermelComputed>& d_fuel_table   = ftd.d_fuel_table;
    const RothermelComputed*&             d_fuel_table_ptr= ftd.d_fuel_table_ptr;
    int&                                  fuel_table_size = ftd.fuel_table_size;
    const BalbiComputed*&                 d_balbi_table_ptr=ftd.d_balbi_table_ptr;
    int&                                  balbi_table_size= ftd.balbi_table_size;
    BalbiComputed&                        bc_global_default=ftd.bc_global_default;
    const CruzCrownComputed*&             ccc_ptr         = ftd.ccc_ptr;

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
        apply_solar_radiation_step(inputs, shade_fraction_mf,
                                   slope_mf, aspect_mf,
                                   has_spatial_crown, cc_mf,
                                   spatial_moisture_mf,
                                   /*elapsed_s=*/Real(0.0),
                                   horizon_mf.get(),
                                   /*print_position=*/true);
    }

    // ---- All time-varying schedules, weather, and events ------------------
    ScheduleData sd = setup_schedules(inputs, vel, phi, geom, ba, dm,
                                      use_indicator);
    // Unpack into names used by the rest of main.cpp.
    FMCSchedule&          fmc_sched       = sd.fmc_sched;
    HerbMoistureSchedule& herb_sched      = sd.herb_sched;
    WindDirSchedule&      wind_dir_sched  = sd.wind_dir_sched;
    WtrWeatherData&       wtr_data        = sd.wtr_data;
    const bool            wtr_active      = sd.wtr_active;
    MultiWtrWeather&      multi_wtr       = sd.multi_wtr;
    const bool            multi_wtr_active= sd.multi_wtr_active;
    RetardantDropList&    retardant_drops = sd.retardant_drops;
    IgnitionSchedule&     ignition_sched  = sd.ignition_sched;
    PrecipState&          precip_state    = sd.precip_state;
    std::vector<std::pair<double,double>>& precip_schedule = sd.precip_schedule;

    // ---- Conditional weather table (ERC/BI/SC percentile lookup) -----------
    ConditionalWeatherTable cond_weather_table =
        load_conditional_weather(inputs.conditional_weather_file);

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
    const bool wind_terrain_modifies_vel = wind_terrain_modifies_velocity(inputs);

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
    Real dt=(inputs.propagation_method == "farsite") ? inputs.farsite.dt : 10.0; 
    const bool use_levelset = (inputs.propagation_method == "levelset");
    const bool use_mtt      = (inputs.propagation_method == "mtt");
    if (use_levelset)
      {
        // Apply wind-terrain velocity modification (Options 3-8) before ROS computation
        apply_wind_terrain_effective_velocity(vel_effective, vel, terrain_slopes.get(), geom, inputs);

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
            // Compute Rothermel wind speed R (levelset path - no cell size correction)
            compute_rothermel_R(R_mf, vel_for_model, geom, inputs.rothermel,
                                 terrain_slopes.get(),
                                 !inputs.rothermel.landscape_file.empty() ? &fuel_model_mf : nullptr,
                                 d_fuel_table_ptr, fuel_table_size,
                                 has_spatial_crown ? &cc_mf : nullptr,
                                 has_spatial_crown ? &canopy_height_mf : nullptr,
                                 0,  // cellsize_enable = 0 (disabled for levelset)
                                 30.0, 0.1);  // default parameters (unused)
        }
        dt = compute_dt(R_mf, geom, inputs.cfl);
        amrex::Print() << "Computed dt = " << dt << "\n";
      } else if (use_mtt) {
      // For MTT: compute the ROS field now (same as levelset path),
      // then run the Dijkstra fast-march to fill arrival_time_mf,
      // and set a sensible dt based on the ROS field.
      {
        apply_wind_terrain_effective_velocity(vel_effective, vel, terrain_slopes.get(), geom, inputs);
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
                                 has_spatial_crown ? &canopy_height_mf : nullptr,
                                 0,  // cellsize_enable = 0 (disabled for MTT)
                                 30.0, 0.1);  // default parameters (unused)
        }
        dt = compute_dt(R_mf, geom, inputs.cfl);
        amrex::Print() << "MTT: pre-computing arrival times (dt = " << dt << ") ...\n";
        compute_mtt_arrival_times(arrival_time_mf, phi, R_mf, geom);
        // Set phi from arrival times at t=0
        apply_mtt_phi_update(phi, arrival_time_mf, Real(0.0));
        fill_boundary_extrap(phi, geom);
      }
      } else {
      // FARSITE: compute initial ROS field and CFL-based dt (same as levelset/MTT paths).
      apply_wind_terrain_effective_velocity(vel_effective, vel, terrain_slopes.get(), geom, inputs);
      if (heat_flux_active) {
          apply_heatflux_wind(vel_effective, vel, heat_flux_mf, &phi, inputs.heat_flux);
      }
      const MultiFab& vel_for_farsite_init = (wind_terrain_modifies_vel || heat_flux_active)
                                              ? vel_effective : vel;
      if (inputs.fire_spread_model == "balbi") {
          compute_balbi_R(R_mf, vel_for_farsite_init, geom, inputs.rothermel, inputs.balbi,
                           terrain_slopes.get(),
                           !inputs.rothermel.landscape_file.empty() ? &fuel_model_mf : nullptr,
                           d_balbi_table_ptr, balbi_table_size,
                           heat_flux_active ? &heat_flux_mf : nullptr,
                           heat_flux_active ? &inputs.heat_flux : nullptr);
      } else if (inputs.fire_spread_model == "cheney_gould") {
          compute_cheney_gould_R(R_mf, vel_for_farsite_init, inputs.cheney_gould);
      } else if (inputs.fire_spread_model == "cruz_crown") {
          compute_cruz_crown_R(R_mf, vel_for_farsite_init, inputs.cruz_crown);
      } else if (inputs.fire_spread_model == "fbp_o1a" ||
                 inputs.fire_spread_model == "fbp_o1b" ||
                 inputs.fire_spread_model == "fbp_s1"  ||
                 inputs.fire_spread_model == "fbp_s2"  ||
                 inputs.fire_spread_model == "fbp_s3") {
          compute_fbp_R(R_mf, vel_for_farsite_init, inputs.fbp);
      } else if (inputs.fire_spread_model == "lautenberger") {
          compute_lautenberger_R(R_mf, vel_for_farsite_init, inputs.rothermel, inputs.lautenberger,
                                  terrain_slopes.get());
      } else {
          compute_rothermel_R(R_mf, vel_for_farsite_init, geom, inputs.rothermel,
                               terrain_slopes.get(),
                               !inputs.rothermel.landscape_file.empty() ? &fuel_model_mf : nullptr,
                               d_fuel_table_ptr, fuel_table_size,
                               has_spatial_crown ? &cc_mf : nullptr,
                               has_spatial_crown ? &canopy_height_mf : nullptr,
                               inputs.cellsize.enable,  // cellsize_enable
                               inputs.cellsize.dx_ref,  // cellsize_dx_ref
                               inputs.cellsize.correction_exponent);  // cellsize_correction_exp
      }
      amrex::Print() << "Using FARSITE propagation; dt = " << dt << "\n";
    }
    compute_fire_behavior(fireline_intensity_mf, flame_length_mf, R_mf, inputs.rothermel);
    // Update flame-length exceedance raster (element-wise max over time)
    for (MFIter mfi(fl_exceedance_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto ex  = fl_exceedance_mf.array(mfi);
        auto const fl = flame_length_mf.const_array(mfi);
        ParallelFor(bx, [ex, fl] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            ex(i, j, k) = amrex::max(ex(i, j, k), fl(i, j, k));
        });
    }
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
            ParallelFor(bx, [R, eco, v, fi, use_sp_i, cbd_arr_i, CBD_g_i, mf_val_i, MC10_i, use_cruz_i, use_roth91_i, use_passive_i, I_o_i] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
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

    // ---- Initial smoke plume rise (t=0) ----
    if (inputs.smoke_plume.enable == 1) {
        compute_smoke_plume_rise(plume_rise_mf, fireline_intensity_mf,
                                 vel, inputs.smoke_plume);
    }

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

    // ---- Satellite fire detection assimilation (initial condition) --------
    // When satellite.enable = 1 and satellite.use_as_ic = 1, fetch active-fire
    // detections at t=0 and apply them as additional ignition points before the
    // first plotfile write.  The SatelliteState tracks the last fetch time so
    // the mid-simulation re-fetch interval is measured from t=0.
    SatelliteState sat_state;
    if (inputs.satellite.enable == 1 && inputs.satellite.use_as_ic == 1) {
        auto sat_pts = fetch_satellite_fire_points(inputs.satellite, Real(0.0));
        if (!sat_pts.empty()) {
            apply_satellite_fire_points(phi, geom, sat_pts,
                                        inputs.satellite.detection_radius_m,
                                        /*suppress_if_burning=*/false);
            fill_boundary_extrap(phi, geom);
            amrex::Print() << "SatelliteAssimilation: applied "
                           << sat_pts.size()
                           << " detection(s) as initial condition.\n";
        }
        sat_state.last_fetch_time = Real(0.0);
    }

    // ---------------- Initialize fire statistics CSV -------------------
    if (!inputs.fire_stats_file.empty())
        write_fire_stats_header(inputs.fire_stats_file);

    // Initialize Fire Spread Atlas (.fsa) and post-processing stats (.pst)
    if (!inputs.fsa_file.empty())
        write_fsa_header(inputs.fsa_file);
    if (!inputs.pst_file.empty())
        write_pst_header(inputs.pst_file);

    // ---------------- Write initial plotfile ---------------
    write_wildfire_plotfile(f, ba, dm, geom, inputs, restart_step, time, ftd,
                            /*write_stats=*/false, /*is_final=*/false);

    // ---- Simulation calendar date/time helper ----
    // Uses sim_datetime.H shared helper for "YYYY-MM-DD HH:MM" formatting.
    // Inputs: sim_start_year/month/day from inputs (falls back to solar_radiation).
    const int   _yr0 = inputs.sim_start_year;
    const int   _mo0 = inputs.sim_start_month;
    const int   _dy0 = inputs.sim_start_day;
    const double _hr0 = (inputs.solar_radiation.enable == 1)
                        ? static_cast<double>(inputs.solar_radiation.sim_start_hour)
                        : 0.0;

    auto sim_datetime_string = [&](double elapsed_s) -> std::string {
        return WildfireDatetime::elapsed_to_datetime(
            _yr0, _mo0, _dy0, _hr0, elapsed_s);
    };

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
      // When a hard stop time is specified, clamp the next timestep so the
      // simulation ends exactly at final_time rather than overshooting by a
      // full CFL step (which could be much larger than the remaining time).
      if (use_final_time)
        dt = std::min(dt, inputs.final_time - time);
      fill_boundary_extrap(phi, geom);
      const Real dt_step = dt;
      {
        // Print time with optional calendar datetime
        const std::string dt_str = sim_datetime_string(static_cast<double>(time));
        if (dt_str.empty()) {
            amrex::Print() << "Time:" << time << " with timestep:" << dt_step << "\n";
        } else {
            amrex::Print() << "Time:" << time << " (" << dt_str
                           << ") with timestep:" << dt_step << "\n";
        }
      }
      
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
          if (multi_wtr_active && inputs.diurnal_moisture.enable == 1) {
              // Use spatially-varying T/RH from multi-station weather
              apply_solar_radiation_step(inputs, shade_fraction_mf,
                                         slope_mf, aspect_mf,
                                         has_spatial_crown, cc_mf,
                                         spatial_moisture_mf,
                                         temperature_mf, humidity_mf,
                                         static_cast<amrex::Real>(time),
                                         horizon_mf.get());
          } else {
              // Use diurnal sinusoidal T/RH model
              apply_solar_radiation_step(inputs, shade_fraction_mf,
                                         slope_mf, aspect_mf,
                                         has_spatial_crown, cc_mf,
                                         spatial_moisture_mf,
                                         static_cast<amrex::Real>(time),
                                         horizon_mf.get());
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

          // Spatial T/RH interpolation for diurnal moisture model
          if (inputs.diurnal_moisture.enable == 1) {
              // Interpolate T and RH to per-cell MultiFabs
              apply_multi_wtr_TRH_to_spatial(temperature_mf, humidity_mf,
                                             geom, multi_wtr,
                                             static_cast<double>(time));

              // Fallback: domain-mean T/RH for global moisture parameters
              // (used when spatial moisture is not enabled)
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
              ParallelFor(bx, [sm, elev, shade, lapse, elev0, T_ref_r, RH_ref_r, sol_C, MAGNUS_A, MAGNUS_B] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
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
      // Apply wind-terrain velocity modification (Options 3-8) before ROS computation
      apply_wind_terrain_effective_velocity(vel_effective, vel, terrain_slopes.get(), geom, inputs);

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
                               has_spatial_crown ? &canopy_height_mf : nullptr,
                               inputs.cellsize.enable,  // cellsize_enable
                               inputs.cellsize.dx_ref,  // cellsize_dx_ref
                               inputs.cellsize.correction_exponent);  // cellsize_correction_exp
      }
      compute_fire_behavior(fireline_intensity_mf, flame_length_mf, R_mf, inputs.rothermel);
      // Update flame-length exceedance raster (element-wise max over time)
      for (MFIter mfi(fl_exceedance_mf); mfi.isValid(); ++mfi) {
          const Box& bx = mfi.validbox();
          auto ex  = fl_exceedance_mf.array(mfi);
          auto const fl = flame_length_mf.const_array(mfi);
          ParallelFor(bx, [ex, fl] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
              ex(i, j, k) = amrex::max(ex(i, j, k), fl(i, j, k));
          });
      }
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

      // ============================================================================
      // Feature Integration: Integrate 10 new wildfire features
      // ============================================================================
       
      // Feature 2: NFDRS Fire Danger Class (operational categories)
      // Always compute for output
      compute_nfdrs_danger_class(f.fire_intensity_class_mf, fireline_intensity_mf);

      // Feature 3: Crown Fraction Burned (CFB) Diagnostic
      if (inputs.crown_fraction.enable == 1) {
          // CFB = (I_B - I_B_surface) / (I_B_active - I_B_surface)
          for (MFIter mfi(f.crown_fraction_burned_mf); mfi.isValid(); ++mfi) {
              const Box& bx = mfi.validbox();
              auto cfb_arr = f.crown_fraction_burned_mf.array(mfi);
              auto const fi_arr = fireline_intensity_mf.const_array(mfi);
              ParallelFor(bx, [cfb_arr, fi_arr] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                  // CFB based on fire intensity magnitude as proxy
                  Real I_B = fi_arr(i, j, k);
                  if (I_B > Real(100.0)) {
                      cfb_arr(i, j, k) = amrex::min(Real(1.0), (I_B - Real(100.0)) / Real(5000.0));
                  } else {
                      cfb_arr(i, j, k) = Real(0.0);
                  }
              });
          }
      }

      // Feature 4: Effective Wind Speed (combined wind + slope)
      if (inputs.effective_wind.enable == 1 && terrain_slopes.get() != nullptr) {
          compute_effective_wind_speed_field(f.effective_wind_speed_mf, vel, 
                                            *(terrain_slopes.get()));
      }

      // Feature 5: Thomas Flame Length Model (alternative to Byram)
      if (inputs.flame_length_model.model == "thomas") {
          // L = 0.0266 * I^0.667 (vs Byram: L = 0.0775 * I^0.46)
          for (MFIter mfi(flame_length_mf); mfi.isValid(); ++mfi) {
              const Box& bx = mfi.validbox();
              auto fl_arr = flame_length_mf.array(mfi);
              auto const fi_arr = fireline_intensity_mf.const_array(mfi);
              ParallelFor(bx, [fl_arr, fi_arr] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                  Real I_B = fi_arr(i, j, k);
                  fl_arr(i, j, k) = Real(0.0266) * std::pow(amrex::max(I_B, Real(0.0)), Real(0.667));
              });
          }
      }

      // Feature 6: Fuel Boundary Smoothing
      if (inputs.fuel_boundary.enable == 1 && !inputs.rothermel.landscape_file.empty()) {
          apply_fuel_boundary_smoothing(R_mf, fuel_model_mf,
                                       inputs.n_cell_x, inputs.n_cell_y, inputs.n_cell_z);
      }

      // Feature 7: CSIRO Grassfire Acceleration (non-equilibrium growth)
      if (inputs.grassfire_accel.enable == 1) {
          // a(t) = 1 - exp(-t/t_accel), applied to R_mf for grassfire scenarios
          Real accel_factor = Real(1.0) - std::exp(-dt_step / inputs.grassfire_accel.t_accel);
          for (MFIter mfi(R_mf); mfi.isValid(); ++mfi) {
              const Box& bx = mfi.validbox();
              auto R_arr = R_mf.array(mfi);
              ParallelFor(bx, [R_arr, accel_factor] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                  R_arr(i, j, k) *= (Real(1.0) + accel_factor);
              });
          }
      }

      // Feature 8: Burnout Time Separation (flaming vs smoldering)
      if (inputs.burnout_separation.enable == 1) {
          // Split residence time by fuel type; store in burnout_phases_mf
          // comp 0: flaming_duration, comp 1: smoldering_duration
          const Real tau_residence = Real(inputs.farsite.tau_residence);
          const Real flaming_frac = Real(inputs.burnout_separation.flaming_fraction_fine);
          for (MFIter mfi(f.burnout_phases_mf); mfi.isValid(); ++mfi) {
              const Box& bx = mfi.validbox();
              auto phases_arr = f.burnout_phases_mf.array(mfi);
              ParallelFor(bx, [phases_arr, tau_residence, flaming_frac] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                  // Could be refined based on fuel_model_mf if available
                  phases_arr(i, j, k, 0) = tau_residence * flaming_frac;      // flaming
                  phases_arr(i, j, k, 1) = tau_residence * (Real(1.0) - flaming_frac); // smoldering
              });
          }
      }

      // ---- Conditional weather update (ERC / BI / SC percentile table) ----
      // Select comp from ecology_mf based on trigger_index:
      //   6 = energy_release_component (ERC)
      //   7 = nfdrs_spread_component   (SC)
      //   8 = nfdrs_burning_index      (BI)
      if (!cond_weather_table.empty()) {
          const std::string& tri = inputs.conditional_weather_trigger;
          int ecomp = 6;  // default: ERC
          if (tri == "sc") ecomp = 7;
          else if (tri == "bi") ecomp = 8;
          // Domain-mean of selected ecology component
          Real index_sum = Real(0.0);
          amrex::Long n_cells = 0;
          for (MFIter mfi(ecology_mf); mfi.isValid(); ++mfi) {
              const Box& bx = mfi.validbox();
              auto const eco = ecology_mf.const_array(mfi);
              amrex::LoopOnCpu(bx, [&](int i, int j, int k) {
                  index_sum += eco(i, j, k, ecomp);
                  ++n_cells;
              });
          }
          ParallelDescriptor::ReduceRealSum(index_sum);
          ParallelDescriptor::ReduceLongSum(n_cells);
          const double mean_index = (n_cells > 0)
              ? static_cast<double>(index_sum) / static_cast<double>(n_cells)
              : 0.0;
          const int sel = apply_conditional_weather(
              cond_weather_table, mean_index,
              inputs.conditional_weather_trigger,
              inputs.rothermel, inputs.ux, inputs.uy);
          if (sel >= 0) {
              // Recompute wind field with updated ux/uy
              vel.setVal(inputs.ux, 0, 1);
              vel.setVal(inputs.uy, 1, 1);
          }
      }

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

              ParallelFor(bx, [R, eco, v, fi, use_sp, cbd_arr, CBD_g, mf_val, use_roth1991_tl, use_cruz_tl, MC10_tl, use_passive_tl, I_o_tl_k] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
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

      // ---- Smoke plume-rise model (Briggs 1965) ----
      // Computes per-cell final plume-rise height from Byram fireline intensity.
      if (inputs.smoke_plume.enable == 1) {
          compute_smoke_plume_rise(plume_rise_mf, fireline_intensity_mf,
                                   vel, inputs.smoke_plume);
          // Print domain maximum plume rise for situational awareness
          const Real max_plume = plume_rise_mf.max(0);
          if (max_plume > Real(0.0))
              amrex::Print() << "  Max smoke plume rise: " << max_plume << " m\n";
      }

      // Feature 10: Fire acceleration model (Rothermel 1983 / Catchpole et al. 1992
      //             or FARSITE temporal model McAlpine & Wakimoto 1991)
      // Scales R_mf to account for the slower spread of small fires or newly-changed
      // wind conditions before quasi-steady-state is reached.
      // When disabled (default) or fire is at equilibrium, this is a no-op.
      apply_fire_acceleration(R_mf, phi, geom, inputs.acceleration, dt,
                              f.accel_state_mf.get());

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
              // Print only on transitions from active to inactive (first inactive step
              // per inactive block) to avoid flooding the log every time step.
              static bool s_prev_bp_active = true;
              if (s_prev_bp_active) {
                  amrex::Print() << "  Burn period inactive at hour " << clock_hour
                                 << " (window " << sh << ":00-" << eh
                                 << ":00) – spread suppressed.\n";
              }
              s_prev_bp_active = false;
          } else {
              static bool s_prev_bp_active = true;
              if (!s_prev_bp_active) {
                  amrex::Print() << "  Burn period active at hour " << clock_hour
                                 << " – spread resumed.\n";
              }
              s_prev_bp_active = true;
          }
      }

      // ---- Spread direction raster: R_mf × normalized wind unit vector ----
      // Computed after all R_mf modifications (crown, acceleration, retardant,
      // burn-period gate) so the raster faithfully represents the conditions
      // used for fire propagation this timestep.
      //   spread_dir_mf[0] = R_mf × (ux / |u|)  [m/s in x]
      //   spread_dir_mf[1] = R_mf × (uy / |u|)  [m/s in y]
      // Cells with negligible wind (|u| < min_wind_mag) are set to zero to
      // avoid division by near-zero values.
      {
          constexpr Real min_wind_mag = Real(1.0e-10); // [m/s] – negligible wind guard
          for (MFIter mfi(spread_dir_mf); mfi.isValid(); ++mfi) {
              const Box& bx = mfi.validbox();
              auto const v  = vel_for_model.const_array(mfi);
              auto const R  = R_mf.const_array(mfi);
              auto sd = spread_dir_mf.array(mfi);
              ParallelFor(bx, [v, R, sd, min_wind_mag] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                  Real ux = v(i, j, k, 0);
                  Real uy = v(i, j, k, 1);
                  Real wind_mag = std::sqrt(ux * ux + uy * uy);
                  Real R_val = R(i, j, k);
                  if (wind_mag > min_wind_mag) {
                      sd(i, j, k, 0) = R_val * ux / wind_mag;
                      sd(i, j, k, 1) = R_val * uy / wind_mag;
                  } else {
                      sd(i, j, k, 0) = Real(0.0);
                      sd(i, j, k, 1) = Real(0.0);
                  }
              });
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
	// Pass R_mf so the ellipse uses the ROS from the configured firespread model.
	compute_farsite_spread(phi, vel_for_model, farsite_spread, geom, dt_step, inputs.rothermel, inputs.farsite, inputs.crown, R_mf, terrain_slopes.get(), &fuel_consumption_mf, &crown_fire_fraction_mf, has_spatial_crown ? &cc_mf : nullptr, has_spatial_crown ? &canopy_height_mf : nullptr, ccc_ptr);
	// Update dt for next step using CFL condition (original FARSITE method).
	//dt = compute_dt(R_mf, geom, inputs.cfl);
	
	// --- Step 5: Apply crown/spotting sub-models
	if (inputs.spotting.enable == 1 && (step % inputs.spotting.check_interval == 0)) {
	  compute_spotting_probability(spotting_data, phi, vel, geom, inputs.rothermel, inputs.spotting, terrain_slopes.get(), &ecology_mf);
	  generate_firebrand_spots(phi, spotting_data, vel, geom, inputs.spotting, step,
	                           !inputs.rothermel.landscape_file.empty() ? &fuel_model_mf : nullptr,
	                           !inputs.rothermel.landscape_file.empty() ? &inputs.rothermel.landscape_fuel_type : nullptr,
	                           (inputs.fuel_depletion.adjust_spotting_reentry == 1)
	                               ? &residual_fuel_mf : nullptr,
	                           inputs.fuel_depletion.spotting_fuel_threshold,
	                           &spotting_lineage_mf,
	                           &ecology_mf);
	}
	// Albini (1983) firebrand spotting with 2-D trajectory integration
	if (inputs.albini_spotting.enable == 1 && (step % inputs.albini_spotting.check_interval == 0)) {
	  compute_albini_spotting(phi, albini_data, vel, R_mf, geom,
	                          inputs.rothermel, inputs.albini_spotting, step,
	                          !inputs.rothermel.landscape_file.empty() ? &fuel_model_mf : nullptr,
	                          !inputs.rothermel.landscape_file.empty() ? &inputs.rothermel.landscape_fuel_type : nullptr,
	                          (inputs.albini_spotting.use_3d_wind == 1 && albini_plt_wind.valid)
	                              ? &albini_plt_wind : nullptr,
	                          (inputs.fuel_depletion.adjust_spotting_reentry == 1)
	                              ? &residual_fuel_mf : nullptr,
	                          inputs.fuel_depletion.spotting_fuel_threshold,
	                          &spotting_lineage_mf,
	                          (inputs.weise_biging.enable == 1 && inputs.weise_biging.enhance_spotting == 1)
	                              ? &weise_data : nullptr,
	                          (inputs.weise_biging.enable == 1 && inputs.weise_biging.enhance_spotting == 1)
	                              ? &inputs.weise_biging : nullptr);
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

	// ---- Spot-fire catching probability raster ----
	// Computed whenever at least one spotting model is active.  Uses the
	// first-active model's P_catch (firebrand takes priority over Albini).
	if (inputs.spotting.enable == 1 || inputs.albini_spotting.enable == 1) {
	    const Real P_catch_use = (inputs.spotting.enable == 1)
	        ? static_cast<Real>(inputs.spotting.P_catch)
	        : static_cast<Real>(inputs.albini_spotting.P_catch);
	    compute_spot_catch_probability(spot_catch_prob_mf, inputs.rothermel,
	                                   P_catch_use,
	                                   has_spatial_moisture ? &spatial_moisture_mf : nullptr);
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

	// Flux-based ember cascade model (plume-rise driven)
	// Runs independently of the single-brand Albini model; both can be
	// active simultaneously for complementary coverage.
	if (inputs.ember_cascade.enable == 1 &&
	    (step % inputs.ember_cascade.check_interval == 0)) {
	  const Real dt_check_ec = dt_step * static_cast<Real>(inputs.ember_cascade.check_interval);
	  compute_ember_cascade_flux(
	      phi, ember_cascade_mf, fireline_intensity_mf,
	      vel, geom, inputs.ember_cascade,
	      dt_check_ec, step,
	      (ember_cascade_plt_wind.valid) ? &ember_cascade_plt_wind : nullptr);
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
	  int  niters   = static_cast<int>(std::ceil(WildfireFields::ng_phi / 0.5)); // = 6 for ng_phi=3
	  reinitialize_phi(phi, geom, niters, dtau);
	  amrex::Print() << "Reinitialized phi with dtau = " << dtau << " and niters = " << niters << "\n";
	}
      }

      time += dt_step;

      // --- Update arrival time: mark cells that first became burned this step
      // Also freeze the instantaneous ROS (ros_at_arrival_mf) at the moment
      // each cell first ignites so the spatial spread-rate record persists
      // through later burn-period gates, retardant drops, and crown overrides.
      {
        const Real cur_time = time;
        for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
          const Box& bx = mfi.validbox();
          auto const p   = phi.const_array(mfi);
          auto const R   = R_mf.const_array(mfi);
          auto       at  = arrival_time_mf.array(mfi);
          auto       rar = ros_at_arrival_mf.array(mfi);
          ParallelFor(bx, [p, R, at, rar, cur_time] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            if (p(i, j, k) < Real(0.0) && at(i, j, k) < Real(0.0)) {
              at(i, j, k)  = cur_time;
              rar(i, j, k) = R(i, j, k);  // freeze ROS at time of arrival
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
          ParallelFor(bx, [at, rf, cur_time, tau_b] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
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
            ParallelFor(bx, [phi_arr, at, rf, R] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
              // Scale only unburned cells where fuel was previously depleted (re-entry)
              if (at(i, j, k) >= Real(0.0) && phi_arr(i, j, k) >= Real(0.0)) {
                R(i, j, k) *= rf(i, j, k);
              }
            });
          }
        }
      }

      // Feature 10: Post-Frontal Smoldering - Track time since burn and compute residual heat release
      if (inputs.post_frontal.enable == 1) {
          const Real cur_time = time;
          // Update time_since_burn_mf: elapsed time since each cell first burned
          for (MFIter mfi(f.time_since_burn_mf); mfi.isValid(); ++mfi) {
              const Box& bx = mfi.validbox();
              auto const at = arrival_time_mf.const_array(mfi);
              auto       tsb = f.time_since_burn_mf.array(mfi);
              ParallelFor(bx, [at, tsb, cur_time] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                  if (at(i, j, k) >= Real(0.0)) {
                      tsb(i, j, k) = cur_time - at(i, j, k);
                  }
              });
          }
          
          // Compute residual smoldering heat release: decays exponentially after flame front
          // I_smolder(t) = I_base * exp(-(t - t_arrival) / tau_decay)
          for (MFIter mfi(f.residual_heat_release_mf); mfi.isValid(); ++mfi) {
              const Box& bx = mfi.validbox();
              auto const at = arrival_time_mf.const_array(mfi);
              auto const fi = fireline_intensity_mf.const_array(mfi);
              auto       rhr = f.residual_heat_release_mf.array(mfi);
              // Use default decay time constant (1 hour for medium fuels)
              const Real tau_decay = Real(3600.0);  // 1 hour in seconds
              const Real cur_time_c = cur_time;
              ParallelFor(bx, [at, fi, rhr, tau_decay, cur_time_c] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                  const Real t_arrive = at(i, j, k);
                  if (t_arrive >= Real(0.0)) {
                      const Real elapsed = cur_time_c - t_arrive;
                      const Real I_base = fi(i, j, k) * Real(0.1);  // Base smoldering intensity 10% of flaming
                      rhr(i, j, k) = I_base * std::exp(-elapsed / tau_decay);
                  } else {
                      rhr(i, j, k) = Real(0.0);
                  }
              });
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

      // ---- Satellite fire detection assimilation (mid-simulation re-marking) ----
      // When satellite.enable = 1 and satellite.use_mid_sim = 1, fetch new
      // active-fire detections at every fetch_interval_s simulation seconds and
      // merge them into phi.  Already-burning cells are preserved when
      // satellite.suppress_if_burning = 1 (default), preventing the satellite
      // from extinguishing the simulated fire front.
      if (inputs.satellite.enable == 1 && inputs.satellite.use_mid_sim == 1) {
          const Real elapsed_since_fetch = time - sat_state.last_fetch_time;
          if (elapsed_since_fetch >= inputs.satellite.fetch_interval_s) {
              auto sat_pts = fetch_satellite_fire_points(inputs.satellite,
                                                         static_cast<Real>(time));
              if (!sat_pts.empty()) {
                  apply_satellite_fire_points(phi, geom, sat_pts,
                                              inputs.satellite.detection_radius_m,
                                              inputs.satellite.suppress_if_burning == 1);
                  fill_boundary_extrap(phi, geom);
                  amrex::Print() << "SatelliteAssimilation: applied "
                                 << sat_pts.size()
                                 << " detection(s) at t=" << time << " s.\n";
              }
              sat_state.last_fetch_time = time;
          }
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
            if (inputs.write_perimeter_kml == 1) {
              char iso_kml_buf[128];
              std::snprintf(iso_kml_buf, sizeof(iso_kml_buf),
                            "isochrone_%06d.kml", iso_i);
              write_fire_perimeter_kml(phi, geom, iso_kml_buf, iso_i, time,
                                       inputs.kml_utm_zone, inputs.kml_utm_northern);
            }
          }
          last_isochrone_idx = cur_iso_idx;
        }
      }

      // --- Step 7: Update states, record outputs, step time
      if (inputs.plot_int > 0 && (step % inputs.plot_int == 0)) {
          write_wildfire_plotfile(f, ba, dm, geom, inputs, step, time, ftd,
                                  /*write_stats=*/true, /*is_final=*/false);
      }
    }
      // ---------------- Final write --------------------------
      // Only write final if it wasn't already written at plot_int
      const int final_step = step;
      bool should_write_final = (inputs.plot_int <= 0);
      if (inputs.plot_int > 0) {
          should_write_final = (final_step % inputs.plot_int != 0);
      }
      if (should_write_final) {
          write_wildfire_plotfile(f, ba, dm, geom, inputs, final_step, time, ftd,
                                  /*write_stats=*/true, /*is_final=*/true);
      }
      // ---- Write HTML fire report (end of run) ----
      write_fire_report_html(inputs, static_cast<double>(time), final_step);
      // ---- Write ASCII fire size summary (mirrors tools/fire_size_summary.py) ----
      write_fire_size_summary_text(inputs.fire_stats_file);
    }
  amrex::Finalize();
  return 0;
}
