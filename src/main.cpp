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

    // Fire ecology diagnostics (ECOLOGY_NCOMP = 4 components):
    //   0 – scorch_height   [m]  (Van Wagner 1973)
    //   1 – prob_ignition   [-]  (Anderson 1970 / Rothermel 1983)
    //   2 – tree_mortality  [-]  (Ryan-Reinhardt 1988 style logistic)
    //   3 – crown_activity  [-]  (0=surface, 1=passive, 2=active crown fire)
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

    // Heat per unit area [BTU/ft²]: I_R × residence_time [min] for burned cells.
    // Zero for unburned cells. Recomputed before each plotfile write.
    MultiFab heat_per_unit_area_mf(ba, dm, 1, 0);
    heat_per_unit_area_mf.setVal(0.0);

    // Vertical vorticity ω_z = ∂v/∂x − ∂u/∂y [s⁻¹].
    // Recomputed from the current velocity field before each plotfile write.
    MultiFab vorticity_mf(ba, dm, 1, 0);
    vorticity_mf.setVal(0.0);

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
    MultiFab viegas_data(ba, dm, VIEGAS_NCOMP, 0);
    viegas_data.setVal(0.0);

    // Heat flux field [W/m²] (1 component, no ghost cells).
    // Initialized from heat_flux.value (uniform) or heat_flux.file (spatially varying).
    MultiFab heat_flux_mf(ba, dm, 1, 0);
    heat_flux_mf.setVal(0.0);


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
                 static_cast<float>(inputs.rothermel.M_lh),
                 static_cast<float>(inputs.rothermel.M_lw)});
            moisture_updated = true;
        }

        if (!moisture_updated) return;

        inputs.rothermel.M_d1   = static_cast<amrex::Real>(m.M_d1);
        inputs.rothermel.M_d10  = static_cast<amrex::Real>(m.M_d10);
        inputs.rothermel.M_d100 = static_cast<amrex::Real>(m.M_d100);
        inputs.rothermel.M_lh   = static_cast<amrex::Real>(m.M_lh);
        inputs.rothermel.M_lw   = static_cast<amrex::Real>(m.M_lw);
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
         inputs.fire_spread_model  == "cruz_crown");

    // Helper flag for Viegas+Balbi coupling
    const bool use_balbi_for_viegas = (inputs.fire_spread_model == "balbi");

    // ---------------- dt from CFL --------------------------
    Real dt=10;
    const bool use_levelset = (inputs.propagation_method == "levelset");
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
        } else {
            // Compute Rothermel wind speed R
            compute_rothermel_R(R_mf, vel_for_model, geom, inputs.rothermel,
                                 terrain_slopes.get(),
                                 !inputs.rothermel.landscape_file.empty() ? &fuel_model_mf : nullptr,
                                 d_fuel_table_ptr, fuel_table_size,
                                 has_spatial_crown ? &cc_mf : nullptr);
        }
        dt = compute_dt(R_mf, geom, inputs.cfl);
        amrex::Print() << "Computed dt = " << dt << "\n";
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
    // Fire emissions (CO₂, CO, PM₂.₅) from fuel load × consumption fraction
    compute_fire_emissions(emissions_mf, phi, fuel_consumption_mf,
                           inputs.rothermel, inputs.emissions);


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
      compute_burnout_time(burnout_time_mf, arrival_time_mf, inputs.farsite.tau_residence);
      compute_vorticity(vorticity_mf, vel, geom);
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
				   "co2_emissions", "co_emissions", "pm25_emissions",
				   "arrival_time", "heat_per_unit_area", "vorticity_z",
				   "cbh", "cbd", "canopy_cover", "canopy_height", "burnout_time"
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
				   "co2_emissions", "co_emissions", "pm25_emissions",
				   "arrival_time", "heat_per_unit_area", "vorticity_z",
				   "cbh", "cbd", "canopy_cover", "canopy_height", "burnout_time"
#endif
      };
      MultiFab plotmf(ba, dm, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 1 + 5 + WEISE_NCOMP + VIEGAS_NCOMP + ECOLOGY_NCOMP + EMISSIONS_NCOMP + 8, 0);
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
      {
        char buf[64];
        std::snprintf(buf, sizeof(buf), "plt%04d", restart_step);
        WriteSingleLevelPlotfile(buf, plotmf, names, geom, time, restart_step);
      }
      // Print burned area and perimeter statistics
      {
        long n_burned = 0, n_perim = 0;
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
      } else {
          compute_rothermel_R(R_mf, vel_for_model, geom, inputs.rothermel,
                               terrain_slopes.get(),
                               !inputs.rothermel.landscape_file.empty() ? &fuel_model_mf : nullptr,
                               d_fuel_table_ptr, fuel_table_size,
                               has_spatial_crown ? &cc_mf : nullptr);
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
      // Fire emissions (CO₂, CO, PM₂.₅)
      compute_fire_emissions(emissions_mf, phi, fuel_consumption_mf,
                             inputs.rothermel, inputs.emissions);
      if (use_levelset) {
	// Traditional level set advection.
	// Pass pre-computed R_mf when a wind-terrain model or non-Rothermel spread
	// model is active (see use_precomp_R_for_advection defined above).
	advect_levelset_weno5z_rk3(phi, vel, geom, dt_step, inputs.rothermel,
                                   terrain_slopes.get(),
                                   use_precomp_R_for_advection ? &R_mf : nullptr);
	dt = compute_dt(R_mf, geom, inputs.cfl);
      } else {
	// --- Step 3: FARSITE elliptical wavelet propagation (Richards 1990)
	// --- Step 4: Merge to new perimeter
	// For wind-terrain models, pass vel_for_model so FARSITE ellipse
	// orientation and ROS also reflect the terrain-corrected wind.
	compute_farsite_spread(phi, vel_for_model, farsite_spread, geom, dt_step, inputs.rothermel, inputs.farsite, inputs.crown, terrain_slopes.get(), &fuel_consumption_mf, &crown_fire_fraction_mf, has_spatial_crown ? &cc_mf : nullptr);
	
	// --- Step 5: Apply crown/spotting sub-models
	if (inputs.spotting.enable == 1 && (step % inputs.spotting.check_interval == 0)) {
	  compute_spotting_probability(spotting_data, phi, vel, geom, inputs.rothermel, inputs.spotting, terrain_slopes.get());
	  generate_firebrand_spots(phi, spotting_data, vel, geom, inputs.spotting, step);
	}
	// Albini (1983) firebrand spotting with 2-D trajectory integration
	if (inputs.albini_spotting.enable == 1 && (step % inputs.albini_spotting.check_interval == 0)) {
	  compute_albini_spotting(phi, albini_data, vel, R_mf, geom,
	                          inputs.rothermel, inputs.albini_spotting, step);
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

      // --- Dynamic fire points: check if a new ignition file has appeared
      if (!inputs.dynamic_fire_points_file.empty()) {
        apply_dynamic_fire_points(phi, geom,
                                  inputs.dynamic_fire_points_file,
                                  inputs.fire_gaussian_sigma);
      }

      // --- Write checkpoint if requested
      if (inputs.chk_int > 0 && ((step - restart_step) % inputs.chk_int == 0)) {
        char chk_buf[64];
        std::snprintf(chk_buf, sizeof(chk_buf), "chk%04d", step);
        write_checkpoint(chk_buf, phi, geom, step, time);
      }

      // --- Step 7: Update states, record outputs, step time
      if (inputs.plot_int > 0 && (step % inputs.plot_int == 0)) {
	// Compute per-plotfile diagnostics
	{
	  const RothermelComputed rc_plt = compute_rothermel_params(inputs.rothermel);
	  compute_heat_per_unit_area(heat_per_unit_area_mf, phi, rc_plt.I_R,
	                             inputs.farsite.tau_residence);
	}
	compute_burnout_time(burnout_time_mf, arrival_time_mf, inputs.farsite.tau_residence);
	compute_vorticity(vorticity_mf, vel, geom);
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
				     "co2_emissions", "co_emissions", "pm25_emissions",
				   "arrival_time", "heat_per_unit_area", "vorticity_z",
				   "cbh", "cbd", "canopy_cover", "canopy_height", "burnout_time"
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
				     "co2_emissions", "co_emissions", "pm25_emissions",
				   "arrival_time", "heat_per_unit_area", "vorticity_z",
				   "cbh", "cbd", "canopy_cover", "canopy_height", "burnout_time"
#endif
	};
	MultiFab plotmf(ba, dm, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 1 + 5 + WEISE_NCOMP + VIEGAS_NCOMP + ECOLOGY_NCOMP + EMISSIONS_NCOMP + 8, 0);
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
	WriteSingleLevelPlotfile(buf, plotmf, names, geom, time, step);
	amrex::Print() << "Wrote " << buf << "\n";
	{
	  long n_burned = 0, n_perim = 0;
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
	if (!inputs.fire_stats_file.empty()) append_fire_stats(phi, geom, &emissions_mf, step, time, inputs.fire_stats_file);
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
	compute_burnout_time(burnout_time_mf, arrival_time_mf, inputs.farsite.tau_residence);
	compute_vorticity(vorticity_mf, vel, geom);
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
				     "co2_emissions", "co_emissions", "pm25_emissions",
				   "arrival_time", "heat_per_unit_area", "vorticity_z",
				   "cbh", "cbd", "canopy_cover", "canopy_height", "burnout_time"
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
				     "co2_emissions", "co_emissions", "pm25_emissions",
				   "arrival_time", "heat_per_unit_area", "vorticity_z",
				   "cbh", "cbd", "canopy_cover", "canopy_height", "burnout_time"
#endif
	};
	MultiFab plotmf(ba, dm, 1 + AMREX_SPACEDIM + AMREX_SPACEDIM + 1 + 4 + 1 + 1 + 4 + 1 + 5 + WEISE_NCOMP + VIEGAS_NCOMP + ECOLOGY_NCOMP + EMISSIONS_NCOMP + 8, 0);
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
	WriteSingleLevelPlotfile(buf, plotmf, names, geom, time, final_step);
	amrex::Print() << "Wrote final " << buf << "\n";
	{
	  long n_burned = 0, n_perim = 0;
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
	if (!inputs.fire_stats_file.empty()) append_fire_stats(phi, geom, &emissions_mf, final_step, time, inputs.fire_stats_file);
	if (inputs.write_perimeter_csv == 1) { char csv_buf[64]; std::snprintf(csv_buf, sizeof(csv_buf), "perimeter_%04d.csv", final_step); write_fire_perimeter_csv(phi, geom, csv_buf); }
	if (inputs.write_perimeter_geojson == 1) { char gjson_buf[64]; std::snprintf(gjson_buf, sizeof(gjson_buf), "perimeter_%04d.geojson", final_step); write_fire_perimeter_geojson(phi, geom, gjson_buf, final_step, time); }
      }
    }
  amrex::Finalize();
  return 0;
}
