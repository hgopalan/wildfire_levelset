// plot_results.H
// Consolidates all wildfire plotfile output into a single function so that
// main.cpp can call it for the initial, per-step, and final outputs without
// duplicating ~140 lines of setup, name-list, MultiFab::Copy, and stat-print
// code three times.
//
// All required headers are assumed to be available via wildfire_includes.H.
//
// Usage in main.cpp:
//   // initial (step 0) – no fire-stats or perimeter side-outputs:
//   write_wildfire_plotfile(f, ba, dm, geom, inputs, restart_step, time, ftd,
//                           /*write_stats=*/false, /*is_final=*/false);
//   // per-step:
//   write_wildfire_plotfile(f, ba, dm, geom, inputs, step, time, ftd,
//                           /*write_stats=*/true,  /*is_final=*/false);
//   // final:
//   write_wildfire_plotfile(f, ba, dm, geom, inputs, final_step, time, ftd,
//                           /*write_stats=*/true,  /*is_final=*/true);
//
// Plotfile variable filtering (inputs.plot_vars):
//   When inputs.plot_vars is non-empty only the listed variable names are
//   written.  Unknown names trigger a warning and are skipped.  If all names
//   are unknown the full plotfile is written as a fallback.
//   Set in the inputs file as a space-separated list, e.g.:
//     plot_vars = phi R arrival_time fireline_intensity flame_length
//   See parse_inputs.H for the complete list of valid variable name strings.

#ifndef PLOT_RESULTS_H
#define PLOT_RESULTS_H

#include "fire_intensity_class.H"
#include "byram_convective_number.H"
#include "flame_tilt_angle.H"
#include "packing_ratio.H"
#include "flame_front_depth.H"
#include "mcarthur_ffdi.H"
#include "chandler_burning_index.H"
#include <unordered_map>
// ============================================================================
// write_wildfire_plotfile
//
// Computes per-plotfile diagnostics (heat-per-unit-area, burnout time,
// vorticity, reaction intensity), assembles a composite plotmf containing
// all output fields, calls WriteSingleLevelPlotfile, prints burned-area and
// perimeter statistics, writes phi_negative / phi_envelope xy files, and
// (when write_stats is true) appends fire-statistics and perimeter CSV/GeoJSON.
//
// Parameters:
//   f           – WildfireFields struct containing all simulation MultiFabs
//   ba          – BoxArray for the computational domain
//   dm          – DistributionMapping for the computational domain
//   geom        – AMReX Geometry
//   inputs      – InputParameters (rothermel, farsite, output file names, etc.)
//   step        – current step number (used for file name)
//   time        – current simulation time [s]
//   ftd         – FuelTableData (GPU lookup tables for burnout/reaction intensity)
//   write_stats – when true, append fire-stats CSV and write perimeter files
//   is_final    – when true, print "Wrote final plt####" instead of "Wrote plt####"
// ============================================================================
inline void write_wildfire_plotfile(
    WildfireFields&        f,
    const BoxArray&        ba,
    const DistributionMapping& dm,
    const Geometry&        geom,
    InputParameters&       inputs,
    int                    step,
    Real                   time,
    const FuelTableData&   ftd,
    bool                   write_stats,
    bool                   is_final)
{
    // Convenience references
    MultiFab& phi                   = f.phi;
    MultiFab& vel                   = f.vel;
    MultiFab& farsite_spread        = f.farsite_spread;
    MultiFab& R_mf                  = f.R_mf;
    MultiFab& spotting_data         = f.spotting_data;
    MultiFab& fuel_consumption_mf   = f.fuel_consumption_mf;
    MultiFab& crown_fire_fraction_mf= f.crown_fire_fraction_mf;
    MultiFab& albini_data           = f.albini_data;
    MultiFab& elevation_mf          = f.elevation_mf;
    MultiFab& slope_mf              = f.slope_mf;
    MultiFab& aspect_mf             = f.aspect_mf;
    MultiFab& fuel_model_mf         = f.fuel_model_mf;
    MultiFab& fireline_intensity_mf = f.fireline_intensity_mf;
    MultiFab& flame_length_mf       = f.flame_length_mf;
    MultiFab& fire_intensity_class_mf = f.fire_intensity_class_mf;
    MultiFab& weise_data            = f.weise_data;
    MultiFab& viegas_data           = f.viegas_data;
    MultiFab& ecology_mf            = f.ecology_mf;
    MultiFab& emissions_mf          = f.emissions_mf;
    MultiFab& arrival_time_mf       = f.arrival_time_mf;
    MultiFab& heat_per_unit_area_mf = f.heat_per_unit_area_mf;
    MultiFab& vorticity_mf          = f.vorticity_mf;
    MultiFab& cbh_mf                = f.cbh_mf;
    MultiFab& cbd_mf                = f.cbd_mf;
    MultiFab& cc_mf                 = f.cc_mf;
    MultiFab& canopy_height_mf      = f.canopy_height_mf;
    MultiFab& burnout_time_mf       = f.burnout_time_mf;
    MultiFab& reaction_intensity_mf = f.reaction_intensity_mf;
    MultiFab& residual_fuel_mf      = f.residual_fuel_mf;
    MultiFab& shade_fraction_mf     = f.shade_fraction_mf;
    MultiFab& ti_full_mf            = f.ti_full_mf;
    MultiFab& ci_full_mf            = f.ci_full_mf;
    MultiFab& spatial_moisture_mf   = f.spatial_moisture_mf;
    MultiFab& plume_rise_mf         = f.plume_rise_mf;
    MultiFab& spread_dir_mf         = f.spread_dir_mf;
    MultiFab& spot_catch_prob_mf    = f.spot_catch_prob_mf;
    MultiFab& spotting_lineage_mf   = f.spotting_lineage_mf;
    MultiFab& ros_at_arrival_mf     = f.ros_at_arrival_mf;
    MultiFab& ember_cascade_mf      = f.ember_cascade_mf;
    MultiFab& fl_exceedance_mf      = f.fl_exceedance_mf;
    
    // New diagnostic fields for 10 feature integration
    MultiFab& crown_fraction_burned_mf = f.crown_fraction_burned_mf;
    MultiFab& effective_wind_speed_mf = f.effective_wind_speed_mf;
    MultiFab& burnout_phases_mf = f.burnout_phases_mf;
    MultiFab& residual_heat_release_mf = f.residual_heat_release_mf;
    MultiFab& time_since_burn_mf = f.time_since_burn_mf;
    MultiFab& temperature_mf        = f.temperature_mf;
    MultiFab& humidity_mf           = f.humidity_mf;

    const bool has_landscape = !inputs.rothermel.landscape_file.empty();

    // ---- Compute per-plotfile diagnostics ----
    {
        const RothermelComputed rc_plt = compute_rothermel_params(inputs.rothermel);
        compute_heat_per_unit_area(heat_per_unit_area_mf, phi, rc_plt.I_R,
                                   inputs.farsite.tau_residence);
    }
    compute_burnout_time(burnout_time_mf, arrival_time_mf,
                         inputs.farsite.tau_residence,
                         has_landscape ? &fuel_model_mf : nullptr,
                         ftd.d_tau_table_ptr, ftd.tau_table_size);
    compute_vorticity(vorticity_mf, vel, geom);
    compute_reaction_intensity(reaction_intensity_mf,
                               has_landscape ? &fuel_model_mf : nullptr,
                               inputs.rothermel,
                               ftd.d_fuel_table_ptr, ftd.fuel_table_size);
    
    // Compute fire intensity class (Byram 1959 6-class system)
    compute_fire_intensity_class(fire_intensity_class_mf, fireline_intensity_mf);

    // ---- Compute new fire behavior diagnostics ----
    // Create temporary residence time MultiFab for diagnostics
    MultiFab residence_time_temp(ba, dm, 1, 0);
    residence_time_temp.setVal(inputs.farsite.tau_residence);

    // Byram convective number
    compute_byram_convective_number(f.convective_number_mf,
                                    fireline_intensity_mf,
                                    residence_time_temp,
                                    vel);

    // Flame tilt angle
    compute_flame_tilt_angle(f.flame_tilt_mf,
                             fireline_intensity_mf,
                             vel);

    // Packing ratio diagnostics
    compute_packing_ratio(f.packing_ratio_mf,
                          f.relative_packing_mf,
                          inputs.rothermel);

    // Flame front depth
    compute_flame_front_depth(f.flame_depth_mf,
                              R_mf,
                              residence_time_temp);

    // McArthur FFDI
    // Default parameters: T=25°C, RH=50%, drought_factor=5.0
    // TODO: Add these to InputParameters when multi-station weather is fully integrated
    const Real default_temp_C = 25.0;
    const Real default_RH_pct = 50.0;
    const Real default_drought = 5.0;
    
    // Compute wind magnitude for FFDI
    Real wind_ms = 0.0;
    for (MFIter mfi(vel); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const u = vel.const_array(mfi);
        Real local_max = 0.0;
        ReduceOps<ReduceOpMax> reduce_op;
        ReduceData<Real> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;
        
        reduce_op.eval(bx, reduce_data,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple {
            const Real u_mag = std::sqrt(u(i,j,k,0)*u(i,j,k,0) + u(i,j,k,1)*u(i,j,k,1));
            return u_mag;
        });
        
        ReduceTuple hv = reduce_data.value();
        local_max = amrex::get<0>(hv);
        wind_ms = amrex::max(wind_ms, local_max);
    }
    ParallelDescriptor::ReduceRealMax(wind_ms);
    
    compute_mcarthur_ffdi_uniform(f.mcarthur_ffdi_mf,
                                   default_temp_C, default_RH_pct, wind_ms,
                                   default_drought);

    // Compute spatially-varying FFWI, SFWI, and FMI diagnostics
    for (MFIter mfi(f.ffwi_mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto ffwi_arr = f.ffwi_mf.array(mfi);
        auto sfwi_arr = f.sfwi_mf.array(mfi);
        auto fmi_arr  = f.fmi_mf.array(mfi);
        auto const T_arr = temperature_mf.const_array(mfi);
        auto const RH_arr = humidity_mf.const_array(mfi);
        auto const vel_arr = vel.const_array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const Real T = T_arr(i, j, k);
            const Real RH = RH_arr(i, j, k);
            const Real ux = vel_arr(i, j, k, 0);
            const Real uy = vel_arr(i, j, k, 1);
            const Real wind_speed_ms = std::sqrt(ux*ux + uy*uy);

            ffwi_arr(i, j, k, 0) = compute_fosberg_fire_weather_index(T, RH, wind_speed_ms);
            
            // Equilibrium moisture content
            Real RH_pct = amrex::max(Real(1.0), amrex::min(RH, Real(100.0)));
            Real EMC = Real(0.03229) + Real(0.281073) * RH_pct - Real(0.000578) * RH_pct * T;
            EMC = amrex::max(Real(1.0), amrex::min(EMC * Real(100.0), Real(50.0)));
            Real eta = Real(1.0) - Real(2.0) * (EMC / Real(30.0)) 
                       + Real(1.5) * std::pow(EMC / Real(30.0), 2.0) 
                       - Real(0.5) * std::pow(EMC / Real(30.0), 3.0);
            eta = amrex::max(Real(0.0), eta);
            ffwi_arr(i, j, k, 1) = EMC;
            ffwi_arr(i, j, k, 2) = eta;

            sfwi_arr(i, j, k) = compute_sharples_fire_weather_index(RH, wind_speed_ms);
            fmi_arr(i, j, k)  = compute_fuel_moisture_index(T, RH);
        });
    }

    // ---- Build output variable name list ----
    Vector<std::string> names = {"phi", "velx", "vely"
#if (AMREX_SPACEDIM == 3)
        , "velz", "farsite_dx", "farsite_dy", "farsite_dz", "R",
        "spot_prob", "spot_count", "spot_dist", "spot_active",
        "fuel_consumption", "crown_fraction",
        "albini_Hz", "albini_count", "albini_dist", "albini_active",
        "elevation", "slope", "aspect", "fuel_model",
        "fireline_intensity", "flame_length", "fire_intensity_class",
        "convective_number", "flame_tilt", "packing_ratio", "relative_packing",
        "flame_depth", "mcarthur_ffdi",
        "weise_flame_height", "weise_flame_tilt",
        "weise_whirl_height", "weise_whirl_radius",
        "weise_angular_velocity", "weise_max_tang_vel",
        "viegas_ROS", "viegas_eruptive_flag",
        "viegas_ROS_excess", "viegas_flame_tilt", "viegas_slope_factor",
        "scorch_height", "prob_ignition", "tree_mortality", "crown_activity",
        "torching_ratio", "crowning_ratio", "energy_release_component",
        "nfdrs_spread_component", "nfdrs_burning_index",
        "co2_emissions", "co_emissions", "pm25_emissions", "ch4_emissions", "n2o_emissions",
        "arrival_time", "heat_per_unit_area", "vorticity_z",
        "cbh", "cbd", "canopy_cover", "canopy_height", "burnout_time",
        "reaction_intensity", "residual_fuel", "shade_fraction",
        "torching_index_kmh", "crowning_index_kmh",
        "moisture_d1", "moisture_d10", "moisture_d100", "moisture_lh", "moisture_lw",
        "plume_rise_m",
        "spread_dir_x", "spread_dir_y",
        "spot_catch_prob",
        "spot_lineage_x", "spot_lineage_y", "albini_lineage_x", "albini_lineage_y",
        "ros_at_arrival",
        "ember_cascade_flux", "ember_cascade_ignition",
        "fl_exceedance",
        "canadian_fwi_ffmc", "canadian_fwi_dmc", "canadian_fwi_dc", "canadian_fwi_isi", "canadian_fwi_bui", "canadian_fwi_fwi",
        "ffwi", "ffwi_emc", "ffwi_eta",
        "sfwi", "fmi",
        // New Feature 3-10 diagnostics
        "crown_fraction_burned", "effective_wind_speed",
        "burnout_flaming_time", "burnout_smoldering_time",
        "residual_heat_release", "time_since_burn"
#else
        , "farsite_dx", "farsite_dy", "R",
        "spot_prob", "spot_count", "spot_dist", "spot_active",
        "fuel_consumption", "crown_fraction",
        "albini_Hz", "albini_count", "albini_dist", "albini_active",
        "elevation", "slope", "aspect", "fuel_model",
        "fireline_intensity", "flame_length", "fire_intensity_class",
        "convective_number", "flame_tilt", "packing_ratio", "relative_packing",
        "flame_depth", "mcarthur_ffdi",
        "weise_flame_height", "weise_flame_tilt",
        "weise_whirl_height", "weise_whirl_radius",
        "weise_angular_velocity", "weise_max_tang_vel",
        "viegas_ROS", "viegas_eruptive_flag",
        "viegas_ROS_excess", "viegas_flame_tilt", "viegas_slope_factor",
        "scorch_height", "prob_ignition", "tree_mortality", "crown_activity",
        "torching_ratio", "crowning_ratio", "energy_release_component",
        "nfdrs_spread_component", "nfdrs_burning_index",
        "co2_emissions", "co_emissions", "pm25_emissions", "ch4_emissions", "n2o_emissions",
        "arrival_time", "heat_per_unit_area", "vorticity_z",
        "cbh", "cbd", "canopy_cover", "canopy_height", "burnout_time",
        "reaction_intensity", "residual_fuel", "shade_fraction",
        "torching_index_kmh", "crowning_index_kmh",
        "moisture_d1", "moisture_d10", "moisture_d100", "moisture_lh", "moisture_lw",
        "plume_rise_m",
        "spread_dir_x", "spread_dir_y",
        "spot_catch_prob",
        "spot_lineage_x", "spot_lineage_y", "albini_lineage_x", "albini_lineage_y",
        "ros_at_arrival",
        "ember_cascade_flux", "ember_cascade_ignition",
        "fl_exceedance",
        "canadian_fwi_ffmc", "canadian_fwi_dmc", "canadian_fwi_dc", "canadian_fwi_isi", "canadian_fwi_bui", "canadian_fwi_fwi",
        "ffwi", "ffwi_emc", "ffwi_eta",
        "sfwi", "fmi",
        // New Feature 3-10 diagnostics
        "crown_fraction_burned", "effective_wind_speed",
        "burnout_flaming_time", "burnout_smoldering_time",
        "residual_heat_release", "time_since_burn"
#endif
    };

    // ---- Assemble composite plotfile MultiFab ----
    // Component layout (offsets are cumulative):
    //   0                                        : phi (1)
    //   1                                        : vel (AMREX_SPACEDIM)
    //   1+DIM                                    : farsite_spread (AMREX_SPACEDIM)
    //   1+2*DIM                                  : R (1)
    //   1+2*DIM+1                                : spotting_data (4)
    //   1+2*DIM+5                                : fuel_consumption (1)
    //   1+2*DIM+6                                : crown_fraction (1)
    //   1+2*DIM+7                                : albini_data (4)
    //   1+2*DIM+11                               : elevation, slope, aspect, fuel_model (1 each = 4)
    //   1+2*DIM+15                               : fireline_intensity, flame_length, fire_intensity_class (3)
    //   1+2*DIM+18                               : convective_number, flame_tilt, packing_ratio, relative_packing, flame_depth, mcarthur_ffdi (6)
    //   1+2*DIM+24                               : weise_data (WEISE_NCOMP)
    //   1+2*DIM+24+WEISE                         : viegas_data (VIEGAS_NCOMP)
    //   1+2*DIM+24+WEISE+VIEGAS                  : ecology_mf (ECOLOGY_NCOMP)
    //   1+2*DIM+24+WEISE+VIEGAS+ECOLOGY          : emissions_mf (EMISSIONS_NCOMP)
    //   ...+EMISSIONS                            : arrival_time (1)
    //   ...+EMISSIONS+1                          : heat_per_unit_area (1)
    //   ...+EMISSIONS+2                          : vorticity_z (1)
    //   ...+EMISSIONS+3 .. +6                    : cbh, cbd, canopy_cover, canopy_height (4)
    //   ...+EMISSIONS+7                          : burnout_time (1)
    //   ...+EMISSIONS+8                          : reaction_intensity (1)
    //   ...+EMISSIONS+9                          : residual_fuel (1)
    //   ...+EMISSIONS+10                         : shade_fraction (1)
    //   ...+EMISSIONS+11                         : torching_index (1)
    //   ...+EMISSIONS+12                         : crowning_index (1)
    //   ...+EMISSIONS+13                         : spatial_moisture (FMS_MOISTURE_NCOMP)
    //   ...+EMISSIONS+13+FMS_MOISTURE_NCOMP      : plume_rise_m (1)
    //   ...+EMISSIONS+14+FMS_MOISTURE_NCOMP      : spread_dir (2)
    //   ...+EMISSIONS+16+FMS_MOISTURE_NCOMP      : spot_catch_prob (1)
    //   ...+EMISSIONS+17+FMS_MOISTURE_NCOMP      : spotting_lineage (SPOT_LINEAGE_NCOMP)
    //   ...+EMISSIONS+17+FMS_MOISTURE_NCOMP+SPOT_LINEAGE_NCOMP : ros_at_arrival (1)
    //   ...+EMISSIONS+18+FMS_MOISTURE_NCOMP+SPOT_LINEAGE_NCOMP : ember_cascade (2)
    //   ...+EMISSIONS+20+FMS_MOISTURE_NCOMP+SPOT_LINEAGE_NCOMP : fl_exceedance (1)
    constexpr int D = AMREX_SPACEDIM;
    const int ncomp = 1 + D + D + 1 + 4 + 1 + 1 + 4 + 4 + 3 + 6  // +6 for new diagnostics
                    + WEISE_NCOMP + VIEGAS_NCOMP + ECOLOGY_NCOMP + EMISSIONS_NCOMP
                    + 13 + FMS_MOISTURE_NCOMP + 1  // +1 for plume_rise_m
                    + 2 + 1 + SPOT_LINEAGE_NCOMP   // spread_dir(2), spot_catch_prob(1), lineage(4)
                    + 1                             // ros_at_arrival (1)
                    + 2                             // ember_cascade (flux + ignition flag)
                    + 1                             // fl_exceedance (1)
                    + 6 + 3 + 1 + 1                 // Canadian FWI(6), FFWI(3), SFWI(1), FMI(1)
                    + 1 + 1 + 2 + 1 + 1;            // new Features: cfb(1), eff_wind(1), burnout_phases(2), residual_heat(1), time_since_burn(1)
    MultiFab plotmf(ba, dm, ncomp, 0);

    int c = 0;
    MultiFab::Copy(plotmf, phi,                   0, c, 1,           0); c += 1;
    MultiFab::Copy(plotmf, vel,                   0, c, D,           0); c += D;
    MultiFab::Copy(plotmf, farsite_spread,        0, c, D,           0); c += D;
    MultiFab::Copy(plotmf, R_mf,                  0, c, 1,           0); c += 1;
    MultiFab::Copy(plotmf, spotting_data,         0, c, 4,           0); c += 4;
    MultiFab::Copy(plotmf, fuel_consumption_mf,   0, c, 1,           0); c += 1;
    MultiFab::Copy(plotmf, crown_fire_fraction_mf,0, c, 1,           0); c += 1;
    MultiFab::Copy(plotmf, albini_data,           0, c, 4,           0); c += 4;
    MultiFab::Copy(plotmf, elevation_mf,          0, c, 1,           0); c += 1;
    MultiFab::Copy(plotmf, slope_mf,              0, c, 1,           0); c += 1;
    MultiFab::Copy(plotmf, aspect_mf,             0, c, 1,           0); c += 1;
    MultiFab::Copy(plotmf, fuel_model_mf,         0, c, 1,           0); c += 1;
    MultiFab::Copy(plotmf, fireline_intensity_mf, 0, c, 1,           0); c += 1;
    MultiFab::Copy(plotmf, flame_length_mf,       0, c, 1,           0); c += 1;
    MultiFab::Copy(plotmf, fire_intensity_class_mf, 0, c, 1,         0); c += 1;
    MultiFab::Copy(plotmf, f.convective_number_mf, 0, c, 1,          0); c += 1;
    MultiFab::Copy(plotmf, f.flame_tilt_mf,       0, c, 1,           0); c += 1;
    MultiFab::Copy(plotmf, f.packing_ratio_mf,    0, c, 1,           0); c += 1;
    MultiFab::Copy(plotmf, f.relative_packing_mf, 0, c, 1,           0); c += 1;
    MultiFab::Copy(plotmf, f.flame_depth_mf,      0, c, 1,           0); c += 1;
    MultiFab::Copy(plotmf, f.mcarthur_ffdi_mf,    0, c, 1,           0); c += 1;
    MultiFab::Copy(plotmf, weise_data,            0, c, WEISE_NCOMP,    0); c += WEISE_NCOMP;
    MultiFab::Copy(plotmf, viegas_data,           0, c, VIEGAS_NCOMP,   0); c += VIEGAS_NCOMP;
    MultiFab::Copy(plotmf, ecology_mf,            0, c, ECOLOGY_NCOMP,  0); c += ECOLOGY_NCOMP;
    MultiFab::Copy(plotmf, emissions_mf,          0, c, EMISSIONS_NCOMP,0); c += EMISSIONS_NCOMP;
    MultiFab::Copy(plotmf, arrival_time_mf,       0, c, 1,           0); c += 1;
    MultiFab::Copy(plotmf, heat_per_unit_area_mf, 0, c, 1,           0); c += 1;
    MultiFab::Copy(plotmf, vorticity_mf,          0, c, 1,           0); c += 1;
    MultiFab::Copy(plotmf, cbh_mf,                0, c, 1,           0); c += 1;
    MultiFab::Copy(plotmf, cbd_mf,                0, c, 1,           0); c += 1;
    MultiFab::Copy(plotmf, cc_mf,                 0, c, 1,           0); c += 1;
    MultiFab::Copy(plotmf, canopy_height_mf,      0, c, 1,           0); c += 1;
    MultiFab::Copy(plotmf, burnout_time_mf,       0, c, 1,           0); c += 1;
    MultiFab::Copy(plotmf, reaction_intensity_mf, 0, c, 1,           0); c += 1;
    MultiFab::Copy(plotmf, residual_fuel_mf,      0, c, 1,           0); c += 1;
    MultiFab::Copy(plotmf, shade_fraction_mf,     0, c, 1,           0); c += 1;
    MultiFab::Copy(plotmf, ti_full_mf,            0, c, 1,           0); c += 1;
    MultiFab::Copy(plotmf, ci_full_mf,            0, c, 1,           0); c += 1;
    MultiFab::Copy(plotmf, spatial_moisture_mf,   0, c, FMS_MOISTURE_NCOMP, 0); c += FMS_MOISTURE_NCOMP;
    MultiFab::Copy(plotmf, plume_rise_mf,         0, c, 1,           0); c += 1;
    MultiFab::Copy(plotmf, spread_dir_mf,         0, c, 2,           0); c += 2;
    MultiFab::Copy(plotmf, spot_catch_prob_mf,    0, c, 1,           0); c += 1;
    MultiFab::Copy(plotmf, spotting_lineage_mf,   0, c, SPOT_LINEAGE_NCOMP, 0); c += SPOT_LINEAGE_NCOMP;
    MultiFab::Copy(plotmf, ros_at_arrival_mf,     0, c, 1,           0); c += 1;
    MultiFab::Copy(plotmf, ember_cascade_mf,      0, c, 2,           0); c += 2;
    MultiFab::Copy(plotmf, fl_exceedance_mf,      0, c, 1,           0); c += 1;
    // FWI, Fosberg, Sharples, and FMI
    MultiFab::Copy(plotmf, f.fwi_mf,               0, c, 6,           0); c += 6;
    MultiFab::Copy(plotmf, f.ffwi_mf,              0, c, 3,           0); c += 3;
    MultiFab::Copy(plotmf, f.sfwi_mf,              0, c, 1,           0); c += 1;
    MultiFab::Copy(plotmf, f.fmi_mf,               0, c, 1,           0); c += 1;
    // New Feature diagnostics
    MultiFab::Copy(plotmf, crown_fraction_burned_mf, 0, c, 1,        0); c += 1;
    MultiFab::Copy(plotmf, effective_wind_speed_mf,  0, c, 1,        0); c += 1;
    MultiFab::Copy(plotmf, burnout_phases_mf,       0, c, 2,        0); c += 2;
    MultiFab::Copy(plotmf, residual_heat_release_mf, 0, c, 1,        0); c += 1;
    MultiFab::Copy(plotmf, time_since_burn_mf,      0, c, 1,        0);

    // ---- Write plotfile ----
    char buf[64];
    std::snprintf(buf, sizeof(buf), "plt%04d", step);

    if (inputs.plot_vars.empty()) {
        // Default: write all variables
        WriteSingleLevelPlotfile(buf, plotmf, names, geom, time, step);
    } else {
        // Filtered: write only the user-requested subset of variables.
        // Build a name→component-index map for O(1) lookup.
        std::unordered_map<std::string,int> name_to_comp;
        name_to_comp.reserve(names.size());
        for (int i = 0; i < (int)names.size(); ++i)
            name_to_comp[names[i]] = i;

        Vector<std::string> filtered_names;
        Vector<int>         filtered_comps;
        for (const auto& req : inputs.plot_vars) {
            auto it = name_to_comp.find(req);
            if (it != name_to_comp.end()) {
                filtered_names.push_back(req);
                filtered_comps.push_back(it->second);
            } else {
                amrex::Print() << "WARNING: plot_vars: unknown variable '" << req
                               << "' – skipped\n";
            }
        }
        if (filtered_names.empty()) {
            amrex::Print() << "WARNING: plot_vars list produced no valid variables; "
                              "writing full plotfile instead.\n";
            WriteSingleLevelPlotfile(buf, plotmf, names, geom, time, step);
        } else {
            MultiFab filtered_mf(ba, dm, (int)filtered_names.size(), 0);
            for (int i = 0; i < (int)filtered_comps.size(); ++i)
                MultiFab::Copy(filtered_mf, plotmf, filtered_comps[i], i, 1, 0);
            WriteSingleLevelPlotfile(buf, filtered_mf, filtered_names, geom, time, step);
        }
    }
    if (write_stats) {
        if (is_final)
            amrex::Print() << "Wrote final " << buf << "\n";
        else
            amrex::Print() << "Wrote " << buf << "\n";
    }

    // ---- Burned area and perimeter statistics ----
    // We use arrival_time_mf >= 0 (set by init_phi_source for the initial
    // ignition and updated each step as new cells burn) rather than phi < 0.
    // In the FARSITE propagation the phi field is reset to 0 after every step
    // so that phi < 0 only marks the current thin fire-front ring; the burned
    // interior is phi = 0 and would be missed.  arrival_time accumulates every
    // cell that has ever been ignited and is correct for all propagation modes.
    //
    // Perimeter uses the edge-crossing method on the arrival_time boundary:
    // every face separating a burned cell (at >= 0) from an unburned cell
    // (at < 0) contributes the face length to the perimeter total.
    {
        amrex::Long n_burned = 0;
        Real perim_m = Real(0.0);
        const auto dxx = geom.CellSize();
        for (MFIter mfi(arrival_time_mf); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            auto const at = arrival_time_mf.const_array(mfi);
            const IntVect bxhi = bx.bigEnd();
            amrex::LoopOnCpu(bx, [&](int i, int j, int k) {
                bool burned = (at(i, j, k, 0) >= Real(0.0));
                if (burned) ++n_burned;
                // x-face between (i,j) and (i+1,j): face length = dxx[1]
                if (i < bxhi[0]) {
                    bool burned_n = (at(i+1, j, k, 0) >= Real(0.0));
                    if (burned != burned_n) perim_m += dxx[1];
                }
                // y-face between (i,j) and (i,j+1): face length = dxx[0]
                if (j < bxhi[1]) {
                    bool burned_n = (at(i, j+1, k, 0) >= Real(0.0));
                    if (burned != burned_n) perim_m += dxx[0];
                }
            });
        }
        ParallelDescriptor::ReduceLongSum(n_burned);
        ParallelDescriptor::ReduceRealSum(perim_m);
        Real cell_area = dxx[0] * dxx[1];
        amrex::Print() << "  Burned area: " << Real(n_burned)*cell_area/1.0e4 << " ha"
                       << "  Perimeter: "   << perim_m/Real(1000.0) << " km\n";
    }

    // ---- Write negative-phi x-y data files ----
    {
        char xy_buf[64];
        std::snprintf(xy_buf, sizeof(xy_buf), "phi_negative_%04d.dat", step);
        write_negative_phi_xy(phi, geom, xy_buf);
        std::snprintf(xy_buf, sizeof(xy_buf), "phi_envelope_%04d.dat", step);
        write_negative_phi_convex_hull(phi, geom, xy_buf);
    }

    // ---- Optional side outputs (fire stats + perimeter files) ----
    if (write_stats) {
        if (!inputs.fire_stats_file.empty())
            append_fire_stats(phi, geom, &emissions_mf, step, time,
                              inputs.fire_stats_file, arrival_time_mf, &R_mf, &flame_length_mf);
        if (!inputs.fsa_file.empty())
            append_fire_spread_atlas(phi, geom, step, time, inputs.fsa_file);
        if (!inputs.pst_file.empty())
            append_pst_stats(phi, geom, arrival_time_mf,
                             &fireline_intensity_mf, &flame_length_mf,
                             &spotting_data,
                             step, time, inputs.pst_file);
        if (inputs.write_perimeter_csv == 1) {
            char csv_buf[64];
            std::snprintf(csv_buf, sizeof(csv_buf), "perimeter_%04d.csv", step);
            write_fire_perimeter_csv(phi, geom, csv_buf);
        }
        if (inputs.write_perimeter_geojson == 1) {
            char gjson_buf[64];
            std::snprintf(gjson_buf, sizeof(gjson_buf), "perimeter_%04d.geojson", step);
            write_fire_perimeter_geojson(phi, geom, gjson_buf, step, time);
        }
        if (inputs.write_perimeter_kml == 1) {
            char kml_buf[64];
            std::snprintf(kml_buf, sizeof(kml_buf), "perimeter_%04d.kml", step);
            write_fire_perimeter_kml(phi, geom, kml_buf, step, time,
                                     inputs.kml_utm_zone, inputs.kml_utm_northern);
        }
    }
}

#endif // PLOT_RESULTS_H