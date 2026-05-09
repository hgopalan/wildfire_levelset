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

    // ---- Build output variable name list ----
    Vector<std::string> names = {"phi", "velx", "vely"
#if (AMREX_SPACEDIM == 3)
        , "velz", "farsite_dx", "farsite_dy", "farsite_dz", "R",
        "spot_prob", "spot_count", "spot_dist", "spot_active",
        "fuel_consumption", "crown_fraction",
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
        "nfdrs_spread_component", "nfdrs_burning_index",
        "co2_emissions", "co_emissions", "pm25_emissions",
        "arrival_time", "heat_per_unit_area", "vorticity_z",
        "cbh", "cbd", "canopy_cover", "canopy_height", "burnout_time",
        "reaction_intensity", "residual_fuel", "shade_fraction",
        "torching_index_kmh", "crowning_index_kmh",
        "moisture_d1", "moisture_d10", "moisture_d100", "moisture_lh", "moisture_lw"
#else
        , "farsite_dx", "farsite_dy", "R",
        "spot_prob", "spot_count", "spot_dist", "spot_active",
        "fuel_consumption", "crown_fraction",
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
        "nfdrs_spread_component", "nfdrs_burning_index",
        "co2_emissions", "co_emissions", "pm25_emissions",
        "arrival_time", "heat_per_unit_area", "vorticity_z",
        "cbh", "cbd", "canopy_cover", "canopy_height", "burnout_time",
        "reaction_intensity", "residual_fuel", "shade_fraction",
        "torching_index_kmh", "crowning_index_kmh",
        "moisture_d1", "moisture_d10", "moisture_d100", "moisture_lh", "moisture_lw"
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
    //   1+2*DIM+15                               : fireline_intensity, flame_length (2)
    //   1+2*DIM+17                               : weise_data (WEISE_NCOMP)
    //   1+2*DIM+17+WEISE                         : viegas_data (VIEGAS_NCOMP)
    //   1+2*DIM+17+WEISE+VIEGAS                  : ecology_mf (ECOLOGY_NCOMP)
    //   1+2*DIM+17+WEISE+VIEGAS+ECOLOGY          : emissions_mf (EMISSIONS_NCOMP)
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
    constexpr int D = AMREX_SPACEDIM;
    const int ncomp = 1 + D + D + 1 + 4 + 1 + 1 + 4 + 1 + 5
                    + WEISE_NCOMP + VIEGAS_NCOMP + ECOLOGY_NCOMP + EMISSIONS_NCOMP
                    + 13 + FMS_MOISTURE_NCOMP;
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
    MultiFab::Copy(plotmf, spatial_moisture_mf,   0, c, FMS_MOISTURE_NCOMP, 0);

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
    // Perimeter is computed using the FARSITE/FlamMap edge-crossing method:
    // for every cell face that separates a burned cell (phi < 0) from an
    // unburned cell (phi >= 0), the face length is added to the perimeter.
    // This correctly accounts for corner cells and diagonal fire fronts,
    // unlike the simpler "perimeter-cell count × sqrt(cell_area)" approach.
    {
        amrex::Long n_burned = 0;
        Real perim_m = Real(0.0);
        const auto dxx = geom.CellSize();
        for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            auto const p = phi.const_array(mfi);
            const IntVect bxhi = bx.bigEnd();
            amrex::LoopOnCpu(bx, [&](int i, int j, int k) {
                Real pv = p(i, j, k, 0);
                if (pv < Real(0.0)) ++n_burned;
                // x-face between (i,j) and (i+1,j): face length = dxx[1]
                if (i < bxhi[0]) {
                    Real pn = p(i+1, j, k, 0);
                    if ((pv < Real(0.0)) != (pn < Real(0.0))) perim_m += dxx[1];
                }
                // y-face between (i,j) and (i,j+1): face length = dxx[0]
                if (j < bxhi[1]) {
                    Real pn = p(i, j+1, k, 0);
                    if ((pv < Real(0.0)) != (pn < Real(0.0))) perim_m += dxx[0];
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
                              inputs.fire_stats_file, &R_mf);
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
    }
}

#endif // PLOT_RESULTS_H