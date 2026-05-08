// terrain_setup.H
// Populates terrain, landscape, crown-layer, and spatial-moisture MultiFabs
// from terrain/LCP files.  Called once at startup after the WildfireFields
// struct is constructed.
//
// Priority rules:
//   elevation / slope / aspect : terrain_file > landscape_file
//   fuel model                  : landscape_file only (terrain file has none)
//   crown layers (CBH/CBD/CC)   : landscape LCP when use_spatial_crown = 1
//   spatial moisture            : FMS file when fms_file is set
#ifndef TERRAIN_SETUP_H
#define TERRAIN_SETUP_H

// All required headers provided by wildfire_includes.H (the umbrella header).
// terrain_setup.H additionally requires lcp_binary_reader.H (has include guard).
#include "lcp_binary_reader.H"

using namespace amrex;

// Populate all terrain / landscape / crown / moisture MultiFabs.
//
// Outputs:
//   terrain_slopes     – 2-comp (slope_x, slope_y) MultiFab; allocated here
//                        when a terrain or landscape file is present.
//   elevation_mf       – ground elevation [m]
//   slope_mf           – terrain slope [degrees]
//   aspect_mf          – terrain aspect [degrees]
//   fuel_model_mf      – Scott & Burgan fuel code per cell
//   cbh_mf             – crown base height [m]   (seeded from inputs.crown.CBH)
//   cbd_mf             – crown bulk density [kg/m³] (seeded from inputs.crown.CBD)
//   cc_mf              – canopy cover [%]
//   canopy_height_mf   – stand height [m]
//   has_spatial_crown  – true when per-cell crown data were loaded from LCP
//   spatial_moisture_mf– per-cell dead/live moisture (FMS_MOISTURE_NCOMP comps)
//   has_spatial_moisture– true when FMS moisture data were loaded
inline void setup_terrain(
    std::unique_ptr<MultiFab>& terrain_slopes,
    MultiFab& elevation_mf,
    MultiFab& slope_mf,
    MultiFab& aspect_mf,
    MultiFab& fuel_model_mf,
    MultiFab& cbh_mf,
    MultiFab& cbd_mf,
    MultiFab& cc_mf,
    MultiFab& canopy_height_mf,
    bool& has_spatial_crown,
    MultiFab& spatial_moisture_mf,
    bool& has_spatial_moisture,
    const BoxArray& ba,
    const DistributionMapping& dm,
    const Geometry& geom,
    const InputParameters& inputs)
{
    // ---- Terrain slopes (slope_x, slope_y) ----
    if (!inputs.rothermel.terrain_file.empty()) {
        terrain_slopes = std::make_unique<MultiFab>(ba, dm, 2, 0);
        compute_slopes_from_terrain(*terrain_slopes, geom,
                                    inputs.rothermel.terrain_file);
        amrex::Print() << "Initialized terrain slopes from terrain file: "
                       << inputs.rothermel.terrain_file << "\n";
        if (!inputs.rothermel.landscape_file.empty())
            amrex::Print() << "NOTE: terrain_file takes precedence for "
                              "elevation/slope/aspect; landscape_file used "
                              "for fuel model only\n";
    } else if (!inputs.rothermel.landscape_file.empty()) {
        terrain_slopes = std::make_unique<MultiFab>(ba, dm, 2, 0);
        compute_slopes_from_landscape(*terrain_slopes, geom,
                                      inputs.rothermel.landscape_file);
        amrex::Print() << "Initialized terrain slopes from landscape file: "
                       << inputs.rothermel.landscape_file << "\n";
    }

    // ---- Elevation ----
    if (!inputs.rothermel.terrain_file.empty())
        compute_elevation_from_terrain(elevation_mf, geom,
                                       inputs.rothermel.terrain_file);
    else if (!inputs.rothermel.landscape_file.empty())
        compute_elevation_from_landscape(elevation_mf, geom,
                                         inputs.rothermel.landscape_file);

    // ---- Slope, aspect, and fuel model ----
    if (!inputs.rothermel.terrain_file.empty() && terrain_slopes) {
        compute_slope_aspect_from_slopes(*terrain_slopes, slope_mf, aspect_mf);
        if (!inputs.rothermel.landscape_file.empty())
            compute_fuel_model_from_landscape(fuel_model_mf, geom,
                                              inputs.rothermel.landscape_file);
    } else if (!inputs.rothermel.landscape_file.empty()) {
        compute_landscape_slope_aspect_fuel(slope_mf, aspect_mf, fuel_model_mf,
                                            geom, inputs.rothermel.landscape_file);
    }

    // ---- Crown spatial layers from binary LCP ----
    has_spatial_crown = false;
    if (inputs.use_spatial_crown == 1 &&
        !inputs.rothermel.landscape_file.empty()) {
        has_spatial_crown = compute_crown_layers_from_lcp(
            &cbh_mf, &cbd_mf, &cc_mf,
            geom, inputs.rothermel.landscape_file,
            &canopy_height_mf);
    }

    // ---- Per-cell spatial moisture from FMS file ----
    has_spatial_moisture = false;
    if (!inputs.fms_file.empty() &&
        !inputs.rothermel.landscape_file.empty()) {
        load_fms_spatial_moisture(spatial_moisture_mf, fuel_model_mf,
                                  geom, inputs.fms_file, inputs.rothermel);
        has_spatial_moisture = true;
    }
}

#endif // TERRAIN_SETUP_H
