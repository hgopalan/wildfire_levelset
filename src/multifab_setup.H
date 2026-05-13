// multifab_setup.H
// Declares WildfireFields – a plain struct that owns every AMReX MultiFab
// used by the wildfire level-set solver.  Constructing one object allocates
// and zero-initialises (or sentinel-initialises) all fields in a single call,
// replacing the ~130-line block of declarations and setVal() calls in main.cpp.
//
// Usage in main.cpp:
//   WildfireFields f(ba, dm, inputs);
//   // Convenience references so the rest of main.cpp is unchanged:
//   MultiFab& phi = f.phi;
//   MultiFab& vel = f.vel;
//   ...
//
// The terrain_slopes unique_ptr starts as nullptr; main.cpp populates it
// conditionally after construction (same logic as before).
#ifndef MULTIFAB_SETUP_H
#define MULTIFAB_SETUP_H

#include <AMReX_MultiFab.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <memory>
#include "parse_inputs.H"
#include "fire_ecology_model.H"    // ECOLOGY_NCOMP
#include "fire_emissions_model.H"  // EMISSIONS_NCOMP
#include "weise_biging_whirl.H"    // WEISE_NCOMP
#include "fms_moisture.H"          // FMS_MOISTURE_NCOMP

using namespace amrex;

// Number of Viegas (2004) diagnostic components stored per cell.
// (Defined here because it is not exported from viegas_model.H.)
static constexpr int VIEGAS_NCOMP = 5;

// Number of spotfire-lineage components:
//   0 firebrand_src_x, 1 firebrand_src_y  (from generate_firebrand_spots)
//   2 albini_src_x,    3 albini_src_y      (from compute_albini_spotting)
// Sentinel -1.0 = no brand has landed in this cell.
static constexpr int SPOT_LINEAGE_NCOMP = 4;

struct WildfireFields {

    // ---- Ghost-cell depth for the level-set field ----
    // WENO5-Z flux divergence stencil reaches ±3 cells.
    static constexpr int ng_phi = 3;

    // ---- Core level-set and wind fields ----
    MultiFab phi;            // Level-set signed-distance function (1 comp, ng_phi ghosts)
    MultiFab vel;            // Ambient wind field (3 comps, 1 ghost)
    MultiFab vel_effective;  // Terrain/heat-flux modified wind (3 comps, 1 ghost)
    MultiFab farsite_spread; // FARSITE ellipse displacement (AMREX_SPACEDIM comps)

    // ---- Spotting diagnostics ----
    // spotting_data: 4 comps (probability, count, distance, active flag)
    MultiFab spotting_data;
    // albini_data: 4 comps (lofting H_z, firebrand count, max distance, active flag)
    MultiFab albini_data;

    // ---- Fire-spread rate ----
    MultiFab R_mf;                  // Rate of spread [m/s] (1 comp)
    MultiFab fuel_consumption_mf;   // Bulk fuel consumption fraction [-] (1 comp)
    MultiFab crown_fire_fraction_mf;// Crown fire fraction [-] (1 comp)

    // ---- Fire behaviour diagnostics ----
    MultiFab fireline_intensity_mf; // Byram fireline intensity [kW/m] (1 comp)
    MultiFab flame_length_mf;       // Byram flame length [m] (1 comp)

    // ---- Fire ecology (ECOLOGY_NCOMP components) ----
    //   0 scorch_height, 1 prob_ignition, 2 tree_mortality, 3 crown_activity,
    //   4 torching_ratio, 5 crowning_ratio, 6 energy_release_component
    MultiFab ecology_mf;

    // ---- Fire emissions (EMISSIONS_NCOMP = 3 components) ----
    //   0 co2, 1 co, 2 pm2.5  [kg/m²]
    MultiFab emissions_mf;

    // ---- Temporal diagnostics ----
    MultiFab arrival_time_mf;  // Simulation time of first ignition [s]; -1 = unburned
    MultiFab burnout_time_mf;  // arrival_time + tau_residence [s]; -1 = unburned
    MultiFab residual_fuel_mf; // Residual fuel fraction [-]; 1 = fully loaded

    // ---- Energy diagnostics ----
    MultiFab heat_per_unit_area_mf; // I_R × tau_residence [BTU/ft²]
    MultiFab heat_flux_mf;          // External heat flux [W/m²]
    MultiFab reaction_intensity_mf; // Rothermel I_R [kW/m²]

    // ---- Fluid diagnostics ----
    MultiFab vorticity_mf; // Vertical vorticity ω_z = ∂v/∂x − ∂u/∂y [s⁻¹]

    // ---- Crown fire thresholds (Scott & Reinhardt 2001) ----
    MultiFab ti_full_mf; // Torching index [km/h]; -1 = disabled
    MultiFab ci_full_mf; // Crowning index [km/h]; -1 = disabled

    // ---- Fire whirl (Weise & Biging 1996, WEISE_NCOMP = 6 components) ----
    //   0 flame_height, 1 flame_tilt, 2 whirl_height,
    //   3 whirl_radius, 4 angular_velocity, 5 max_tang_vel
    MultiFab weise_data;

    // ---- Eruptive fire (Viegas 2004, VIEGAS_NCOMP = 5 components) ----
    //   0 ROS, 1 eruptive_flag, 2 ROS_excess, 3 flame_tilt, 4 slope_factor
    MultiFab viegas_data;

    // ---- Solar radiation ----
    MultiFab shade_fraction_mf; // Per-cell shade fraction [-] (0=sun, 1=shade)

    // ---- Terrain and landscape fields ----
    MultiFab elevation_mf;  // Ground elevation [m]
    MultiFab slope_mf;      // Terrain slope [degrees]
    MultiFab aspect_mf;     // Terrain aspect [degrees]
    MultiFab fuel_model_mf; // Scott & Burgan fuel model code (integer stored as Real)

    // ---- Crown spatial layers (from binary LCP file) ----
    MultiFab cbh_mf;          // Crown base height [m]
    MultiFab cbd_mf;          // Crown bulk density [kg/m³]
    MultiFab cc_mf;           // Canopy cover [%]
    MultiFab canopy_height_mf;// Stand (canopy) height [m]

    // ---- Spatially-varying fuel moisture (FMS_MOISTURE_NCOMP = 5 components) ----
    //   0 M_d1, 1 M_d10, 2 M_d100, 3 M_lh, 4 M_lw  (fractions)
    MultiFab spatial_moisture_mf;

    // ---- Smoke plume rise (Briggs 1965) ----
    MultiFab plume_rise_mf;     // Plume-rise height [m] per fire-front cell (1 comp)

    // ---- Spread direction raster (2 comps: spread_vel_x, spread_vel_y) ----
    // Stores the fire spread velocity vector at each cell:
    //   spread_vel = R_mf × (wind_unit_vector)  [m/s in x and y]
    // Computed from effective wind and finalized R_mf each timestep.
    MultiFab spread_dir_mf;

    // ---- Spot-fire catching (ignition) probability raster (1 comp) ----
    // P_catch_local = P_catch × max(0, (M_x − M_f_local) / M_x)
    // Spatially varies when spatial fuel-moisture data are available.
    MultiFab spot_catch_prob_mf;

    // ---- Spotfire lineage (SPOT_LINEAGE_NCOMP = 4 comps) ----
    //   0 firebrand_src_x  – x-coord [m] of the source fire cell (firebrand model)
    //   1 firebrand_src_y  – y-coord [m] of the source fire cell (firebrand model)
    //   2 albini_src_x     – x-coord [m] of the source fire cell (Albini model)
    //   3 albini_src_y     – y-coord [m] of the source fire cell (Albini model)
    // Sentinel -1.0 indicates no brand has landed in that cell for that model.
    MultiFab spotting_lineage_mf;

    // ---- Terrain slopes helper (2 comps: slope_x, slope_y) ----
    // nullptr when no terrain or landscape file is configured.
    std::unique_ptr<MultiFab> terrain_slopes;

    // -----------------------------------------------------------------------
    // Constructor: allocate and initialise every field.
    // ba/dm define the box decomposition; inputs provides default values for
    // fields that are seeded from fuel and crown parameters.
    // -----------------------------------------------------------------------
    WildfireFields(const BoxArray& ba,
                   const DistributionMapping& dm,
                   const InputParameters& inputs)
        // Core fields
        : phi                  (ba, dm, 1,                 ng_phi)
        , vel                  (ba, dm, 3,                 1)
        , vel_effective        (ba, dm, 3,                 1)
        , farsite_spread       (ba, dm, AMREX_SPACEDIM,    0)
        // Spotting
        , spotting_data        (ba, dm, 4,                 0)
        , albini_data          (ba, dm, 4,                 0)
        // Fire spread
        , R_mf                 (ba, dm, 1,                 0)
        , fuel_consumption_mf  (ba, dm, 1,                 0)
        , crown_fire_fraction_mf(ba, dm, 1,                0)
        // Fire behaviour
        , fireline_intensity_mf(ba, dm, 1,                 0)
        , flame_length_mf      (ba, dm, 1,                 0)
        // Ecology / emissions
        , ecology_mf           (ba, dm, ECOLOGY_NCOMP,     0)
        , emissions_mf         (ba, dm, EMISSIONS_NCOMP,   0)
        // Temporal
        , arrival_time_mf      (ba, dm, 1,                 0)
        , burnout_time_mf      (ba, dm, 1,                 0)
        , residual_fuel_mf     (ba, dm, 1,                 0)
        // Energy
        , heat_per_unit_area_mf(ba, dm, 1,                 0)
        , heat_flux_mf         (ba, dm, 1,                 0)
        , reaction_intensity_mf(ba, dm, 1,                 0)
        // Fluid
        , vorticity_mf         (ba, dm, 1,                 0)
        // Crown thresholds
        , ti_full_mf           (ba, dm, 1,                 0)
        , ci_full_mf           (ba, dm, 1,                 0)
        // Wind / fire whirl
        , weise_data           (ba, dm, WEISE_NCOMP,       0)
        , viegas_data          (ba, dm, VIEGAS_NCOMP,      0)
        // Solar
        , shade_fraction_mf    (ba, dm, 1,                 0)
        // Terrain / landscape
        , elevation_mf         (ba, dm, 1,                 0)
        , slope_mf             (ba, dm, 1,                 0)
        , aspect_mf            (ba, dm, 1,                 0)
        , fuel_model_mf        (ba, dm, 1,                 0)
        // Crown spatial layers
        , cbh_mf               (ba, dm, 1,                 0)
        , cbd_mf               (ba, dm, 1,                 0)
        , cc_mf                (ba, dm, 1,                 0)
        , canopy_height_mf     (ba, dm, 1,                 0)
        // Spatial moisture
        , spatial_moisture_mf  (ba, dm, FMS_MOISTURE_NCOMP, 0)
        // Plume rise
        , plume_rise_mf        (ba, dm, 1,                 0)
        // Spread direction raster
        , spread_dir_mf        (ba, dm, 2,                 0)
        // Spot catching probability
        , spot_catch_prob_mf   (ba, dm, 1,                 0)
        // Spotfire lineage
        , spotting_lineage_mf  (ba, dm, SPOT_LINEAGE_NCOMP, 0)
        // terrain_slopes starts as nullptr
        , terrain_slopes       (nullptr)
    {
        // ---- Zero / sentinel initialisation ----
        fuel_consumption_mf  .setVal(Real(0.0));
        crown_fire_fraction_mf.setVal(Real(0.0));
        albini_data          .setVal(Real(0.0));
        fireline_intensity_mf.setVal(Real(0.0));
        flame_length_mf      .setVal(Real(0.0));
        ecology_mf           .setVal(Real(0.0));
        emissions_mf         .setVal(Real(0.0));
        arrival_time_mf      .setVal(Real(-1.0)); // -1 = not yet burned
        burnout_time_mf      .setVal(Real(-1.0)); // -1 = not yet burned
        residual_fuel_mf     .setVal(Real(1.0));  //  1 = fully loaded
        heat_per_unit_area_mf.setVal(Real(0.0));
        heat_flux_mf         .setVal(Real(0.0));
        reaction_intensity_mf.setVal(Real(0.0));
        vorticity_mf         .setVal(Real(0.0));
        ti_full_mf           .setVal(Real(-1.0)); // -1 = not computed / disabled
        ci_full_mf           .setVal(Real(-1.0)); // -1 = not computed / disabled
        weise_data           .setVal(Real(0.0));
        viegas_data          .setVal(Real(0.0));
        shade_fraction_mf    .setVal(Real(0.0));
        elevation_mf         .setVal(Real(0.0));
        slope_mf             .setVal(Real(0.0));
        aspect_mf            .setVal(Real(0.0));
        fuel_model_mf        .setVal(Real(0.0));
        // Crown layers seeded from global crown parameters
        cbh_mf          .setVal(Real(inputs.crown.CBH));
        cbd_mf          .setVal(Real(inputs.crown.CBD));
        cc_mf           .setVal(Real(0.0));
        canopy_height_mf.setVal(Real(0.0));
        // Spatial moisture seeded from global Rothermel moisture values
        spatial_moisture_mf.setVal(static_cast<Real>(inputs.rothermel.M_d1),   0, 1, 0);
        spatial_moisture_mf.setVal(static_cast<Real>(inputs.rothermel.M_d10),  1, 1, 0);
        spatial_moisture_mf.setVal(static_cast<Real>(inputs.rothermel.M_d100), 2, 1, 0);
        spatial_moisture_mf.setVal(static_cast<Real>(inputs.rothermel.M_lh),   3, 1, 0);
        spatial_moisture_mf.setVal(static_cast<Real>(inputs.rothermel.M_lw),   4, 1, 0);
        plume_rise_mf.setVal(Real(0.0));
        // Spread direction and spot catching probability
        spread_dir_mf      .setVal(Real(0.0));
        spot_catch_prob_mf .setVal(Real(0.0));
        // Spotfire lineage: sentinel -1 = no brand has landed
        spotting_lineage_mf.setVal(Real(-1.0));
    }
};

#endif // MULTIFAB_SETUP_H
