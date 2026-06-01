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
    // ember_cascade_mf: 2 comps (landing flux density [embers/m²/s], ignition flag)
    MultiFab ember_cascade_mf;

    // ---- Fire-spread rate ----
    MultiFab R_mf;                  // Rate of spread [m/s] (1 comp)
    MultiFab fuel_consumption_mf;   // Bulk fuel consumption fraction [-] (1 comp)
    MultiFab crown_fire_fraction_mf;// Crown fire fraction [-] (1 comp)

    // ---- Fire behaviour diagnostics ----
    MultiFab fireline_intensity_mf;      // Byram fireline intensity [kW/m] (1 comp)
    MultiFab flame_length_mf;            // Byram flame length [m] (1 comp)
    MultiFab fire_intensity_class_mf;    // Fire intensity class (1-6, Byram 1959) (1 comp)
    MultiFab convective_number_mf;       // Byram convective number Nc [-] (1 comp)
    MultiFab flame_tilt_mf;              // Flame tilt angle [degrees from vertical] (1 comp)
    MultiFab packing_ratio_mf;           // Packing ratio β [-] (1 comp)
    MultiFab relative_packing_mf;        // Relative packing ratio β/β_opt [-] (1 comp)
    MultiFab flame_depth_mf;             // Flame front depth D_f = R × τ_res [m] (1 comp)
    MultiFab mcarthur_ffdi_mf;           // McArthur Forest Fire Danger Index [-] (1 comp)

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
    // ROS at the moment each cell first ignited [m/s]; -1 = unburned.
    // Frozen at ignition time – not modified by later burn-period gates,
    // retardant, or crown-fire overrides, providing a persistent spatial
    // record of how fast the fire was spreading when it reached each cell.
    MultiFab ros_at_arrival_mf;

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

    // ---- Weather fields (multi-station interpolation) ----
    MultiFab temperature_mf;    // Per-cell air temperature [°C] (1 comp)
    MultiFab humidity_mf;       // Per-cell relative humidity [%] (1 comp)

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

    // ---- Flame length exceedance raster (1 comp) ----
    // Maximum observed Byram flame length [m] per cell over the simulation.
    // Initialised to 0; updated each timestep via element-wise max with
    // flame_length_mf so that it records the worst-case flame length seen
    // in each cell, regardless of whether the cell is actively burning.
    MultiFab fl_exceedance_mf;

    // ---- Feature 3: Crown Fraction Burned diagnostic ----
    // CFB = (I_B - I_B_surface) / (I_B_active_crown - I_B_surface) for crown fire analysis
    MultiFab crown_fraction_burned_mf;

    // ---- Feature 4: Effective Wind Speed (Wind + Slope combined) ----
    // Scalar effective wind speed combining ambient wind and slope effects [m/s]
    MultiFab effective_wind_speed_mf;

    // ---- Feature 8: Burnout Time Separation ----
    // Tracks flaming phase duration and smoldering phase separately
    // comp 0: flaming_duration [s], comp 1: smoldering_duration [s]
    MultiFab burnout_phases_mf;

    // ---- Feature 10: Post-Frontal Residual Heat Release ----
    // Tracks residual smoldering heat release post-flame-front [kW/m²]
    // Decays exponentially after fire front passage
    MultiFab residual_heat_release_mf;

    // ---- Feature 10: Time-since-burnfront (for post-frontal decay) ----
    // Records elapsed time since each cell was burned, for post-frontal decay
    // -1 = not yet burned; >= 0 = seconds since arrival
    MultiFab time_since_burn_mf;

    // ---- Terrain slopes helper (2 comps: slope_x, slope_y) ----
    // nullptr when no terrain or landscape file is configured.
    std::unique_ptr<MultiFab> terrain_slopes;

    // ---- Fire acceleration state (3 components for FARSITE temporal model) ----
    //   0 ROS_current    – current instantaneous ROS [m/s]
    //   1 elapsed_time   – time elapsed since entering current acceleration state [s]
    //   2 ROS_equilibrium – equilibrium (quasi-steady) ROS target [m/s]
    // Only allocated when acceleration.use_temporal = 1; nullptr otherwise.
    std::unique_ptr<MultiFab> accel_state_mf;

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
        , ember_cascade_mf     (ba, dm, 2,                 0)
        // Fire spread
        , R_mf                 (ba, dm, 1,                 0)
        , fuel_consumption_mf  (ba, dm, 1,                 0)
        , crown_fire_fraction_mf(ba, dm, 1,                0)
        // Fire behaviour
        , fireline_intensity_mf(ba, dm, 1,                 0)
        , flame_length_mf      (ba, dm, 1,                 0)
        , fire_intensity_class_mf(ba, dm, 1,               0)
        , convective_number_mf (ba, dm, 1,                 0)
        , flame_tilt_mf        (ba, dm, 1,                 0)
        , packing_ratio_mf     (ba, dm, 1,                 0)
        , relative_packing_mf  (ba, dm, 1,                 0)
        , flame_depth_mf       (ba, dm, 1,                 0)
        , mcarthur_ffdi_mf     (ba, dm, 1,                 0)
        // Ecology / emissions
        , ecology_mf           (ba, dm, ECOLOGY_NCOMP,     0)
        , emissions_mf         (ba, dm, EMISSIONS_NCOMP,   0)
        // Temporal
        , arrival_time_mf      (ba, dm, 1,                 0)
        , burnout_time_mf      (ba, dm, 1,                 0)
        , residual_fuel_mf     (ba, dm, 1,                 0)
        , ros_at_arrival_mf    (ba, dm, 1,                 0)
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
        , spotting_lineage_mf(ba, dm, SPOT_LINEAGE_NCOMP, 0)
        // Weather fields
        , temperature_mf(ba, dm, 1, 0)
        , humidity_mf(ba, dm, 1, 0)
        , fl_exceedance_mf     (ba, dm, 1,                 0)
        // Feature 3: Crown fraction burned
        , crown_fraction_burned_mf(ba, dm, 1,              0)
        // Feature 4: Effective wind speed
        , effective_wind_speed_mf(ba, dm, 1,              0)
        // Feature 8: Burnout phases (flaming/smoldering)
        , burnout_phases_mf    (ba, dm, 2,                 0)
        // Feature 10: Residual heat release and time tracking
        , residual_heat_release_mf(ba, dm, 1,             0)
        , time_since_burn_mf   (ba, dm, 1,                 0)
        // terrain_slopes starts as nullptr
        , terrain_slopes       (nullptr)
        // accel_state_mf starts as nullptr; allocated if use_temporal = 1
        , accel_state_mf       (nullptr)
    {
        // ---- Zero / sentinel initialisation ----
        fuel_consumption_mf  .setVal(Real(0.0));
        crown_fire_fraction_mf.setVal(Real(0.0));
        albini_data          .setVal(Real(0.0));
        ember_cascade_mf     .setVal(Real(0.0));
        fireline_intensity_mf.setVal(Real(0.0));
        flame_length_mf      .setVal(Real(0.0));
        fire_intensity_class_mf.setVal(Real(0.0));
        convective_number_mf .setVal(Real(0.0));
        flame_tilt_mf        .setVal(Real(0.0));
        packing_ratio_mf     .setVal(Real(0.0));
        relative_packing_mf  .setVal(Real(1.0));  // β/β_opt default ~1
        flame_depth_mf       .setVal(Real(0.0));
        mcarthur_ffdi_mf     .setVal(Real(0.0));
        ecology_mf           .setVal(Real(0.0));
        emissions_mf         .setVal(Real(0.0));
        arrival_time_mf      .setVal(Real(-1.0)); // -1 = not yet burned
        burnout_time_mf      .setVal(Real(-1.0)); // -1 = not yet burned
        residual_fuel_mf     .setVal(Real(1.0));  //  1 = fully loaded
        ros_at_arrival_mf    .setVal(Real(-1.0)); // -1 = not yet burned
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

        // Default T/RH values (updated each timestep by multi_wtr or diurnal_moisture)
        temperature_mf.setVal(Real(25.0));  // 25°C default
        humidity_mf.setVal(Real(30.0));     // 30% RH default
        plume_rise_mf.setVal(Real(0.0));
        // Spread direction and spot catching probability
        spread_dir_mf      .setVal(Real(0.0));
        spot_catch_prob_mf .setVal(Real(0.0));
        // Spotfire lineage: sentinel -1 = no brand has landed
        spotting_lineage_mf.setVal(Real(-1.0));
        // Flame-length exceedance: 0 = no fire observed yet
        fl_exceedance_mf.setVal(Real(0.0));
        // Feature 3: Crown fraction burned (0-1, or -1 if disabled/not computed)
        crown_fraction_burned_mf.setVal(Real(-1.0));
        // Feature 4: Effective wind speed [m/s]
        effective_wind_speed_mf.setVal(Real(0.0));
        // Feature 8: Burnout phases
        burnout_phases_mf.setVal(Real(0.0));  // comp 0: flaming, comp 1: smoldering
        // Feature 10: Residual heat release and time tracking
        residual_heat_release_mf.setVal(Real(0.0));
        time_since_burn_mf.setVal(Real(-1.0));  // -1 = not yet burned
    }
};

#endif // MULTIFAB_SETUP_H
