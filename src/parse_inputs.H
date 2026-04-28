#ifndef PARSE_INPUTS_H
#define PARSE_INPUTS_H

#include <AMReX_ParmParse.H>
#include <string>

struct RothermelParams {
    // --- Fuel (NFFL Model 4 / Southern California chaparral defaults) ---
    amrex::Real w0;           // oven-dry total fuel load [lb/ft²]        default: 0.230
    amrex::Real sigma;        // surface-area-to-volume ratio [ft⁻¹]      default: 1739.0
    amrex::Real delta;        // fuel bed depth [ft]                       default: 2.0
    amrex::Real M_f;          // fuel moisture content (fraction)          default: 0.08
    amrex::Real M_x;          // moisture of extinction (fraction)         default: 0.30
    amrex::Real h_heat;       // heat content [BTU/lb]                     default: 8000.0
    amrex::Real S_T;          // total mineral content (fraction)          default: 0.0555
    amrex::Real S_e;          // effective mineral content (fraction)      default: 0.010
    amrex::Real rho_p;        // particle density [lb/ft³]                 default: 32.0
    // --- Terrain slope (constant over domain) ---
    amrex::Real slope_x;      // ∂z/∂x = tan(slope angle in x-direction)  default: 0.0
    amrex::Real slope_y;      // ∂z/∂y = tan(slope angle in y-direction)  default: 0.0
    std::string terrain_file; // Path to terrain file (X Y Z format) for spatially-varying slopes
    std::string landscape_file; // Path to FARSITE landscape file (.lcp or ASCII) with elevation, slope, aspect, fuel data
    // --- Unit conversions ---
    // Velocity field is in m/s; Rothermel needs ft/min
    amrex::Real wind_conv;    // ft/min per (sim velocity unit)            default: 196.85
    // Rothermel ROS is in ft/min; convert back to sim length/time units
    amrex::Real ros_conv;     // sim length/time per (ft/min)              default: 0.00508
};

struct InputParameters {
    // Grid & domain
    int n_cell_x, n_cell_y, n_cell_z;
    int max_grid;
    amrex::Real plo_x, plo_y, plo_z;
    amrex::Real phi_x, phi_y, phi_z;

    // Time & output
    int reinit_int;
    int nsteps;
    amrex::Real cfl;
    int plot_int;

    // Velocity
    amrex::Real ux, uy, uz;
    std::string velocity_file;
    int use_time_dependent_wind;      // 1 to use time-dependent wind fields, 0 for static (default: 0)
    amrex::Real wind_time_spacing;    // Time spacing in seconds between wind field files (default: 60.0)

    // Source selection
    std::string source_type;

    // Sphere parameters
    amrex::Real cx, cy, cz, radius;
    // Box Extents 
    amrex::Real xmin, xmax, ymin, ymax, zmin, zmax;
    // Ellipse parameters
    amrex::Real ellipse_center_x, ellipse_center_y, ellipse_center_z;
    amrex::Real ellipse_radius_x, ellipse_radius_y, ellipse_radius_z;
    
    // EB (Embedded Boundary) implicit function parameters
    std::string eb_type;  // "plane", "cylinder", "sphere", "ellipsoid"
    amrex::Real eb_param1, eb_param2, eb_param3;
    amrex::Real eb_param4, eb_param5, eb_param6;

    // Rothermel fire spread model parameters
    RothermelParams rothermel;
    
    // FARSITE ellipse model parameters (Richards 1990)
    struct FARSITEParams {
        amrex::Real length_to_width_ratio;  // L/W ratio of ellipse (default: 3.0)
        amrex::Real phi_threshold;           // threshold for identifying fire front (default: 0.1)
        int enable;                          // 1 to enable FARSITE, 0 to disable
        int use_anderson_LW;                 // 1 to use Anderson (1983) L/W ratio based on wind speed, 0 to use fixed coefficients
        // Richards' (1990) elliptical expansion coefficients
        amrex::Real coeff_a;                 // head fire coefficient (default: 1.0)
        amrex::Real coeff_b;                 // flank fire coefficient (default: 0.5)
        amrex::Real coeff_c;                 // backing fire coefficient (default: 0.2)
        // Bulk Fuel Consumption Fraction Model parameters
        int use_bulk_fuel_consumption;       // 1 to enable bulk fuel consumption model, 0 to disable (default: 0)
        amrex::Real tau_residence;           // Residence time for fuel consumption in seconds (default: 60.0)
        amrex::Real f_consumed_max;          // Maximum fuel consumption fraction (0.0-1.0) (default: 0.9)
        amrex::Real f_consumed_min;          // Minimum fuel consumption fraction (0.0-1.0) (default: 0.5)
    } farsite;
    
    // Firebrand spotting model parameters
    struct SpottingParams {
        int enable;                          // 1 to enable spotting, 0 to disable (default: 0)
        amrex::Real P_base;                  // Base spotting probability (0.0-1.0) (default: 0.02)
        amrex::Real k_wind;                  // Wind speed coefficient for probability (default: 0.3)
        amrex::Real I_critical;              // Critical fire intensity for spotting (default: 1000.0)
        amrex::Real d_mean;                  // Mean spotting distance (default: 0.1)
        amrex::Real d_sigma;                 // Spotting distance std deviation for lognormal (default: 0.5)
        amrex::Real d_lambda;                // Spotting distance decay rate for exponential (default: 10.0)
        std::string distance_model;          // "lognormal" or "exponential" (default: "lognormal")
        amrex::Real lateral_spread_angle;    // Angular spread perpendicular to wind in degrees (default: 15.0)
        amrex::Real spot_radius;             // Radius of new spot fires (default: 0.02)
        int random_seed;                     // Seed for random number generator (0 = use time) (default: 0)
        int check_interval;                  // Check for spotting every N timesteps (default: 5)
    } spotting;
    
    // Van Wagner crown fire initiation model parameters (Van Wagner 1977)
    struct CrownInitiationParams {
        int enable;                          // 1 to enable crown initiation model, 0 to disable (default: 0)
        amrex::Real CBH;                     // Canopy base height [m] (default: 4.0)
        amrex::Real CBD;                     // Canopy bulk density [kg/m³] (default: 0.15)
        amrex::Real FMC;                     // Foliar moisture content [%] (default: 100.0)
        amrex::Real crown_fraction_weight;   // Crown fire weighting factor (0-2) (default: 1.0)
        int use_metric_units;                // 1 for metric (m, kW), 0 for imperial (ft, BTU) (default: 1)
    } crown;
    
    // Albini (1983) firebrand spotting model with 2-D trajectory integration
    struct AlbiniSpottingParams {
        int enable;                          // 1 to enable Albini spotting, 0 to disable (default: 0)
        amrex::Real terminal_velocity;       // Firebrand terminal descent velocity [m/s] (default: 1.0)
        amrex::Real P_base;                  // Max probability of launching a firebrand per fire-front cell (default: 0.01)
        amrex::Real I_B_min;                 // Minimum Byram fire line intensity to allow spotting [kW/m] (default: 10.0)
        amrex::Real spot_radius;             // Radius of new spot-fire ignition zone [m] (default: 5.0)
        int random_seed;                     // Seed for RNG (0 = use system time) (default: 0)
        int check_interval;                  // Run spotting check every N timesteps (default: 5)
        int n_traj_steps;                    // Number of forward-Euler steps for 2-D trajectory (default: 100)
    } albini_spotting;

    // Skip level set evolution option
    int skip_levelset;                       // 1 to skip level set advection and use initial phi, 0 for normal operation (default: 0)

    // CSV fire points initialization (takes precedence over source_type when non-empty)
    std::string fire_points_file;            // Path to CSV file with X Y [Z] ignition point coordinates
    amrex::Real fire_gaussian_sigma;         // Gaussian sigma for each ignition point (physical units); <=0 means auto (3 cells)
};

void parse_inputs(InputParameters& inputs);

#endif
