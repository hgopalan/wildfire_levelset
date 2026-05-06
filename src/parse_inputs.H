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
    // --- Per-class fuel moistures (Rothermel 1972 multi-class path) ---
    // Default to M_f for dead classes; live class defaults are set in parse_inputs.cpp.
    // When all are equal to M_f, behaviour is identical to the single-class path.
    amrex::Real M_d1;         // 1-hr  dead moisture (fraction)            default: M_f
    amrex::Real M_d10;        // 10-hr dead moisture (fraction)            default: M_f
    amrex::Real M_d100;       // 100-hr dead moisture (fraction)           default: M_f
    amrex::Real M_lh;         // live herbaceous moisture (fraction)       default: 0.90
    amrex::Real M_lw;         // live woody moisture (fraction)            default: 1.20
    // --- Per-class fuel loads and SAVs (populated from FuelModel database) ---
    // When all loads are zero the single-class path (w0/sigma/M_f/M_x) is used.
    amrex::Real w_d1;         // 1-hr  dead load [lb/ft²]                  default: 0.0
    amrex::Real sigma_d1;     // 1-hr  dead SAV  [ft⁻¹]                    default: 0.0
    amrex::Real w_d10;        // 10-hr dead load [lb/ft²]                  default: 0.0
    amrex::Real w_d100;       // 100-hr dead load [lb/ft²]                 default: 0.0
    amrex::Real w_lh;         // live herbaceous load [lb/ft²]             default: 0.0
    amrex::Real sigma_lh;     // live herbaceous SAV [ft⁻¹]                default: 0.0
    amrex::Real w_lw;         // live woody load [lb/ft²]                  default: 0.0
    amrex::Real sigma_lw;     // live woody SAV [ft⁻¹]                     default: 0.0
    // --- Terrain slope (constant over domain) ---
    amrex::Real slope_x;      // ∂z/∂x = tan(slope angle in x-direction)  default: 0.0
    amrex::Real slope_y;      // ∂z/∂y = tan(slope angle in y-direction)  default: 0.0
    std::string terrain_file; // Path to terrain file (X Y Z format) for spatially-varying slopes
    std::string landscape_file; // Path to FARSITE landscape file (.lcp or ASCII) with elevation, slope, aspect, fuel data
    std::string landscape_fuel_type; // Fuel model system in landscape file: "13" (FBFM13) or "40" (FBFM40)
    // --- Andrews (2018) wind adjustments ---
    int use_waf;              // 1: apply Wind Adjustment Factor (converts 20-ft→midflame wind)  default: 0
    int use_wind_limit;       // 1: apply maximum effective wind speed (MEWS) limit              default: 0
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
    amrex::Real final_time;   // simulation stop time [s]; if > 0 overrides nsteps
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
        int use_anderson_LW;                 // 1 to use Anderson (1983) L/W ratio based on wind speed, 0 to use fixed coefficients
        // Richards' (1990) elliptical expansion coefficients
        amrex::Real coeff_a;                 // head fire coefficient (default: 1.0)
        amrex::Real coeff_b;                 // flank fire coefficient (default: 0.5)
        amrex::Real coeff_c;                 // backing fire coefficient (default: 0.2)
        // Fire shape model selection (only used with propagation_method = farsite)
        // "richards"           – Richards (1990) ellipse (default; linear head/flank/backing combination)
        // "catchpole_demestre" – Catchpole & de Mestre (1986) true double-ellipse from fire-front origin
        // "wilson"             – Wilson (1988) single ellipse with fire origin at rear focus
        // "lemniscate"         – Alexander et al. Limaçon / lemniscate-shaped wavelet
        std::string fire_shape_model;
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
    
    // Cheney & Gould (1995) / Cheney, Gould & Catchpole (1998) grassland fire
    // spread model parameters
    struct CheneyGouldParams {
        amrex::Real moisture;  // dead fine fuel moisture content [%]  default: 10.0
        amrex::Real curing;    // degree of curing [0–1; 1=fully cured] default: 1.0
    } cheney_gould;

    // Cruz, Alexander & Wakimoto (2005) crown fire spread model parameters
    // Calibrated for western North American conifer forests (continental US).
    // Simple algebraic formula: R = 11.02 × U₁₀^0.90 × CBD^0.19 × exp(-0.17 × MC₁₀)
    struct CruzCrownParams {
        amrex::Real CBD;   // canopy bulk density [kg/m³]              default: 0.10
        amrex::Real MC10;  // 10-h timelag fine fuel moisture content [%] default: 10.0
    } cruz_crown;

    // Balbi (2009) physical fire spread model parameters
    struct BalbiParams {
        amrex::Real T_a;         // Ambient temperature [K] (default: 300.0)
        amrex::Real T_f;         // Mean flame temperature [K] (default: 1000.0)
        amrex::Real T_i;         // Ignition temperature [K] (default: 600.0)
        amrex::Real delta_H;     // Latent heat of vaporization of water [J/kg] (default: 2.26e6)
        amrex::Real C_pf;        // Specific heat of dry fuel [J/(kg·K)] (default: 1800.0)
        amrex::Real r_00;        // Radiation factor [m] (default: 2.5e-4)
        amrex::Real tau_0;       // Residence-time coefficient [s/m] (default: 75591.0)
    } balbi;

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

    // Weise & Biging (1996) fire whirl model parameters
    struct WeiseBigingParams {
        int enable;                          // 1 to enable, 0 to disable (default: 0)
        amrex::Real c_r;                     // whirl core radius-to-height ratio (default: 0.1)
        amrex::Real I_B_min;                 // minimum Byram fireline intensity to compute whirl [kW/m] (default: 1.0)
    } weise_biging;

    // Viegas (2004) eruptive fire model parameters
    struct ViegasParams {
        int enable;                          // 1 to enable Viegas diagnostics, 0 to disable (default: 0)
        amrex::Real a_V;                     // Viegas slope coefficient (default: 1.83)
        amrex::Real tan_phi_c;               // Critical slope for eruptive fire tan(phi_c) (default: 0.4, ~22°)
        amrex::Real T_a;                     // Ambient temperature [K] for flame-tilt calculation (default: 300.0)
        amrex::Real T_f;                     // Flame temperature [K] for flame-tilt calculation (default: 1000.0)
    } viegas;

    // Turbulent wind perturbation model.
    // Adds stochastic velocity perturbations to the ambient wind field at each
    // timestep.  The unperturbed base wind is always preserved internally so
    // the perturbation is computed relative to the clean background.
    //
    //  "none"           (default) – no perturbation
    //  "ou_process"     – Ornstein-Uhlenbeck temporally-correlated noise.
    //                     When L_c = 0: domain-uniform; all cells receive the
    //                     same (u', v') increment each step.
    //                     When L_c > 0: per-cell spatially correlated noise
    //                     generated by driving the OU update with white noise
    //                     smoothed via a separable Gaussian kernel of length
    //                     scale L_c [m].
    //  "spectral_noise" – Random Fourier Feature spectral noise with OU
    //                     temporal evolution.  N_modes random wavenumbers are
    //                     drawn at init from the 2-D isotropic Gaussian PSD
    //                     (length scale L_c); scalar OU amplitudes per mode
    //                     evolve each step; the perturbation field is
    //                     reconstructed on the GPU as a superposition of
    //                     cosine waves.
    //  "direction_walk" – bounded cumulative random walk of wind direction.
    //                     Wind speed is preserved; only direction fluctuates.
    //
    // Parameters (under the "turb_wind." prefix):
    //   model       – model name (default: "none")
    //   theta       – OU reversion rate [1/s] (default: 0.1)
    //   sigma       – OU stationary std dev [m/s] (default: 0.5)
    //   L_c         – spatial correlation length [m]; 0 = domain-uniform (default: 0.0)
    //   N_modes     – number of Fourier modes for spectral_noise (default: 32)
    //   sigma_theta – direction_walk angular std dev [rad/step] (default: 0.1)
    //   theta_max   – direction_walk max cumulative deviation [rad] (default: 0.5236 ≈ 30°)
    //   random_seed – RNG seed; 0 = use system clock (default: 0)
    struct TurbWindParams {
        std::string model;           // model name (default: "none")
        amrex::Real theta;           // OU reversion rate [1/s] (default: 0.1)
        amrex::Real sigma;           // OU stationary std dev [m/s] (default: 0.5)
        amrex::Real L_c;             // spatial correlation length [m] (default: 0.0)
        int N_modes;                 // Fourier modes for spectral_noise (default: 32)
        amrex::Real sigma_theta;     // direction_walk angular std dev [rad/step] (default: 0.1)
        amrex::Real theta_max;       // direction_walk max cumulative angle [rad] (default: 0.5236)
        int random_seed;             // RNG seed; 0 = system clock (default: 0)
    } turb_wind;

    // Wind–terrain feedback model for Rothermel ROS enhancement.
    // Selects how terrain-induced or fire-feedback winds modify the Rothermel
    // wind field before (or after) the ROS computation.
    //
    //  "none"                  (Option 1) – default Rothermel, no modification
    //  "viegas_ros"            (Option 2) – replace R_mf with max(R_rothermel, R_viegas)
    //                                       in eruptive cells (requires terrain slopes)
    //  "viegas_wind"           (Option 3) – add Viegas buoyancy-induced upslope wind only
    //                                       in eruptive cells (tan_phi > tan_phi_c)
    //  "canyon_wind"           (Option 4) – Rothermel (1983) canyon wind: scale ambient
    //                                       wind by (1 + k_canyon * tan_phi)
    //  "viegas_neto"           (Option 5) – Viegas & Neto (1994): add buoyancy-induced
    //                                       upslope wind at every cell (no eruptive threshold)
    //  "pimont"                (Option 6) – Pimont et al. (2009): scale ambient wind by
    //                                       exp(k_pimont * tan_phi) (CFD-calibrated)
    //  "windninja_ridge_canyon"(Option 7) – WindNinja empirical ridge speed-up and
    //                                       canyon channeling based on wind-slope alignment:
    //                                       ridge  (alignment > 0): f = 1 + k_ridge * tan_phi * alignment
    //                                       canyon (alignment < 0): f = 1 + k_canyon_wn * tan_phi * |alignment|
    //
    // For "viegas_ros", "viegas_wind", and "viegas_neto" the Viegas diagnostic
    // model is automatically enabled (viegas.enable is set to 1).
    // Buoyancy velocity parameters (T_a, T_f, delta) are shared with the
    // ViegasParams and RothermelParams structs respectively.
    struct WindTerrainParams {
        std::string model;        // wind-terrain model name (default: "none")
        amrex::Real k_canyon;     // Option 4: terrain channeling coefficient (default: 1.0)
        amrex::Real k_pimont;     // Option 6: exponential slope correction coefficient (default: 0.5)
        amrex::Real k_ridge;      // Option 7: ridge speed-up coefficient (default: 1.0)
        amrex::Real k_canyon_wn;  // Option 7: canyon channeling coefficient (default: 0.5)
    } wind_terrain;

    // Heat flux MultiFab parameters for fire-induced wind corrections.
    // Provides upward velocity (fire plume buoyancy) and induced horizontal
    // inflow (fire entrainment) based on WindNinja-style parameterisations.
    //
    //  heat_flux.value   – uniform heat flux [W/m²] (0 = disabled)
    //  heat_flux.file    – path to ASCII file (X Y Q columns; overrides .value)
    //  heat_flux.rho_air – air density [kg/m³] (default: 1.2)
    //  heat_flux.Cp_air  – specific heat of air [J/(kg·K)] (default: 1005.0)
    //  heat_flux.T_a     – ambient temperature [K] (default: 300.0)
    //  heat_flux.ref_height – reference height for buoyancy [m] (default: 10.0)
    //  heat_flux.k_upward   – upward velocity coefficient (default: 1.0)
    //  heat_flux.k_induced  – induced inflow coefficient (default: 0.5)
    //  heat_flux.enable_upward  – 1 to add upward velocity term (default: 0)
    //  heat_flux.enable_induced – 1 to add induced horizontal inflow (default: 0)
    struct HeatFluxParams {
        amrex::Real heat_flux_value;  // uniform heat flux [W/m²] (default: 0.0)
        std::string heat_flux_file;   // path to X Y Q ASCII file (default: "")
        amrex::Real rho_air;          // air density [kg/m³] (default: 1.2)
        amrex::Real Cp_air;           // specific heat of air [J/(kg·K)] (default: 1005.0)
        amrex::Real T_a;              // ambient temperature [K] (default: 300.0)
        amrex::Real ref_height;       // reference height for buoyancy [m] (default: 10.0)
        amrex::Real k_upward;         // upward velocity coefficient (default: 1.0)
        amrex::Real k_induced;        // induced inflow coefficient (default: 0.5)
        int enable_upward;            // 1 to enable upward velocity term (default: 0)
        int enable_induced;           // 1 to enable induced horizontal inflow (default: 0)
    } heat_flux;

    // Fire ecology diagnostics (always computed, written to every plotfile)
    // Parameters for scorch height, probability of ignition, tree mortality,
    // and passive/active crown fire classification.
    struct FireEcologyParams {
        amrex::Real T_a_C;           // ambient temperature [°C] for scorch height / P_i (default: 25.0)
        amrex::Real solar_heating_F; // solar temperature increment [°F] added to fuel temp for P_i (default: 25.0)
        amrex::Real tree_height;     // mean stand height [m] for tree mortality (default: 10.0)
    } fire_ecology;

    // Fire emissions (CO₂, CO, PM₂.₅) computed from fuel consumption
    // Follows WRF-Fire / QUIC-Fire emission factor convention.
    struct EmissionsParams {
        amrex::Real EF_CO2;              // CO₂ emission factor [kg CO₂/kg fuel] (default: 1.57)
        amrex::Real EF_CO;               // CO  emission factor [kg CO /kg fuel] (default: 0.102)
        amrex::Real EF_PM25;             // PM₂.₅ emission factor [kg PM₂.₅/kg fuel] (default: 0.0162)
        amrex::Real default_consumed_frac; // default consumption fraction for burned cells without
                                           // bulk fuel model (default: 0.7)
    } emissions;

    // Skip level set evolution option
    // Kept for backward compatibility; prefer propagation_method = "farsite" instead.
    // When set to 1, overrides propagation_method to "farsite".
    int skip_levelset;                       // deprecated: use propagation_method = farsite

    // Fire spread model selection: "rothermel" (default), "balbi", "cheney_gould",
    // or "cruz_crown"
    std::string fire_spread_model;

    // Propagation method: "levelset" (default), "farsite", or "mtt"
    //
    //  "levelset"  – WENO5-Z/RK3 level-set advection (default).
    //  "farsite"   – FARSITE Huygens elliptical wavelet propagation.
    //  "mtt"       – Minimum Travel Time (Finney 2002): arrival times are
    //                pre-computed once via a Dijkstra fast-marching sweep using
    //                the selected ROS model; phi is updated each step as
    //                phi = arrival_time - current_time.
    std::string propagation_method;

    // ---- FARSITE barrier polygon / firebreak files ----
    // List of CSV files, each containing X Y [Z] coordinates of barrier
    // vertices (polyline or closed polygon).  At every time step any barrier
    // cell that is burning (phi < 0) is extinguished (phi set to +half_cell).
    // Specify multiple files as a space-separated list:
    //   barrier_files = road.csv ridge_break.csv
    // Empty list = disabled (default).
    std::vector<std::string> barrier_files;

    // CSV fire points initialization (takes precedence over source_type when non-empty)
    std::string fire_points_file;            // Path to CSV file with X Y [Z] ignition point coordinates
    amrex::Real fire_gaussian_sigma;         // Gaussian sigma for each ignition point (physical units); <=0 means auto (3 cells)

    // Dynamic fire points file (checked at every time step)
    std::string dynamic_fire_points_file;    // Path to a CSV file polled each step; when the file is present its
                                             // ignition points are merged into phi and the file is deleted/ignored
                                             // until it reappears (default: "" = disabled)

    // Checkpoint / restart options
    int chk_int;                             // Write a checkpoint file every chk_int steps; <=0 disables (default: -1)
    std::string restart_chkfile;             // Path to a checkpoint directory to restart from (default: "" = cold start)

    // ---- FARSITE fuel adjustment file (.adj) ----
    // Per-fuel-model rate-of-spread multipliers.  Applied once at startup to
    // the Rothermel lookup table (and to RothermelParams when no landscape file
    // is present).  Empty string = disabled.
    std::string fuel_adj_file;               // path to FARSITE .adj file (default: "")
    int         fuel_adj_model;              // global fuel model code for single-model runs (default: 0)

    // ---- Time-varying fuel moisture schedule (.fmd) ----
    // When set, moisture values are updated at every timestep by interpolating
    // between entries in the schedule file.  Start date anchors the schedule to
    // simulation time 0.  When all start fields are 0 the first record's
    // timestamp defines t = 0.
    std::string fmd_file;                    // path to FMD-format schedule file (default: "")
    int fmd_start_year;                      // reference year  (default: 0 → auto)
    int fmd_start_month;                     // reference month (default: 0 → auto)
    int fmd_start_day;                       // reference day   (default: 0 → auto)
    int fmd_start_hour;                      // reference hour  (0-2300, default: 0)
    int fmd_fuel_model;                      // which fuel model to look up (0 = global/all, default: 0)

    // ---- Fire statistics time series ----
    // When non-empty, a CSV file is created at startup and a row is appended at
    // every plotfile write with: step, time_s, burned_area_ha, perimeter_km,
    // active_front_cells, mean_co2, mean_co, mean_pm25.
    std::string fire_stats_file;             // output CSV path (default: "fire_stats.csv")

    // ---- Fire perimeter output ----
    // Controls whether CSV and GeoJSON perimeter files are written alongside
    // the existing .dat (phi_negative / phi_envelope) files at every plotfile step.
    int write_perimeter_csv;                 // 1 to write perimeter_NNNN.csv (default: 1)
    int write_perimeter_geojson;             // 1 to write perimeter_NNNN.geojson (default: 0)

    // ---- Crown spatial layers from binary LCP ----
    // When a binary LCP with crown fuel layers is loaded, these flags control
    // whether per-cell CBH, CBD, and canopy cover MultiFabs are populated and
    // written to plotfiles.  The global crown.* scalars are used as fallbacks
    // when the landscape does not carry crown layers.
    int use_spatial_crown;                   // 1 to use per-cell CBH/CBD from LCP (default: 1)

    // ---- Diurnal fuel moisture model (Nelson 2000 / FARSITE-style) ----
    // When enabled and no FMD schedule is provided (fmd_file is empty), dead
    // fuel moisture is updated each timestep using a sinusoidal diurnal cycle
    // of temperature and relative humidity, driving the Nelson (2000) EMC
    // equations.
    //
    // The diurnal cycle is:
    //   T(t)  = T_mean + 0.5*(T_max-T_min)*sin(2π*(t+t_T_peak-t_ref)/86400)
    //   RH(t) = RH_mean - 0.5*(RH_max-RH_min)*sin(2π*(t+t_T_peak-t_ref)/86400)
    //
    // where t_T_peak [s from midnight] is when T reaches its daily maximum
    // (default: 50400 = 14:00).
    //
    // EMC per size class (Nelson 2000 / Simard 1968):
    //   M_1hr   = EMC(T, RH)
    //   M_10hr  ≈ EMC(T, RH) * 1.10  (10-hr timelag correction)
    //   M_100hr ≈ EMC(T, RH) * 1.30  (100-hr timelag correction)
    //
    // Input prefix: "diurnal_moisture."
    struct DiurnalMoistureParams {
        int enable;             // 1 to enable diurnal moisture cycle (default: 0)
        amrex::Real T_min;      // minimum daily temperature [°C]  (default: 15.0)
        amrex::Real T_max;      // maximum daily temperature [°C]  (default: 35.0)
        amrex::Real RH_min;     // minimum daily relative humidity [%] (default: 10.0)
        amrex::Real RH_max;     // maximum daily relative humidity [%] (default: 60.0)
        amrex::Real t_start_s;  // seconds from midnight at t=0 (default: 36000 = 10:00)
        amrex::Real t_T_peak_s; // seconds from midnight when T is maximum (default: 50400 = 14:00)
    } diurnal_moisture;

};

void parse_inputs(InputParameters& inputs);

#endif
