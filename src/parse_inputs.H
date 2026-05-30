#ifndef PARSE_INPUTS_H
#define PARSE_INPUTS_H

#include <AMReX_ParmParse.H>
#include <string>

// SatelliteParams is defined in satellite_assimilation.H, included here so
// InputParameters can embed it as a member.
#include "satellite_assimilation.H"

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
    amrex::Real M_d1000;      // 1000-hr dead moisture (fraction)          default: M_f
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
    // --- Fuel continuity (horizontal fuel coverage) ---
    // Accounts for heterogeneous fuel distributions where fuels are discontinuous.
    // Applied as a multiplicative factor to the computed ROS.
    //   continuity = 1.0 → continuous fuel bed (default)
    //   continuity = 0.5 → 50% fuel coverage (patchy fuels)
    //   continuity = 0.0 → no fuel (gaps)
    // References:
    //   Finney, M.A. (2006). "An overview of FlamMap fire modeling capabilities."
    //     USDA Forest Service Proceedings RMRS-P-41, pp. 213-220.
    //   Parsons, R.A., et al. (2017). "Linking 3D spatial models of fuels and fire."
    //     Ecological Modelling, 222(3):679-691.
    amrex::Real fuel_continuity; // horizontal fuel continuity factor [0-1]  default: 1.0
    // --- Terrain slope (constant over domain) ---
    amrex::Real slope_x;      // ∂z/∂x = tan(slope angle in x-direction)  default: 0.0
    amrex::Real slope_y;      // ∂z/∂y = tan(slope angle in y-direction)  default: 0.0
    std::string terrain_file; // Path to terrain file (X Y Z format) for spatially-varying slopes
    std::string landscape_file; // Path to FARSITE landscape file (.lcp or ASCII) with elevation, slope, aspect, fuel data
    std::string landscape_fuel_type; // Fuel model system in landscape file: "13" (FBFM13) or "40" (FBFM40)
    // --- Andrews (2018) wind adjustments ---
    int use_waf;              // 1: apply Wind Adjustment Factor (converts 20-ft→midflame wind)  default: 0
    int use_wind_limit;       // 1: apply maximum effective wind speed (MEWS) limit              default: 0
    // WAF formula selection (only relevant when use_waf = 1):
    //   "andrews"      (default) – logarithmic Albini & Baughman (1979) formula for open fuels;
    //                              FARSITE sheltered formula for closed-canopy fuels.
    //   "behaviorplus"           – BehavePlus linear formula (WAF = 0.36 + 0.004 × h_in) for
    //                              open/shrub fuels; exponential canopy attenuation for forest fuels.
    std::string waf_formula;  // WAF formula: "andrews" (default) or "behaviorplus"
    // Canopy attenuation coefficient for the BehavePlus exponential canopy WAF.
    // Larger values produce stronger wind sheltering under a closed canopy.
    // Only used when waf_formula = "behaviorplus" and canopy data are available.
    amrex::Real waf_canopy_alpha; // canopy attenuation coefficient α_c (default: 1.5)
    // --- Live herbaceous fuel curing transfer ---
    // When M_lh < herb_curing_threshold, a fraction of the live herbaceous
    // load is treated as 1-hr dead fuel (FARSITE curing model).
    // Set to 0.0 to disable (default).  Typical FARSITE value: 0.30 (30%).
    amrex::Real herb_curing_threshold; // M_lh fraction below which curing begins  default: 0.0 (disabled)
    // --- Unit conversions ---
    // Velocity field is in m/s; Rothermel needs ft/min
    amrex::Real wind_conv;    // ft/min per (sim velocity unit)            default: 196.85
    // Rothermel ROS is in ft/min; convert back to sim length/time units
    amrex::Real ros_conv;     // sim length/time per (ft/min)              default: 0.00508
    // --- Slope/wind interaction cross-term (Rothermel 1983 / Anderson 1982) ---
    // When enabled, replaces the additive R = R₀(1 + φ_w + φ_s) formula with
    // R = R₀(1 + φ_w + φ_s + k_sw × φ_w × φ_s).
    // This captures the positive interaction when wind blows upslope:
    // combined spread exceeds the sum of individual effects.
    //   k_slope_wind_cross = 0.0 → additive only (default, backward-compatible)
    //   k_slope_wind_cross = 1.0 → fully multiplicative (1+φ_w)(1+φ_s) form
    // Reference: Rothermel (1983) GTR INT-143; Anderson (1982) Res. Note INT-328.
    int         use_slope_wind_cross;   // 1 to enable cross-term (default: 0)
    amrex::Real k_slope_wind_cross;     // cross-term coupling coefficient (default: 1.0)
    // --- FARSITE-style vectorial slope/wind combination (new option) ---
    // When use_slope_wind_vectors = 1 (and terrain is available), the wind and
    // slope effects are added as vectors before computing the scalar ROS:
    //   φ_w_vec = φ_w × wind_unit_vector
    //   φ_s_vec = φ_s × upslope_unit_vector
    //   φ_combined = |φ_w_vec + φ_s_vec|
    //   R = R₀ × (1 + φ_combined)
    // This replaces the scalar sum R = R₀(1 + φ_w + φ_s) and naturally captures
    // the effect of wind/slope alignment: aligned vectors add constructively,
    // opposed vectors partially cancel.  The cross-term (use_slope_wind_cross)
    // is not applied when this option is active (they are mutually exclusive).
    // When no terrain data is present, falls back to scalar additive.
    // Reference: Finney (1998) FARSITE; Rothermel (1983) GTR INT-143.
    int         use_slope_wind_vectors; // 1 = FARSITE vectorial mode (default: 0)
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
        // Feature 12: ellipse shape scaling for active crown fire (backward compatible; default off)
        // When enabled (=1) and active crown fire is detected, the L/W ratio is multiplied by
        // crown_lw_scale producing a more elongated (wind-driven) ellipse consistent with
        // active crown fire behaviour (Finney 1998 / Van Wagner 1977 active-crown criteria).
        int scale_ellipse_with_crown;        // 1 to scale L/W for active crown fire (default: 0)
        amrex::Real crown_lw_scale;         // L/W multiplier for active crown fire (default: 1.5)
        amrex::Real dt;                      // Farsite time step is tricky. Using fixed time step (default:10)
        // Gaussian smoothing of FARSITE spread-point stamping (signed-distance disk around each
        // propagated position, analogous to fire_points_file / fire_gaussian_sigma).
        //   < 0  : disabled – single-cell marking, preserving the original behaviour (default: -1)
        //   = 0  : auto – sigma is set to 3 * min(dx, dy [, dz])
        //   > 0  : user-specified sigma [m]
        amrex::Real gaussian_sigma;          // sigma for SDF smoothing of FARSITE spread points [m]
    } farsite;
    
    // Cell size effects correction parameters (FARSITE only)
    // Empirical correction for fire spread rate based on cell/grid resolution.
    // Accounts for sub-grid fire behavior: finer grids capture more of the fire
    // front geometry, leading to higher effective spread rates.
    // FARSITE Reference: Finney (1998) Section 3.4; based on empirical observations
    // that fire spread speed depends on grid resolution.
    //
    // Correction formula: R_corrected = R_base × (dx_ref / dx_geom)^correction_exponent
    // where dx_geom = min(dx, dy [, dz]) is the actual grid cell size.
    //
    // - Smaller grid → R_corrected > R_base (fires spread faster at high resolution)
    // - Larger grid  → R_corrected < R_base (fires spread slower at low resolution)
    //
    // Applied only when propagation_method = "farsite" and enable = 1.
    // Not applied to levelset propagation (scale-independent).
    struct CellSizeParams {
        int enable;                         // 1 to enable cell size correction, 0 to disable (default: 0)
        amrex::Real dx_ref;                 // reference cell size [m] (default: 30.0)
        amrex::Real correction_exponent;    // power-law exponent (default: 0.1)
    } cellsize;
    
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
        // Spot-fire catching (ignition) probability
        // When P_catch < 1.0, a landed firebrand ignites the receiving fuel with
        // probability P_catch × f_moisture, where f_moisture = max(0,(M_x−M_f)/M_x)
        // (Anderson 1970 / Rothermel 1972 moisture-damping factor).
        // Set to 1.0 (default) to always ignite on contact (original behaviour).
        amrex::Real P_catch;                 // max catching probability [0-1] (default: 1.0)
        
        // Ignition delay parameters (FARSITE-compatible)
        // Models delay between firebrand landing and actual ignition based on fuel moisture.
        // τ_delay = τ_base × (M_x / (M_x - M_f))
        int enable_delay;                    // 1 = enable ignition delay, 0 = instant ignition (default: 0)
        amrex::Real tau_base;                // base delay time at zero moisture [s] (default: 120.0)
    } spotting;
    
    // Van Wagner crown fire initiation model parameters (Van Wagner 1977)
    struct CrownInitiationParams {
        int enable;                          // 1 to enable crown initiation model, 0 to disable (default: 0)
        amrex::Real CBH;                     // Canopy base height [m] (default: 4.0)
        amrex::Real CBD;                     // Canopy bulk density [kg/m³] (default: 0.15)
        amrex::Real FMC;                     // Foliar moisture content [%] (default: 100.0)
        amrex::Real crown_fraction_weight;   // Crown fire weighting factor (0-2) (default: 1.0)
        int use_metric_units;                // 1 for metric (m, kW), 0 for imperial (ft, BTU) (default: 1)
        // Active crown fire ROS model selection:
        //   0 – Van Wagner (1977) proxy: R_crown = 3.0 / CBD * moisture_factor [m/min] (default)
        //   1 – Cruz, Alexander & Wakimoto (2005): R = 11.02 × U₁₀^0.90 × CBD^0.19 × exp(-0.17 × MC₁₀)
        //       Uses cruz_crown.CBD and cruz_crown.MC10 for canopy/moisture inputs.
        int use_cruz_crown;                  // 1 = Cruz et al. (2005) active crown ROS (default: 0)
        // Rothermel (1991) active crown fire ROS option:
        //   0 - do not use (keep existing crown model selection)  [default]
        //   1 - use Rothermel (1991): R_crown = 3.34 x R_surface
        int use_rothermel1991_crown;         // 1 = Rothermel (1991) (default: 0)
        // Van Wagner (1977) passive crown fire interpolation blending:
        //   0 – abrupt max(R_surface, R_crown)  [current behaviour, default]
        //   1 – smooth blend via CF = (I_B/I_o)^(2/3) following Van Wagner (1977)
        int use_passive_blend;               // 1 = Van Wagner passive blend (default: 0)
        // Passive crown fire ellipse coefficients (used when crown fire is detected
        // and farsite.scale_ellipse_with_crown = 1). These produce a more elongated
        // fire shape consistent with active crown fire.  Blended against surface
        // coefficients proportional to the crowning fraction CF.
        amrex::Real a_crown;                 // head fire ellipse coefficient for crown fire (default: 1.5)
        amrex::Real b_crown;                 // flank fire ellipse coefficient for crown fire (default: 0.4)
        amrex::Real c_crown;                 // backing fire ellipse coefficient for crown fire (default: 0.1)
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
        // ---- 3-D wind from massconsistent_amr ----
        // When use_3d_wind = 1 and plt_wind_file is non-empty, the horizontal
        // wind used for firebrand trajectory integration is read from an AMReX
        // plotfile produced by https://github.com/hgopalan/massconsistent_amr.
        // The plotfile is expected to contain variables "u", "v", and "w" on a
        // 3-D Cartesian grid.  The horizontal wind is obtained by column-averaging
        // u and v over all vertical levels at each (x, y) location, which is then
        // used in place of the 2-D velocity MultiFab for trajectory integration.
        // (default: disabled)
        int use_3d_wind;                     // 1 to read 3-D wind from plt file (default: 0)
        std::string plt_wind_file;           // path to massconsistent_amr plt directory (default: "")
        // Spot-fire catching (ignition) probability for Albini model.
        // Applied identically to SpottingParams.P_catch: a landed firebrand
        // ignites with probability P_catch × f_moisture before placing a new
        // ignition zone in phi.  Default 1.0 = always ignite.
        amrex::Real P_catch;                 // max catching probability [0-1] (default: 1.0)
    } albini_spotting;

    // Weise & Biging (1996) fire whirl model parameters
    struct WeiseBigingParams {
        int enable;                          // 1 to enable, 0 to disable (default: 0)
        amrex::Real c_r;                     // whirl core radius-to-height ratio (default: 0.1)
        amrex::Real I_B_min;                 // minimum Byram fireline intensity to compute whirl [kW/m] (default: 1.0)
        int enhance_spotting;                // 1 to enhance Albini spotting with whirl effects, 0 to disable (default: 0)
        amrex::Real alpha;                   // spotting enhancement factor coefficient (default: 1.0)
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
    //  "farsite_wind"          (Option 8) – FARSITE empirical wind stream simulation using
    //                                       terrain curvature to identify ridges/valleys and
    //                                       modify both wind speed and direction
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
        // Option 8 (FARSITE wind stream simulation) parameters:
        amrex::Real k_ridge_farsite;  // FARSITE ridge speed-up coefficient (default: 1.5)
        amrex::Real k_shelter;        // lee-side sheltering coefficient (default: 0.6)
        amrex::Real k_valley;         // valley channeling coefficient (default: 0.8)
        amrex::Real k_deflection;     // wind direction deflection coefficient (default: 0.3)
        amrex::Real min_curvature;    // minimum curvature threshold (default: 0.0001)
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
        // When couple_to_ros = 1, the computed probability-of-ignition (P_i) is used to
        // scale R_mf before advection: R_eff = R * max(P_i, p_i_floor).
        // This makes fire spread slower through wet or cold fuel where sustained ignition
        // is less likely (Anderson 1970 / Rothermel 1983 coupling).
        int couple_to_ros;           // 1 to scale R_mf by P_ignition (default: 0)
        amrex::Real p_ignition_floor;// minimum P_i used as a floor when scaling ROS (default: 0.05)
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

    // ---- Multiple scheduled ignitions ----
    // CSV file of (time_s, cx, cy, radius [, source_type]) rows.
    // Each row schedules one ignition event that fires when the simulation
    // time enters (time_s - dt, time_s].  Allows multi-start scenarios such
    // as lightning strikes, spotfire confirmations, or drip-torch runs.
    // (default: "" = disabled)
    std::string ignition_schedule_file;      // path to ignition schedule CSV (default: "")

    // ---- FARSITE .wtr single-file hourly weather ----
    // When non-empty, a FARSITE .wtr file is parsed and its wind, temperature,
    // RH, and precipitation records populate the WindDirSchedule and override
    // the diurnal_moisture parameters, replacing wind_dir_schedule_file and
    // manual diurnal_moisture.T_min/T_max/RH_min/RH_max settings.
    // (default: "" = disabled)
    std::string wtr_file;                    // path to FARSITE .wtr weather file (default: "")
    int         wtr_start_year;              // calendar year  of first .wtr record (0=auto)
    int         wtr_start_month;             // calendar month of first .wtr record (0=auto)
    int         wtr_start_day;              // calendar day   of first .wtr record (0=auto)
    int         wtr_start_hour;             // HHMM hour      of first .wtr record (0=auto)

    // ---- Multiple weather stations with spatial interpolation ----
    // When non-empty, a CSV station-list file is parsed to build spatially-
    // varying wind and weather fields via inverse-distance weighting (IDW) of
    // per-station .wtr records.  Each row in the station list is:
    //   station_id, x_m, y_m, wtr_file
    // The IDW velocity field replaces the constant ux/uy wind at each timestep.
    // Domain-mean T/RH from IDW centroid are used for the diurnal moisture model.
    // (default: "" = disabled; multi_wtr_file takes precedence over wtr_file)
    std::string multi_wtr_file;              // path to station list CSV (default: "")
    amrex::Real multi_wtr_idw_power;         // IDW exponent p (default: 2.0)

    // ---- Aerial retardant suppression model ----
    // CSV file of drop zones: (drop_time_s, effectiveness, tau_decay_s, n_vertices)
    // followed by n_vertices rows of (x y) polygon vertices.
    // At each timestep, active drops scale R_mf by
    //   f = 1 - effectiveness * exp(-(t - drop_time) / tau_decay).
    // (default: "" = disabled)
    std::string retardant_file;              // path to retardant schedule CSV (default: "")

    // ---- Fuel moisture conditioning period ----
    // Before the main time loop, run the Nelson (2000) diurnal EMC model for
    // N days of weather to arrive at realistic starting moisture values.
    // When conditioning.wtr_file is non-empty, the .wtr weather data provides
    // the conditioning T/RH/precipitation.  Otherwise, diurnal_moisture.* are used.
    // conditioning.n_days = 0 (default) disables conditioning.
    struct ConditioningParams {
        int         n_days;      // number of conditioning days (default: 0 = disabled)
        std::string wtr_file;    // optional separate .wtr file for conditioning weather
    } conditioning;

    // ---- Elevation lapse-rate T/RH correction ----
    // When enabled (= 1), the diurnal temperature and relative humidity used to
    // compute per-cell EMC are adjusted for cell elevation via a standard
    // atmospheric lapse rate before feeding into apply_solar_emc_to_spatial_moisture.
    //   T_cell  = T_air - lapse_rate_C_per_m * (elev_cell - elev_ref_m)
    //   RH_cell corrected via Clausius–Clapeyron approximation.
    // Requires: a landscape or terrain file providing elevation data.
    // (default: disabled)
    int    use_elevation_lapse;              // 1 to enable lapse-rate correction (default: 0)
    amrex::Real lapse_rate_C_per_m;         // temperature lapse rate [°C/m] (default: 0.0065)
    amrex::Real lapse_ref_elevation_m;      // reference elevation for lapse [m] (default: 0.0)

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

    // ---- FARSITE Fire Spread Atlas (.fsa) output ----
    // When non-empty, a cumulative fire perimeter archive is written to this
    // file.  Each plotfile step appends a SNAPSHOT block containing the
    // angle-sorted perimeter points (same cells as perimeter_NNNN.csv) with
    // a header giving the step, simulation time, and point count.  The
    // resulting single file collects all perimeter snapshots across the full
    // simulation, mirroring FARSITE's Fire Spread Atlas concept.
    // (default: "" = disabled)
    std::string fsa_file;                    // cumulative .fsa archive path

    // ---- FARSITE Post-processing Statistics (.pst) output ----
    // When non-empty, a FARSITE-compatible statistics file is written with one
    // row per plotfile step.  Columns: step, time_h, burned_ha, perimeter_km,
    // active_cells, max_fireline_intensity_kW_m, max_flame_length_m,
    // max_spotting_dist_m.  Extends the fire_stats_file CSV with maximum
    // observed fire behaviour metrics.
    // (default: "" = disabled)
    std::string pst_file;                    // .pst statistics path

    // ---- Automatic HTML report generation ----
    // When non-empty, a self-contained HTML report is written at end-of-run,
    // summarising simulation parameters, per-step fire statistics (read from
    // fire_stats_file), and SVG time-series charts of burned area, perimeter,
    // head ROS, and cumulative emissions.  The file has no external dependencies
    // and can be opened in any browser without an internet connection.
    // Requires fire_stats_file to be set so that data are available.
    std::string fire_report_file;            // output HTML path (default: "" = disabled)

    // ---- Timed isochrone output ----
    // When isochrone_interval_s > 0, a fire perimeter file is written at each
    // simulation clock-time that is an exact multiple of isochrone_interval_s.
    // Files are named isochrone_NNNNNN.csv (or .geojson) where NNNNNN is the
    // isochrone index (0, 1, 2, ...).  The file format follows the same
    // convention as perimeter_NNNN.csv / .geojson from write_perimeter_csv /
    // write_perimeter_geojson.
    // Set to 0 or negative to disable (default: 0.0 = disabled).
    amrex::Real isochrone_interval_s;        // isochrone write interval [s] (default: 0.0)

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

    // ---- Cap 10: Custom fuel model (text-based, inline in inputs file) ----
    // Allows users to define a fully custom fuel model without modifying the
    // compiled fuel database.  When ``custom_fuel.enable = 1``, all
    // ``custom_fuel.*`` parameters override the Rothermel parameter defaults
    // (or any database fuel model selected via ``rothermel.fuel_model``).
    //
    // This makes it straightforward to test a one-off fuel model — such as
    // exotic grass species, post-fire shrubs, or experimental chaparral blends
    // — directly from the inputs file without recompilation.
    //
    // Supported parameters (all optional; unset fields inherit the active
    // RothermelParams values at the time custom_fuel is applied):
    //
    //   custom_fuel.enable   = 1            # must be 1 to activate
    //   custom_fuel.name     = "my_fuel"    # label printed at startup
    //   custom_fuel.w0       = 0.150        # oven-dry total fuel load [lb/ft²]
    //   custom_fuel.sigma    = 2000.0       # surface-area-to-volume ratio [ft⁻¹]
    //   custom_fuel.delta    = 1.5          # fuel bed depth [ft]
    //   custom_fuel.M_x      = 0.20         # moisture of extinction [fraction]
    //   custom_fuel.h_heat   = 8000.0       # heat content [BTU/lb]
    //   custom_fuel.S_T      = 0.0555       # total mineral content [fraction]
    //   custom_fuel.S_e      = 0.010        # effective mineral content [fraction]
    //   custom_fuel.rho_p    = 32.0         # particle density [lb/ft³]
    struct CustomFuelParams {
        int enable;              // 1 to activate custom fuel override (default: 0)
        std::string name;        // descriptive label (default: "custom")
        amrex::Real w0;          // oven-dry total fuel load [lb/ft²]        (-1 = inherit)
        amrex::Real sigma;       // surface-area-to-volume ratio [ft⁻¹]      (-1 = inherit)
        amrex::Real delta;       // fuel bed depth [ft]                       (-1 = inherit)
        amrex::Real M_x;         // moisture of extinction [fraction]         (-1 = inherit)
        amrex::Real h_heat;      // heat content [BTU/lb]                     (-1 = inherit)
        amrex::Real S_T;         // total mineral content [fraction]          (-1 = inherit)
        amrex::Real S_e;         // effective mineral content [fraction]      (-1 = inherit)
        amrex::Real rho_p;       // particle density [lb/ft³]                 (-1 = inherit)
    } custom_fuel;

    // ---- Flux-based ember cascade model (plume-rise driven) ----
    // Models dense firebrand showers released by a rising smoke plume using a
    // continuous landing-flux approach.  Rather than tracking individual brands
    // (as in the Albini trajectory model), every fire-front cell emits a source
    // flux of embers proportional to its Byram fireline intensity.  Each source
    // contributes a 2-D Gaussian footprint downwind (centred at the wind-advected
    // mean landing position, spread by turbulent diffusion from the plume), and
    // the total landing flux density [embers/m²/s] is accumulated across all
    // active sources.  Cells where the integrated landing density exceeds a
    // threshold ignition density are marked as spot fires.
    //
    // Plume lofting height (Albini 1983):
    //   H_z = 12.2 × I_B^(1/3)   [I_B in kW/m, H_z in m]
    //
    // Source flux per unit fire-line length [embers/m/s]:
    //   F_src = k_flux × (I_B / I_B_ref)^flux_exp
    //
    // Flight time and mean landing displacement:
    //   t_f = H_z / terminal_velocity   [s]
    //   dx  = u × t_f,   dy = v × t_f  [m]
    //
    // Gaussian spread (turbulent diffusion from plume):
    //   sigma = sigma_base + k_sigma × H_z   [m]
    //
    // Landing flux density at a cell (x, y):
    //   F_land += F_src × (cell_area / (2π σ²)) × exp(-(Δx² + Δy²) / (2σ²))
    //
    // Ignition:
    //   cells where F_land × dt_check / N_min_density > uniform[0,1]
    //   are ignited (stochastic Poisson arrival with P = 1−exp(−λ), λ =
    //   F_land × dt_check / N_min_density).
    //
    // 3-D wind: when require_3d_wind = 1 the height-averaged horizontal wind
    // from an AMReX massconsistent_amr plt file is used for the mean transport
    // displacement.  An assertion aborts the run if no valid plt file is
    // provided (plt_wind_file must be non-empty).
    //
    // Input prefix: "ember_cascade."
    struct EmberCascadeParams {
        int enable;                   // 1 to enable flux-based cascade (default: 0)
        amrex::Real I_B_min;          // minimum Byram fireline intensity [kW/m] to emit embers (default: 10.0)
        amrex::Real terminal_velocity;// firebrand terminal descent velocity [m/s] (default: 1.0)
        amrex::Real k_flux;           // ember source flux coefficient [embers/m/s] (default: 1.0)
        amrex::Real I_B_ref;          // reference intensity [kW/m] for flux scaling (default: 100.0)
        amrex::Real flux_exp;         // intensity exponent for source flux (default: 1.0)
        amrex::Real sigma_base;       // minimum Gaussian spread radius [m] (default: 50.0)
        amrex::Real k_sigma;          // spread scale factor × H_z [m/m] (default: 0.1)
        amrex::Real n_sigma_cutoff;   // Gaussian kernel truncation: contributions beyond
                                      // n_sigma_cutoff × sigma from the mean are ignored (default: 4.0)
        amrex::Real N_min_density;    // minimum ember density for ignition [embers/m²/s] (default: 1.0e-3)
        amrex::Real spot_radius;      // radius of new spot-fire ignition zone [m] (default: 5.0)
        int check_interval;           // run cascade check every N timesteps (default: 5)
        int random_seed;              // RNG seed; 0 = use system time (default: 0)
        // ---- 3-D wind from massconsistent_amr ----
        // When require_3d_wind = 1 the height-averaged horizontal wind from a
        // massconsistent_amr plt file is used.  If plt_wind_file is empty or
        // the file cannot be read the run aborts with an error message.
        int require_3d_wind;          // 1 = abort if no 3-D plt wind provided (default: 0)
        int use_3d_wind;              // 1 to load 3-D plt wind for transport (default: 0)
        std::string plt_wind_file;    // path to massconsistent_amr plt directory (default: "")
    } ember_cascade;

    // ---- Feature 8: Albini (1979) torching-tree spotting model ----
    // Triggered by passive or active crown fire (crown_activity >= min_crown_activity).
    // Uses flame length from Byram's formula and canopy height to estimate lofting
    // height, then places spot fires at a wind-driven distance.
    // Input prefix: "torching_spotting."
    struct TorchingSpottingParams {
        int enable;              // 1 to enable torching-tree spotting (default: 0)
        amrex::Real k_torch;    // spotting distance constant [m^0.5 / (m/s)]; d = k*U*sqrt(H_eff) (default: 4.24)
        amrex::Real I_B_min;    // minimum Byram fireline intensity [kW/m] for torching (default: 100.0)
        amrex::Real spot_radius;// radius of new spot-fire ignition zone [m] (default: 5.0)
        amrex::Real P_base;     // probability of spot fire per torching cell per check (default: 0.05)
        int random_seed;        // RNG seed (0 = system time, default: 0)
        int check_interval;     // run check every N steps (default: 5)
        int min_crown_activity; // minimum crown_activity to trigger: 1=passive+active, 2=active-only (default: 1)
    } torching_spotting;

    // ---- Feature 9: Persistent per-cell fuel load depletion ----
    // Tracks residual (unconsumed) fuel fraction as a per-cell MultiFab that
    // decays exponentially after the fire front passes.  Written to every
    // plotfile as "residual_fuel".  Optionally feeds back into ROS.
    // Input prefix: "fuel_depletion."
    struct FuelDepletionParams {
        int enable;              // 1 to track per-cell fuel depletion (default: 0)
        amrex::Real tau_burnout; // fine-fuel burnout time constant [s] (default: 3600.0)
        int couple_to_ros;       // 1 to scale ROS by residual fuel in re-entry cells (default: 0)
        // Post-fire fuel adjustment for re-entry spot fires.
        // When adjust_spotting_reentry = 1, firebrand spot fires landing in cells
        // that have been previously burned (arrival_time >= 0) have their
        // catching probability (P_catch) scaled by residual_fuel_mf so that spots
        // landing in fuel-depleted areas are less likely to ignite.
        // Applies to both the basic firebrand spotting and the Albini spotting model.
        int adjust_spotting_reentry; // 1 to scale P_catch by residual fuel (default: 0)
        amrex::Real spotting_fuel_threshold; // residual fuel fraction below which spots cannot ignite (default: 0.05)
    } fuel_depletion;

    // ---- Feature 10: Fire acceleration model (Rothermel 1983 / Catchpole et al. 1992) ----
    // When a fire is small it spreads slower than the quasi-steady-state ROS.
    // The acceleration factor alpha = 1 - exp(-r_fire / L_acc) where r_fire is
    // the effective fire radius (from burned area) and L_acc is the minimum
    // characteristic fire dimension for full ROS.  When fire is large
    // (r_fire >> L_acc) alpha → 1 and there is no modification.
    // Input prefix: "acceleration."
    struct AccelerationParams {
        int enable;              // 1 to enable fire acceleration model (default: 0)
        amrex::Real L_acc;      // acceleration length scale [m]; fire reaches full ROS when r ≈ L_acc (default: 50.0)
        // FARSITE-style temporal acceleration model (McAlpine & Wakimoto 1991)
        int use_temporal;        // 1 to use FARSITE temporal model (default: 0 = use size-based)
        amrex::Real A_point;     // point ignition acceleration constant [1/min] (default: 0.115)
        amrex::Real A_line;      // line ignition acceleration constant [1/min] (default: 0.300)
        amrex::Real perim_limit; // perimeter threshold for point→line switch [m] (default: 402.3 = 20 chains)
        // Wind-onset time-lag model (applies to temporal model only)
        int enable_wind_lag;     // 1 to track wind changes and apply time lag (default: 0)
        amrex::Real tau_wind;    // time constant for wind response [s] (default: 180.0 = 3 min)
    } acceleration;

    // ---- FMC seasonal / phenological schedule ----
    // When fmc_schedule_file is non-empty, the foliar moisture content used by
    // the Van Wagner (1977) crown fire initiation model is updated each timestep
    // from the schedule.  When fmc_use_farsite_curve = 1, a built-in parametric
    // phenological curve calibrated to Southern California chaparral is used
    // instead of a CSV file.
    // Input prefix: "fmc_schedule."
    struct FMCScheduleParams {
        int         enable;              // 1 to enable (default: 0)
        std::string file;                // path to two-column CSV (doy, fmc_pct); empty → use curve
        int         use_farsite_curve;   // 1 to use built-in parametric curve (default: 0)
        int         start_doy;           // day-of-year at t=0 (default: 1 = Jan 1)
        // Parametric curve breakpoints (used when use_farsite_curve = 1)
        int    spring_start;   // doy when green-up begins  (default: 90)
        int    summer_peak;    // doy when FMC reaches max  (default: 150)
        int    fall_start;     // doy when curing begins    (default: 240)
        int    fall_end;       // doy when curing ends      (default: 300)
        amrex::Real fmc_min;   // dormant / cured FMC [%]   (default: 85.0)
        amrex::Real fmc_max;   // peak green FMC [%]        (default: 140.0)
    } fmc_schedule;

    // ---- Precipitation wetting of dead fuels ----
    // When diurnal_moisture.enable = 1 AND precip_rain_rate_mm_hr > 0, dead
    // fuel moisture is additionally driven by a wetting model.  Rain rate can
    // be either a constant (precip_rain_rate_mm_hr) or read from a time series
    // file (precip_schedule_file: two-column CSV of time_s, rain_mm_hr).
    // Input prefix: "diurnal_moisture." (extends the existing struct).
    // New fields:
    amrex::Real precip_rain_rate_mm_hr;   // constant rain rate [mm/hr] (default: 0.0)
    std::string precip_schedule_file;      // time-series CSV of (time_s, rain_mm_hr); overrides constant
    amrex::Real precip_threshold_mm_hr;   // min rain rate to trigger wetting (default: 0.25)
    amrex::Real M_sat;                     // saturation moisture content [fraction] (default: 1.20)

    // ---- Polygon / polyline fire ignition ----
    // When source_type = "polygon" or "polyline", the initial fire perimeter is
    // rasterised from a CSV vertex file.
    // Input prefix: none (uses existing source_type / fire_polygon_file).
    std::string fire_polygon_file;         // path to CSV vertex file (X, Y [, Z])
    amrex::Real polyline_width;            // half-width [m] for polyline ignition (default: 10.0)
    amrex::Real fire_polygon_z_level;      // Z centre of ignition layer in 3-D (default: 0.5)

    // ---- Per-cell live canopy moisture from FARSITE .fms scenario file ----
    // When fms_file is non-empty, a per-cell moisture MultiFab is populated
    // from the .fms scenario and used to override global rothermel moistures
    // on a per-cell basis.  Requires a landscape file to provide per-cell fuel
    // model codes.
    std::string fms_file;                  // path to FARSITE .fms file (default: "")

    // ---- Compact wind direction schedule ----
    // When wind_dir_schedule_file is non-empty, the ambient wind is updated
    // each timestep from a compact (time_s, speed_ms, dir_deg) CSV schedule.
    // This overrides the constant ux/uy and use_time_dependent_wind modes.
    std::string wind_dir_schedule_file;    // path to compact wind schedule CSV (default: "")

    // ---- Solar radiation shading and per-cell shade-adjusted EMC ----
    // Computes per-cell terrain and canopy shading from the solar position
    // (derived from date/time/lat/lon) and the terrain slope/aspect fields.
    // The shade fraction reduces the solar heating of fine fuels, yielding
    // higher EMC in shaded cells and lower EMC in insolated cells.
    //
    // Requires: diurnal_moisture.enable = 1 to supply T_air and RH for the
    //           Simard (1968) EMC equations.  When diurnal_moisture is off,
    //           a warning is issued and the model is disabled.
    //
    // Slope and aspect must come from terrain_file or landscape_file.
    // Canopy shading optionally uses the canopy_cover layer from LCP.
    //
    // Input prefix: "solar_radiation."
    struct SolarRadiationParams {
        int enable;                  // 1 to enable solar shading (default: 0)
        amrex::Real latitude;        // site latitude [deg N; negative = S] (default: 40.0)
        amrex::Real longitude;       // site longitude [deg E; negative = W] (default: -120.0)
        int year;                    // calendar year at simulation t=0 (default: 2024)
        int month;                   // calendar month at t=0 (1-12, default: 7)
        int day;                     // calendar day at t=0 (1-31, default: 1)
        amrex::Real sim_start_hour;  // local clock time [decimal hours] at t=0 (default: 8.0)
        amrex::Real timezone_offset; // UTC offset [hours] (e.g., -8 = PST; default: -8.0)
        amrex::Real solar_heating_C; // max solar heating of fuel above air temp [°C] (default: 17.0)
        int use_canopy_shading;      // 1 to include canopy cover in shade fraction (default: 1)
        amrex::Real cloud_cover;     // domain-uniform cloud fraction [0, 1] (0 = clear sky; default: 0.0)
        // FARSITE full topographic horizon scan (expensive – off by default).
        // For each cell, 8 compass-direction ray marches find the maximum terrain-
        // elevation angle (horizon angle); the cell is shaded when the solar
        // elevation is below that angle.  Requires a global MPI gather of the
        // elevation field – can be costly on large domains or many ranks.
        // Set use_topographic_horizon = 0 to skip this step entirely.
        int         use_topographic_horizon; // 1 = FARSITE 8-dir horizon scan (default: 0)
        amrex::Real horizon_scan_max_dist_m; // max ray-march distance [m]; 0 = full domain (default: 0.0)
    } solar_radiation;

    // ---- Canadian FBP System fire spread model parameters ----
    // Used when fire_spread_model = "fbp_o1a", "fbp_o1b", "fbp_s1", "fbp_s2", or "fbp_s3"
    // Reference: Forestry Canada (1992). FCFDG Information Report ST-X-3.
    //
    // fuel_type : one of "o1a", "o1b", "s1", "s2", "s3"
    // moisture  : dead fine fuel moisture content [%]  (default: 10.0)
    // curing    : degree of curing [0-100 %] for O1a/O1b only  (default: 80.0)
    struct FBPParams {
        std::string fuel_type;  // "o1a", "o1b", "s1", "s2", "s3"  (default: "o1a")
        amrex::Real moisture;   // dead fine fuel moisture [%]   (default: 10.0)
        amrex::Real curing;     // curing percentage [0-100]     (default: 80.0)
    } fbp;

    // ---- Lautenberger (2013) physics-based cellular spread model parameters ----
    // Used when fire_spread_model = "lautenberger"
    // Reference: Lautenberger, C. (2013). Fire Safety Journal 62:289-298.
    //
    // A_L : pre-factor [m^2/s]        (default: 1.05e-5)
    // B_L : moisture sensitivity [-] (default: 2.5)
    // C_L : wind correction [(m/s)^-1] (default: 0.40)
    // D_L : slope correction [-]     (default: 0.50)
    struct LautenbergerParams {
        amrex::Real A_L;   // pre-factor [m^2/s]         (default: 1.05e-5)
        amrex::Real B_L;   // moisture sensitivity [-]   (default: 2.5)
        amrex::Real C_L;   // wind coefficient [(m/s)^-1] (default: 0.40)
        amrex::Real D_L;   // slope coefficient [-]      (default: 0.50)
    } lautenberger;

    // ---- Scott & Reinhardt (2001) full TI/CI computation ----
    // Computes the full Torching Index (TI) and Crowning Index (CI) as minimum
    // 20-ft open wind speeds [km/h] using bisection over the Rothermel formula.
    // Written to plotfiles as "torching_index_kmh" and "crowning_index_kmh".
    struct ScottReinhardtFullParams {
        int enable;              // 1 to enable bisection-based TI/CI (default: 0)
        amrex::Real U_max_kmh;  // maximum wind speed for bisection [km/h] (default: 200.0)
    } scott_reinhardt_full;

    // ---- Live herbaceous moisture / curing schedule ----
    // When enabled, the live herbaceous moisture (M_lh) is updated each
    // timestep from a seasonal schedule.  Two modes:
    //   1. CSV file (herb_moisture_schedule.file): (doy, m_lh_pct) rows.
    //   2. Built-in parametric curing curve (use_curing_curve = 1).
    // Input prefix: "herb_moisture_schedule."
    struct HerbMoistureScheduleParams {
        int         enable;          // 1 to enable (default: 0)
        std::string file;            // path to CSV file (empty → use built-in curve)
        int         use_curing_curve;// 1 to use built-in parametric curve (default: 0)
        int         start_doy;       // day-of-year at t=0 (default: 1 = Jan 1)
        // Parametric curing curve breakpoints
        int    spring_start;   // doy when green-up begins  (default: 90)
        int    summer_peak;    // doy when M_lh reaches max (default: 150)
        int    fall_start;     // doy when curing begins    (default: 200)
        int    fall_end;       // doy when curing ends      (default: 270)
        amrex::Real m_lh_min;  // dormant / cured M_lh [%]  (default: 30.0)
        amrex::Real m_lh_max;  // peak green M_lh [%]       (default: 180.0)
    } herb_moisture_schedule;

    // ---- Live fuel moisture FMC seasonal link ----
    // When fmc_schedule is enabled, live herbaceous and woody moistures are
    // linearly scaled between their summer-peak and winter-dormant values in
    // proportion to the seasonal FMC fraction.  This links the same phenological
    // cycle used for crown fire initiation to live fuel loading in Rothermel.
    // Input prefix: "live_fuel_seasonal."
    struct LiveFuelSeasonalParams {
        int enable;              // 1 to enable seasonal live fuel moisture (default: 0)
        amrex::Real M_lh_summer; // peak-season live herbaceous moisture [fraction] (default: 1.20)
        amrex::Real M_lh_winter; // dormant-season live herbaceous moisture [fraction] (default: 0.30)
        amrex::Real M_lw_summer; // peak-season live woody moisture [fraction] (default: 1.50)
        amrex::Real M_lw_winter; // dormant-season live woody moisture [fraction] (default: 0.60)
    } live_fuel_seasonal;

    // ---- Ground suppression lines (dozer / hand crew) ----
    // CSV file defining time-scheduled suppression lines: polylines that become
    // active firebreaks once their construction is complete.
    // File format (whitespace or comma delimited; # or ! lines are comments):
    //   LINE_START_S  LINE_END_S  LINE_TYPE  N_VERTICES
    //   X1  Y1
    //   X2  Y2
    //   ...                   (polyline vertices, one per line)
    //   [next line]
    //
    // LINE_TYPE: 0 = dozer line, 1 = hand crew line (currently both treated as total firebreak)
    // Cells along the line become firebreaks once current_time >= LINE_END_S.
    // (default: "" = disabled)
    std::string suppression_file;          // path to suppression schedule CSV (default: "")

    // ---- Post-fire fuel model replacement ----
    // CSV file mapping original fuel model codes to post-fire replacement codes.
    // After the fire front passes through a cell (phi flips negative) the
    // fuel model code is replaced using this lookup table.
    // File format: original_code  replacement_code
    // (default: "" = disabled)
    std::string post_fire_fuel_map_file;   // path to post-fire fuel map CSV (default: "")

    // ---- CCFR (Crown fuel consumption rate) active crown fire criterion ----
    // Scott & Reinhardt (2001): active crown fire requires
    //   R_active >= R'_SA × CCFR   where CCFR = (CBD - CBD_min) / CBD_range
    // This replaces the simpler R >= R'_SA criterion from Van Wagner (1977).
    // Input prefix: "crown." (extends existing CrownInitiationParams).
    // New fields:
    amrex::Real ccfr_cbd_min;    // minimum CBD for full passivity [kg/m³] (default: 0.05)
    amrex::Real ccfr_cbd_range;  // CBD range over which CCFR goes from 0 to 1 (default: 0.20)
    int         use_ccfr;        // 1 to enable CCFR active-crown criterion (default: 0)

    // ---- Conditional weather / ERC percentile lookup table ----
    // When enabled, simulation conditions are mapped to wind/moisture scenarios
    // from a user-supplied percentile table.
    // File format (CSV, header row): erc_percentile, wind_speed_ms, wind_dir_deg,
    //                                M_d1_pct, M_d10_pct, M_d100_pct
    // (default: "" = disabled)
    //
    // conditional_weather.trigger_index selects which fire-danger index is used
    // to look up the scenario:
    //   "erc"  – NFDRS Energy Release Component (default; existing behaviour)
    //   "bi"   – NFDRS Burning Index  (ecology_mf component 8)
    //   "sc"   – NFDRS Spread Component (ecology_mf component 7)
    std::string conditional_weather_file;       // path to percentile-weather CSV
    std::string conditional_weather_trigger;    // trigger index: "erc" | "bi" | "sc"

    // ---- Burn-period controls (diurnal active-spread window) ----
    // When enabled, fire spread (R_mf) is gated to a daily active-burn time
    // window [start_hour, end_hour) in local clock hours.  Outside this window
    // R_mf is zeroed so the level-set front does not advance; moisture and all
    // other diagnostics continue to evolve normally.
    //
    // This mimics FARSITE's and FSPro's operational "burn period" concept and
    // is particularly useful for multi-day simulations where overnight moisture
    // recovery dramatically slows or halts spread.
    //
    // The local clock time is computed as:
    //   clock_hour = fmod(sim_start_hour + time_s / 3600.0, 24.0)
    // where sim_start_hour comes from solar_radiation.sim_start_hour (default 8.0)
    // or from the standalone burn_period.sim_start_hour override.
    //
    // Wrapping is handled correctly when the window crosses midnight
    // (e.g., start_hour = 22.0, end_hour = 6.0 = active 22:00–06:00).
    //
    // Input prefix: "burn_period."
    struct BurnPeriodParams {
        int         enable;          // 1 to enable burn-period gating (default: 0)
        amrex::Real start_hour;      // local hour when fire becomes active (default: 10.0 = 10:00 AM)
        amrex::Real end_hour;        // local hour when fire becomes inactive (default: 20.0 = 8:00 PM)
        amrex::Real sim_start_hour;  // local clock hour at simulation t=0 (default: inherit from
                                     //   solar_radiation.sim_start_hour if set, else 0.0)
    } burn_period;

    // ---- Plotfile variable filter ----
    // When non-empty, only the listed variable names are written to each AMReX
    // plotfile.  If the list is empty (default) every field is written.
    //
    // Usage in inputs file (space-separated):
    //   plot_vars = phi R arrival_time fireline_intensity flame_length
    //
    // Available variable names (use these exact strings):
    //   phi                  – level-set function (fire front where phi = 0)
    //   velx, vely [, velz]  – wind velocity components [m/s]
    //   farsite_dx, farsite_dy [, farsite_dz] – FARSITE spread displacements
    //   R                    – rate of spread [m/s]
    //   spot_prob            – firebrand spotting probability
    //   spot_count           – firebrand count per cell
    //   spot_dist            – firebrand landing distance [m]
    //   spot_active          – spotting ignition flag
    //   fuel_consumption     – bulk fuel consumption fraction
    //   crown_fraction       – crown fire fraction
    //   albini_Hz            – Albini lofting height [m]
    //   albini_count         – firebrands launched per cell
    //   albini_dist          – max firebrand landing distance [m]
    //   albini_active        – Albini spot ignition flag
    //   elevation            – terrain elevation [m]
    //   slope                – terrain slope [degrees]
    //   aspect               – terrain aspect [degrees]
    //   fuel_model           – fuel model number (from landscape)
    //   fireline_intensity   – Byram fireline intensity [kW/m]
    //   flame_length         – Byram flame length [m]
    //   fire_intensity_class – Fire intensity class (1-6, Byram 1959)
    //   weise_flame_height   – Weise & Biging flame height [m]
    //   weise_flame_tilt     – Weise & Biging flame tilt angle [rad]
    //   weise_whirl_height   – fire whirl height [m]
    //   weise_whirl_radius   – fire whirl radius [m]
    //   weise_angular_velocity – fire whirl angular velocity [rad/s]
    //   weise_max_tang_vel   – fire whirl max tangential velocity [m/s]
    //   viegas_ROS           – Viegas eruptive rate of spread [m/s]
    //   viegas_eruptive_flag – Viegas eruptive fire flag (0/1)
    //   viegas_ROS_excess    – Viegas ROS excess above threshold
    //   viegas_flame_tilt    – Viegas flame tilt angle [rad]
    //   viegas_slope_factor  – Viegas slope amplification factor
    //   scorch_height        – Van Wagner scorch height [m]
    //   prob_ignition        – Anderson probability of ignition [0-1]
    //   tree_mortality       – Ryan-Reinhardt tree mortality [0-1]
    //   crown_activity       – crown activity class (0=surface,1=passive,2=active)
    //   torching_ratio       – Scott-Reinhardt torching ratio
    //   crowning_ratio       – Scott-Reinhardt crowning ratio
    //   energy_release_component – NFDRS energy release component
    //   nfdrs_spread_component   – NFDRS spread component
    //   nfdrs_burning_index      – NFDRS burning index
    //   co2_emissions        – CO₂ emission rate [kg/m²/s]
    //   co_emissions         – CO emission rate [kg/m²/s]
    //   pm25_emissions       – PM₂.₅ emission rate [kg/m²/s]
    //   arrival_time         – fire arrival time [s]
    //   heat_per_unit_area   – heat per unit area [kJ/m²]
    //   vorticity_z          – vertical vorticity [1/s]
    //   cbh                  – canopy base height [m]
    //   cbd                  – canopy bulk density [kg/m³]
    //   canopy_cover         – canopy cover fraction [0-1]
    //   canopy_height        – canopy height [m]
    //   burnout_time         – cell burnout time [s]
    //   reaction_intensity   – Rothermel reaction intensity [BTU/ft²/min]
    //   residual_fuel        – residual (unconsumed) fuel fraction
    //   shade_fraction       – terrain/canopy shade fraction [0-1]
    //   torching_index_kmh   – Scott-Reinhardt torching index [km/h]
    //   crowning_index_kmh   – Scott-Reinhardt crowning index [km/h]
    //   moisture_d1          – 1-hr dead fuel moisture [fraction]
    //   moisture_d10         – 10-hr dead fuel moisture [fraction]
    //   moisture_d100        – 100-hr dead fuel moisture [fraction]
    //   moisture_lh          – live herbaceous moisture [fraction]
    //   moisture_lw          – live woody moisture [fraction]
    //   plume_rise_m         – Briggs (1965) smoke plume-rise height [m]
    //   fl_exceedance        – maximum observed flame length per cell [m]
    std::vector<std::string> plot_vars;  // subset of variable names to write; empty = write all

    // ---- Smoke plume-rise model (Briggs 1965 / 1969 with Pasquill-Gifford stability) ----
    // Computes the effective smoke plume-rise height Δh [m] above the fire
    // front using the Briggs (1965/1969) buoyancy-dominated formula.
    // Optionally applies Pasquill-Gifford atmospheric stability class corrections
    // to adjust plume rise based on atmospheric conditions (FARSITE-compatible).
    // The result is written to every plotfile as "plume_rise_m".
    // Input prefix: "smoke_plume."
    struct SmokePlumeParams {
        int  enable;                    // 1 to enable plume-rise computation (default: 0)
        amrex::Real T_a;                // ambient temperature [K]  (default: 303.15 = 30 °C)
        amrex::Real rho_a;              // air density [kg/m³]      (default: 1.20)
        amrex::Real Cp_a;               // specific heat [J/(kg·K)] (default: 1005.0)
        
        // Pasquill-Gifford atmospheric stability class (FARSITE-compatible)
        // Classes A-F representing atmospheric stability conditions:
        //   A = extremely unstable (strong solar radiation, light winds)
        //   B = moderately unstable
        //   C = slightly unstable
        //   D = neutral (overcast or moderate wind)
        //   E = slightly stable (clear night, light winds)
        //   F = moderately stable (clear night, very light winds)
        // Default "D" = neutral conditions
        std::string stability_class;    // stability class: "A", "B", "C", "D", "E", or "F" (default: "D")
        int use_stability_correction;   // 1 = apply stability correction to plume rise, 0 = use base Briggs (default: 0)
    } smoke_plume;

    // ---- KML perimeter export ----
    // When write_perimeter_kml = 1, a .kml file is written alongside the CSV
    // and GeoJSON perimeter files at every plotfile step and for isochrones.
    // Coordinates are converted from simulation UTM [m] to WGS-84 lon/lat
    // using the UTM zone number below.  When kml_utm_zone = 0, no conversion
    // is applied and the raw UTM coordinates are written (not valid WGS-84).
    // The KML file can be opened directly in Google Earth.
    int write_perimeter_kml;   // 1 to write perimeter_NNNN.kml (default: 0)
    int kml_utm_zone;          // UTM zone number 1-60; 0 = raw coords (default: 0)
    int kml_utm_northern;      // 1 = northern hemisphere, 0 = southern (default: 1)

    // ---- Simulation start date/time (for log and report display) ----
    // When non-zero, these provide the calendar start date that is combined
    // with solar_radiation.sim_start_hour to display human-readable datetime
    // strings in per-step log output and the HTML fire report.
    // If all three are 0 (default), the values from solar_radiation.year/
    // month/day are used when solar_radiation.enable = 1; otherwise
    // datetime display in logs and report is omitted.
    // Input prefix: "sim_datetime."
    int sim_start_year;        // simulation start calendar year  (default: 0 → auto)
    int sim_start_month;       // simulation start calendar month (default: 0 → auto)
    int sim_start_day;         // simulation start calendar day   (default: 0 → auto)

    // ============================================================================
    // McArthur-style Moisture Response Time Scaling
    // ============================================================================
    // Temperature and humidity-dependent time constant scaling for dead fuel
    // moisture drying. Higher temperatures and lower humidity accelerate drying.
    // Input prefix: "mcarthur_moisture."
    // References: McArthur (1967), Noble et al. (1980)
    struct McArthurMoistureParams {
        int enable;                  // 1 to enable McArthur scaling (default: 0)
        amrex::Real T_ref;           // reference temperature [°C] (default: 20.0)
        amrex::Real k_T;             // temperature sensitivity coefficient (default: 0.05)
        amrex::Real k_RH;            // humidity sensitivity coefficient (default: 0.3)
    } mcarthur_moisture;

    // ============================================================================
    // Enhanced FMC Phenology Models
    // ============================================================================
    // Additional phenology models beyond the FARSITE curve in FMCScheduleParams.
    // Input prefix: "fmc_phenology."
    // References: Jolly et al. (2005), Dennison et al. (2008)
    struct FMCPhenologyParams {
        std::string model;           // "none", "sinusoidal", "gdd" (default: "none")
        amrex::Real fmc_min;         // minimum FMC [%] (default: 85.0)
        amrex::Real fmc_max;         // maximum FMC [%] (default: 140.0)
        amrex::Real doy_offset;      // day-of-year offset for sinusoidal (default: 0.0)
        // GDD parameters
        amrex::Real T_base;          // base temperature for GDD [°C] (default: 5.0)
        amrex::Real GDD_mid;         // midpoint GDD for green-up [°C·day] (default: 500.0)
        amrex::Real k_gdd;           // logistic growth rate (default: 0.01)
    } fmc_phenology;

    // ============================================================================
    // Ember Accumulation Model
    // ============================================================================
    // Tracks accumulated ember density and computes probabilistic ignition.
    // Extends ember_cascade and Albini spotting with persistence/burnout.
    // Input prefix: "ember_accumulation."
    // References: Sardoy et al. (2007), Koo et al. (2010)
    struct EmberAccumulationParams {
        int enable;                  // 1 to enable ember accumulation (default: 0)
        amrex::Real k_decay;         // decay rate [1/s] (default: 1/600 = 10 min burnout)
        amrex::Real rho_threshold;   // minimum density for ignition [embers/m²] (default: 10.0)
        amrex::Real k_ignition;      // ignition rate coefficient (default: 1.0)
        amrex::Real spot_radius;     // radius of ignition zone [m] (default: 5.0)
        int use_moisture_damping;    // 1 to apply fuel moisture damping (default: 1)
    } ember_accumulation;

    // ============================================================================
    // Periodic Wind Gust Model
    // ============================================================================
    // Simple deterministic sinusoidal gust factor for wind variability.
    // Input prefix: "wind_gust."
    // References: Finney et al. (2015), Linn et al. (2012)
    struct WindGustParams {
        int enable;                  // 1 to enable periodic gusts (default: 0)
        amrex::Real amplitude;       // fractional amplitude ±X% (default: 0.2 = ±20%)
        amrex::Real period;          // period of gust cycle [s] (default: 300 = 5 min)
    } wind_gust;

    // ============================================================================
    // Slope-Dependent Flame Tilt Enhancement
    // ============================================================================
    // Enhanced radiation preheating with explicit slope effects on flame tilt.
    // Input prefix: "flame_tilt."
    // References: Anderson (1968), Butler et al. (2004)
    struct FlameTiltParams {
        int enable;                  // 1 to enable slope-dependent tilt (default: 0)
        amrex::Real k_slope;         // slope tilt coefficient (default: 0.5)
        amrex::Real view_factor;     // radiative view factor (default: 0.4)
    } flame_tilt;

    // ---- Real-time satellite fire detection assimilation ----
    // Full configuration is documented in satellite_assimilation.H.
    // Input prefix: "satellite."
    SatelliteParams satellite;

};

void parse_inputs(InputParameters& inputs);

#endif
