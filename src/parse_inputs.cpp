// parse_inputs.cpp
#include "parse_inputs.H"
#include "fuel_database.H"
#include <AMReX_ParmParse.H>

using namespace amrex;

void parse_inputs(InputParameters& p)
{
    ParmParse pp;

    // ---------------- Inputs: grid & domain ----------------
    // You can specify either n_cell (cubic) or n_cell_x/y/z
    int n_cell = 64;            
    pp.query("n_cell", n_cell);
    p.n_cell_x = n_cell;         pp.query("n_cell_x", p.n_cell_x);
    p.n_cell_y = n_cell;         pp.query("n_cell_y", p.n_cell_y);
    p.n_cell_z = n_cell;         pp.query("n_cell_z", p.n_cell_z);

    p.max_grid = 32;             pp.query("max_grid_size", p.max_grid);

    p.plo_x = 0.0;               pp.query("prob_lo_x", p.plo_x);
    p.plo_y = 0.0;               pp.query("prob_lo_y", p.plo_y);
    p.plo_z = 0.0;               pp.query("prob_lo_z", p.plo_z);
    p.phi_x = 1.0;               pp.query("prob_hi_x", p.phi_x);
    p.phi_y = 1.0;               pp.query("prob_hi_y", p.phi_y);
    p.phi_z = 1.0;               pp.query("prob_hi_z", p.phi_z);

    // ---------------- Inputs: time & output ----------------
    p.reinit_int   = 20;         pp.query("reinit_int",   p.reinit_int);
    p.nsteps       = 300;        pp.query("nsteps",       p.nsteps);
    p.final_time   = -1.0;       pp.query("final_time",   p.final_time);
    p.cfl          = 0.5;        pp.query("cfl",          p.cfl);
    p.plot_int     = 50;         pp.query("plot_int",     p.plot_int);

    // ---------------- Inputs: velocity ---------------------
    p.ux = 0.25;                 pp.query("u_x", p.ux);
    p.uy = 0.0;                  pp.query("u_y", p.uy);
    p.uz = 0.0;                  pp.query("u_z", p.uz);
    p.velocity_file = "";        pp.query("velocity_file", p.velocity_file);
    p.use_time_dependent_wind = 0;     pp.query("use_time_dependent_wind", p.use_time_dependent_wind);
    p.wind_time_spacing = 60.0;        pp.query("wind_time_spacing", p.wind_time_spacing);

    // ---------------- Inputs: source selection --------------
    p.source_type = "sphere";    pp.query("source_type", p.source_type);

    // ---------------- Inputs: sphere ------------------------
    p.cx = 0.5;                  pp.query("center_x", p.cx);
    p.cy = 0.5;                  pp.query("center_y", p.cy);
    p.cz = 0.5;                  pp.query("center_z", p.cz);
    p.radius = 0.25;             pp.query("sphere_radius",   p.radius);

    p.xmin = 32;                 pp.query("box_xmin", p.xmin);
    p.ymin = 32;                 pp.query("box_ymin", p.ymin);
    p.zmin = 32;                 pp.query("box_zmin", p.zmin);
    p.xmax = 64;                 pp.query("box_xmax", p.xmax);
    p.ymax = 64;                 pp.query("box_ymax", p.ymax);
    p.zmax = 64;                 pp.query("box_zmax", p.zmax);

    // ---------------- Inputs: ellipse ----------------------
    p.ellipse_center_x = 0.5;    pp.query("ellipse_center_x", p.ellipse_center_x);
    p.ellipse_center_y = 0.5;    pp.query("ellipse_center_y", p.ellipse_center_y);
    p.ellipse_center_z = 0.5;    pp.query("ellipse_center_z", p.ellipse_center_z);
    p.ellipse_radius_x = 0.3;    pp.query("ellipse_radius_x", p.ellipse_radius_x);
    p.ellipse_radius_y = 0.2;    pp.query("ellipse_radius_y", p.ellipse_radius_y);
    p.ellipse_radius_z = 0.15;   pp.query("ellipse_radius_z", p.ellipse_radius_z);

    // ---------------- Inputs: EB implicit function ----------
    p.eb_type = "sphere";        pp.query("eb_type", p.eb_type);
    p.eb_param1 = 0.5;           pp.query("eb_param1", p.eb_param1);
    p.eb_param2 = 0.5;           pp.query("eb_param2", p.eb_param2);
    p.eb_param3 = 0.5;           pp.query("eb_param3", p.eb_param3);
    p.eb_param4 = 0.25;          pp.query("eb_param4", p.eb_param4);
    p.eb_param5 = 0.0;           pp.query("eb_param5", p.eb_param5);
    p.eb_param6 = 0.0;           pp.query("eb_param6", p.eb_param6);

    // -------- Rothermel fire spread model --------
    // Check if user wants to use a fuel model from the database
    std::string fuel_model_name = "";
    pp.query("rothermel.fuel_model", fuel_model_name);
    
    // Set defaults (custom Southern California chaparral from original code)
    p.rothermel.w0        = 0.230;
    p.rothermel.sigma     = 1739.0;
    p.rothermel.delta     = 2.0;
    p.rothermel.M_f       = 0.08;
    p.rothermel.M_x       = 0.30;
    p.rothermel.h_heat    = 8000.0;
    p.rothermel.S_T       = 0.0555;
    p.rothermel.S_e       = 0.010;
    p.rothermel.rho_p     = 32.0;
    p.rothermel.slope_x   = 0.0;
    p.rothermel.slope_y   = 0.0;
    p.rothermel.terrain_file = "";
    p.rothermel.landscape_file = "";
    p.rothermel.wind_conv = 196.85;
    p.rothermel.ros_conv  = 0.00508;
    // Per-class fuel loads initialised to 0 (triggers single-class fallback
    // until a fuel model or explicit overrides provide non-zero values)
    p.rothermel.w_d1    = 0.0;
    p.rothermel.sigma_d1= 0.0;
    p.rothermel.w_d10   = 0.0;
    p.rothermel.w_d100  = 0.0;
    p.rothermel.w_lh    = 0.0;
    p.rothermel.sigma_lh= 0.0;
    p.rothermel.w_lw    = 0.0;
    p.rothermel.sigma_lw= 0.0;
    // Per-class moistures are set after M_f has its final value (see below)
    
    // Apply fuel model from database if specified
    if (!fuel_model_name.empty()) {
        FuelModel model;
        if (lookup_fuel_model(fuel_model_name, model)) {
            Print() << "Using fuel model: " << model.name << " - " << model.description << "\n";
            apply_fuel_model(model, p.rothermel);
            print_fuel_model_info(model);
        } else {
            Print() << "WARNING: Fuel model '" << fuel_model_name << "' not found in database.\n";
            Print() << "         Using default values (custom Southern California chaparral).\n";
            print_available_fuel_models();
        }
    }
    
    // Allow individual parameter overrides (these take precedence over fuel model)
    pp.query("rothermel.w0",        p.rothermel.w0);
    pp.query("rothermel.sigma",     p.rothermel.sigma);
    pp.query("rothermel.delta",     p.rothermel.delta);
    pp.query("rothermel.M_f",       p.rothermel.M_f);
    pp.query("rothermel.M_x",       p.rothermel.M_x);
    pp.query("rothermel.h_heat",    p.rothermel.h_heat);
    pp.query("rothermel.S_T",       p.rothermel.S_T);
    pp.query("rothermel.S_e",       p.rothermel.S_e);
    pp.query("rothermel.rho_p",     p.rothermel.rho_p);
    pp.query("rothermel.slope_x",   p.rothermel.slope_x);
    pp.query("rothermel.slope_y",   p.rothermel.slope_y);
    pp.query("rothermel.terrain_file", p.rothermel.terrain_file);
    pp.query("rothermel.landscape_file", p.rothermel.landscape_file);
    p.rothermel.landscape_fuel_type = "13";
    pp.query("rothermel.landscape_fuel_type", p.rothermel.landscape_fuel_type);
    if (p.rothermel.landscape_fuel_type != "13" && p.rothermel.landscape_fuel_type != "40") {
        amrex::Abort("rothermel.landscape_fuel_type must be '13' (FBFM13) or '40' (FBFM40)");
    }
    if (!p.rothermel.landscape_file.empty()) {
        Print() << "Landscape fuel model system: "
                << (p.rothermel.landscape_fuel_type == "40"
                    ? "FBFM40 (Scott & Burgan 40)"
                    : "FBFM13 (Anderson 13)")
                << " (rothermel.landscape_fuel_type = "
                << p.rothermel.landscape_fuel_type << ")\n";
    }
    pp.query("rothermel.wind_conv", p.rothermel.wind_conv);
    pp.query("rothermel.ros_conv",  p.rothermel.ros_conv);
    p.rothermel.use_waf         = 0;    pp.query("rothermel.use_waf",         p.rothermel.use_waf);
    p.rothermel.use_wind_limit  = 0;    pp.query("rothermel.use_wind_limit",  p.rothermel.use_wind_limit);

    // Per-class fuel load overrides (take precedence over fuel model database)
    pp.query("rothermel.w_d1",    p.rothermel.w_d1);
    pp.query("rothermel.sigma_d1", p.rothermel.sigma_d1);
    pp.query("rothermel.w_d10",   p.rothermel.w_d10);
    pp.query("rothermel.w_d100",  p.rothermel.w_d100);
    pp.query("rothermel.w_lh",    p.rothermel.w_lh);
    pp.query("rothermel.sigma_lh", p.rothermel.sigma_lh);
    pp.query("rothermel.w_lw",    p.rothermel.w_lw);
    pp.query("rothermel.sigma_lw",p.rothermel.sigma_lw);

    // Per-class fuel moistures — defaults derived from M_f (which now has its
    // final value) so that the multi-class path degrades gracefully to the
    // single-class result when only rothermel.M_f is specified.
    p.rothermel.M_d1   = p.rothermel.M_f;
    p.rothermel.M_d10  = p.rothermel.M_f;
    p.rothermel.M_d100 = p.rothermel.M_f;
    p.rothermel.M_lh   = 0.90;   // live herbaceous moisture default
    p.rothermel.M_lw   = 1.20;   // live woody moisture default
    pp.query("rothermel.M_d1",   p.rothermel.M_d1);
    pp.query("rothermel.M_d10",  p.rothermel.M_d10);
    pp.query("rothermel.M_d100", p.rothermel.M_d100);
    pp.query("rothermel.M_lh",   p.rothermel.M_lh);
    pp.query("rothermel.M_lw",   p.rothermel.M_lw);

    // -------- FARSITE ellipse model parameters (Richards 1990) --------
    p.farsite.use_anderson_LW = 0;               pp.query("farsite.use_anderson_LW", p.farsite.use_anderson_LW);
    p.farsite.length_to_width_ratio = 3.0;       pp.query("farsite.length_to_width_ratio", p.farsite.length_to_width_ratio);
    p.farsite.phi_threshold = 0.1;               pp.query("farsite.phi_threshold", p.farsite.phi_threshold);
    p.farsite.coeff_a = 1.0;                     pp.query("farsite.coeff_a", p.farsite.coeff_a);
    p.farsite.coeff_b = 0.5;                     pp.query("farsite.coeff_b", p.farsite.coeff_b);
    p.farsite.coeff_c = 0.2;                     pp.query("farsite.coeff_c", p.farsite.coeff_c);

    // Fire shape model (only used when propagation_method = farsite)
    p.farsite.fire_shape_model = "richards";     pp.query("farsite.fire_shape_model", p.farsite.fire_shape_model);
    if (p.farsite.fire_shape_model != "richards" &&
        p.farsite.fire_shape_model != "catchpole_demestre" &&
        p.farsite.fire_shape_model != "wilson" &&
        p.farsite.fire_shape_model != "lemniscate") {
        amrex::Abort("farsite.fire_shape_model must be 'richards', 'catchpole_demestre', "
                     "'wilson', or 'lemniscate'");
    }
    Print() << "FARSITE fire shape model: " << p.farsite.fire_shape_model << "\n";
    
    // -------- Bulk Fuel Consumption Fraction Model parameters --------
    p.farsite.use_bulk_fuel_consumption = 0;     pp.query("farsite.use_bulk_fuel_consumption", p.farsite.use_bulk_fuel_consumption);
    p.farsite.tau_residence = 60.0;              pp.query("farsite.tau_residence", p.farsite.tau_residence);
    p.farsite.f_consumed_max = 0.9;              pp.query("farsite.f_consumed_max", p.farsite.f_consumed_max);
    p.farsite.f_consumed_min = 0.5;              pp.query("farsite.f_consumed_min", p.farsite.f_consumed_min);
    
    // Validate bulk fuel consumption parameters
    if (p.farsite.use_bulk_fuel_consumption == 1) {
        if (p.farsite.f_consumed_max < 0.0 || p.farsite.f_consumed_max > 1.0) {
            amrex::Abort("farsite.f_consumed_max must be between 0.0 and 1.0");
        }
        if (p.farsite.f_consumed_min < 0.0 || p.farsite.f_consumed_min > 1.0) {
            amrex::Abort("farsite.f_consumed_min must be between 0.0 and 1.0");
        }
        if (p.farsite.f_consumed_min > p.farsite.f_consumed_max) {
            amrex::Abort("farsite.f_consumed_min must be less than or equal to farsite.f_consumed_max");
        }
        if (p.farsite.tau_residence <= 0.0) {
            amrex::Abort("farsite.tau_residence must be greater than 0");
        }
    }

    // -------- Firebrand spotting model parameters --------
    p.spotting.enable = 0;                       pp.query("spotting.enable", p.spotting.enable);
    p.spotting.P_base = 0.02;                    pp.query("spotting.P_base", p.spotting.P_base);
    p.spotting.k_wind = 0.3;                     pp.query("spotting.k_wind", p.spotting.k_wind);
    p.spotting.I_critical = 1000.0;              pp.query("spotting.I_critical", p.spotting.I_critical);
    p.spotting.d_mean = 0.1;                     pp.query("spotting.d_mean", p.spotting.d_mean);
    p.spotting.d_sigma = 0.5;                    pp.query("spotting.d_sigma", p.spotting.d_sigma);
    p.spotting.d_lambda = 10.0;                  pp.query("spotting.d_lambda", p.spotting.d_lambda);
    p.spotting.distance_model = "lognormal";     pp.query("spotting.distance_model", p.spotting.distance_model);
    p.spotting.lateral_spread_angle = 15.0;      pp.query("spotting.lateral_spread_angle", p.spotting.lateral_spread_angle);
    p.spotting.spot_radius = 0.02;               pp.query("spotting.spot_radius", p.spotting.spot_radius);
    p.spotting.random_seed = 0;                  pp.query("spotting.random_seed", p.spotting.random_seed);
    p.spotting.check_interval = 5;               pp.query("spotting.check_interval", p.spotting.check_interval);
    
    // Validate spotting parameters
    if (p.spotting.enable == 1) {
        if (p.spotting.P_base < 0.0 || p.spotting.P_base > 1.0) {
            amrex::Abort("spotting.P_base must be between 0.0 and 1.0");
        }
        if (p.spotting.check_interval < 1) {
            amrex::Abort("spotting.check_interval must be at least 1");
        }
        if (p.spotting.distance_model != "lognormal" && p.spotting.distance_model != "exponential") {
            amrex::Abort("spotting.distance_model must be either 'lognormal' or 'exponential'");
        }
    }

    // -------- Van Wagner crown fire initiation model parameters --------
    p.crown.enable = 0;                          pp.query("crown.enable", p.crown.enable);
    p.crown.CBH = 4.0;                           pp.query("crown.CBH", p.crown.CBH);
    p.crown.CBD = 0.15;                          pp.query("crown.CBD", p.crown.CBD);
    p.crown.FMC = 100.0;                         pp.query("crown.FMC", p.crown.FMC);
    p.crown.crown_fraction_weight = 1.0;         pp.query("crown.crown_fraction_weight", p.crown.crown_fraction_weight);
    p.crown.use_metric_units = 1;                pp.query("crown.use_metric_units", p.crown.use_metric_units);
    
    // Validate crown fire parameters
    if (p.crown.enable == 1) {
        if (p.crown.CBH <= 0.0) {
            amrex::Abort("crown.CBH (canopy base height) must be greater than 0");
        }
        if (p.crown.CBD <= 0.0) {
            amrex::Abort("crown.CBD (canopy bulk density) must be greater than 0");
        }
        if (p.crown.FMC < 50.0 || p.crown.FMC > 300.0) {
            amrex::Abort("crown.FMC (foliar moisture content) must be between 50% and 300%");
        }
        if (p.crown.crown_fraction_weight < 0.0 || p.crown.crown_fraction_weight > 2.0) {
            amrex::Abort("crown.crown_fraction_weight must be between 0.0 and 2.0");
        }
    }

    // -------- Fire spread model selection --------
    // Reads "fire_spread_model" key ("rothermel" or "balbi").
    // Legacy key "balbi.enable = 1" is also accepted for backward compatibility.
    p.fire_spread_model = "rothermel";
    pp.query("fire_spread_model", p.fire_spread_model);
    {
        int balbi_enable_legacy = 0;
        pp.query("balbi.enable", balbi_enable_legacy);
        if (balbi_enable_legacy == 1) {
            p.fire_spread_model = "balbi";
            Print() << "NOTE: balbi.enable is deprecated; use fire_spread_model = balbi\n";
        }
    }
    if (p.fire_spread_model != "rothermel" && p.fire_spread_model != "balbi"
        && p.fire_spread_model != "cheney_gould"
        && p.fire_spread_model != "cruz_crown") {
        amrex::Abort("fire_spread_model must be 'rothermel', 'balbi', 'cheney_gould', or 'cruz_crown'");
    }

    // -------- Propagation method selection --------
    // Reads "propagation_method" key ("levelset" or "farsite").
    // Legacy keys "skip_levelset" and "farsite.enable" are also accepted.
    p.propagation_method = "levelset";
    pp.query("propagation_method", p.propagation_method);
    {
        int skip_ls_legacy = 0;
        int farsite_en_legacy = -1;
        pp.query("skip_levelset", skip_ls_legacy);
        pp.query("farsite.enable", farsite_en_legacy);
        if (skip_ls_legacy == 1) {
            p.propagation_method = "farsite";
            Print() << "NOTE: skip_levelset is deprecated; use propagation_method = farsite\n";
        } else if (skip_ls_legacy == 0 && farsite_en_legacy == 0) {
            // Explicitly disabled FARSITE → levelset
            p.propagation_method = "levelset";
        }
        p.skip_levelset = (p.propagation_method == "farsite" || p.propagation_method == "mtt") ? 1 : 0;
    }
    if (p.propagation_method != "levelset" &&
        p.propagation_method != "farsite"  &&
        p.propagation_method != "mtt") {
        amrex::Abort("propagation_method must be 'levelset', 'farsite', or 'mtt'");
    }
    Print() << "Propagation method: " << p.propagation_method << "\n";

    // -------- Cheney & Gould (1995) grassland fire spread model --------
    p.cheney_gould.moisture = 10.0;  pp.query("cheney_gould.moisture", p.cheney_gould.moisture);
    p.cheney_gould.curing   = 1.0;   pp.query("cheney_gould.curing",   p.cheney_gould.curing);

    // -------- Cruz, Alexander & Wakimoto (2005) crown fire spread model --------
    p.cruz_crown.CBD  = 0.10;   pp.query("cruz_crown.CBD",  p.cruz_crown.CBD);
    p.cruz_crown.MC10 = 10.0;   pp.query("cruz_crown.MC10", p.cruz_crown.MC10);

    if (p.fire_spread_model == "cruz_crown") {
        if (p.cruz_crown.CBD <= 0.0)
            amrex::Abort("cruz_crown.CBD (canopy bulk density) must be > 0 kg/m3");
        if (p.cruz_crown.MC10 < 0.0)
            amrex::Abort("cruz_crown.MC10 (10-h fuel moisture) must be >= 0%");
        Print() << "Fire spread model: Cruz, Alexander & Wakimoto (2005) crown fire\n";
        Print() << "  CBD=" << p.cruz_crown.CBD << " kg/m3"
                << "  MC10=" << p.cruz_crown.MC10 << " %\n";
    }

    // -------- Balbi (2009) physical fire spread model --------
    p.balbi.T_a      = 300.0;      pp.query("balbi.T_a",      p.balbi.T_a);
    p.balbi.T_f      = 1000.0;     pp.query("balbi.T_f",      p.balbi.T_f);
    p.balbi.T_i      = 600.0;      pp.query("balbi.T_i",      p.balbi.T_i);
    p.balbi.delta_H  = 2.26e6;     pp.query("balbi.delta_H",  p.balbi.delta_H);
    p.balbi.C_pf     = 1800.0;     pp.query("balbi.C_pf",     p.balbi.C_pf);
    p.balbi.r_00     = 2.5e-4;     pp.query("balbi.r_00",     p.balbi.r_00);
    p.balbi.tau_0    = 75591.0;    pp.query("balbi.tau_0",    p.balbi.tau_0);

    if (p.fire_spread_model == "balbi") {
        if (p.balbi.T_a <= 0.0)
            amrex::Abort("balbi.T_a must be > 0 K");
        if (p.balbi.T_f <= p.balbi.T_a)
            amrex::Abort("balbi.T_f must be > balbi.T_a");
        if (p.balbi.T_i <= p.balbi.T_a)
            amrex::Abort("balbi.T_i must be > balbi.T_a");
        if (p.balbi.tau_0 <= 0.0)
            amrex::Abort("balbi.tau_0 must be > 0");
        if (p.balbi.r_00 <= 0.0)
            amrex::Abort("balbi.r_00 must be > 0");
        Print() << "Fire spread model: Balbi (2009)\n";
        Print() << "  T_a=" << p.balbi.T_a << " K  T_f=" << p.balbi.T_f
                << " K  T_i=" << p.balbi.T_i << " K\n";
        Print() << "  C_pf=" << p.balbi.C_pf << " J/(kg·K)"
                << "  r_00=" << p.balbi.r_00 << " m"
                << "  tau_0=" << p.balbi.tau_0 << " s/m\n";
    } else if (p.fire_spread_model == "cheney_gould") {
        Print() << "Fire spread model: Cheney & Gould (1995) grassland\n";
        Print() << "  moisture=" << p.cheney_gould.moisture
                << " %  curing=" << p.cheney_gould.curing << "\n";
    } else {
        Print() << "Fire spread model: Rothermel (1972)\n";
    }

    // Print Andrews (2018) wind adjustment settings
    if (p.rothermel.use_waf == 1 || p.rothermel.use_wind_limit == 1) {
        Print() << "Andrews (2018) wind adjustments enabled:";
        if (p.rothermel.use_waf == 1)
            Print() << " WAF (20-ft→midflame)";
        if (p.rothermel.use_wind_limit == 1)
            Print() << " MEWS-limit";
        Print() << "\n";
    }

    // -------- Albini (1983) firebrand spotting with 2-D trajectory --------
    p.albini_spotting.enable             = 0;        pp.query("albini_spotting.enable",             p.albini_spotting.enable);
    p.albini_spotting.terminal_velocity  = 1.0;      pp.query("albini_spotting.terminal_velocity",  p.albini_spotting.terminal_velocity);
    p.albini_spotting.P_base             = 0.01;     pp.query("albini_spotting.P_base",             p.albini_spotting.P_base);
    p.albini_spotting.I_B_min            = 10.0;     pp.query("albini_spotting.I_B_min",            p.albini_spotting.I_B_min);
    p.albini_spotting.spot_radius        = 5.0;      pp.query("albini_spotting.spot_radius",        p.albini_spotting.spot_radius);
    p.albini_spotting.random_seed        = 0;        pp.query("albini_spotting.random_seed",        p.albini_spotting.random_seed);
    p.albini_spotting.check_interval     = 5;        pp.query("albini_spotting.check_interval",     p.albini_spotting.check_interval);
    p.albini_spotting.n_traj_steps       = 100;      pp.query("albini_spotting.n_traj_steps",       p.albini_spotting.n_traj_steps);

    // Validate Albini spotting parameters
    if (p.albini_spotting.enable == 1) {
        if (p.albini_spotting.terminal_velocity <= 0.0) {
            amrex::Abort("albini_spotting.terminal_velocity must be greater than 0");
        }
        if (p.albini_spotting.P_base < 0.0 || p.albini_spotting.P_base > 1.0) {
            amrex::Abort("albini_spotting.P_base must be between 0.0 and 1.0");
        }
        if (p.albini_spotting.check_interval < 1) {
            amrex::Abort("albini_spotting.check_interval must be at least 1");
        }
        if (p.albini_spotting.n_traj_steps < 1) {
            amrex::Abort("albini_spotting.n_traj_steps must be at least 1");
        }
    }

    // -------- Weise & Biging (1996) fire whirl model --------
    p.weise_biging.enable  = 0;     pp.query("weise_biging.enable",  p.weise_biging.enable);
    p.weise_biging.c_r     = 0.1;   pp.query("weise_biging.c_r",     p.weise_biging.c_r);
    p.weise_biging.I_B_min = 1.0;   pp.query("weise_biging.I_B_min", p.weise_biging.I_B_min);

    if (p.weise_biging.enable == 1) {
        if (p.weise_biging.c_r <= 0.0)
            amrex::Abort("weise_biging.c_r must be > 0");
        if (p.weise_biging.I_B_min < 0.0)
            amrex::Abort("weise_biging.I_B_min must be >= 0");
        Print() << "Weise & Biging (1996) fire whirl model enabled:\n";
        Print() << "  c_r=" << p.weise_biging.c_r
                << "  I_B_min=" << p.weise_biging.I_B_min << " kW/m\n";
    }

    // -------- Viegas (2004) eruptive fire model --------
    p.viegas.enable     = 0;      pp.query("viegas.enable",     p.viegas.enable);
    p.viegas.a_V        = 1.83;   pp.query("viegas.a_V",        p.viegas.a_V);
    p.viegas.tan_phi_c  = 0.4;    pp.query("viegas.tan_phi_c",  p.viegas.tan_phi_c);
    p.viegas.T_a        = 300.0;  pp.query("viegas.T_a",        p.viegas.T_a);
    p.viegas.T_f        = 1000.0; pp.query("viegas.T_f",        p.viegas.T_f);

    // -------- Turbulent wind perturbation model --------
    p.turb_wind.model       = "none"; pp.query("turb_wind.model",       p.turb_wind.model);
    p.turb_wind.theta       = 0.1;    pp.query("turb_wind.theta",       p.turb_wind.theta);
    p.turb_wind.sigma       = 0.5;    pp.query("turb_wind.sigma",       p.turb_wind.sigma);
    p.turb_wind.L_c         = 0.0;    pp.query("turb_wind.L_c",         p.turb_wind.L_c);
    p.turb_wind.N_modes     = 32;     pp.query("turb_wind.N_modes",     p.turb_wind.N_modes);
    p.turb_wind.sigma_theta = 0.1;    pp.query("turb_wind.sigma_theta", p.turb_wind.sigma_theta);
    p.turb_wind.theta_max   = 0.5236; pp.query("turb_wind.theta_max",   p.turb_wind.theta_max);
    p.turb_wind.random_seed = 0;      pp.query("turb_wind.random_seed", p.turb_wind.random_seed);

    // Validate
    if (p.turb_wind.model != "none"             &&
        p.turb_wind.model != "ou_process"       &&
        p.turb_wind.model != "spectral_noise"   &&
        p.turb_wind.model != "direction_walk") {
        amrex::Abort("turb_wind.model must be one of: "
                     "none, ou_process, spectral_noise, direction_walk");
    }
    if (p.turb_wind.model == "ou_process") {
        if (p.turb_wind.theta <= 0.0)
            amrex::Abort("turb_wind.theta must be > 0");
        if (p.turb_wind.sigma <= 0.0)
            amrex::Abort("turb_wind.sigma must be > 0");
        if (p.turb_wind.L_c < 0.0)
            amrex::Abort("turb_wind.L_c must be >= 0");
        Print() << "Turbulent wind model: ou_process\n";
        Print() << "  theta=" << p.turb_wind.theta << " s^-1"
                << "  sigma=" << p.turb_wind.sigma << " m/s";
        if (p.turb_wind.L_c > 0.0)
            Print() << "  L_c=" << p.turb_wind.L_c << " m (spatially correlated)";
        else
            Print() << "  L_c=0 (domain-uniform)";
        Print() << "\n";
    }
    if (p.turb_wind.model == "spectral_noise") {
        if (p.turb_wind.theta <= 0.0)
            amrex::Abort("turb_wind.theta must be > 0");
        if (p.turb_wind.sigma <= 0.0)
            amrex::Abort("turb_wind.sigma must be > 0");
        if (p.turb_wind.L_c <= 0.0)
            amrex::Abort("turb_wind.L_c must be > 0 for spectral_noise");
        if (p.turb_wind.N_modes < 1)
            amrex::Abort("turb_wind.N_modes must be >= 1");
        Print() << "Turbulent wind model: spectral_noise (Random Fourier Features)\n";
        Print() << "  theta=" << p.turb_wind.theta << " s^-1"
                << "  sigma=" << p.turb_wind.sigma << " m/s"
                << "  L_c="   << p.turb_wind.L_c   << " m"
                << "  N_modes=" << p.turb_wind.N_modes << "\n";
    }
    if (p.turb_wind.model == "direction_walk") {
        if (p.turb_wind.sigma_theta <= 0.0)
            amrex::Abort("turb_wind.sigma_theta must be > 0");
        if (p.turb_wind.theta_max <= 0.0)
            amrex::Abort("turb_wind.theta_max must be > 0");
        Print() << "Turbulent wind model: direction_walk\n";
        Print() << "  sigma_theta=" << p.turb_wind.sigma_theta << " rad/step"
                << "  theta_max=" << p.turb_wind.theta_max << " rad ("
                << p.turb_wind.theta_max * (180.0 / 3.14159265) << " deg)\n";
    }

    // -------- Wind-terrain feedback model for Rothermel wind enhancement --------
    // Parsed before the viegas validation block so auto-enable can set viegas.enable.
    p.wind_terrain.model       = "none";  pp.query("wind_terrain.model",       p.wind_terrain.model);
    p.wind_terrain.k_canyon    = 1.0;    pp.query("wind_terrain.k_canyon",    p.wind_terrain.k_canyon);
    p.wind_terrain.k_pimont    = 0.5;    pp.query("wind_terrain.k_pimont",    p.wind_terrain.k_pimont);
    p.wind_terrain.k_ridge     = 1.0;    pp.query("wind_terrain.k_ridge",     p.wind_terrain.k_ridge);
    p.wind_terrain.k_canyon_wn = 0.5;    pp.query("wind_terrain.k_canyon_wn", p.wind_terrain.k_canyon_wn);

    // Validate model name
    if (p.wind_terrain.model != "none"                   &&
        p.wind_terrain.model != "viegas_ros"             &&
        p.wind_terrain.model != "viegas_wind"            &&
        p.wind_terrain.model != "canyon_wind"            &&
        p.wind_terrain.model != "viegas_neto"            &&
        p.wind_terrain.model != "pimont"                 &&
        p.wind_terrain.model != "windninja_ridge_canyon") {
        amrex::Abort("wind_terrain.model must be one of: "
                     "none, viegas_ros, viegas_wind, canyon_wind, viegas_neto, pimont, "
                     "windninja_ridge_canyon");
    }

    // Validate model-specific parameters
    if (p.wind_terrain.model == "canyon_wind" && p.wind_terrain.k_canyon <= 0.0) {
        amrex::Abort("wind_terrain.k_canyon must be > 0");
    }
    if (p.wind_terrain.model == "pimont" && p.wind_terrain.k_pimont <= 0.0) {
        amrex::Abort("wind_terrain.k_pimont must be > 0");
    }
    if (p.wind_terrain.model == "windninja_ridge_canyon") {
        if (p.wind_terrain.k_ridge <= 0.0)
            amrex::Abort("wind_terrain.k_ridge must be > 0");
        if (p.wind_terrain.k_canyon_wn <= 0.0)
            amrex::Abort("wind_terrain.k_canyon_wn must be > 0");
    }

    // Auto-enable Viegas diagnostics for Viegas-based wind-terrain models
    if (p.wind_terrain.model == "viegas_ros"  ||
        p.wind_terrain.model == "viegas_wind" ||
        p.wind_terrain.model == "viegas_neto") {
        p.viegas.enable = 1;
    }

    if (p.viegas.enable == 1) {
        if (p.viegas.a_V <= 0.0)
            amrex::Abort("viegas.a_V must be > 0");
        if (p.viegas.tan_phi_c < 0.0)
            amrex::Abort("viegas.tan_phi_c must be >= 0");
        if (p.viegas.T_a <= 0.0)
            amrex::Abort("viegas.T_a must be > 0 K");
        if (p.viegas.T_f <= p.viegas.T_a)
            amrex::Abort("viegas.T_f must be > viegas.T_a");
        Print() << "Viegas (2004) eruptive fire diagnostics enabled:\n";
        Print() << "  a_V=" << p.viegas.a_V << " (dimensionless)"
                << "  tan_phi_c=" << p.viegas.tan_phi_c << " (dimensionless, ~"
                << static_cast<int>(std::atan(p.viegas.tan_phi_c) * 180.0 / 3.14159265) << " deg)"
                << "  T_a=" << p.viegas.T_a << " K"
                << "  T_f=" << p.viegas.T_f << " K\n";
    }

    // Print wind-terrain model info
    if (p.wind_terrain.model == "none") {
        Print() << "Wind-terrain model: none (Option 1 – default Rothermel)\n";
    } else if (p.wind_terrain.model == "viegas_ros") {
        Print() << "Wind-terrain model: viegas_ros (Option 2 – Viegas ROS as spread rate)\n";
        Print() << "  R_final = max(R_rothermel, R_viegas) in eruptive cells"
                << " (tan_phi > " << p.viegas.tan_phi_c << ")\n";
    } else if (p.wind_terrain.model == "viegas_wind") {
        Print() << "Wind-terrain model: viegas_wind (Option 3 – Viegas-induced upslope wind, eruptive cells)\n";
        Print() << "  v_b = sqrt(g * delta * (T_f - T_a) / T_a) added to ambient wind"
                << " where tan_phi > " << p.viegas.tan_phi_c << "\n";
    } else if (p.wind_terrain.model == "canyon_wind") {
        Print() << "Wind-terrain model: canyon_wind (Option 4 – Rothermel 1983 canyon wind)\n";
        Print() << "  U_eff = U * (1 + k_canyon * tan_phi),  k_canyon = "
                << p.wind_terrain.k_canyon << "\n";
    } else if (p.wind_terrain.model == "viegas_neto") {
        Print() << "Wind-terrain model: viegas_neto (Option 5 – Viegas & Neto 1994 buoyancy wind)\n";
        Print() << "  U_ind = v_b * tan_phi added upslope at all cells\n";
    } else if (p.wind_terrain.model == "pimont") {
        Print() << "Wind-terrain model: pimont (Option 6 – Pimont et al. 2009 slope correction)\n";
        Print() << "  U_eff = U * exp(k_pimont * tan_phi),  k_pimont = "
                << p.wind_terrain.k_pimont << "\n";
    } else if (p.wind_terrain.model == "windninja_ridge_canyon") {
        Print() << "Wind-terrain model: windninja_ridge_canyon (Option 7 – WindNinja ridge/canyon)\n";
        Print() << "  Ridge  (wind upslope):   f = 1 + k_ridge * tan_phi * alignment,"
                << "  k_ridge = " << p.wind_terrain.k_ridge << "\n";
        Print() << "  Canyon (wind downslope): f = 1 + k_canyon_wn * tan_phi * |alignment|,"
                << "  k_canyon_wn = " << p.wind_terrain.k_canyon_wn << "\n";
    }

    // -------- Heat flux MultiFab parameters --------
    p.heat_flux.heat_flux_value = 0.0;    pp.query("heat_flux.value",          p.heat_flux.heat_flux_value);
    p.heat_flux.heat_flux_file  = "";     pp.query("heat_flux.file",           p.heat_flux.heat_flux_file);
    p.heat_flux.rho_air         = 1.2;    pp.query("heat_flux.rho_air",        p.heat_flux.rho_air);
    p.heat_flux.Cp_air          = 1005.0; pp.query("heat_flux.Cp_air",         p.heat_flux.Cp_air);
    p.heat_flux.T_a             = 300.0;  pp.query("heat_flux.T_a",            p.heat_flux.T_a);
    p.heat_flux.ref_height      = 10.0;   pp.query("heat_flux.ref_height",     p.heat_flux.ref_height);
    p.heat_flux.k_upward        = 1.0;    pp.query("heat_flux.k_upward",       p.heat_flux.k_upward);
    p.heat_flux.k_induced       = 0.5;    pp.query("heat_flux.k_induced",      p.heat_flux.k_induced);
    p.heat_flux.enable_upward   = 0;      pp.query("heat_flux.enable_upward",  p.heat_flux.enable_upward);
    p.heat_flux.enable_induced  = 0;      pp.query("heat_flux.enable_induced", p.heat_flux.enable_induced);

    // Validate heat flux parameters
    if (p.heat_flux.enable_upward == 1 || p.heat_flux.enable_induced == 1) {
        if (p.heat_flux.rho_air <= 0.0)
            amrex::Abort("heat_flux.rho_air must be > 0");
        if (p.heat_flux.Cp_air <= 0.0)
            amrex::Abort("heat_flux.Cp_air must be > 0");
        if (p.heat_flux.T_a <= 0.0)
            amrex::Abort("heat_flux.T_a must be > 0 K");
        if (p.heat_flux.ref_height <= 0.0)
            amrex::Abort("heat_flux.ref_height must be > 0");
        const bool has_hf = (p.heat_flux.heat_flux_value > 0.0 ||
                             !p.heat_flux.heat_flux_file.empty());
        if (!has_hf)
            amrex::Abort("heat_flux.enable_upward/enable_induced requires "
                         "heat_flux.value > 0 or heat_flux.file to be set");
        Print() << "Heat flux wind model enabled:\n";
        if (!p.heat_flux.heat_flux_file.empty())
            Print() << "  Q from file: " << p.heat_flux.heat_flux_file << "\n";
        else
            Print() << "  Q = " << p.heat_flux.heat_flux_value << " W/m²  (uniform)\n";
        if (p.heat_flux.enable_upward == 1)
            Print() << "  Upward velocity term: k_upward = " << p.heat_flux.k_upward
                    << "  ref_height = " << p.heat_flux.ref_height << " m\n";
        if (p.heat_flux.enable_induced == 1)
            Print() << "  Induced inflow term: k_induced = " << p.heat_flux.k_induced << "\n";
    }

    // -------- CSV fire points initialization --------
    p.fire_points_file  = "";                    pp.query("fire_points_file",     p.fire_points_file);
    p.fire_gaussian_sigma = -1.0;                pp.query("fire_gaussian_sigma",  p.fire_gaussian_sigma);

    // -------- Dynamic fire points file (polled each time step) --------
    p.dynamic_fire_points_file = "";             pp.query("dynamic_fire_points_file", p.dynamic_fire_points_file);

    // -------- Checkpoint / restart options --------
    p.chk_int        = -1;                       pp.query("chk_int",         p.chk_int);
    p.restart_chkfile = "";                      pp.query("restart_chkfile", p.restart_chkfile);

    // -------- FARSITE fuel adjustment file (.adj) --------
    p.fuel_adj_file  = "";  pp.query("fuel_adj_file",  p.fuel_adj_file);
    p.fuel_adj_model = 0;   pp.query("fuel_adj_model", p.fuel_adj_model);
    if (!p.fuel_adj_file.empty())
        Print() << "Fuel adjustment file: " << p.fuel_adj_file << "\n";

    // -------- Time-varying fuel moisture schedule (.fmd) --------
    p.fmd_file        = "";  pp.query("fmd_file",        p.fmd_file);
    p.fmd_start_year  = 0;   pp.query("fmd_start_year",  p.fmd_start_year);
    p.fmd_start_month = 0;   pp.query("fmd_start_month", p.fmd_start_month);
    p.fmd_start_day   = 0;   pp.query("fmd_start_day",   p.fmd_start_day);
    p.fmd_start_hour  = 0;   pp.query("fmd_start_hour",  p.fmd_start_hour);
    p.fmd_fuel_model  = 0;   pp.query("fmd_fuel_model",  p.fmd_fuel_model);
    if (!p.fmd_file.empty())
        Print() << "Fuel moisture schedule file: " << p.fmd_file << "\n";

    // -------- Fire statistics time series --------
    p.fire_stats_file = "fire_stats.csv";
    pp.query("fire_stats_file", p.fire_stats_file);

    // -------- Fire perimeter output --------
    p.write_perimeter_csv     = 1;  pp.query("write_perimeter_csv",     p.write_perimeter_csv);
    p.write_perimeter_geojson = 0;  pp.query("write_perimeter_geojson", p.write_perimeter_geojson);

    // -------- Crown spatial layers from binary LCP --------
    p.use_spatial_crown = 1;  pp.query("use_spatial_crown", p.use_spatial_crown);

    // -------- Diurnal fuel moisture model --------
    p.diurnal_moisture.enable      = 0;       pp.query("diurnal_moisture.enable",      p.diurnal_moisture.enable);
    p.diurnal_moisture.T_min       = 15.0;    pp.query("diurnal_moisture.T_min",       p.diurnal_moisture.T_min);
    p.diurnal_moisture.T_max       = 35.0;    pp.query("diurnal_moisture.T_max",       p.diurnal_moisture.T_max);
    p.diurnal_moisture.RH_min      = 10.0;    pp.query("diurnal_moisture.RH_min",      p.diurnal_moisture.RH_min);
    p.diurnal_moisture.RH_max      = 60.0;    pp.query("diurnal_moisture.RH_max",      p.diurnal_moisture.RH_max);
    p.diurnal_moisture.t_start_s   = 36000.0; pp.query("diurnal_moisture.t_start_s",   p.diurnal_moisture.t_start_s);
    p.diurnal_moisture.t_T_peak_s  = 50400.0; pp.query("diurnal_moisture.t_T_peak_s",  p.diurnal_moisture.t_T_peak_s);

    if (p.diurnal_moisture.enable == 1) {
        if (p.diurnal_moisture.T_max <= p.diurnal_moisture.T_min)
            amrex::Abort("diurnal_moisture.T_max must be > T_min");
        if (p.diurnal_moisture.RH_min < 0.0 || p.diurnal_moisture.RH_max > 100.0 ||
            p.diurnal_moisture.RH_max <= p.diurnal_moisture.RH_min)
            amrex::Abort("diurnal_moisture: RH_min/RH_max must satisfy 0 <= RH_min < RH_max <= 100");
        Print() << "Diurnal moisture model enabled (Nelson 2000 EMC):\n";
        Print() << "  T: [" << p.diurnal_moisture.T_min << ", " << p.diurnal_moisture.T_max << "] °C"
                << "  RH: [" << p.diurnal_moisture.RH_min << ", " << p.diurnal_moisture.RH_max << "] %\n";
        Print() << "  t_start=" << p.diurnal_moisture.t_start_s << " s from midnight"
                << "  t_T_peak=" << p.diurnal_moisture.t_T_peak_s << " s from midnight\n";
        if (!p.fmd_file.empty())
            Print() << "  WARNING: diurnal_moisture is overridden by fmd_file when both are set\n";
    }

    // -------- Fire ecology diagnostics --------
    // Scorch height (Van Wagner 1973), probability of ignition (Anderson 1970),
    // tree mortality (Ryan-Reinhardt 1988), and crown activity classification
    // (Van Wagner 1977 / Scott-Reinhardt 2001).  Always computed; uses crown.*
    // for CBH/CBD/FMC regardless of whether crown.enable is set.
    p.fire_ecology.T_a_C           = 25.0;   pp.query("fire_ecology.T_a_C",           p.fire_ecology.T_a_C);
    p.fire_ecology.solar_heating_F = 25.0;   pp.query("fire_ecology.solar_heating_F", p.fire_ecology.solar_heating_F);
    p.fire_ecology.tree_height     = 10.0;   pp.query("fire_ecology.tree_height",     p.fire_ecology.tree_height);
    if (p.fire_ecology.tree_height <= 0.0)
        amrex::Abort("fire_ecology.tree_height must be > 0 m");

    // -------- Fire emissions (CO2, CO, PM2.5) --------
    // Emission factors from Seiler & Crutzen (1980) / WRF-Fire defaults.
    p.emissions.EF_CO2              = 1.570;  pp.query("emissions.EF_CO2",              p.emissions.EF_CO2);
    p.emissions.EF_CO               = 0.102;  pp.query("emissions.EF_CO",               p.emissions.EF_CO);
    p.emissions.EF_PM25             = 0.0162; pp.query("emissions.EF_PM25",             p.emissions.EF_PM25);
    p.emissions.default_consumed_frac = 0.7;  pp.query("emissions.default_consumed_frac", p.emissions.default_consumed_frac);
    if (p.emissions.EF_CO2  < 0.0) amrex::Abort("emissions.EF_CO2 must be >= 0");
    if (p.emissions.EF_CO   < 0.0) amrex::Abort("emissions.EF_CO must be >= 0");
    if (p.emissions.EF_PM25 < 0.0) amrex::Abort("emissions.EF_PM25 must be >= 0");
    if (p.emissions.default_consumed_frac < 0.0 || p.emissions.default_consumed_frac > 1.0)
        amrex::Abort("emissions.default_consumed_frac must be in [0,1]");

    // -------- FARSITE barrier polygon / firebreak files --------
    // Accept a space-separated list: barrier_files = file1.csv file2.csv
    {
        std::vector<std::string> bfiles;
        pp.queryarr("barrier_files", bfiles);
        p.barrier_files = bfiles;
        if (!p.barrier_files.empty()) {
            Print() << "Barrier polygon files (" << p.barrier_files.size() << "):\n";
            for (const auto& f : p.barrier_files)
                Print() << "  " << f << "\n";
        }
    }

}
