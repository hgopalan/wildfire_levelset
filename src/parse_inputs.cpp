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
    p.cfl          = 0.5;        pp.query("cfl",          p.cfl);
    p.plot_int     = 50;         pp.query("plot_int",     p.plot_int);

    // ---------------- Inputs: velocity ---------------------
    p.ux = 0.25;                 pp.query("u_x", p.ux);
    p.uy = 0.0;                  pp.query("u_y", p.uy);
    p.uz = 0.0;                  pp.query("u_z", p.uz);
    p.velocity_file = "";        pp.query("velocity_file", p.velocity_file);

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
    p.rothermel.wind_conv = 196.85;
    p.rothermel.ros_conv  = 0.00508;
    
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
    pp.query("rothermel.wind_conv", p.rothermel.wind_conv);
    pp.query("rothermel.ros_conv",  p.rothermel.ros_conv);

    // -------- FARSITE ellipse model parameters (Richards 1990) --------
    p.farsite.enable = 1;                        pp.query("farsite.enable", p.farsite.enable);
    p.farsite.use_anderson_LW = 0;               pp.query("farsite.use_anderson_LW", p.farsite.use_anderson_LW);
    p.farsite.length_to_width_ratio = 3.0;       pp.query("farsite.length_to_width_ratio", p.farsite.length_to_width_ratio);
    p.farsite.phi_threshold = 0.1;               pp.query("farsite.phi_threshold", p.farsite.phi_threshold);
    p.farsite.coeff_a = 1.0;                     pp.query("farsite.coeff_a", p.farsite.coeff_a);
    p.farsite.coeff_b = 0.5;                     pp.query("farsite.coeff_b", p.farsite.coeff_b);
    p.farsite.coeff_c = 0.2;                     pp.query("farsite.coeff_c", p.farsite.coeff_c);
    
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

    // -------- Skip level set option --------
    p.skip_levelset = 0;                         pp.query("skip_levelset", p.skip_levelset);
}
