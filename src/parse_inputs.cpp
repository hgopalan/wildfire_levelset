// parse_inputs.cpp
#include "parse_inputs.H"
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

    p.bx = 32;                 pp.query("box_size_x", p.bx);
    p.by = 32;                 pp.query("box_size_y", p.by);
    p.bz = 32;                 pp.query("box_size_z", p.bz);

    // ---------------- Inputs: dynamic AMR for negative phi -----
    p.amr_enable_negative_phi_refine = 1;
    pp.query("amr_enable_negative_phi_refine", p.amr_enable_negative_phi_refine);
    p.amr_regrid_int = 10;
    pp.query("amr_regrid_int", p.amr_regrid_int);
    p.amr_refine_ratio = 2;
    pp.query("amr_refine_ratio", p.amr_refine_ratio);
    p.amr_max_refinements = 1;
    pp.query("amr_max_refinements", p.amr_max_refinements);
    p.amr_tag_phi_threshold = 0.0;
    pp.query("amr_tag_phi_threshold", p.amr_tag_phi_threshold);

    // -------- Rothermel fire spread model (Southern California chaparral defaults) --------
    p.rothermel.w0        = 0.230;    pp.query("rothermel.w0",        p.rothermel.w0);
    p.rothermel.sigma     = 1739.0;   pp.query("rothermel.sigma",     p.rothermel.sigma);
    p.rothermel.delta     = 2.0;      pp.query("rothermel.delta",     p.rothermel.delta);
    p.rothermel.M_f       = 0.08;     pp.query("rothermel.M_f",       p.rothermel.M_f);
    p.rothermel.M_x       = 0.30;     pp.query("rothermel.M_x",       p.rothermel.M_x);
    p.rothermel.h_heat    = 8000.0;   pp.query("rothermel.h_heat",    p.rothermel.h_heat);
    p.rothermel.S_T       = 0.0555;   pp.query("rothermel.S_T",       p.rothermel.S_T);
    p.rothermel.S_e       = 0.010;    pp.query("rothermel.S_e",       p.rothermel.S_e);
    p.rothermel.rho_p     = 32.0;     pp.query("rothermel.rho_p",     p.rothermel.rho_p);
    p.rothermel.slope_x   = 0.0;      pp.query("rothermel.slope_x",   p.rothermel.slope_x);
    p.rothermel.slope_y   = 0.0;      pp.query("rothermel.slope_y",   p.rothermel.slope_y);
    p.rothermel.wind_conv = 196.85;   pp.query("rothermel.wind_conv", p.rothermel.wind_conv);
    p.rothermel.ros_conv  = 0.00508;  pp.query("rothermel.ros_conv",  p.rothermel.ros_conv);

    // -------- FARSITE ellipse model parameters --------
    p.farsite.enable = 1;                        pp.query("farsite.enable", p.farsite.enable);
    p.farsite.length_to_width_ratio = 3.0;       pp.query("farsite.length_to_width_ratio", p.farsite.length_to_width_ratio);
    p.farsite.phi_threshold = 0.1;               pp.query("farsite.phi_threshold", p.farsite.phi_threshold);
}
