#ifndef PARSE_INPUTS_H
#define PARSE_INPUTS_H

#include <AMReX_ParmParse.H>
#include <string>

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

    // Source selection
    std::string source_type;

    // Sphere parameters
    amrex::Real cx, cy, cz, radius;
    // Box Size 
    amrex::Real bx, by, bz;

    // Dynamic AMR-style refinement for negative phi
    int amr_enable_negative_phi_refine;
    int amr_regrid_int;
    int amr_refine_ratio;
    int amr_max_refinements;
    amrex::Real amr_tag_phi_threshold;
};

void parse_inputs(InputParameters& inputs);

#endif
