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
    int reinit_iters;
    amrex::Real reinit_dtau;
    int nsteps;
    amrex::Real cfl;
    int plot_int;

    // Velocity
    amrex::Real ux, uy, uz;

    // Source selection
    std::string source_type;

    // Sphere parameters
    amrex::Real cx, cy, cz, radius;
};

void parse_inputs(InputParameters& inputs);

#endif
