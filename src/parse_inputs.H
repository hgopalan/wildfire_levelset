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

    // Source selection
    std::string source_type;

    // Sphere parameters
    amrex::Real cx, cy, cz, radius;
    // Box Extents 
    amrex::Real xmin, xmax, ymin, ymax, zmin, zmax;

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
    } farsite;
    
    // Skip level set evolution option
    int skip_levelset;                       // 1 to skip level set advection and use initial phi, 0 for normal operation (default: 0)
};

void parse_inputs(InputParameters& inputs);

#endif
