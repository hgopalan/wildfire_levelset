// compute_ros_dispatch.H
// Dispatches to the appropriate fire-spread ROS model based on
// inputs.fire_spread_model.  Called once before the time loop (for the
// initial CFL dt) and once per time step.
//
// Extracted from main.cpp to avoid repeating the identical if-else chain
// three times.
#ifndef COMPUTE_ROS_DISPATCH_H
#define COMPUTE_ROS_DISPATCH_H

#include "parse_inputs.H"
#include "compute_rothermel_R.H"
#include "compute_balbi_R.H"
#include "compute_cheney_gould_R.H"
#include "compute_cruz_crown_R.H"
#include "compute_fbp_R.H"
#include "compute_lautenberger_R.H"

using namespace amrex;

// Dispatch to the appropriate fire-spread ROS kernel.
//
// Parameters:
//   R_mf                 – output: ROS field [m/s]
//   vel_for_model        – effective wind field (ambient or terrain-modified)
//   geom                 – AMReX Geometry
//   inputs               – InputParameters (provides model name and all sub-params)
//   terrain_slopes       – terrain slope MultiFab (nullptr when no terrain file)
//   fuel_model_mf        – per-cell fuel model MultiFab (nullptr when no landscape)
//   d_balbi_table_ptr    – device Balbi lookup table (nullptr when not used)
//   balbi_table_size     – number of entries in the Balbi table (0 if not used)
//   d_fuel_table_ptr     – device Rothermel lookup table (nullptr when not used)
//   fuel_table_size      – number of entries in the Rothermel table (0 if not used)
//   has_spatial_crown    – whether per-cell crown data (cc_mf, canopy_height_mf) exist
//   cc_mf                – canopy cover MultiFab (used when has_spatial_crown)
//   canopy_height_mf     – canopy height MultiFab (used when has_spatial_crown)
//   heat_flux_active     – whether heat-flux wind correction is enabled
//   heat_flux_mf         – heat flux MultiFab (nullptr when heat_flux_active is false)
inline void dispatch_compute_ros(
    MultiFab& R_mf,
    const MultiFab& vel_for_model,
    const Geometry& geom,
    const InputParameters& inputs,
    const MultiFab* terrain_slopes,
    const MultiFab* fuel_model_mf,
    const BalbiComputed* d_balbi_table_ptr,
    int balbi_table_size,
    const RothermelComputed* d_fuel_table_ptr,
    int fuel_table_size,
    bool has_spatial_crown,
    const MultiFab* cc_mf,
    const MultiFab* canopy_height_mf,
    bool heat_flux_active,
    const MultiFab* heat_flux_mf,
    int cellsize_enable = 0,
    amrex::Real cellsize_dx_ref = 30.0,
    amrex::Real cellsize_correction_exp = 0.1)
{
    if (inputs.fire_spread_model == "balbi") {
        compute_balbi_R(R_mf, vel_for_model, geom, inputs.rothermel, inputs.balbi,
                        terrain_slopes,
                        fuel_model_mf,
                        d_balbi_table_ptr, balbi_table_size,
                        heat_flux_active ? heat_flux_mf : nullptr,
                        heat_flux_active ? &inputs.heat_flux : nullptr);
    } else if (inputs.fire_spread_model == "cheney_gould") {
        compute_cheney_gould_R(R_mf, vel_for_model, inputs.cheney_gould);
    } else if (inputs.fire_spread_model == "cruz_crown") {
        compute_cruz_crown_R(R_mf, vel_for_model, inputs.cruz_crown);
    } else if (inputs.fire_spread_model == "fbp_o1a" ||
               inputs.fire_spread_model == "fbp_o1b" ||
               inputs.fire_spread_model == "fbp_s1"  ||
               inputs.fire_spread_model == "fbp_s2"  ||
               inputs.fire_spread_model == "fbp_s3") {
        compute_fbp_R(R_mf, vel_for_model, inputs.fbp);
    } else if (inputs.fire_spread_model == "lautenberger") {
        compute_lautenberger_R(R_mf, vel_for_model, inputs.rothermel,
                               inputs.lautenberger, terrain_slopes);
    } else {
        // Default: Rothermel (1972)
        compute_rothermel_R(R_mf, vel_for_model, geom, inputs.rothermel,
                            terrain_slopes,
                            fuel_model_mf,
                            d_fuel_table_ptr, fuel_table_size,
                            has_spatial_crown ? cc_mf : nullptr,
                            has_spatial_crown ? canopy_height_mf : nullptr,
                            cellsize_enable, cellsize_dx_ref, cellsize_correction_exp);
    }

}

#endif // COMPUTE_ROS_DISPATCH_H
