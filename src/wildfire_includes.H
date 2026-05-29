// wildfire_includes.H
// Single convenience header that pulls in every AMReX, standard-library, and
// project header needed by main.cpp.  Include this one file instead of the
// original ~55-line include block.
#ifndef WILDFIRE_INCLUDES_H
#define WILDFIRE_INCLUDES_H

// ---- AMReX headers ----
#include <AMReX.H>
#include <AMReX_Array4.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_BoxArray.H>
#include <AMReX_VisMF.H>

// ---- Standard-library headers ----
#include <algorithm>
#include <cmath>
#include <vector>
#include <string>
#include <memory>
#include <fstream>
#include <sstream>

using namespace amrex;

// ---- Project headers: core numerics and geometry ----
#include "numerical_schemes.H"
#include "initial_conditions.H"
#include "compute_dt.H"
#include "boundary_conditions.H"
#include "reinitialization.H"
#include "advection.H"
#include "velocity_field.H"
#include "parse_inputs.H"

// ---- Project headers: terrain and landscape ----
#include "terrain_slope.H"
#include "landscape_file.H"

// ---- Project headers: fire-spread models ----
#include "farsite_ellipse.H"
#include "compute_rothermel_R.H"
#include "balbi_model.H"
#include "compute_balbi_R.H"
#include "andrews_model.H"
#include "cheney_gould_model.H"
#include "compute_cheney_gould_R.H"
#include "cruz_crown_model.H"
#include "compute_cruz_crown_R.H"
#include "fbp_model.H"
#include "compute_fbp_R.H"
#include "lautenberger_model.H"
#include "compute_lautenberger_R.H"

// ---- Project headers: ROS dispatch and crown override (extracted helpers) ----
#include "compute_ros_dispatch.H"
#include "crown_ros_override.H"

// ---- Project headers: spotting ----
#include "firebrand_spotting.H"
#include "albini_spotting.H"
#include "albini_torching_spotting.H"
#include "ember_cascade_flux.H"
#include "scott_spotting_table.H"

// ---- Project headers: fire diagnostics and ecology ----
#include "compute_fire_behavior.H"
#include "fire_ecology_model.H"
#include "fire_emissions_model.H"
#include "heat_per_unit_area.H"
#include "burnout_time.H"
#include "vorticity_model.H"
#include "compute_reaction_intensity.H"
#include "fire_acceleration.H"

// ---- Project headers: wind models ----
#include "weise_biging_whirl.H"
#include "viegas_model.H"
#include "wind_terrain_models.H"
#include "heat_flux_model.H"
#include "turb_wind.H"

// ---- Project headers: fuel and moisture ----
#include "fuel_adjustment.H"
#include "fuel_moisture_scheduler.H"
#include "fmc_schedule.H"
#include "herb_moisture_schedule.H"
#include "precipitation_moisture.H"
#include "fms_moisture.H"
#include "solar_radiation.H"

// ---- Project headers: propagation methods ----
#include "mtt_propagation.H"

// ---- Project headers: events and scheduling ----
#include "ignition_schedule.H"
#include "polygon_ignition.H"
#include "wind_dir_schedule.H"
#include "barrier_polygons.H"
#include "retardant_drop.H"
#include "satellite_assimilation.H"

// ---- Project headers: weather inputs ----
#include "wtr_weather.H"
#include "multi_wtr_weather.H"
#include "conditional_weather.H"

// ---- Project headers: field setup ----
#include "multifab_setup.H"

// ---- Project headers: extracted initialisation helpers ----
#include "init_phi_source.H"
#include "terrain_setup.H"
#include "fuel_table_setup.H"
#include "schedule_setup.H"

// ---- Project headers: output ----
#include "write_xy_data.H"
#include "plot_results.H"
#include "fire_report.H"
#include "smoke_plume_rise.H"
#include "sim_datetime.H"

// ---- Project headers: new fire behavior diagnostics ----
#include "byram_convective_number.H"
#include "flame_tilt_angle.H"
#include "packing_ratio.H"
#include "flame_front_depth.H"
#include "mcarthur_ffdi.H"

// ---- Project headers: new moisture and fuel models ----
#include "canadian_fwi.H"
#include "two_fuel_blending.H"
#include "duff_moisture_smoldering.H"
#include "moisture_of_extinction.H"

#endif // WILDFIRE_INCLUDES_H
