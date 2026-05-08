// init_phi_source.H
// Dispatches phi initialisation to the appropriate geometry primitive or file
// source based on inputs.source_type / inputs.restart_chkfile /
// inputs.fire_points_file, then marks cells that are already burned in the
// initial state (arrival_time = 0).
//
// Extracted from main.cpp to keep the main function concise.
#ifndef INIT_PHI_SOURCE_H
#define INIT_PHI_SOURCE_H

// All required headers provided by wildfire_includes.H (the umbrella header).

using namespace amrex;

// Initialise phi from the source specified in inputs and fill ghost cells.
// restart_step and time are set when restarting from a checkpoint; otherwise
// they keep their default values of 0 / 0.0.
inline void init_phi_source(
    MultiFab& phi,
    MultiFab& arrival_time_mf,
    const Geometry& geom,
    const InputParameters& inputs,
    bool use_indicator,
    int& restart_step,
    Real& time)
{
    if (!inputs.restart_chkfile.empty()) {
        read_checkpoint(inputs.restart_chkfile, phi, restart_step, time);
        fill_boundary_extrap(phi, geom);
    } else if (!inputs.fire_points_file.empty()) {
        init_phi_from_fire_points_csv(phi, geom, inputs.fire_points_file,
                                      inputs.fire_gaussian_sigma);
        fill_boundary_extrap(phi, geom);
    } else if (inputs.source_type == "sphere") {
        if (use_indicator)
            init_phi_sphere_indicator(phi, geom,
                inputs.cx, inputs.cy, inputs.cz, inputs.radius);
        else
            init_phi_sphere(phi, geom,
                inputs.cx, inputs.cy, inputs.cz, inputs.radius);
        fill_boundary_extrap(phi, geom);
    } else if (inputs.source_type == "box") {
        if (use_indicator)
            init_phi_box_indicator(phi, geom,
                inputs.xmin, inputs.ymin, inputs.zmin,
                inputs.xmax, inputs.ymax, inputs.zmax);
        else
            init_phi_box(phi, geom,
                inputs.xmin, inputs.ymin, inputs.zmin,
                inputs.xmax, inputs.ymax, inputs.zmax);
        fill_boundary_extrap(phi, geom);
    } else if (inputs.source_type == "ellipse") {
        if (use_indicator)
            init_phi_ellipse_indicator(phi, geom,
                inputs.ellipse_center_x, inputs.ellipse_center_y, inputs.ellipse_center_z,
                inputs.ellipse_radius_x, inputs.ellipse_radius_y, inputs.ellipse_radius_z);
        else
            init_phi_ellipse(phi, geom,
                inputs.ellipse_center_x, inputs.ellipse_center_y, inputs.ellipse_center_z,
                inputs.ellipse_radius_x, inputs.ellipse_radius_y, inputs.ellipse_radius_z);
        fill_boundary_extrap(phi, geom);
    } else if (inputs.source_type == "eb") {
        if (use_indicator)
            init_phi_from_eb_implicit_indicator(phi, geom, inputs.eb_type,
                inputs.eb_param1, inputs.eb_param2, inputs.eb_param3,
                inputs.eb_param4, inputs.eb_param5, inputs.eb_param6);
        else
            init_phi_from_eb_implicit(phi, geom, inputs.eb_type,
                inputs.eb_param1, inputs.eb_param2, inputs.eb_param3,
                inputs.eb_param4, inputs.eb_param5, inputs.eb_param6);
        fill_boundary_extrap(phi, geom);
    } else if (inputs.source_type == "polygon") {
        std::vector<Real> poly_xs, poly_ys;
        read_polygon_vertices(inputs.fire_polygon_file, poly_xs, poly_ys);
        init_phi_from_polygon(phi, geom, poly_xs, poly_ys, inputs.fire_polygon_z_level);
        fill_boundary_extrap(phi, geom);
    } else if (inputs.source_type == "polyline") {
        std::vector<Real> poly_xs, poly_ys;
        read_polygon_vertices(inputs.fire_polygon_file, poly_xs, poly_ys);
        init_phi_from_polyline(phi, geom, poly_xs, poly_ys,
                               inputs.polyline_width, inputs.fire_polygon_z_level);
        fill_boundary_extrap(phi, geom);
    } else {
        amrex::Abort("Invalid source_type: " + inputs.source_type);
    }

    // Mark cells already burned in the initial state (arrival_time = 0.0)
    const Real t0 = Real(0.0);
    for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const p  = phi.const_array(mfi);
        auto       at = arrival_time_mf.array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            if (p(i, j, k) < Real(0.0) && at(i, j, k) < Real(0.0))
                at(i, j, k) = t0;
        });
    }
}

#endif // INIT_PHI_SOURCE_H
