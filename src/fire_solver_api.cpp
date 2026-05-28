// ============================================================================
// fire_solver_api.cpp
// Implementation of Python API for wildfire_levelset solver
// ============================================================================

#include "fire_solver_api.H"
#include "init_phi_source.H"
#include "initial_conditions.H"
#include "terrain_setup.H"
#include "fuel_table_setup.H"
#include "velocity_field.H"
#include "compute_ros_dispatch.H"
#include "compute_dt.H"
#include "compute_fire_behavior.H"
#include "advection.H"
#include "plot_results.H"
#include "boundary_conditions.H"

using namespace amrex;

// Global singleton instance
std::unique_ptr<FireSolverState> g_fire_solver_state = nullptr;

// ============================================================================
// Helper function to ensure AMReX is initialized
// ============================================================================
static bool ensure_amrex_initialized() {
    static bool amrex_init_done = false;
    if (!amrex_init_done) {
        if (!amrex::Initialized()) {
            int argc = 0;
            char** argv = nullptr;
            amrex::Initialize(argc, argv, false);
        }
        amrex_init_done = true;
    }
    return true;
}

// ============================================================================
// fire_solver_initialize
// ============================================================================
bool fire_solver_initialize(const std::string& inputs_file) {
    ensure_amrex_initialized();
    
    // Clean up any existing state
    if (g_fire_solver_state && g_fire_solver_state->initialized) {
        fire_solver_finalize();
    }
    
    // Create new state
    g_fire_solver_state = std::make_unique<FireSolverState>();
    FireSolverState& state = *g_fire_solver_state;
    
    try {
        // Parse inputs from file
        state.inputs = std::make_unique<InputParameters>();
        
        // Read the inputs file
        ParmParse::Initialize(0, nullptr, inputs_file.c_str());
        parse_inputs(*state.inputs);
        
        // Setup geometry
        IntVect dom_lo(0, 0);
        IntVect dom_hi(state.inputs->n_cell_x-1, state.inputs->n_cell_y-1);
        Box domain(dom_lo, dom_hi);
        
        RealBox rb({state.inputs->plo_x, state.inputs->plo_y},
                   {state.inputs->phi_x, state.inputs->phi_y});
        
        Array<int,AMREX_SPACEDIM> is_periodic{0, 0};
        state.geom = std::make_unique<Geometry>(domain, &rb, 0, is_periodic.data());
        
        // Setup grids
        state.ba = std::make_unique<BoxArray>(domain);
        state.ba->maxSize(state.inputs->max_grid);
        state.dm = std::make_unique<DistributionMapping>(*state.ba);
        
        // Allocate all MultiFabs
        state.fields = std::make_unique<WildfireFields>(*state.ba, *state.dm, *state.inputs);
        
        // Initialize phi (fire front level set)
        state.use_levelset = (state.inputs->propagation_method == "levelset");
        state.use_mtt = (state.inputs->propagation_method == "mtt");
        state.use_farsite = (state.inputs->propagation_method == "farsite");
        
        bool use_indicator = state.use_farsite;
        init_phi_source(state.fields->phi, state.fields->arrival_time_mf, *state.geom,
                       *state.inputs, use_indicator, state.step, state.time);
        
        // Initialize velocity field
        if (!state.inputs->velocity_file.empty()) {
            if (state.inputs->use_time_dependent_wind == 1) {
                amrex::Print() << "Using time-dependent wind fields with spacing = " 
                              << state.inputs->wind_time_spacing << " seconds\n";
                update_time_dependent_velocity(state.fields->vel, *state.geom,
                                              state.inputs->velocity_file, 0.0,
                                              state.inputs->wind_time_spacing,
                                              state.wind_x_data1, state.wind_y_data1,
                                              state.wind_u_data1, state.wind_v_data1,
                                              state.wind_x_data2, state.wind_y_data2,
                                              state.wind_u_data2, state.wind_v_data2,
                                              state.current_wind_field_index,
                                              state.next_wind_field_index);
            } else {
                init_velocity_from_file(state.fields->vel, *state.geom, state.inputs->velocity_file);
            }
        } else {
            init_velocity_constant(state.fields->vel, *state.geom,
                                 state.inputs->ux, state.inputs->uy, state.inputs->uz);
        }
        
        // Setup turbulent wind if needed
        state.turb_wind_active = (state.inputs->turb_wind.model != "none");
        if (state.turb_wind_active) {
            state.vel_base = std::make_unique<MultiFab>(*state.ba, *state.dm, 3, 1);
            MultiFab::Copy(*state.vel_base, state.fields->vel, 0, 0, 3, 1);
            init_turb_wind_state(state.turb_state, state.inputs->turb_wind);
            
            if (state.inputs->turb_wind.model == "ou_process" &&
                state.inputs->turb_wind.L_c > Real(0.0)) {
                state.ou_state_mf = std::make_unique<MultiFab>(*state.ba, *state.dm, 2, 0);
                state.ou_state_mf->setVal(Real(0.0));
            }
        }
        
        // Load 3D wind for spotting if specified
        if (state.inputs->albini_spotting.enable == 1 &&
            state.inputs->albini_spotting.use_3d_wind == 1 &&
            !state.inputs->albini_spotting.plt_wind_file.empty()) {
            if (!read_plt_wind_file(state.inputs->albini_spotting.plt_wind_file,
                                   state.albini_plt_wind)) {
                amrex::Print() << "Warning: Failed to read 3D wind plt file for Albini spotting\n";
            }
        }
        
        if (state.inputs->ember_cascade.enable == 1 &&
            (state.inputs->ember_cascade.use_3d_wind == 1 ||
             state.inputs->ember_cascade.require_3d_wind == 1) &&
            !state.inputs->ember_cascade.plt_wind_file.empty()) {
            if (!read_plt_wind_file(state.inputs->ember_cascade.plt_wind_file,
                                   state.ember_cascade_plt_wind)) {
                amrex::Print() << "Warning: Failed to read 3D wind plt file for ember cascade\n";
            }
        }
        
        // Initialize terrain if specified
        if (!state.inputs->rothermel.landscape_file.empty()) {
            auto landscape_result = load_landscape_file(
                state.inputs->rothermel.landscape_file,
                state.inputs->rothermel.landscape_fuel_type,
                *state.geom);
            
            state.has_spatial_crown = landscape_result.has_crown;
            state.has_spatial_moisture = landscape_result.has_moisture;
            
            // Copy landscape data to MultiFabs
            if (landscape_result.has_elevation) {
                copy_to_multifab(state.fields->elevation_mf, landscape_result.elevation, *state.geom);
                copy_to_multifab(state.fields->slope_mf, landscape_result.slope, *state.geom);
                copy_to_multifab(state.fields->aspect_mf, landscape_result.aspect, *state.geom);
            }
            
            if (landscape_result.has_fuel) {
                copy_to_multifab(state.fields->fuel_model_mf, landscape_result.fuel_model, *state.geom);
            }
            
            if (state.has_spatial_crown) {
                copy_to_multifab(state.fields->cbh_mf, landscape_result.cbh, *state.geom);
                copy_to_multifab(state.fields->cbd_mf, landscape_result.cbd, *state.geom);
                copy_to_multifab(state.fields->cc_mf, landscape_result.cc, *state.geom);
                copy_to_multifab(state.fields->canopy_height_mf, landscape_result.canopy_height, *state.geom);
            }
            
            if (state.has_spatial_moisture) {
                for (int i = 0; i < 6; ++i) {
                    copy_to_multifab_comp(state.fields->spatial_moisture_mf, i,
                                         landscape_result.moisture[i], *state.geom);
                }
            }
        }
        
        // Setup terrain slopes if needed
        if (!state.inputs->rothermel.terrain_file.empty() ||
            !state.inputs->rothermel.landscape_file.empty()) {
            state.fields->terrain_slopes = std::make_unique<MultiFab>(*state.ba, *state.dm, 2, 1);
            setup_terrain_slopes(*state.fields->terrain_slopes, state.fields->slope_mf,
                               state.fields->aspect_mf, *state.geom, *state.inputs);
        }
        
        // Build fuel lookup tables
        if (!state.inputs->rothermel.landscape_file.empty()) {
            state.h_fuel_table = build_fuel_rothermel_table(
                *state.inputs,
                state.inputs->rothermel.landscape_fuel_type);
            state.fuel_table_size = static_cast<int>(state.h_fuel_table.size());
            
            if (state.fuel_table_size > 0) {
                state.d_fuel_table.resize(state.fuel_table_size);
                Gpu::copy(Gpu::hostToDevice,
                         state.h_fuel_table.begin(),
                         state.h_fuel_table.end(),
                         state.d_fuel_table.begin());
                state.d_fuel_table_ptr = state.d_fuel_table.data();
            }
        }
        
        // Compute initial ROS and dt
        if (state.use_levelset || state.use_farsite) {
            compute_rothermel_R(state.fields->R_mf, state.fields->vel, *state.geom,
                              *state.inputs,
                              state.fields->terrain_slopes.get(),
                              !state.inputs->rothermel.landscape_file.empty() ?
                                  &state.fields->fuel_model_mf : nullptr,
                              state.d_fuel_table_ptr, state.fuel_table_size,
                              state.has_spatial_crown ? &state.fields->cc_mf : nullptr,
                              state.has_spatial_crown ? &state.fields->canopy_height_mf : nullptr,
                              0, 30.0, 0.1);
            
            state.dt = compute_dt(state.fields->R_mf, *state.geom, state.inputs->cfl);
            
            compute_fire_behavior(state.fields->fireline_intensity_mf,
                                state.fields->flame_length_mf,
                                state.fields->R_mf,
                                *state.inputs);
        }
        
        state.initialized = true;
        amrex::Print() << "Fire solver initialized successfully\n";
        amrex::Print() << "  Grid: " << state.inputs->n_cell_x << " × "
                      << state.inputs->n_cell_y << "\n";
        amrex::Print() << "  Domain: [" << state.inputs->plo_x << ", "
                      << state.inputs->phi_x << "] × ["
                      << state.inputs->plo_y << ", " << state.inputs->phi_y << "]\n";
        amrex::Print() << "  Initial dt: " << state.dt << " s\n";
        
        return true;
        
    } catch (const std::exception& e) {
        amrex::Print() << "Error initializing fire solver: " << e.what() << "\n";
        g_fire_solver_state.reset();
        return false;
    }
}

// ============================================================================
// fire_solver_advance
// ============================================================================
Real fire_solver_advance() {
    if (!g_fire_solver_state || !g_fire_solver_state->initialized) {
        amrex::Print() << "Error: Fire solver not initialized\n";
        return -1.0;
    }
    
    FireSolverState& state = *g_fire_solver_state;
    state.step++;
    
    try {
        // Fill boundary conditions
        fill_boundary_extrap(state.fields->phi, *state.geom);
        
        const Real dt_step = state.dt;
        
        amrex::Print() << "Step " << state.step << ": Time = " << state.time
                      << ", dt = " << dt_step << "\n";
        
        // Update time-dependent wind if needed
        if (!state.inputs->velocity_file.empty() &&
            state.inputs->use_time_dependent_wind == 1) {
            MultiFab& wind_target = state.turb_wind_active ? *state.vel_base : state.fields->vel;
            update_time_dependent_velocity(wind_target, *state.geom,
                                          state.inputs->velocity_file,
                                          state.time, state.inputs->wind_time_spacing,
                                          state.wind_x_data1, state.wind_y_data1,
                                          state.wind_u_data1, state.wind_v_data1,
                                          state.wind_x_data2, state.wind_y_data2,
                                          state.wind_u_data2, state.wind_v_data2,
                                          state.current_wind_field_index,
                                          state.next_wind_field_index);
        }
        
        // Apply turbulent wind if active
        if (state.turb_wind_active) {
            apply_turb_wind(state.fields->vel, *state.vel_base,
                          state.ou_state_mf.get(), state.turb_state,
                          dt_step, state.inputs->turb_wind, *state.geom);
        }
        
        // Compute ROS for this timestep
        compute_rothermel_R(state.fields->R_mf, state.fields->vel, *state.geom,
                          *state.inputs,
                          state.fields->terrain_slopes.get(),
                          !state.inputs->rothermel.landscape_file.empty() ?
                              &state.fields->fuel_model_mf : nullptr,
                          state.d_fuel_table_ptr, state.fuel_table_size,
                          state.has_spatial_crown ? &state.fields->cc_mf : nullptr,
                          state.has_spatial_crown ? &state.fields->canopy_height_mf : nullptr,
                          0, 30.0, 0.1);
        
        // Advect phi using level set method
        if (state.use_levelset) {
            advect_phi_rk2(state.fields->phi, state.fields->R_mf,
                          *state.geom, dt_step);
        } else if (state.use_farsite) {
            // FARSITE elliptical spread would go here
            // For now, use simple level set as fallback
            advect_phi_rk2(state.fields->phi, state.fields->R_mf,
                          *state.geom, dt_step);
        }
        
        // Update fire behavior fields
        compute_fire_behavior(state.fields->fireline_intensity_mf,
                            state.fields->flame_length_mf,
                            state.fields->R_mf,
                            *state.inputs);
        
        // Update arrival time for newly burned cells
        for (MFIter mfi(state.fields->arrival_time_mf); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            auto const& phi_arr = state.fields->phi.const_array(mfi);
            auto arrival_arr = state.fields->arrival_time_mf.array(mfi);
            const Real current_time = state.time + dt_step;
            
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                // If cell just burned (phi crossed zero), record arrival time
                if (phi_arr(i,j,k) <= Real(0.0) && arrival_arr(i,j,k) < Real(0.0)) {
                    arrival_arr(i,j,k) = current_time;
                }
            });
        }
        
        // Compute next timestep
        state.dt = compute_dt(state.fields->R_mf, *state.geom, state.inputs->cfl);
        state.time += dt_step;
        
        return dt_step;
        
    } catch (const std::exception& e) {
        amrex::Print() << "Error advancing fire solver: " << e.what() << "\n";
        return -1.0;
    }
}

// ============================================================================
// fire_solver_get_time
// ============================================================================
void fire_solver_get_time(Real& time, int& step) {
    if (g_fire_solver_state && g_fire_solver_state->initialized) {
        time = g_fire_solver_state->time;
        step = g_fire_solver_state->step;
    } else {
        time = 0.0;
        step = 0;
    }
}

// ============================================================================
// fire_solver_get_geometry
// ============================================================================
void fire_solver_get_geometry(
    int& nx, int& ny,
    double& xmin, double& xmax,
    double& ymin, double& ymax,
    double& dx, double& dy) {
    
    if (!g_fire_solver_state || !g_fire_solver_state->initialized) {
        nx = ny = 0;
        xmin = xmax = ymin = ymax = dx = dy = 0.0;
        return;
    }
    
    FireSolverState& state = *g_fire_solver_state;
    nx = state.inputs->n_cell_x;
    ny = state.inputs->n_cell_y;
    xmin = static_cast<double>(state.inputs->plo_x);
    xmax = static_cast<double>(state.inputs->phi_x);
    ymin = static_cast<double>(state.inputs->plo_y);
    ymax = static_cast<double>(state.inputs->phi_y);
    dx = static_cast<double>(state.geom->CellSize(0));
    dy = static_cast<double>(state.geom->CellSize(1));
}

// ============================================================================
// fire_solver_update_wind
// ============================================================================
bool fire_solver_update_wind(int nx, int ny,
                             const double* u_data,
                             const double* v_data) {
    if (!g_fire_solver_state || !g_fire_solver_state->initialized) {
        amrex::Print() << "Error: Fire solver not initialized\n";
        return false;
    }
    
    if (nx != g_fire_solver_state->inputs->n_cell_x ||
        ny != g_fire_solver_state->inputs->n_cell_y) {
        amrex::Print() << "Error: Wind grid size mismatch\n";
        return false;
    }
    
    FireSolverState& state = *g_fire_solver_state;
    MultiFab& vel_target = state.turb_wind_active ? *state.vel_base : state.fields->vel;
    
    // Copy wind data into velocity MultiFab
    for (MFIter mfi(vel_target); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto vel_arr = vel_target.array(mfi);
        
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const int idx = j * nx + i;
            vel_arr(i, j, k, 0) = static_cast<Real>(u_data[idx]);
            vel_arr(i, j, k, 1) = static_cast<Real>(v_data[idx]);
            vel_arr(i, j, k, 2) = Real(0.0);  // w-component = 0 for 2D
        });
    }
    
    return true;
}

// ============================================================================
// fire_solver_update_wind_3d
// ============================================================================
bool fire_solver_update_wind_3d(
    int nx, int ny, int nz,
    double xmin, double xmax,
    double ymin, double ymax,
    double zmin, double zmax,
    const double* u_data,
    const double* v_data,
    const double* w_data) {
    
    if (!g_fire_solver_state || !g_fire_solver_state->initialized) {
        amrex::Print() << "Error: Fire solver not initialized\n";
        return false;
    }
    
    // Update 3D wind for spotting calculations
    FireSolverState& state = *g_fire_solver_state;
    
    // Load into albini_plt_wind structure
    if (state.inputs->albini_spotting.enable == 1 &&
        state.inputs->albini_spotting.use_3d_wind == 1) {
        bool success = load_plt_wind_from_arrays(
            nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax,
            u_data, v_data, w_data,
            state.albini_plt_wind);
        if (!success) {
            amrex::Print() << "Warning: Failed to load 3D wind for Albini spotting\n";
        }
    }
    
    // Also load into ember_cascade_plt_wind if needed
    if (state.inputs->ember_cascade.enable == 1 &&
        state.inputs->ember_cascade.use_3d_wind == 1) {
        bool success = load_plt_wind_from_arrays(
            nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax,
            u_data, v_data, w_data,
            state.ember_cascade_plt_wind);
        if (!success) {
            amrex::Print() << "Warning: Failed to load 3D wind for ember cascade\n";
        }
    }
    
    // Extract 2D column-averaged wind and update velocity field
    if (state.albini_plt_wind.valid && state.albini_plt_wind.nx_2d == nx &&
        state.albini_plt_wind.ny_2d == ny) {
        return fire_solver_update_wind(nx, ny,
                                       state.albini_plt_wind.u2d.data(),
                                       state.albini_plt_wind.v2d.data());
    }
    
    return true;
}

// ============================================================================
// Helper function to extract MultiFab as flattened vector
// ============================================================================
static std::vector<double> extract_multifab(const MultiFab& mf, int comp = 0) {
    const int nx = mf.nGrow() > 0 ? mf.nGrowVect()[0] : 0;
    const int ny = mf.nGrow() > 0 ? mf.nGrowVect()[1] : 0;
    
    std::vector<double> result;
    
    // Get total size
    long total_cells = 0;
    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
        total_cells += mfi.validbox().numPts();
    }
    result.reserve(total_cells);
    
    // Extract data in Fortran order (column-major)
    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const& arr = mf.const_array(mfi);
        
        for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
            for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
                result.push_back(static_cast<double>(arr(i, j, 0, comp)));
            }
        }
    }
    
    return result;
}

// ============================================================================
// fire_solver_get_phi
// ============================================================================
std::vector<double> fire_solver_get_phi() {
    if (!g_fire_solver_state || !g_fire_solver_state->initialized) {
        return std::vector<double>();
    }
    return extract_multifab(g_fire_solver_state->fields->phi);
}

// ============================================================================
// fire_solver_get_ros
// ============================================================================
std::vector<double> fire_solver_get_ros() {
    if (!g_fire_solver_state || !g_fire_solver_state->initialized) {
        return std::vector<double>();
    }
    return extract_multifab(g_fire_solver_state->fields->R_mf);
}

// ============================================================================
// fire_solver_get_intensity
// ============================================================================
std::vector<double> fire_solver_get_intensity() {
    if (!g_fire_solver_state || !g_fire_solver_state->initialized) {
        return std::vector<double>();
    }
    return extract_multifab(g_fire_solver_state->fields->fireline_intensity_mf);
}

// ============================================================================
// fire_solver_get_flame_length
// ============================================================================
std::vector<double> fire_solver_get_flame_length() {
    if (!g_fire_solver_state || !g_fire_solver_state->initialized) {
        return std::vector<double>();
    }
    return extract_multifab(g_fire_solver_state->fields->flame_length_mf);
}

// ============================================================================
// fire_solver_get_wind
// ============================================================================
void fire_solver_get_wind(std::vector<double>& u_data, std::vector<double>& v_data) {
    if (!g_fire_solver_state || !g_fire_solver_state->initialized) {
        u_data.clear();
        v_data.clear();
        return;
    }
    u_data = extract_multifab(g_fire_solver_state->fields->vel, 0);
    v_data = extract_multifab(g_fire_solver_state->fields->vel, 1);
}

// ============================================================================
// fire_solver_get_arrival_time
// ============================================================================
std::vector<double> fire_solver_get_arrival_time() {
    if (!g_fire_solver_state || !g_fire_solver_state->initialized) {
        return std::vector<double>();
    }
    return extract_multifab(g_fire_solver_state->fields->arrival_time_mf);
}

// ============================================================================
// fire_solver_write_plotfile
// ============================================================================
bool fire_solver_write_plotfile(const std::string& plotfile_name) {
    if (!g_fire_solver_state || !g_fire_solver_state->initialized) {
        amrex::Print() << "Error: Fire solver not initialized\n";
        return false;
    }
    
    try {
        FireSolverState& state = *g_fire_solver_state;
        write_wildfire_plotfile(*state.fields, *state.ba, *state.dm, *state.geom,
                               *state.inputs, state.step, state.time, state.ftd,
                               false, false, plotfile_name);
        return true;
    } catch (const std::exception& e) {
        amrex::Print() << "Error writing plotfile: " << e.what() << "\n";
        return false;
    }
}

// ============================================================================
// fire_solver_finalize
// ============================================================================
void fire_solver_finalize() {
    if (g_fire_solver_state) {
        amrex::Print() << "Finalizing fire solver\n";
        g_fire_solver_state.reset();
    }
}

// ============================================================================
// fire_solver_is_initialized
// ============================================================================
bool fire_solver_is_initialized() {
    return (g_fire_solver_state && g_fire_solver_state->initialized);
}
