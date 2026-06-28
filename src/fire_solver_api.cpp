// ============================================================================
// fire_solver_api.cpp
// Stub implementation of Python API for wildfire_levelset solver
// ============================================================================

#include "fire_solver_api.H"

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
    
    try {
        g_fire_solver_state = std::make_unique<FireSolverState>();
        FireSolverState& state = *g_fire_solver_state;
        
        // Parse inputs
        state.inputs = std::make_unique<InputParameters>();
        parse_inputs(*state.inputs);
        
        // Create geometry
        IntVect dom_lo(0, 0);
        IntVect dom_hi(state.inputs->n_cell_x - 1, state.inputs->n_cell_y - 1);
        Box domain(dom_lo, dom_hi);
        RealBox rb({state.inputs->plo_x, state.inputs->plo_y},
                   {state.inputs->phi_x, state.inputs->phi_y});
        Array<int, AMREX_SPACEDIM> is_periodic{0, 0};
        
        state.geom = std::make_unique<Geometry>(domain, &rb, 0, is_periodic.data());
        
        // Create grids and distribution
        state.ba = std::make_unique<BoxArray>(domain);
        state.ba->maxSize(state.inputs->max_grid);
        state.dm = std::make_unique<DistributionMapping>(*state.ba);
        
        // Create fields
        state.fields = std::make_unique<WildfireFields>(*state.ba, *state.dm, *state.inputs);
        
        state.initialized = true;
        amrex::Print() << "Fire solver initialized successfully (stub implementation)\n";
        
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
        const Real dt_step = state.dt;
        state.time += dt_step;
        
        amrex::Print() << "Step " << state.step << ": Time = " << state.time
                       << ", dt = " << dt_step << " (stub)\n";
        
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
    double& dx, double& dy)
{
    if (!g_fire_solver_state || !g_fire_solver_state->initialized) {
        nx = ny = 0;
        xmin = ymin = 0.0;
        xmax = ymax = 0.0;
        dx = dy = 0.0;
        return;
    }
    
    nx = g_fire_solver_state->inputs->n_cell_x;
    ny = g_fire_solver_state->inputs->n_cell_y;
    xmin = g_fire_solver_state->inputs->plo_x;
    xmax = g_fire_solver_state->inputs->phi_x;
    ymin = g_fire_solver_state->inputs->plo_y;
    ymax = g_fire_solver_state->inputs->phi_y;
    
    dx = (xmax - xmin) / nx;
    dy = (ymax - ymin) / ny;
}

// ============================================================================
// fire_solver_update_wind
// ============================================================================
bool fire_solver_update_wind(
    int nx, int ny,
    const double* u_data,
    const double* v_data)
{
    if (!g_fire_solver_state || !g_fire_solver_state->initialized) {
        return false;
    }
    
    try {
        // Stub: just store the wind data
        amrex::Print() << "Wind updated: " << nx << " x " << ny << " (stub)\n";
        return true;
    } catch (...) {
        return false;
    }
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
    const double* w_data)
{
    if (!g_fire_solver_state || !g_fire_solver_state->initialized) {
        return false;
    }
    
    try {
        // Stub: just store the 3D wind data
        amrex::Print() << "3D Wind updated: " << nx << " x " << ny << " x " << nz << " (stub)\n";
        return true;
    } catch (...) {
        return false;
    }
}

// ============================================================================
// fire_solver_get_phi
// ============================================================================
std::vector<double> fire_solver_get_phi() {
    std::vector<double> result;
    
    if (!g_fire_solver_state || !g_fire_solver_state->initialized) {
        return result;
    }
    
    try {
        const MultiFab& phi = g_fire_solver_state->fields->phi;
        int nx = g_fire_solver_state->inputs->n_cell_x;
        int ny = g_fire_solver_state->inputs->n_cell_y;
        
        result.resize(nx * ny);
        
        // Copy phi data from MultiFab to vector (stub: just fill with zeros)
        std::fill(result.begin(), result.end(), 1.0);
        
        return result;
    } catch (...) {
        return result;
    }
}

// ============================================================================
// fire_solver_get_ros
// ============================================================================
std::vector<double> fire_solver_get_ros() {
    std::vector<double> result;
    
    if (!g_fire_solver_state || !g_fire_solver_state->initialized) {
        return result;
    }
    
    try {
        int nx = g_fire_solver_state->inputs->n_cell_x;
        int ny = g_fire_solver_state->inputs->n_cell_y;
        
        result.resize(nx * ny);
        std::fill(result.begin(), result.end(), 0.0);
        
        return result;
    } catch (...) {
        return result;
    }
}

// ============================================================================
// fire_solver_get_intensity
// ============================================================================
std::vector<double> fire_solver_get_intensity() {
    std::vector<double> result;
    
    if (!g_fire_solver_state || !g_fire_solver_state->initialized) {
        return result;
    }
    
    try {
        int nx = g_fire_solver_state->inputs->n_cell_x;
        int ny = g_fire_solver_state->inputs->n_cell_y;
        
        result.resize(nx * ny);
        std::fill(result.begin(), result.end(), 0.0);
        
        return result;
    } catch (...) {
        return result;
    }
}

// ============================================================================
// fire_solver_get_flame_length
// ============================================================================
std::vector<double> fire_solver_get_flame_length() {
    std::vector<double> result;
    
    if (!g_fire_solver_state || !g_fire_solver_state->initialized) {
        return result;
    }
    
    try {
        int nx = g_fire_solver_state->inputs->n_cell_x;
        int ny = g_fire_solver_state->inputs->n_cell_y;
        
        result.resize(nx * ny);
        std::fill(result.begin(), result.end(), 0.0);
        
        return result;
    } catch (...) {
        return result;
    }
}

// ============================================================================
// fire_solver_get_wind
// ============================================================================
void fire_solver_get_wind(std::vector<double>& u_data, std::vector<double>& v_data) {
    if (!g_fire_solver_state || !g_fire_solver_state->initialized) {
        return;
    }
    
    try {
        int nx = g_fire_solver_state->inputs->n_cell_x;
        int ny = g_fire_solver_state->inputs->n_cell_y;
        
        u_data.resize(nx * ny);
        v_data.resize(nx * ny);
        
        std::fill(u_data.begin(), u_data.end(), 0.0);
        std::fill(v_data.begin(), v_data.end(), 0.0);
    } catch (...) {
        // Silent failure
    }
}

// ============================================================================
// fire_solver_get_arrival_time
// ============================================================================
std::vector<double> fire_solver_get_arrival_time() {
    std::vector<double> result;
    
    if (!g_fire_solver_state || !g_fire_solver_state->initialized) {
        return result;
    }
    
    try {
        int nx = g_fire_solver_state->inputs->n_cell_x;
        int ny = g_fire_solver_state->inputs->n_cell_y;
        
        result.resize(nx * ny);
        std::fill(result.begin(), result.end(), -1.0);
        
        return result;
    } catch (...) {
        return result;
    }
}

// ============================================================================
// fire_solver_write_plotfile
// ============================================================================
bool fire_solver_write_plotfile(const std::string& plotfile_name) {
    if (!g_fire_solver_state || !g_fire_solver_state->initialized) {
        return false;
    }
    
    try {
        amrex::Print() << "Writing plotfile: " << plotfile_name << " (stub)\n";
        return true;
    } catch (...) {
        return false;
    }
}

// ============================================================================
// fire_solver_finalize
// ============================================================================
void fire_solver_finalize() {
    if (g_fire_solver_state) {
        g_fire_solver_state.reset();
        amrex::Print() << "Fire solver finalized\n";
    }
}

// ============================================================================
// fire_solver_is_initialized
// ============================================================================
bool fire_solver_is_initialized() {
    return g_fire_solver_state && g_fire_solver_state->initialized;
}

// ============================================================================
// PHASE 1: Core Configuration & Properties Functions
// ============================================================================

// ============================================================================
// fire_solver_get_config
// ============================================================================
std::map<std::string, double> fire_solver_get_config() {
    std::map<std::string, double> config;
    
    if (!g_fire_solver_state || !g_fire_solver_state->initialized) {
        return config;
    }
    
    auto& inputs = *g_fire_solver_state->inputs;
    config["nx"] = inputs.n_cell_x;
    config["ny"] = inputs.n_cell_y;
    config["xmin"] = inputs.plo_x;
    config["xmax"] = inputs.phi_x;
    config["ymin"] = inputs.plo_y;
    config["ymax"] = inputs.phi_y;
    config["dx"] = (inputs.phi_x - inputs.plo_x) / inputs.n_cell_x;
    config["dy"] = (inputs.phi_y - inputs.plo_y) / inputs.n_cell_y;
    config["time"] = g_fire_solver_state->time;
    config["step"] = g_fire_solver_state->step;
    config["dt"] = g_fire_solver_state->dt;
    
    return config;
}

// ============================================================================
// fire_solver_get_vertical_domain
// ============================================================================
void fire_solver_get_vertical_domain(double& zmin, double& zmax) {
    zmin = 0.0;
    zmax = 100.0;  // Default vertical domain
    
    if (g_fire_solver_state && g_fire_solver_state->initialized) {
        // Could read from inputs if vertical domain is configured
        zmin = 0.0;
        zmax = 100.0;  // Placeholder
    }
}

// ============================================================================
// fire_solver_get_rothermel_properties
// ============================================================================
std::map<std::string, double> fire_solver_get_rothermel_properties() {
    std::map<std::string, double> props;
    
    if (!g_fire_solver_state || !g_fire_solver_state->initialized) {
        return props;
    }
    
    // Placeholder Rothermel model properties (NFFL fuel model 1)
    props["model_number"] = 1;
    props["fuel_load_dead_1h"] = 3.0;  // tons/acre
    props["fuel_load_dead_10h"] = 2.0;
    props["fuel_load_dead_100h"] = 1.5;
    props["fuel_load_live_herbaceous"] = 0.0;
    props["fuel_load_live_woody"] = 0.5;
    props["fuel_moisture_dead_1h"] = 5.0;  // percent
    props["fuel_moisture_dead_10h"] = 6.0;
    props["fuel_moisture_dead_100h"] = 7.0;
    props["fuel_moisture_live_herbaceous"] = 50.0;
    props["surface_area_to_volume"] = 1500.0;  // ft^2/ft^3
    props["heat_content"] = 8000.0;  // BTU/lb
    props["mineral_content"] = 5.0;  // percent
    props["effective_mineral_content"] = 1.0;
    
    return props;
}

// ============================================================================
// fire_solver_get_wind_ros_relationship
// ============================================================================
std::map<std::string, double> fire_solver_get_wind_ros_relationship() {
    std::map<std::string, double> relationship;
    
    if (!g_fire_solver_state || !g_fire_solver_state->initialized) {
        return relationship;
    }
    
    // Placeholder wind-ROS relationship parameters
    relationship["B_coefficient"] = 0.0185;
    relationship["C_coefficient"] = 0.34;
    relationship["E_coefficient"] = 0.460;
    relationship["wind_reduction_factor"] = 0.4;  // Fraction of ref wind
    relationship["midflame_wind_speed"] = 0.0;
    relationship["wind_direction_effect"] = 1.0;
    
    return relationship;
}

// ============================================================================
// fire_solver_get_spread_parameters
// ============================================================================
std::map<std::string, double> fire_solver_get_spread_parameters() {
    std::map<std::string, double> params;
    
    if (!g_fire_solver_state || !g_fire_solver_state->initialized) {
        return params;
    }
    
    // Placeholder spread model parameters
    params["reaction_intensity"] = 0.0;  // BTU/ft^2/min
    params["propagation_index"] = 0.0;
    params["richardson_number"] = 0.0;
    params["level_set_thickness"] = 1.0;  // cells
    params["normal_speed"] = 0.0;  // m/s
    
    return params;
}

// ============================================================================
// fire_solver_update_rothermel_fuel_load
// ============================================================================
bool fire_solver_update_rothermel_fuel_load(
    const std::vector<double>& dead_load,
    const std::vector<double>& live_load)
{
    if (!g_fire_solver_state || !g_fire_solver_state->initialized) {
        return false;
    }
    
    try {
        // Placeholder: Update fuel load in solver state
        // In a real implementation, this would update the fuel model
        amrex::Print() << "Updated fuel loads: " << dead_load.size() 
                       << " dead categories, " << live_load.size() 
                       << " live categories\n";
        return true;
    } catch (...) {
        return false;
    }
}

// ============================================================================
// fire_solver_validate_domain_compatibility
// ============================================================================
std::map<std::string, bool> fire_solver_validate_domain_compatibility(
    int wind_nx, int wind_ny,
    double wind_xmin, double wind_xmax,
    double wind_ymin, double wind_ymax)
{
    std::map<std::string, bool> compatibility;
    compatibility["compatible"] = false;
    compatibility["x_bounds_match"] = false;
    compatibility["y_bounds_match"] = false;
    compatibility["resolution_match"] = false;
    
    if (!g_fire_solver_state || !g_fire_solver_state->initialized) {
        return compatibility;
    }
    
    auto& inputs = *g_fire_solver_state->inputs;
    
    // Check domain compatibility
    double tolerance = 1e-6;
    bool x_match = (std::abs(wind_xmin - inputs.plo_x) < tolerance &&
                    std::abs(wind_xmax - inputs.phi_x) < tolerance);
    bool y_match = (std::abs(wind_ymin - inputs.plo_y) < tolerance &&
                    std::abs(wind_ymax - inputs.phi_y) < tolerance);
    bool res_match = (wind_nx == inputs.n_cell_x && wind_ny == inputs.n_cell_y);
    
    compatibility["x_bounds_match"] = x_match;
    compatibility["y_bounds_match"] = y_match;
    compatibility["resolution_match"] = res_match;
    compatibility["compatible"] = (x_match && y_match && res_match);
    
    return compatibility;
}

// ============================================================================
// PHASE 2: Terrain & Spatial Features Functions
// ============================================================================

// ============================================================================
// fire_solver_update_terrain
// ============================================================================
bool fire_solver_update_terrain(
    const double* elevation,
    const double* slope,
    const double* aspect,
    int nx, int ny)
{
    if (!g_fire_solver_state || !g_fire_solver_state->initialized) {
        return false;
    }
    
    try {
        // Placeholder: Store terrain data in solver state
        // In real implementation, would update MultiFab with terrain fields
        amrex::Print() << "Updated terrain: " << nx << " x " << ny 
                       << " elevation/slope/aspect data\n";
        return true;
    } catch (...) {
        return false;
    }
}

// ============================================================================
// fire_solver_get_terrain_info
// ============================================================================
std::map<std::string, std::vector<double>> fire_solver_get_terrain_info() {
    std::map<std::string, std::vector<double>> terrain_info;
    
    if (!g_fire_solver_state || !g_fire_solver_state->initialized) {
        return terrain_info;
    }
    
    int total_cells = g_fire_solver_state->inputs->n_cell_x * 
                      g_fire_solver_state->inputs->n_cell_y;
    
    // Placeholder: Return flat terrain (elevation = 0, slope = 0, aspect = 0)
    terrain_info["elevation"] = std::vector<double>(total_cells, 0.0);
    terrain_info["slope"] = std::vector<double>(total_cells, 0.0);
    terrain_info["aspect"] = std::vector<double>(total_cells, 0.0);
    terrain_info["slope_ros_factor"] = std::vector<double>(total_cells, 1.0);
    
    return terrain_info;
}

// ============================================================================
// fire_solver_get_ros_at_location
// ============================================================================
double fire_solver_get_ros_at_location(double x, double y) {
    if (!g_fire_solver_state || !g_fire_solver_state->initialized) {
        return 0.0;
    }
    
    // Placeholder: Interpolate ROS at location
    // Would require bilinear interpolation from ROS field
    return 0.5;  // m/s - placeholder value
}

// ============================================================================
// fire_solver_interpolate_field
// ============================================================================
double fire_solver_interpolate_field(
    const std::string& field_name,
    double x, double y)
{
    if (!g_fire_solver_state || !g_fire_solver_state->initialized) {
        return 0.0;
    }
    
    // Placeholder: Interpolate requested field at location
    if (field_name == "phi") {
        return 1.0;  // Unburned
    } else if (field_name == "ros") {
        return 0.5;  // m/s
    } else if (field_name == "intensity") {
        return 100.0;  // kW/m
    }
    
    return 0.0;
}

// ============================================================================
// PHASE 3 & 4: Enhanced Fire State Fields
// ============================================================================

// ============================================================================
// fire_solver_get_ros_x
// ============================================================================
std::vector<double> fire_solver_get_ros_x() {
    std::vector<double> ros_x;
    
    if (!g_fire_solver_state || !g_fire_solver_state->initialized) {
        return ros_x;
    }
    
    int total_cells = g_fire_solver_state->inputs->n_cell_x * 
                      g_fire_solver_state->inputs->n_cell_y;
    ros_x = std::vector<double>(total_cells, 0.0);
    
    return ros_x;
}

// ============================================================================
// fire_solver_get_ros_y
// ============================================================================
std::vector<double> fire_solver_get_ros_y() {
    std::vector<double> ros_y;
    
    if (!g_fire_solver_state || !g_fire_solver_state->initialized) {
        return ros_y;
    }
    
    int total_cells = g_fire_solver_state->inputs->n_cell_x * 
                      g_fire_solver_state->inputs->n_cell_y;
    ros_y = std::vector<double>(total_cells, 0.0);
    
    return ros_y;
}

// ============================================================================
// fire_solver_get_ros_wind
// ============================================================================
std::vector<double> fire_solver_get_ros_wind() {
    std::vector<double> ros_wind;
    
    if (!g_fire_solver_state || !g_fire_solver_state->initialized) {
        return ros_wind;
    }
    
    int total_cells = g_fire_solver_state->inputs->n_cell_x * 
                      g_fire_solver_state->inputs->n_cell_y;
    ros_wind = std::vector<double>(total_cells, 0.0);
    
    return ros_wind;
}

// ============================================================================
// fire_solver_get_ros_slope
// ============================================================================
std::vector<double> fire_solver_get_ros_slope() {
    std::vector<double> ros_slope;
    
    if (!g_fire_solver_state || !g_fire_solver_state->initialized) {
        return ros_slope;
    }
    
    int total_cells = g_fire_solver_state->inputs->n_cell_x * 
                      g_fire_solver_state->inputs->n_cell_y;
    ros_slope = std::vector<double>(total_cells, 0.0);
    
    return ros_slope;
}

// ============================================================================
// fire_solver_get_residence_time
// ============================================================================
std::vector<double> fire_solver_get_residence_time() {
    std::vector<double> residence_time;
    
    if (!g_fire_solver_state || !g_fire_solver_state->initialized) {
        return residence_time;
    }
    
    int total_cells = g_fire_solver_state->inputs->n_cell_x * 
                      g_fire_solver_state->inputs->n_cell_y;
    residence_time = std::vector<double>(total_cells, 300.0);  // 5 minutes
    
    return residence_time;
}

// ============================================================================
// fire_solver_get_fuel_consumption
// ============================================================================
std::vector<double> fire_solver_get_fuel_consumption() {
    std::vector<double> consumption;
    
    if (!g_fire_solver_state || !g_fire_solver_state->initialized) {
        return consumption;
    }
    
    int total_cells = g_fire_solver_state->inputs->n_cell_x * 
                      g_fire_solver_state->inputs->n_cell_y;
    consumption = std::vector<double>(total_cells, 5.0);  // kg/m^2
    
    return consumption;
}

// ============================================================================
// fire_solver_get_front_curvature
// ============================================================================
std::vector<double> fire_solver_get_front_curvature() {
    std::vector<double> curvature;
    
    if (!g_fire_solver_state || !g_fire_solver_state->initialized) {
        return curvature;
    }
    
    int total_cells = g_fire_solver_state->inputs->n_cell_x * 
                      g_fire_solver_state->inputs->n_cell_y;
    curvature = std::vector<double>(total_cells, 0.0);
    
    return curvature;
}

// ============================================================================
// fire_solver_get_spread_direction
// ============================================================================
std::vector<double> fire_solver_get_spread_direction() {
    std::vector<double> direction;
    
    if (!g_fire_solver_state || !g_fire_solver_state->initialized) {
        return direction;
    }
    
    int total_cells = g_fire_solver_state->inputs->n_cell_x * 
                      g_fire_solver_state->inputs->n_cell_y;
    direction = std::vector<double>(total_cells, 0.0);  // radians
    
    return direction;
}

// ============================================================================
// PHASE 5: Advanced Ignition & Control Functions
// ============================================================================

// ============================================================================
// fire_solver_set_ignition_region
// ============================================================================
bool fire_solver_set_ignition_region(
    double xmin, double xmax,
    double ymin, double ymax,
    double time)
{
    if (!g_fire_solver_state || !g_fire_solver_state->initialized) {
        return false;
    }
    
    try {
        amrex::Print() << "Set ignition region: X=[" << xmin << ", " << xmax 
                       << "], Y=[" << ymin << ", " << ymax << "] at t=" << time << "\n";
        return true;
    } catch (...) {
        return false;
    }
}

// ============================================================================
// fire_solver_set_ignition_from_array
// ============================================================================
bool fire_solver_set_ignition_from_array(const double* phi_init) {
    if (!g_fire_solver_state || !g_fire_solver_state->initialized) {
        return false;
    }
    
    try {
        amrex::Print() << "Set ignition from custom phi field\n";
        return true;
    } catch (...) {
        return false;
    }
}

// ============================================================================
// fire_solver_set_spread_model
// ============================================================================
bool fire_solver_set_spread_model(const std::string& model_name) {
    if (!g_fire_solver_state || !g_fire_solver_state->initialized) {
        return false;
    }
    
    try {
        if (model_name == "levelset") {
            g_fire_solver_state->use_levelset = true;
            g_fire_solver_state->use_mtt = false;
            g_fire_solver_state->use_farsite = false;
        } else if (model_name == "richards") {
            g_fire_solver_state->use_levelset = false;
            g_fire_solver_state->use_mtt = true;
            g_fire_solver_state->use_farsite = false;
        } else if (model_name == "hybrid") {
            g_fire_solver_state->use_levelset = true;
            g_fire_solver_state->use_mtt = true;
            g_fire_solver_state->use_farsite = false;
        } else {
            return false;
        }
        amrex::Print() << "Spread model set to: " << model_name << "\n";
        return true;
    } catch (...) {
        return false;
    }
}

// ============================================================================
// fire_solver_step_with_subcycles
// ============================================================================
double fire_solver_step_with_subcycles(
    double target_dt,
    int max_subcycles)
{
    if (!g_fire_solver_state || !g_fire_solver_state->initialized) {
        return -1.0;
    }
    
    try {
        double dt_subcycle = target_dt / max_subcycles;
        for (int i = 0; i < max_subcycles; ++i) {
            g_fire_solver_state->time += dt_subcycle;
        }
        g_fire_solver_state->step++;
        amrex::Print() << "Subcycled: " << max_subcycles << " cycles, dt=" << target_dt << "\n";
        return target_dt;
    } catch (...) {
        return -1.0;
    }
}

// ============================================================================
// fire_solver_get_timestep_recommendation
// ============================================================================
double fire_solver_get_timestep_recommendation() {
    if (!g_fire_solver_state || !g_fire_solver_state->initialized) {
        return 0.01;
    }
    
    // Placeholder: Could compute adaptive timestep based on CFL, ROS, etc.
    double cfl = 0.8;
    double dx = g_fire_solver_state->inputs ? 
        (g_fire_solver_state->inputs->phi_x - g_fire_solver_state->inputs->plo_x) / 
        g_fire_solver_state->inputs->n_cell_x : 100.0;
    
    return cfl * dx / 5.0;  // Assuming ~5 m/s max ROS
}

// ============================================================================
// PHASE 6: Surface Fluxes & Emissions Functions
// ============================================================================

// ============================================================================
// fire_solver_get_all_surface_fluxes
// ============================================================================
std::map<std::string, std::vector<double>> fire_solver_get_all_surface_fluxes() {
    std::map<std::string, std::vector<double>> fluxes;
    
    if (!g_fire_solver_state || !g_fire_solver_state->initialized) {
        return fluxes;
    }
    
    int total_cells = g_fire_solver_state->inputs->n_cell_x * 
                      g_fire_solver_state->inputs->n_cell_y;
    
    // Placeholder surface fluxes
    fluxes["heat_flux"] = std::vector<double>(total_cells, 100.0);  // kW/m^2
    fluxes["sensible_heat"] = std::vector<double>(total_cells, 80.0);  // kW/m^2
    fluxes["latent_heat"] = std::vector<double>(total_cells, 20.0);  // kW/m^2
    fluxes["radiation"] = std::vector<double>(total_cells, 30.0);  // kW/m^2
    fluxes["momentum_flux"] = std::vector<double>(total_cells, 0.1);  // N/m^2
    fluxes["co2_flux"] = std::vector<double>(total_cells, 0.0001);  // kg/m^2/s
    fluxes["pm25_flux"] = std::vector<double>(total_cells, 0.00005);  // kg/m^2/s
    fluxes["smoke_height"] = std::vector<double>(total_cells, 500.0);  // m
    
    return fluxes;
}

// ============================================================================
// fire_solver_get_emission_factors
// ============================================================================
std::map<std::string, double> fire_solver_get_emission_factors() {
    std::map<std::string, double> factors;
    
    if (!g_fire_solver_state || !g_fire_solver_state->initialized) {
        return factors;
    }
    
    // Placeholder emission factors (per unit fuel consumed)
    factors["co2"] = 1.634;  // kg CO2 / kg fuel
    factors["co"] = 0.130;   // kg CO / kg fuel
    factors["ch4"] = 0.008;  // kg CH4 / kg fuel
    factors["pm25"] = 0.011; // kg PM2.5 / kg fuel
    factors["nox"] = 0.005;  // kg NOx / kg fuel
    factors["so2"] = 0.001;  // kg SO2 / kg fuel
    
    return factors;
}

// ============================================================================
// PHASE 7: Advanced I/O & Checkpointing Functions
// ============================================================================

// ============================================================================
// fire_solver_write_checkpoint
// ============================================================================
bool fire_solver_write_checkpoint(const std::string& filename) {
    if (!g_fire_solver_state || !g_fire_solver_state->initialized) {
        return false;
    }
    
    try {
        amrex::Print() << "Writing checkpoint to: " << filename << "\n";
        // Placeholder: Would write MultiFab and state to file
        return true;
    } catch (...) {
        return false;
    }
}

// ============================================================================
// fire_solver_read_checkpoint
// ============================================================================
bool fire_solver_read_checkpoint(const std::string& filename) {
    if (!g_fire_solver_state || !g_fire_solver_state->initialized) {
        return false;
    }
    
    try {
        amrex::Print() << "Reading checkpoint from: " << filename << "\n";
        // Placeholder: Would read MultiFab and state from file
        return true;
    } catch (...) {
        return false;
    }
}

// ============================================================================
// fire_solver_get_checkpoint_data
// ============================================================================
std::map<std::string, std::vector<double>> fire_solver_get_checkpoint_data() {
    std::map<std::string, std::vector<double>> checkpoint;
    
    if (!g_fire_solver_state || !g_fire_solver_state->initialized) {
        return checkpoint;
    }
    
    int total_cells = g_fire_solver_state->inputs->n_cell_x * 
                      g_fire_solver_state->inputs->n_cell_y;
    
    // Placeholder: Would include all fields needed to restart
    checkpoint["phi"] = std::vector<double>(total_cells, 1.0);
    checkpoint["ros"] = std::vector<double>(total_cells, 0.0);
    checkpoint["intensity"] = std::vector<double>(total_cells, 0.0);
    checkpoint["time_scalar"] = {g_fire_solver_state->time};
    checkpoint["step_scalar"] = {static_cast<double>(g_fire_solver_state->step)};
    
    return checkpoint;
}

// ============================================================================
// PHASE 8: Atmosphere Coupling & Diagnostics Functions
// ============================================================================

// ============================================================================
// fire_solver_set_fire_atmosphere_feedback_enabled
// ============================================================================
void fire_solver_set_fire_atmosphere_feedback_enabled(bool enabled) {
    if (g_fire_solver_state && g_fire_solver_state->initialized) {
        amrex::Print() << "Fire-atmosphere feedback: " 
                       << (enabled ? "ENABLED" : "DISABLED") << "\n";
        // Would set internal flag to enable/disable coupling
    }
}

// ============================================================================
// fire_solver_get_buoyancy_driven_winds
// ============================================================================
std::map<std::string, std::vector<double>> fire_solver_get_buoyancy_driven_winds() {
    std::map<std::string, std::vector<double>> winds;
    
    if (!g_fire_solver_state || !g_fire_solver_state->initialized) {
        return winds;
    }
    
    int total_cells = g_fire_solver_state->inputs->n_cell_x * 
                      g_fire_solver_state->inputs->n_cell_y;
    
    // Placeholder: Induced wind from fire plume
    winds["u_induced"] = std::vector<double>(total_cells, 0.0);
    winds["v_induced"] = std::vector<double>(total_cells, 0.0);
    winds["w_induced"] = std::vector<double>(total_cells, 0.0);
    
    return winds;
}

// ============================================================================
// fire_solver_get_coupling_statistics
// ============================================================================
std::map<std::string, double> fire_solver_get_coupling_statistics() {
    std::map<std::string, double> stats;
    
    if (!g_fire_solver_state || !g_fire_solver_state->initialized) {
        return stats;
    }
    
    // Placeholder coupling statistics
    stats["total_heat_release"] = 1e12;  // J
    stats["max_flame_height"] = 50.0;    // m
    stats["wind_speed_at_fire"] = 5.0;   // m/s
    stats["fire_to_wind_feedback_strength"] = 0.5;  // 0-1
    stats["ros_wind_coupling_factor"] = 0.8;  // 0-1
    stats["max_ros"] = 2.5;  // m/s
    
    return stats;
}

// ============================================================================
// PHASE 9: GPU & Performance Functions
// ============================================================================

// ============================================================================
// fire_solver_set_accelerated_ros_computation
// ============================================================================
void fire_solver_set_accelerated_ros_computation(bool enabled) {
    if (g_fire_solver_state && g_fire_solver_state->initialized) {
        amrex::Print() << "GPU acceleration: " 
                       << (enabled ? "ENABLED" : "DISABLED") << "\n";
        // Would enable/disable GPU offloading
    }
}

// ============================================================================
// fire_solver_profile_ros_calculation
// ============================================================================
std::map<std::string, double> fire_solver_profile_ros_calculation() {
    std::map<std::string, double> profile;
    
    if (!g_fire_solver_state || !g_fire_solver_state->initialized) {
        return profile;
    }
    
    // Placeholder profiling results (in milliseconds)
    profile["windfield_interpolation"] = 1.2;
    profile["rothermel_calc"] = 3.5;
    profile["terrain_effects"] = 0.8;
    profile["levelset_advection"] = 2.1;
    profile["total_ros_calc"] = 7.6;
    
    return profile;
}

// ============================================================================
// PHASE 10: Enhanced Error Handling & Diagnostics
// ============================================================================

// ============================================================================
// fire_solver_get_wind_at_surface
// ============================================================================
std::map<std::string, std::vector<double>> fire_solver_get_wind_at_surface() {
    std::map<std::string, std::vector<double>> wind_data;
    
    if (!g_fire_solver_state || !g_fire_solver_state->initialized) {
        return wind_data;
    }
    
    int total_cells = g_fire_solver_state->inputs->n_cell_x * 
                      g_fire_solver_state->inputs->n_cell_y;
    
    // Get wind components
    std::vector<double> u_data, v_data;
    fire_solver_get_wind(u_data, v_data);
    
    wind_data["u"] = u_data;
    wind_data["v"] = v_data;
    
    // Compute derived quantities
    std::vector<double> wind_speed(total_cells);
    std::vector<double> wind_direction(total_cells);
    
    for (int i = 0; i < total_cells; ++i) {
        wind_speed[i] = std::sqrt(u_data[i]*u_data[i] + v_data[i]*v_data[i]);
        wind_direction[i] = std::atan2(v_data[i], u_data[i]);  // radians
    }
    
    wind_data["wind_speed"] = wind_speed;
    wind_data["wind_direction"] = wind_direction;
    wind_data["w"] = std::vector<double>(total_cells, 0.0);
    
    return wind_data;
}
