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
