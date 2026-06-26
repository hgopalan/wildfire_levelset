// ============================================================================
// fire_solver_api.H
// API for Python bindings to control wildfire_levelset solver from Python
//
// This header provides functions to:
// - Initialize the fire solver from an inputs file
// - Advance the simulation by one timestep
// - Extract MultiFab data (fire state, wind, etc.) as numpy arrays
// - Update wind field from external sources (e.g., massconsistent_amr)
// - Clean up solver state
//
// The solver state is stored in a global singleton structure to persist
// between Python function calls.
// ============================================================================

#ifndef FIRE_SOLVER_API_H
#define FIRE_SOLVER_API_H

#include "wildfire_includes.H"
#include "parse_inputs.H"
#include "multifab_setup.H"
#include "fmc_schedule.H"
#include "herb_moisture_schedule.H"
#include "wind_dir_schedule.H"
#include "multi_wtr_weather.H"
#include "satellite_assimilation.H"
#include "ignition_schedule.H"
#include "cruz_crown_model.H"
#include <memory>
#include <string>

// ============================================================================
// Forward declarations and placeholder types for Python bindings
// ============================================================================

// Placeholder for fuel depletion tracking state
struct FuelDepletionTracker {
    // Placeholder implementation to allow compilation with pybind11 bindings
    // When fuel depletion is implemented, this will contain residual fuel state
};

// Placeholder for ignition schedule state
struct IgnitionScheduleState {
    // Placeholder implementation to allow compilation with pybind11 bindings
    // When integrated with scheduler, this will track firing status
};

// ============================================================================
// FireSolverState
//
// Singleton structure that holds all solver state between Python calls.
// This allows initialization once, then multiple time-step advances.
// ============================================================================
struct FireSolverState {
    // Initialization flag
    bool initialized = false;
    
    // Input parameters
    std::unique_ptr<InputParameters> inputs;
    
    // Geometry and distribution
    std::unique_ptr<amrex::Geometry> geom;
    std::unique_ptr<amrex::BoxArray> ba;
    std::unique_ptr<amrex::DistributionMapping> dm;
    
    // All field data (phi, velocity, ROS, etc.)
    std::unique_ptr<WildfireFields> fields;
    
    // Time tracking
    amrex::Real time = 0.0;
    amrex::Real dt = 0.01;
    int step = 0;
    
    // Additional state for wind field management
    std::vector<amrex::Real> wind_x_data1, wind_y_data1, wind_u_data1, wind_v_data1;
    std::vector<amrex::Real> wind_x_data2, wind_y_data2, wind_u_data2, wind_v_data2;
    int current_wind_field_index = -1;
    int next_wind_field_index = -1;
    
    // Turbulent wind state
    std::unique_ptr<amrex::MultiFab> vel_base;
    std::unique_ptr<amrex::MultiFab> ou_state_mf;
    TurbWindState turb_state;
    bool turb_wind_active = false;
    
    // 3D wind data for spotting
    PltWindData albini_plt_wind;
    PltWindData ember_cascade_plt_wind;
    
    // Fire acceleration state
    std::unique_ptr<amrex::MultiFab> accel_state_mf;
    
    // Terrain and landscape data
    bool has_spatial_crown = false;
    bool has_spatial_moisture = false;
    
    // Fuel lookup tables
    std::vector<RothermelComputed> h_fuel_table;
    amrex::Gpu::DeviceVector<RothermelComputed> d_fuel_table;
    int fuel_table_size = 0;
    RothermelComputed* d_fuel_table_ptr = nullptr;
    
    // Balbi lookup table
    std::vector<BalbiComputed> h_balbi_table;
    amrex::Gpu::DeviceVector<BalbiComputed> d_balbi_table;
    int balbi_table_size = 0;
    BalbiComputed* d_balbi_table_ptr = nullptr;
    BalbiComputed bc_global_default;
    
    // FMC schedule
    std::vector<FMCRecord> fmc_sched;
    
    // Herb moisture schedule
    std::vector<HerbMoistureRecord> herb_sched;
    
    // Wind direction schedule
    std::vector<WindDirRecord> wind_dir_sched;
    
    // Multi-station weather
    MultiWtrWeather multi_wtr;
    bool multi_wtr_active = false;
    
    // Barrier cells
    std::vector<amrex::IntVect> barrier_cells;
    
    // Satellite assimilation state
    SatelliteState sat_state;
    
    // Horizon data for solar radiation
    std::unique_ptr<amrex::MultiFab> horizon_mf;
    
    // Fuel depletion tracking
    FuelDepletionTracker ftd;
    
    // Propagation method flags
    bool use_levelset = false;
    bool use_mtt = false;
    bool use_farsite = false;
    
    // Heat flux active flag
    bool heat_flux_active = false;
    bool wind_terrain_modifies_vel = false;
    
    // Viegas flags
    bool use_balbi_for_viegas = false;
    
    // Ignition schedule state
    IgnitionScheduleState ign_sched_state;
    
    // Crown fire computed parameters
    CruzCrownComputed ccc_global;
};

// Global singleton instance
extern std::unique_ptr<FireSolverState> g_fire_solver_state;

// ============================================================================
// API Functions
// ============================================================================

// Initialize the fire solver from an inputs file
// Returns true on success, false on failure
bool fire_solver_initialize(const std::string& inputs_file);

// Advance the fire simulation by one timestep
// Returns the timestep used (dt) or negative value on error
amrex::Real fire_solver_advance();

// Get current simulation time and step number
void fire_solver_get_time(amrex::Real& time, int& step);

// Get geometry information (grid dimensions, domain bounds, cell sizes)
void fire_solver_get_geometry(
    int& nx, int& ny,
    double& xmin, double& xmax,
    double& ymin, double& ymax,
    double& dx, double& dy);

// Update wind field from external 2D arrays (column-averaged wind)
// This allows coupling with external wind solvers
bool fire_solver_update_wind(
    int nx, int ny,
    const double* u_data,
    const double* v_data);

// Update wind field from external 3D arrays (full 3D wind)
// This allows coupling with 3D wind solvers like massconsistent_amr
bool fire_solver_update_wind_3d(
    int nx, int ny, int nz,
    double xmin, double xmax,
    double ymin, double ymax,
    double zmin, double zmax,
    const double* u_data,
    const double* v_data,
    const double* w_data);

// Extract fire front level set field (phi) as a flattened array
// phi < 0: burned, phi > 0: unburned, phi ≈ 0: fire front
// Returns data in Fortran (column-major) order for numpy compatibility
std::vector<double> fire_solver_get_phi();

// Extract rate of spread (ROS) field as a flattened array
std::vector<double> fire_solver_get_ros();

// Extract fire line intensity field as a flattened array
std::vector<double> fire_solver_get_intensity();

// Extract flame length field as a flattened array
std::vector<double> fire_solver_get_flame_length();

// Extract current wind field as flattened arrays
void fire_solver_get_wind(std::vector<double>& u_data, std::vector<double>& v_data);

// Extract arrival time field as a flattened array
std::vector<double> fire_solver_get_arrival_time();

// Write plotfile at current time
bool fire_solver_write_plotfile(const std::string& plotfile_name);

// Clean up and finalize the fire solver
void fire_solver_finalize();

// Check if solver is initialized
bool fire_solver_is_initialized();

#endif // FIRE_SOLVER_API_H
