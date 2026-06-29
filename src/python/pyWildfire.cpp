// pyWildfire.cpp - Python bindings for wildfire_levelset wind data interface
// 
// Provides Python interface to load 3D wind data from arrays (e.g., numpy or pyAMReX MultiFab)
// into the PltWindData structure used by wildfire_levelset for spotting calculations.
//
// Build with: cmake -S . -B build -DLEVELSET_BUILD_PYTHON_BINDINGS=ON
//             cmake --build build
//
// Usage from Python:
//   import pyWildfire
//   wind_data = pyWildfire.load_wind_from_arrays(nx, ny, nz, 
//                                                  xmin, xmax, ymin, ymax, zmin, zmax,
//                                                  u_array, v_array, w_array)

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <AMReX.H>
#include <AMReX_Print.H>

// Include the plt_wind_reader header which contains PltWindData and helper functions
#include "plt_wind_reader.H"

// Include the fire solver API for time-stepping and state management
#include "fire_solver_api.H"

namespace py = pybind11;
using namespace pybind11::literals;

// Wrapper function to load wind data from numpy arrays
// Returns a dictionary with the populated PltWindData fields
py::dict load_wind_from_arrays_py(
    int nx, int ny, int nz,
    double xmin, double xmax,
    double ymin, double ymax,
    double zmin, double zmax,
    py::array_t<double> u_array,
    py::array_t<double> v_array,
    py::object w_array_obj = py::none())
{
    // Initialize AMReX if not already initialized
    static bool amrex_initialized = false;
    if (!amrex_initialized) {
        int argc = 0;
        char** argv = nullptr;
        amrex::Initialize(argc, argv);
        amrex_initialized = true;
    }

    // Verify array shapes
    auto u_info = u_array.request();
    auto v_info = v_array.request();
    
    const int expected_size = nx * ny * nz;
    
    if (u_info.size != expected_size) {
        throw std::runtime_error("u_array size mismatch: expected " + 
                                 std::to_string(expected_size) + 
                                 ", got " + std::to_string(u_info.size));
    }
    if (v_info.size != expected_size) {
        throw std::runtime_error("v_array size mismatch: expected " + 
                                 std::to_string(expected_size) + 
                                 ", got " + std::to_string(v_info.size));
    }

    const double* u_ptr = static_cast<const double*>(u_info.ptr);
    const double* v_ptr = static_cast<const double*>(v_info.ptr);
    const double* w_ptr = nullptr;

    // Check if w_array was provided
    if (!w_array_obj.is_none()) {
        try {
            py::array_t<double> w_array = w_array_obj.cast<py::array_t<double>>();
            auto w_info = w_array.request();
            if (w_info.size != expected_size) {
                throw std::runtime_error("w_array size mismatch: expected " + 
                                         std::to_string(expected_size) + 
                                         ", got " + std::to_string(w_info.size));
            }
            w_ptr = static_cast<const double*>(w_info.ptr);
        } catch (const py::cast_error&) {
            throw std::runtime_error("w_array must be a numpy array if provided");
        }
    }
    else {
        // w_array is optional - will be filled with zeros in load function
        w_ptr = nullptr;
    }

    // Create PltWindData structure and populate it
    PltWindData pwd;
    bool success = load_plt_wind_from_arrays(
        nx, ny, nz,
        xmin, xmax, ymin, ymax, zmin, zmax,
        u_ptr, v_ptr, w_ptr,
        pwd
    );

    if (!success) {
        throw std::runtime_error("Failed to load wind data from arrays");
    }

    // Return results as a dictionary
    // Note: We return 2D column-averaged wind which is what wildfire_levelset uses
    py::dict result;
    result["valid"] = pwd.valid;
    result["nx_2d"] = pwd.nx_2d;
    result["ny_2d"] = pwd.ny_2d;
    result["n_points"] = pwd.n_points;
    result["xmin"] = pwd.lo2d[0];
    result["xmax"] = pwd.hi2d[0];
    result["ymin"] = pwd.lo2d[1];
    result["ymax"] = pwd.hi2d[1];
    
    // Convert 2D wind to numpy arrays (column-averaged horizontal wind)
    auto u2d_np = py::array_t<double>({pwd.ny_2d, pwd.nx_2d});
    auto v2d_np = py::array_t<double>({pwd.ny_2d, pwd.nx_2d});
    
    auto u2d_buf = u2d_np.mutable_unchecked<2>();
    auto v2d_buf = v2d_np.mutable_unchecked<2>();
    
    for (int j = 0; j < pwd.ny_2d; ++j) {
        for (int i = 0; i < pwd.nx_2d; ++i) {
            u2d_buf(j, i) = pwd.u2d[j * pwd.nx_2d + i];
            v2d_buf(j, i) = pwd.v2d[j * pwd.nx_2d + i];
        }
    }
    
    result["u2d"] = u2d_np;
    result["v2d"] = v2d_np;
    
    return result;
}

// Wrapper functions for fire solver API

py::dict fire_solver_init_py(const std::string& inputs_file) {
    bool success = fire_solver_initialize(inputs_file);
    
    py::dict result;
    result["success"] = success;
    
    if (success) {
        int nx, ny;
        double xmin, xmax, ymin, ymax, dx, dy;
        fire_solver_get_geometry(nx, ny, xmin, xmax, ymin, ymax, dx, dy);
        
        result["nx"] = nx;
        result["ny"] = ny;
        result["xmin"] = xmin;
        result["xmax"] = xmax;
        result["ymin"] = ymin;
        result["ymax"] = ymax;
        result["dx"] = dx;
        result["dy"] = dy;
    }
    
    return result;
}

py::dict fire_solver_step_py() {
    amrex::Real dt = fire_solver_advance();
    
    amrex::Real time;
    int step;
    fire_solver_get_time(time, step);
    
    py::dict result;
    result["success"] = (dt >= 0.0);
    result["dt"] = dt;
    result["time"] = time;
    result["step"] = step;
    
    return result;
}

py::dict fire_solver_get_state_py() {
    if (!fire_solver_is_initialized()) {
        return py::dict();
    }
    
    amrex::Real time;
    int step;
    fire_solver_get_time(time, step);
    
    int nx, ny;
    double xmin, xmax, ymin, ymax, dx, dy;
    fire_solver_get_geometry(nx, ny, xmin, xmax, ymin, ymax, dx, dy);
    
    // Extract phi (level set)
    std::vector<double> phi_vec = fire_solver_get_phi();
    auto phi_np = py::array_t<double>({ny, nx});
    auto phi_buf = phi_np.mutable_unchecked<2>();
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            phi_buf(j, i) = phi_vec[j * nx + i];
        }
    }
    
    // Extract ROS
    std::vector<double> ros_vec = fire_solver_get_ros();
    auto ros_np = py::array_t<double>({ny, nx});
    auto ros_buf = ros_np.mutable_unchecked<2>();
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            ros_buf(j, i) = ros_vec[j * nx + i];
        }
    }
    
    // Extract intensity
    std::vector<double> intensity_vec = fire_solver_get_intensity();
    auto intensity_np = py::array_t<double>({ny, nx});
    auto intensity_buf = intensity_np.mutable_unchecked<2>();
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            intensity_buf(j, i) = intensity_vec[j * nx + i];
        }
    }
    
    // Extract flame length
    std::vector<double> fl_vec = fire_solver_get_flame_length();
    auto fl_np = py::array_t<double>({ny, nx});
    auto fl_buf = fl_np.mutable_unchecked<2>();
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            fl_buf(j, i) = fl_vec[j * nx + i];
        }
    }
    
    // Extract wind
    std::vector<double> u_vec, v_vec;
    fire_solver_get_wind(u_vec, v_vec);
    auto u_np = py::array_t<double>({ny, nx});
    auto v_np = py::array_t<double>({ny, nx});
    auto u_buf = u_np.mutable_unchecked<2>();
    auto v_buf = v_np.mutable_unchecked<2>();
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            u_buf(j, i) = u_vec[j * nx + i];
            v_buf(j, i) = v_vec[j * nx + i];
        }
    }
    
    // Extract arrival time
    std::vector<double> arrival_vec = fire_solver_get_arrival_time();
    auto arrival_np = py::array_t<double>({ny, nx});
    auto arrival_buf = arrival_np.mutable_unchecked<2>();
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            arrival_buf(j, i) = arrival_vec[j * nx + i];
        }
    }
    
    py::dict result;
    result["time"] = time;
    result["step"] = step;
    result["nx"] = nx;
    result["ny"] = ny;
    result["xmin"] = xmin;
    result["xmax"] = xmax;
    result["ymin"] = ymin;
    result["ymax"] = ymax;
    result["dx"] = dx;
    result["dy"] = dy;
    result["phi"] = phi_np;
    result["ros"] = ros_np;
    result["intensity"] = intensity_np;
    result["flame_length"] = fl_np;
    result["u_wind"] = u_np;
    result["v_wind"] = v_np;
    result["arrival_time"] = arrival_np;
    
    return result;
}

bool fire_solver_update_wind_py(py::array_t<double> u_array, py::array_t<double> v_array) {
    auto u_info = u_array.request();
    auto v_info = v_array.request();
    
    if (u_info.ndim != 2 || v_info.ndim != 2) {
        throw std::runtime_error("Wind arrays must be 2D");
    }
    
    int ny = static_cast<int>(u_info.shape[0]);
    int nx = static_cast<int>(u_info.shape[1]);
    
    const double* u_ptr = static_cast<const double*>(u_info.ptr);
    const double* v_ptr = static_cast<const double*>(v_info.ptr);
    
    return fire_solver_update_wind(nx, ny, u_ptr, v_ptr);
}

bool fire_solver_update_wind_3d_py(
    int nx, int ny, int nz,
    double xmin, double xmax,
    double ymin, double ymax,
    double zmin, double zmax,
    py::array_t<double> u_array,
    py::array_t<double> v_array,
    py::array_t<double> w_array) {
    
    auto u_info = u_array.request();
    auto v_info = v_array.request();
    auto w_info = w_array.request();
    
    const double* u_ptr = static_cast<const double*>(u_info.ptr);
    const double* v_ptr = static_cast<const double*>(v_info.ptr);
    const double* w_ptr = static_cast<const double*>(w_info.ptr);
    
    return fire_solver_update_wind_3d(nx, ny, nz,
                                      xmin, xmax, ymin, ymax, zmin, zmax,
                                      u_ptr, v_ptr, w_ptr);
}

// Module definition
PYBIND11_MODULE(pyWildfire, m) {
    m.doc() = "Python bindings for wildfire_levelset: wind data interface and fire solver control";

    // ========================================================================
    // Wind data loading (existing functionality)
    // ========================================================================
    m.def("load_wind_from_arrays", &load_wind_from_arrays_py,
          py::arg("nx"), py::arg("ny"), py::arg("nz"),
          py::arg("xmin"), py::arg("xmax"),
          py::arg("ymin"), py::arg("ymax"),
          py::arg("zmin"), py::arg("zmax"),
          py::arg("u_array"),
          py::arg("v_array"),
          py::arg("w_array") = py::none(),
          R"pbdoc(
        Load 3D wind field data from numpy arrays and compute 2D column-averaged wind.

        This function takes 3D wind velocity components and computes the 2D column-averaged
        horizontal wind field used by wildfire_levelset for fire spotting calculations.

        Parameters
        ----------
        nx, ny, nz : int
            Grid dimensions in x, y, z directions
        xmin, xmax, ymin, ymax, zmin, zmax : float
            Physical domain bounds in meters (UTM coordinates)
        u_array : numpy.ndarray
            3D array of u-component velocities (m/s), shape (nz, ny, nx) or flattened
        v_array : numpy.ndarray
            3D array of v-component velocities (m/s), shape (nz, ny, nx) or flattened
        w_array : numpy.ndarray, optional
            3D array of w-component velocities (m/s), shape (nz, ny, nx) or flattened
            If not provided, w is assumed to be zero

        Returns
        -------
        dict
            Dictionary containing:
            - valid : bool - Whether data was loaded successfully
            - nx_2d, ny_2d : int - 2D grid dimensions
            - n_points : int - Total number of 3D points
            - xmin, xmax, ymin, ymax : float - 2D domain bounds
            - u2d : numpy.ndarray - Column-averaged u-component, shape (ny, nx)
            - v2d : numpy.ndarray - Column-averaged v-component, shape (ny, nx)

        Notes
        -----
        Arrays should be in Fortran (column-major) order with indexing [k, j, i]
        where k is the fastest-varying index (z-direction).

        Examples
        --------
        >>> import numpy as np
        >>> import pyWildfire
        >>> 
        >>> # Create synthetic 3D wind field
        >>> nx, ny, nz = 8, 8, 4
        >>> u = np.full((nz, ny, nx), 5.0)  # 5 m/s westerly
        >>> v = np.full((nz, ny, nx), 0.5)  # 0.5 m/s southerly
        >>> 
        >>> # Load and average
        >>> result = pyWildfire.load_wind_from_arrays(
        ...     nx, ny, nz,
        ...     329900.0, 330500.0,  # UTM X bounds
        ...     3774900.0, 3775500.0,  # UTM Y bounds
        ...     0.0, 40.0,  # Z bounds (meters AGL)
        ...     u.flatten(), v.flatten()
        ... )
        >>> 
        >>> print(result['u2d'].shape)  # (8, 8)
        >>> print(result['u2d'].mean())  # ~5.0 m/s
        )pbdoc");

    // ========================================================================
    // Fire solver control (new functionality)
    // ========================================================================
    m.def("initialize", &fire_solver_init_py,
          py::arg("inputs_file"),
          R"pbdoc(
        Initialize the wildfire solver from an inputs file.
        
        This sets up the simulation domain, loads landscape data, fuel models,
        and prepares all necessary fields for fire propagation.
        
        Parameters
        ----------
        inputs_file : str
            Path to the inputs file (e.g., "inputs.i")
        
        Returns
        -------
        dict
            Dictionary containing:
            - success : bool - Whether initialization succeeded
            - nx, ny : int - Grid dimensions
            - xmin, xmax, ymin, ymax : float - Domain bounds (meters)
            - dx, dy : float - Cell sizes (meters)
        
        Examples
        --------
        >>> import pyWildfire
        >>> result = pyWildfire.initialize("inputs.i")
        >>> if result['success']:
        ...     print(f"Initialized {result['nx']}x{result['ny']} grid")
        )pbdoc");
    
    m.def("advance", &fire_solver_step_py,
          R"pbdoc(
        Advance the fire simulation by one timestep.
        
        This computes the rate of spread, advects the fire front, updates
        fire behavior fields (intensity, flame length), and computes the
        next timestep based on CFL condition.
        
        Returns
        -------
        dict
            Dictionary containing:
            - success : bool - Whether timestep succeeded
            - dt : float - Timestep used (seconds)
            - time : float - Current simulation time (seconds)
            - step : int - Current timestep number
        
        Examples
        --------
        >>> result = pyWildfire.advance()
        >>> print(f"Stepped to t={result['time']} s with dt={result['dt']} s")
        )pbdoc");
    
    m.def("get_state", &fire_solver_get_state_py,
          R"pbdoc(
        Extract the current state of the fire simulation.
        
        Returns all fire-related fields as numpy arrays for visualization
        and analysis.
        
        Returns
        -------
        dict
            Dictionary containing:
            - time : float - Current simulation time (seconds)
            - step : int - Current timestep number
            - nx, ny : int - Grid dimensions
            - xmin, xmax, ymin, ymax : float - Domain bounds (meters)
            - dx, dy : float - Cell sizes (meters)
            - phi : ndarray - Level set field (ny, nx), < 0 burned, > 0 unburned
            - ros : ndarray - Rate of spread (m/s), shape (ny, nx)
            - intensity : ndarray - Fire line intensity (kW/m), shape (ny, nx)
            - flame_length : ndarray - Flame length (m), shape (ny, nx)
            - u_wind, v_wind : ndarray - Wind components (m/s), shape (ny, nx)
            - arrival_time : ndarray - Fire arrival time (s), shape (ny, nx)
        
        Examples
        --------
        >>> state = pyWildfire.get_state()
        >>> import matplotlib.pyplot as plt
        >>> plt.contourf(state['phi'], levels=[-1, 0, 1])
        >>> plt.title(f"Fire front at t={state['time']} s")
        >>> plt.show()
        )pbdoc");
    
    m.def("update_wind", &fire_solver_update_wind_py,
          py::arg("u_wind"),
          py::arg("v_wind"),
          R"pbdoc(
        Update the wind field from external 2D arrays.
        
        This allows coupling with external wind solvers by passing
        column-averaged wind components.
        
        Parameters
        ----------
        u_wind : ndarray
            U-component of wind (m/s), shape (ny, nx)
        v_wind : ndarray
            V-component of wind (m/s), shape (ny, nx)
        
        Returns
        -------
        bool
            True if wind was updated successfully
        
        Examples
        --------
        >>> import numpy as np
        >>> # Update to uniform 10 m/s easterly wind
        >>> u = np.full((ny, nx), 10.0)
        >>> v = np.zeros((ny, nx))
        >>> pyWildfire.update_wind(u, v)
        )pbdoc");
    
    m.def("update_wind_3d", &fire_solver_update_wind_3d_py,
          py::arg("nx"), py::arg("ny"), py::arg("nz"),
          py::arg("xmin"), py::arg("xmax"),
          py::arg("ymin"), py::arg("ymax"),
          py::arg("zmin"), py::arg("zmax"),
          py::arg("u_array"),
          py::arg("v_array"),
          py::arg("w_array"),
          R"pbdoc(
        Update the wind field from external 3D arrays.
        
        This allows coupling with 3D wind solvers (e.g., massconsistent_amr)
        by passing full 3D wind velocity fields. The wind is used for spotting
        calculations and column-averaged for surface fire spread.
        
        Parameters
        ----------
        nx, ny, nz : int
            Grid dimensions
        xmin, xmax, ymin, ymax, zmin, zmax : float
            Domain bounds (meters)
        u_array, v_array, w_array : ndarray
            3D wind velocity components (m/s), flattened in Fortran order
        
        Returns
        -------
        bool
            True if wind was updated successfully
        
        Examples
        --------
        >>> # Update from massconsistent_amr wind solver
        >>> u_3d = wind_solver.get_u()  # shape (nz, ny, nx)
        >>> v_3d = wind_solver.get_v()
        >>> w_3d = wind_solver.get_w()
        >>> pyWildfire.update_wind_3d(
        ...     nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax,
        ...     u_3d.flatten('F'), v_3d.flatten('F'), w_3d.flatten('F'))
        )pbdoc");
    
    m.def("write_plotfile", &fire_solver_write_plotfile,
          py::arg("plotfile_name"),
          R"pbdoc(
        Write the current state to an AMReX plotfile.
        
        Parameters
        ----------
        plotfile_name : str
            Name/path of the plotfile to write
        
        Returns
        -------
        bool
            True if plotfile was written successfully
        )pbdoc");
    
    m.def("finalize", &fire_solver_finalize,
          R"pbdoc(
        Clean up and finalize the fire solver.
        
        This releases all allocated memory and prepares for shutdown.
        After calling this, you must call initialize() again before
        using the solver.
        )pbdoc");
    
    m.def("is_initialized", &fire_solver_is_initialized,
          R"pbdoc(
        Check if the fire solver is currently initialized.
        
        Returns
        -------
        bool
            True if solver is initialized and ready to use
        )pbdoc");

    // ========================================================================
    // PHASE 1: Core Configuration & Properties (7 functions)
    // ========================================================================
    
    m.def("get_config", &fire_solver_get_config,
          R"pbdoc(
        Get comprehensive solver configuration.
        
        Returns
        -------
        dict
            Configuration with keys: nx, ny, xmin, xmax, ymin, ymax, dx, dy, time, step, dt
        )pbdoc");
    
    m.def("get_vertical_domain", 
        [](){ 
            double zmin, zmax;
            fire_solver_get_vertical_domain(zmin, zmax);
            return py::dict("zmin"_a=zmin, "zmax"_a=zmax);
        },
        R"pbdoc(
        Get vertical domain bounds.
        
        Returns
        -------
        dict
            Dict with 'zmin' and 'zmax' in meters
        )pbdoc");
    
    m.def("get_rothermel_properties", &fire_solver_get_rothermel_properties,
          R"pbdoc(
        Get Rothermel fuel model properties.
        
        Returns
        -------
        dict
            Fuel model properties for all size classes and moisture
        )pbdoc");
    
    m.def("get_wind_ros_relationship", &fire_solver_get_wind_ros_relationship,
          R"pbdoc(
        Get wind-ROS interaction parameters.
        
        Returns
        -------
        dict
            B-coefficient and other Rothermel wind parameters
        )pbdoc");
    
    m.def("get_spread_parameters", &fire_solver_get_spread_parameters,
          R"pbdoc(
        Get spread model parameters.
        
        Returns
        -------
        dict
            Parameters for Richards/Balbi model if enabled
        )pbdoc");
    
    m.def("update_rothermel_fuel_load",
        [](py::array_t<double> dead_load, py::array_t<double> live_load){
            auto dead_info = dead_load.request();
            auto live_info = live_load.request();
            const double* dead_ptr = static_cast<const double*>(dead_info.ptr);
            const double* live_ptr = static_cast<const double*>(live_info.ptr);
            std::vector<double> dead_vec(dead_ptr, dead_ptr + dead_info.size);
            std::vector<double> live_vec(live_ptr, live_ptr + live_info.size);
            return fire_solver_update_rothermel_fuel_load(dead_vec, live_vec);
        },
        py::arg("dead_load"), py::arg("live_load"),
        R"pbdoc(
        Update Rothermel fuel loads.
        
        Parameters
        ----------
        dead_load : ndarray
            Dead fuel loads (tons/acre)
        live_load : ndarray
            Live fuel loads (tons/acre)
        
        Returns
        -------
        bool
            True if update succeeded
        )pbdoc");
    
    m.def("validate_domain_compatibility",
        &fire_solver_validate_domain_compatibility,
        py::arg("wind_nx"), py::arg("wind_ny"),
        py::arg("wind_xmin"), py::arg("wind_xmax"),
        py::arg("wind_ymin"), py::arg("wind_ymax"),
        R"pbdoc(
        Validate compatibility between fire and wind domains.
        
        Parameters
        ----------
        wind_nx, wind_ny : int
            Wind grid dimensions
        wind_xmin, wind_xmax, wind_ymin, wind_ymax : float
            Wind domain bounds
        
        Returns
        -------
        dict
            Compatibility flags
        )pbdoc");

    // ========================================================================
    // PHASE 2: Terrain & Spatial Features (4 functions)
    // ========================================================================
    
    m.def("update_terrain",
        [](py::array_t<double> elevation, 
           py::array_t<double> slope,
           py::array_t<double> aspect,
           int nx, int ny){
            auto elev_info = elevation.request();
            auto slope_info = slope.request();
            auto aspect_info = aspect.request();
            const double* elev_ptr = static_cast<const double*>(elev_info.ptr);
            const double* slope_ptr = static_cast<const double*>(slope_info.ptr);
            const double* aspect_ptr = static_cast<const double*>(aspect_info.ptr);
            return fire_solver_update_terrain(elev_ptr, slope_ptr, aspect_ptr, nx, ny);
        },
        py::arg("elevation"), py::arg("slope"), py::arg("aspect"), 
        py::arg("nx"), py::arg("ny"),
        R"pbdoc(
        Update terrain data (elevation, slope, aspect).
        
        Parameters
        ----------
        elevation : ndarray
            Elevation map (m)
        slope : ndarray
            Slope (degrees)
        aspect : ndarray
            Aspect (degrees)
        nx, ny : int
            Grid dimensions
        
        Returns
        -------
        bool
            True if update succeeded
        )pbdoc");
    
    m.def("get_terrain_info", &fire_solver_get_terrain_info,
          R"pbdoc(
        Get terrain information fields.
        
        Returns
        -------
        dict
            Dict with elevation, slope, aspect arrays
        )pbdoc");
    
    m.def("get_ros_at_location",
        &fire_solver_get_ros_at_location,
        py::arg("x"), py::arg("y"),
        R"pbdoc(
        Query ROS at a specific location.
        
        Parameters
        ----------
        x, y : float
            Location coordinates (meters)
        
        Returns
        -------
        float
            ROS at location (m/s)
        )pbdoc");
    
    m.def("interpolate_field",
        &fire_solver_interpolate_field,
        py::arg("field_name"), py::arg("x"), py::arg("y"),
        R"pbdoc(
        Interpolate any field to a specific location.
        
        Parameters
        ----------
        field_name : str
            Field to interpolate ('phi', 'ros', 'intensity', etc.)
        x, y : float
            Location coordinates (meters)
        
        Returns
        -------
        float
            Field value at location
        )pbdoc");

    // ========================================================================
    // PHASE 3 & 4: Fire State Fields (8 functions)
    // ========================================================================
    
    m.def("get_ros_x", &fire_solver_get_ros_x,
          R"pbdoc(Get ROS x-component field as numpy array)pbdoc");
    
    m.def("get_ros_y", &fire_solver_get_ros_y,
          R"pbdoc(Get ROS y-component field as numpy array)pbdoc");
    
    m.def("get_ros_wind", &fire_solver_get_ros_wind,
          R"pbdoc(Get wind-driven ROS component as numpy array)pbdoc");
    
    m.def("get_ros_slope", &fire_solver_get_ros_slope,
          R"pbdoc(Get slope-driven ROS component as numpy array)pbdoc");
    
    m.def("get_residence_time", &fire_solver_get_residence_time,
          R"pbdoc(Get fuel residence time field as numpy array)pbdoc");
    
    m.def("get_fuel_consumption", &fire_solver_get_fuel_consumption,
          R"pbdoc(Get fuel consumption field as numpy array)pbdoc");
    
    m.def("get_front_curvature", &fire_solver_get_front_curvature,
          R"pbdoc(Get fire front curvature field as numpy array)pbdoc");
    
    m.def("get_spread_direction", &fire_solver_get_spread_direction,
          R"pbdoc(Get primary spread direction field as numpy array (radians))pbdoc");

    // ========================================================================
    // PHASE 5: Advanced Ignition & Control (5 functions)
    // ========================================================================
    
    m.def("set_ignition_region",
        &fire_solver_set_ignition_region,
        py::arg("xmin"), py::arg("xmax"), py::arg("ymin"), py::arg("ymax"), py::arg("time"),
        R"pbdoc(
        Set rectangular ignition region.
        
        Parameters
        ----------
        xmin, xmax, ymin, ymax : float
            Region bounds (meters)
        time : float
            Ignition time (seconds)
        
        Returns
        -------
        bool
            True if ignition was set
        )pbdoc");
    
    m.def("set_ignition_from_array",
        [](py::array_t<double> phi_init){
            auto info = phi_init.request();
            const double* ptr = static_cast<const double*>(info.ptr);
            return fire_solver_set_ignition_from_array(ptr);
        },
        py::arg("phi_init"),
        R"pbdoc(
        Set custom ignition pattern from level-set array.
        
        Parameters
        ----------
        phi_init : ndarray
            Initial level-set field
        
        Returns
        -------
        bool
            True if ignition was set
        )pbdoc");
    
    m.def("set_spread_model",
        &fire_solver_set_spread_model,
        py::arg("model_name"),
        R"pbdoc(
        Set propagation method.
        
        Parameters
        ----------
        model_name : str
            'levelset', 'richards', or 'hybrid'
        
        Returns
        -------
        bool
            True if model was set
        )pbdoc");
    
    m.def("step_with_subcycles",
        &fire_solver_step_with_subcycles,
        py::arg("target_dt"), py::arg("max_subcycles"),
        R"pbdoc(
        Advance with subcycling control.
        
        Parameters
        ----------
        target_dt : float
            Target timestep (seconds)
        max_subcycles : int
            Maximum number of subcycles
        
        Returns
        -------
        float
            Actual timestep used
        )pbdoc");
    
    m.def("get_timestep_recommendation", &fire_solver_get_timestep_recommendation,
          R"pbdoc(
        Get recommended next timestep.
        
        Returns
        -------
        float
            Recommended timestep (seconds)
        )pbdoc");

    // ========================================================================
    // PHASE 6: Surface Fluxes & Emissions (2 functions)
    // ========================================================================
    
    m.def("get_all_surface_fluxes", &fire_solver_get_all_surface_fluxes,
          R"pbdoc(
        Get all surface flux components.
        
        Returns
        -------
        dict
            Dict with heat_flux, sensible_heat, latent_heat, radiation, 
            momentum_flux, co2_flux, pm25_flux, smoke_height arrays
        )pbdoc");
    
    m.def("get_emission_factors", &fire_solver_get_emission_factors,
          R"pbdoc(
        Get emission factors for species (per unit fuel consumed).
        
        Returns
        -------
        dict
            Dict with co2, co, ch4, pm25, nox, so2 factors
        )pbdoc");

    // ========================================================================
    // PHASE 7: Advanced I/O & Checkpointing (3 functions)
    // ========================================================================
    
    m.def("write_checkpoint",
        &fire_solver_write_checkpoint,
        py::arg("filename"),
        R"pbdoc(
        Write checkpoint to file.
        
        Parameters
        ----------
        filename : str
            Checkpoint file path
        
        Returns
        -------
        bool
            True if write succeeded
        )pbdoc");
    
    m.def("read_checkpoint",
        &fire_solver_read_checkpoint,
        py::arg("filename"),
        R"pbdoc(
        Read checkpoint from file.
        
        Parameters
        ----------
        filename : str
            Checkpoint file path
        
        Returns
        -------
        bool
            True if read succeeded
        )pbdoc");
    
    m.def("get_checkpoint_data", &fire_solver_get_checkpoint_data,
          R"pbdoc(
        Get checkpoint-compatible data dictionary.
        
        Returns
        -------
        dict
            Dict with all fields needed to restart
        )pbdoc");

    // ========================================================================
    // PHASE 8: Atmosphere Coupling & Diagnostics (3 functions)
    // ========================================================================
    
    m.def("set_fire_atmosphere_feedback_enabled",
        &fire_solver_set_fire_atmosphere_feedback_enabled,
        py::arg("enabled"),
        R"pbdoc(
        Enable/disable fire-atmosphere feedback.
        
        Parameters
        ----------
        enabled : bool
            True to enable feedback
        )pbdoc");
    
    m.def("get_buoyancy_driven_winds", &fire_solver_get_buoyancy_driven_winds,
          R"pbdoc(
        Get induced wind from fire plume.
        
        Returns
        -------
        dict
            Dict with u_induced, v_induced, w_induced fields
        )pbdoc");
    
    m.def("get_coupling_statistics", &fire_solver_get_coupling_statistics,
          R"pbdoc(
        Get fire-atmosphere coupling statistics.
        
        Returns
        -------
        dict
            Dict with total_heat_release, max_flame_height, wind_speed_at_fire, etc.
        )pbdoc");

    // ========================================================================
    // PHASE 9: GPU & Performance (2 functions)
    // ========================================================================
    
    m.def("set_accelerated_ros_computation",
        &fire_solver_set_accelerated_ros_computation,
        py::arg("enabled"),
        R"pbdoc(
        Enable/disable GPU acceleration for ROS computation.
        
        Parameters
        ----------
        enabled : bool
            True to enable GPU acceleration
        )pbdoc");
    
    m.def("profile_ros_calculation", &fire_solver_profile_ros_calculation,
          R"pbdoc(
        Profile ROS calculation bottlenecks.
        
        Returns
        -------
        dict
            Dict with timing breakdown (milliseconds)
        )pbdoc");

    // ========================================================================
    // PHASE 10: Enhanced Diagnostics (1 function)
    // ========================================================================
    
    m.def("get_wind_at_surface", &fire_solver_get_wind_at_surface,
          R"pbdoc(
        Get wind at surface with derived fields.
        
        Returns
        -------
        dict
            Dict with u, v, w, wind_speed, wind_direction arrays
        )pbdoc");

    // Add version information
    m.attr("__version__") = "0.3.0";
}
