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

    // Add version information
    m.attr("__version__") = "0.2.0";
}
