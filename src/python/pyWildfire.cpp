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
    py::array_t<double> w_array = py::none())
{
    // Initialize AMReX if not already initialized
    static bool amrex_initialized = false;
    if (!amrex_initialized) {
        amrex::Initialize(0, nullptr);
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

    if (!w_array.is_none()) {
        auto w_info = w_array.request();
        if (w_info.size != expected_size) {
            throw std::runtime_error("w_array size mismatch: expected " + 
                                     std::to_string(expected_size) + 
                                     ", got " + std::to_string(w_info.size));
        }
        w_ptr = static_cast<const double*>(w_info.ptr);
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

// Module definition
PYBIND11_MODULE(pyWildfire, m) {
    m.doc() = "Python bindings for wildfire_levelset wind data interface";

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

    // Add version information
    m.attr("__version__") = "0.1.0";
}
