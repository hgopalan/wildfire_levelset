"""
pyWildfire - Python bindings for wildfire_levelset

This module provides Python interface to load 3D wind field data into 
wildfire_levelset's wind reading infrastructure. It enables:

1. Loading wind data from numpy arrays
2. Computing 2D column-averaged wind for fire spotting calculations
3. Future: Direct integration with pyAMReX MultiFab objects

Example usage:
    >>> import numpy as np
    >>> import pyWildfire
    >>> 
    >>> # Create synthetic 3D wind field
    >>> nx, ny, nz = 8, 8, 4
    >>> u = np.full((nz, ny, nx), 5.0)  # 5 m/s westerly
    >>> v = np.full((nz, ny, nx), 0.5)  # 0.5 m/s southerly
    >>> 
    >>> # Load and compute column-averaged wind
    >>> result = pyWildfire.load_wind_from_arrays(
    ...     nx, ny, nz,
    ...     329900.0, 330500.0,  # UTM X bounds (meters)
    ...     3774900.0, 3775500.0,  # UTM Y bounds (meters)
    ...     0.0, 40.0,  # Z bounds (meters AGL)
    ...     u.flatten('F'), v.flatten('F')  # Fortran order
    ... )
    >>> 
    >>> print(f"2D wind grid: {result['nx_2d']} × {result['ny_2d']}")
    >>> print(f"Mean u-wind: {result['u2d'].mean():.2f} m/s")
"""

from .pyWildfire import load_wind_from_arrays, __version__

__all__ = ['load_wind_from_arrays', '__version__']
