#!/usr/bin/env python3
"""
wind_solver_interface.py - Abstract interface for wind solver coupling

Defines the contract that any wind solver implementation must follow
for integration with wildfire_levelset in a two-way coupling framework.
"""

from abc import ABC, abstractmethod
from typing import Dict, Tuple, Optional
import numpy as np


class WindSolverInterface(ABC):
    """
    Abstract base class for wind solvers coupled to wildfire simulations.

    Any wind solver (massconsistent_amr, WRF, etc.) that couples with
    wildfire_levelset should inherit from this class and implement
    the abstract methods.
    """

    @abstractmethod
    def initialize(self, config_file: str) -> bool:
        """
        Initialize the wind solver from a configuration file.

        Parameters:
            config_file: Path to wind solver configuration file

        Returns:
            True if initialization succeeded, False otherwise
        """
        pass

    @abstractmethod
    def solve(self, time: float, heat_source: Optional[Dict] = None) -> bool:
        """
        Compute wind field at the given simulation time.

        Parameters:
            time: Simulation time [seconds]
            heat_source: Optional dictionary containing heat flux data from fire:
                        - 'heat_flux': (ny, nx) array [W/m²]
                        - 'buoyancy': (ny, nx) array [K]
                        - 'fire_footprint': (ny, nx) binary mask
                        When provided, wind solver should account for heat forcing
                        (buoyancy-driven flow).

        Returns:
            True if solve succeeded, False otherwise
        """
        pass

    @abstractmethod
    def get_velocity_arrays(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Extract computed 3D velocity field from wind solver.

        Returns:
            (u_3d, v_3d, w_3d): Tuple of 3D velocity arrays
                - u_3d: (nz, ny, nx) u-component [m/s]
                - v_3d: (nz, ny, nx) v-component [m/s]
                - w_3d: (nz, ny, nx) w-component [m/s]
        """
        pass

    @abstractmethod
    def get_domain_info(self) -> Dict[str, float]:
        """
        Get domain information from wind solver.

        Returns:
            Dictionary with keys:
            - 'nx', 'ny', 'nz': Grid dimensions
            - 'xmin', 'xmax', 'ymin', 'ymax': Horizontal domain bounds [m]
            - 'zmin', 'zmax': Vertical domain bounds [m]
            - 'dx', 'dy': Horizontal grid spacing [m]
            - 'dz': Vertical grid spacing [m]
        """
        pass

    @abstractmethod
    def finalize(self) -> bool:
        """
        Clean up wind solver (deallocate memory, close files, etc.).

        Returns:
            True if finalization succeeded
        """
        pass

    def is_heat_aware(self) -> bool:
        """
        Check if this wind solver supports heat source forcing.

        Override in subclass if solver supports buoyancy forcing.

        Returns:
            True if solver can accept and process heat_source parameter
        """
        return False


class SyntheticWindSolver(WindSolverInterface):
    """
    Synthetic wind solver for testing and demonstration.

    Generates a simple time-varying wind field without reading external data.
    Useful for testing two-way coupling logic before integrating real wind solver.
    """

    def __init__(self):
        """Initialize synthetic wind solver."""
        self.initialized = False
        self.nx = 32
        self.ny = 32
        self.nz = 8
        self.xmin = 330000.0
        self.xmax = 330960.0
        self.ymin = 3775000.0
        self.ymax = 3775960.0
        self.zmin = 0.0
        self.zmax = 100.0
        self.dx = (self.xmax - self.xmin) / self.nx
        self.dy = (self.ymax - self.ymin) / self.ny
        self.dz = (self.zmax - self.zmin) / self.nz
        self.current_time = 0.0
        self.u_3d = None
        self.v_3d = None
        self.w_3d = None

    def initialize(self, config_file: str = None) -> bool:
        """Initialize synthetic wind solver."""
        self.initialized = True
        print("✓ Synthetic wind solver initialized")
        return True

    def solve(self, time: float, heat_source: Optional[Dict] = None) -> bool:
        """
        Generate synthetic wind field.

        Parameters:
            time: Simulation time [seconds]
            heat_source: Optional heat flux data (for demonstration, ignored)

        Returns:
            True if solve succeeded
        """
        if not self.initialized:
            return False

        self.current_time = time

        # Create coordinate arrays
        z = np.linspace(self.zmin, self.zmax, self.nz)
        y = np.linspace(self.ymin, self.ymax, self.ny)
        x = np.linspace(self.xmin, self.xmax, self.nx)

        # Initialize 3D wind
        self.u_3d = np.zeros((self.nz, self.ny, self.nx))
        self.v_3d = np.zeros((self.nz, self.ny, self.nx))
        self.w_3d = np.zeros((self.nz, self.ny, self.nx))

        # Base wind speed (5 m/s westerly with diurnal variation)
        hour = (time / 3600.0) % 24.0
        diurnal_factor = 0.7 + 0.3 * np.sin(2 * np.pi * (hour - 6.0) / 24.0)
        u_ref = 5.0 * diurnal_factor
        v_ref = 1.0 * diurnal_factor

        # Log-law wind profile
        z_ref = 10.0  # Reference height
        z0 = 0.1  # Roughness length
        kappa = 0.41  # von Kármán constant

        for k in range(self.nz):
            z_k = z[k]
            if z_k > z0:
                u_scale = np.log(z_k / z0) / np.log(z_ref / z0)
                v_scale = np.log(z_k / z0) / np.log(z_ref / z0)
            else:
                u_scale = 0.0
                v_scale = 0.0

            # Add spatial variation
            for j in range(self.ny):
                for i in range(self.nx):
                    x_norm = (x[i] - self.xmin) / (self.xmax - self.xmin)
                    y_norm = (y[j] - self.ymin) / (self.ymax - self.ymin)

                    # Terrain channeling effect (synthetic)
                    spatial_mod = 1.0 + 0.2 * np.sin(4 * np.pi * x_norm)

                    self.u_3d[k, j, i] = u_ref * u_scale * spatial_mod
                    self.v_3d[k, j, i] = v_ref * v_scale

            # Small vertical velocity
            self.w_3d[k, :, :] = 0.05 * np.sin(np.pi * z_k / self.zmax)

        # Optional: Modify wind based on heat source
        if heat_source is not None and 'heat_flux' in heat_source:
            self._apply_heat_forcing(heat_source)

        return True

    def _apply_heat_forcing(self, heat_source: Dict):
        """
        Apply simple buoyancy forcing from heat source.

        Reduces horizontal wind and increases vertical velocity
        in high-heat regions.

        Parameters:
            heat_source: Dictionary with 'heat_flux', 'buoyancy', 'fire_footprint'
        """
        heat_flux = heat_source.get('heat_flux', np.zeros((self.ny, self.nx)))
        buoyancy = heat_source.get('buoyancy', np.zeros((self.ny, self.nx)))
        fire_footprint = heat_source.get('fire_footprint', np.zeros((self.ny, self.nx)))

        if np.max(heat_flux) < 1.0:  # No significant heat
            return

        # Normalize heat flux
        heat_norm = np.clip(heat_flux / np.max(heat_flux), 0.0, 1.0)

        # Modify wind field: reduce horizontal speed by ~15% in high-heat areas
        wind_reduction = 0.15 * heat_norm

        for k in range(self.nz):
            self.u_3d[k, :, :] *= (1.0 - wind_reduction)
            self.v_3d[k, :, :] *= (1.0 - wind_reduction)

            # Increase vertical velocity
            self.w_3d[k, :, :] += 0.2 * np.sin(np.pi * (k + 0.5) / self.nz) * heat_norm

    def get_velocity_arrays(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Extract 3D velocity arrays."""
        return self.u_3d, self.v_3d, self.w_3d

    def get_domain_info(self) -> Dict[str, float]:
        """Get wind solver domain information."""
        return {
            'nx': self.nx,
            'ny': self.ny,
            'nz': self.nz,
            'xmin': self.xmin,
            'xmax': self.xmax,
            'ymin': self.ymin,
            'ymax': self.ymax,
            'zmin': self.zmin,
            'zmax': self.zmax,
            'dx': self.dx,
            'dy': self.dy,
            'dz': self.dz,
        }

    def finalize(self) -> bool:
        """Clean up wind solver."""
        self.initialized = False
        return True

    def is_heat_aware(self) -> bool:
        """This solver can respond to heat forcing."""
        return True


class MockWindSolver(WindSolverInterface):
    """
    Mock wind solver that does nothing (for testing state extraction).
    """

    def __init__(self):
        """Initialize mock wind solver."""
        self.initialized = False

    def initialize(self, config_file: str = None) -> bool:
        """Initialize mock solver."""
        self.initialized = True
        return True

    def solve(self, time: float, heat_source: Optional[Dict] = None) -> bool:
        """Mock solve (do nothing)."""
        return self.initialized

    def get_velocity_arrays(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return dummy velocity arrays."""
        return np.zeros((5, 32, 32)), np.zeros((5, 32, 32)), np.zeros((5, 32, 32))

    def get_domain_info(self) -> Dict[str, float]:
        """Return dummy domain info."""
        return {
            'nx': 32, 'ny': 32, 'nz': 5,
            'xmin': 330000.0, 'xmax': 330960.0,
            'ymin': 3775000.0, 'ymax': 3775960.0,
            'zmin': 0.0, 'zmax': 100.0,
            'dx': 30.0, 'dy': 30.0, 'dz': 20.0,
        }

    def finalize(self) -> bool:
        """Clean up."""
        self.initialized = False
        return True


if __name__ == "__main__":
    print("Wind Solver Interface - Unit Test")
    print("=" * 70)

    # Test synthetic wind solver
    wind = SyntheticWindSolver()
    wind.initialize()

    # Solve without heat
    wind.solve(0.0)
    u, v, w = wind.get_velocity_arrays()
    print(f"\nWind field (no heat):")
    print(f"  u: mean={u.mean():.2f} m/s, max={u.max():.2f} m/s")
    print(f"  v: mean={v.mean():.2f} m/s, max={v.max():.2f} m/s")
    print(f"  w: mean={w.mean():.3f} m/s, max={w.max():.3f} m/s")

    # Solve with heat
    heat_source = {
        'heat_flux': np.ones((32, 32)) * 50000.0,  # 50 kW/m²
        'buoyancy': np.ones((32, 32)) * 5.0,  # 5 K
        'fire_footprint': np.ones((32, 32)),
    }
    wind.solve(0.0, heat_source=heat_source)
    u, v, w = wind.get_velocity_arrays()
    print(f"\nWind field (with heat):")
    print(f"  u: mean={u.mean():.2f} m/s, max={u.max():.2f} m/s")
    print(f"  v: mean={v.mean():.2f} m/s, max={v.max():.2f} m/s")
    print(f"  w: mean={w.mean():.3f} m/s, max={w.max():.3f} m/s")

    wind.finalize()
    print("\n✓ Wind solver interface working correctly")
