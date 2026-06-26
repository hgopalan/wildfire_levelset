#!/usr/bin/env python3
"""
surface_flux_extractor.py - Extract surface heat fluxes from fire solver state

Converts fire solver output (intensity, flame length, ROS) into heat flux terms
that can be passed to a wind solver with buoyancy forcing.

Physical relationships:
- Intensity [kW/m] = heat release per unit fire perimeter
- Column-integrated flux [kW/m²] ≈ intensity * flame_front_density
- Sensible heat flux [W/m²] = intensity * conversion_factor
"""

import numpy as np
from typing import Dict, Tuple


class SurfaceFluxExtractor:
    """Extract heat fluxes from wildfire solver state for wind coupling."""

    def __init__(self, fire_grid_shape: Tuple[int, int], dx: float, dy: float):
        """
        Initialize flux extractor.

        Parameters:
            fire_grid_shape: (ny, nx) grid dimensions
            dx, dy: Grid cell sizes (meters)
        """
        self.ny, self.nx = fire_grid_shape
        self.dx = dx
        self.dy = dy
        self.cell_area = dx * dy

    def compute_column_heat_flux(self, state: Dict) -> np.ndarray:
        """
        Compute column-integrated sensible heat flux from fire state.

        This estimates the vertical heat release from fire at each grid point.
        Uses fire line intensity and flame length as proxies for heat release height.

        Parameters:
            state: Dictionary from fire.get_state() with at least:
                   - 'intensity': (ny, nx) array [kW/m]
                   - 'flame_length': (ny, nx) array [m]

        Returns:
            heat_flux: (ny, nx) array [W/m²] - column-integrated sensible heat flux
                       Values in burned/burning areas, 0 elsewhere
        """
        intensity = state.get('intensity', np.zeros((self.ny, self.nx)))  # [kW/m]
        flame_length = state.get('flame_length', np.zeros((self.ny, self.nx)))  # [m]
        phi = state.get('phi', np.ones((self.ny, self.nx)))  # Level set

        # Mark actively burning cells (phi near or below 0)
        # Use a narrow band: -2*max(dx,dy) < phi < 0
        active_band = max(self.dx, self.dy) * 2.0
        burning = (phi > -active_band) & (phi <= 0.0)

        # Convert intensity [kW/m] to heat release rate per unit area [W/m²]
        # Assumption: flame occupies a vertical column with height ~ flame_length
        # Heat release per unit area ≈ intensity / characteristic_width
        # Characteristic width typically 1-2 m for grassfire flame front
        characteristic_flame_width = 1.5  # meters

        # Sensible heat flux [W/m²] = intensity [kW/m] / width [m] * 1000
        intensity_W_per_m = intensity * 1000.0  # Convert kW to W
        heat_flux = np.zeros_like(intensity)

        # Only in burning area
        heat_flux[burning] = intensity_W_per_m[burning] / characteristic_flame_width

        return heat_flux

    def compute_buoyancy_parameter(self, state: Dict, T_ambient: float = 300.0,
                                   T_flame: float = 1200.0) -> np.ndarray:
        """
        Compute buoyancy parameter (Δθ/θ or density anomaly) from heat flux.

        Used by wind solver to compute updraft velocity and buoyancy forcing.

        Parameters:
            state: Fire state dictionary
            T_ambient: Ambient air temperature [K], default 300 K (27°C)
            T_flame: Typical flame temperature [K], default 1200 K

        Returns:
            buoyancy: (ny, nx) array [K] - potential temperature anomaly
        """
        heat_flux = self.compute_column_heat_flux(state)

        # Convert heat flux to temperature anomaly
        # Q = ρ * Cp * w' * Δθ
        # Δθ ≈ Q / (ρ * Cp * w_scale)
        # where w_scale is a characteristic updraft velocity (~1-3 m/s)

        rho_air = 1.2  # Air density [kg/m³] at sea level
        Cp_air = 1005.0  # Specific heat of air [J/(kg·K)]
        w_scale = 2.0  # Characteristic updraft velocity [m/s]

        # Potential temperature anomaly [K]
        buoyancy = heat_flux / (rho_air * Cp_air * w_scale)

        return buoyancy

    def compute_fire_footprint(self, state: Dict) -> np.ndarray:
        """
        Compute binary footprint of active fire area.

        Used to apply heat sources only where fire is spreading.

        Parameters:
            state: Fire state dictionary

        Returns:
            footprint: (ny, nx) binary array, 1.0 = active fire, 0.0 = no fire
        """
        phi = state.get('phi', np.ones((self.ny, self.nx)))
        intensity = state.get('intensity', np.zeros((self.ny, self.nx)))

        # Active fire: intensity > 0 and near the fire front
        active_band = max(self.dx, self.dy) * 2.0
        burning = (phi > -active_band) & (phi <= 0.0) & (intensity > 10.0)  # 10 kW/m threshold

        footprint = burning.astype(float)
        return footprint

    def compute_total_heat_release(self, state: Dict) -> float:
        """
        Compute total heat release rate across domain.

        Useful for monitoring and diagnostic purposes.

        Parameters:
            state: Fire state dictionary

        Returns:
            total_heat: Total heat release [W] (or [MW])
        """
        heat_flux = self.compute_column_heat_flux(state)
        total_heat_W = np.sum(heat_flux) * self.cell_area
        total_heat_MW = total_heat_W / 1.0e6

        return total_heat_MW

    def get_flux_dict(self, state: Dict) -> Dict[str, np.ndarray]:
        """
        Compute all heat flux terms and return as dictionary.

        Convenient interface for passing to wind solver.

        Parameters:
            state: Fire state dictionary

        Returns:
            flux_dict: Dictionary with keys:
                - 'heat_flux': (ny, nx) sensible heat flux [W/m²]
                - 'buoyancy': (ny, nx) potential temperature anomaly [K]
                - 'fire_footprint': (ny, nx) binary active fire mask
                - 'total_heat_MW': scalar total heat release [MW]
                - 'intensity': (ny, nx) original fire line intensity [kW/m]
                - 'flame_length': (ny, nx) original flame length [m]
        """
        return {
            'heat_flux': self.compute_column_heat_flux(state),
            'buoyancy': self.compute_buoyancy_parameter(state),
            'fire_footprint': self.compute_fire_footprint(state),
            'total_heat_MW': self.compute_total_heat_release(state),
            'intensity': state.get('intensity', np.zeros((self.ny, self.nx))),
            'flame_length': state.get('flame_length', np.zeros((self.ny, self.nx))),
        }

    def apply_heat_scaling(self, heat_flux: np.ndarray, scale_factor: float = 1.0) -> np.ndarray:
        """
        Apply scaling factor to heat fluxes (e.g., for sensitivity studies).

        Parameters:
            heat_flux: (ny, nx) heat flux array [W/m²]
            scale_factor: Multiplicative scaling factor

        Returns:
            scaled_heat_flux: (ny, nx) scaled array
        """
        return heat_flux * scale_factor

    def smooth_heat_field(self, heat_flux: np.ndarray, kernel_size: int = 3) -> np.ndarray:
        """
        Apply spatial smoothing to heat flux field (optional).

        Can reduce numerical noise when passing to wind solver.

        Parameters:
            heat_flux: (ny, nx) array
            kernel_size: Size of smoothing kernel (3, 5, 7, etc.)

        Returns:
            smoothed_flux: (ny, nx) smoothed array
        """
        from scipy.ndimage import uniform_filter

        # Only smooth in areas with heat flux
        mask = heat_flux > 0.0

        smoothed = np.copy(heat_flux)
        smoothed[mask] = uniform_filter(heat_flux[mask], size=kernel_size)

        return smoothed


# Convenience functions for standalone use

def extract_heat_flux(fire_state: Dict, grid_info: Dict) -> Dict[str, np.ndarray]:
    """
    Standalone function to extract all heat flux terms from fire state.

    Parameters:
        fire_state: Dictionary from fire.get_state()
        grid_info: Dictionary with keys: 'nx', 'ny', 'dx', 'dy'

    Returns:
        Dictionary with all heat flux fields
    """
    ny, nx = grid_info['ny'], grid_info['nx']
    dx, dy = grid_info['dx'], grid_info['dy']

    extractor = SurfaceFluxExtractor((ny, nx), dx, dy)
    return extractor.get_flux_dict(fire_state)


def compute_heat_weighted_wind_forcing(heat_flux: np.ndarray,
                                       u_wind: np.ndarray,
                                       v_wind: np.ndarray,
                                       max_flux: float = 1e5) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute wind field modifications due to heat flux (simple model).

    This is a simplified wind response to heat: updrafts reduce horizontal wind speed
    in high-heat regions.

    Parameters:
        heat_flux: (ny, nx) heat flux [W/m²]
        u_wind: (ny, nx) base u-component wind [m/s]
        v_wind: (ny, nx) base v-component wind [m/s]
        max_flux: Normalization reference [W/m²]

    Returns:
        (u_modified, v_modified): Wind components with heat effects applied
    """
    # Normalize heat flux to [0, 1]
    heat_normalized = np.clip(heat_flux / max_flux, 0.0, 1.0)

    # Reduce horizontal wind by ~20% where heat is high
    # (physical interpretation: strong vertical updraft reduces horizontal component)
    wind_reduction = 0.2 * heat_normalized

    u_modified = u_wind * (1.0 - wind_reduction)
    v_modified = v_wind * (1.0 - wind_reduction)

    return u_modified, v_modified


if __name__ == "__main__":
    # Example usage
    import sys

    print("Surface Flux Extractor - Unit Test")
    print("=" * 70)

    # Create mock fire state
    ny, nx = 32, 32
    dx, dy = 30.0, 30.0

    mock_state = {
        'intensity': np.random.exponential(100.0, (ny, nx)) * (np.random.rand(ny, nx) > 0.7),
        'flame_length': np.ones((ny, nx)) * 5.0,
        'phi': np.random.randn(ny, nx) * 50.0,
    }

    # Extract fluxes
    extractor = SurfaceFluxExtractor((ny, nx), dx, dy)
    fluxes = extractor.get_flux_dict(mock_state)

    print(f"\nFire grid: {ny} × {nx}")
    print(f"Grid spacing: dx={dx} m, dy={dy} m")
    print(f"\nHeat flux statistics:")
    print(f"  Mean: {fluxes['heat_flux'].mean():.1f} W/m²")
    print(f"  Max:  {fluxes['heat_flux'].max():.1f} W/m²")
    print(f"  Total: {fluxes['total_heat_MW']:.2f} MW")
    print(f"\nBuoyancy statistics:")
    print(f"  Mean: {fluxes['buoyancy'].mean():.3f} K")
    print(f"  Max:  {fluxes['buoyancy'].max():.3f} K")
    print(f"\nFire footprint:")
    print(f"  Active cells: {int(fluxes['fire_footprint'].sum())} / {ny*nx}")
    print(f"  Active fraction: {fluxes['fire_footprint'].sum() / (ny*nx) * 100:.1f}%")

    print("\n✓ Surface flux extraction successful")
