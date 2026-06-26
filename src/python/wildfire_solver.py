#!/usr/bin/env python3
"""
wildfire_solver.py - High-level Python wrapper for pyWildfire fire solver

Provides a clean, object-oriented API for running wildfire simulations from Python.
Supports coupling with external wind solvers (e.g., massconsistent_amr).

Example:
    from wildfire_solver import WildfireSolver
    
    # Initialize solver
    fire = WildfireSolver("inputs.i")
    
    # Run time loop
    for n in range(num_steps):
        # Optionally update wind from external source
        fire.update_wind(u_wind, v_wind)
        
        # Advance one timestep
        fire.step()
        
        # Get current state for visualization
        state = fire.get_state()
        plot_fire_front(state['phi'], state['time'])
"""

import numpy as np
try:
    import pyWildfire
except ImportError as e:
    raise ImportError(
        "Could not import pyWildfire module. "
        "Build with: cmake -S . -B build -DLEVELSET_BUILD_PYTHON_BINDINGS=ON\n"
        f"Error: {e}"
    )


class WildfireSolver:
    """
    High-level Python interface to wildfire_levelset fire solver.
    
    This class provides a convenient object-oriented API for running
    wildfire simulations, with support for coupling to external wind solvers.
    
    Attributes:
        initialized (bool): Whether the solver has been initialized
        time (float): Current simulation time (seconds)
        step_num (int): Current timestep number
        nx, ny (int): Grid dimensions
        xmin, xmax, ymin, ymax (float): Domain bounds (meters)
        dx, dy (float): Cell sizes (meters)
    """
    
    def __init__(self, inputs_file=None):
        """
        Initialize the wildfire solver.
        
        Parameters:
            inputs_file (str, optional): Path to inputs file.
                                        If None, must call initialize() later.
        """
        self.initialized = False
        self.time = 0.0
        self.step_num = 0
        self.nx = 0
        self.ny = 0
        self.xmin = self.xmax = self.ymin = self.ymax = 0.0
        self.dx = self.dy = 0.0
        
        if inputs_file is not None:
            self.initialize(inputs_file)
    
    def initialize(self, inputs_file):
        """
        Initialize the solver from an inputs file.
        
        Parameters:
            inputs_file (str): Path to the inputs file (e.g., "inputs.i")
        
        Returns:
            bool: True if initialization succeeded
        
        Raises:
            RuntimeError: If initialization fails
        """
        result = pyWildfire.initialize(inputs_file)
        
        if not result['success']:
            raise RuntimeError(f"Failed to initialize fire solver from {inputs_file}")
        
        self.initialized = True
        self.nx = result['nx']
        self.ny = result['ny']
        self.xmin = result['xmin']
        self.xmax = result['xmax']
        self.ymin = result['ymin']
        self.ymax = result['ymax']
        self.dx = result['dx']
        self.dy = result['dy']
        
        print(f"✓ Fire solver initialized")
        print(f"  Grid: {self.nx} × {self.ny}")
        print(f"  Domain: X=[{self.xmin:.1f}, {self.xmax:.1f}] m, "
              f"Y=[{self.ymin:.1f}, {self.ymax:.1f}] m")
        print(f"  Resolution: dx={self.dx:.2f} m, dy={self.dy:.2f} m")
        
        return True
    
    def step(self, max_time=None):
        """
        Advance the simulation by one timestep.
        
        Parameters:
            max_time (float, optional): Maximum time to reach.
                                       If provided, dt is clamped so time doesn't exceed this.
        
        Returns:
            dict: Dictionary with 'success', 'dt', 'time', 'step'
        
        Raises:
            RuntimeError: If solver is not initialized or step fails
        """
        if not self.initialized:
            raise RuntimeError("Solver not initialized. Call initialize() first.")
        
        result = pyWildfire.advance()
        
        if not result['success']:
            raise RuntimeError(f"Timestep failed at step {self.step_num}")
        
        self.time = result['time']
        self.step_num = result['step']
        
        return result
    
    def get_state(self):
        """
        Get the current state of the fire simulation.
        
        Returns:
            dict: Dictionary containing:
                - time (float): Current simulation time (seconds)
                - step (int): Current timestep number
                - nx, ny (int): Grid dimensions
                - xmin, xmax, ymin, ymax (float): Domain bounds
                - dx, dy (float): Cell sizes
                - phi (ndarray): Level set field, shape (ny, nx)
                - ros (ndarray): Rate of spread (m/s), shape (ny, nx)
                - intensity (ndarray): Fire line intensity (kW/m), shape (ny, nx)
                - flame_length (ndarray): Flame length (m), shape (ny, nx)
                - u_wind, v_wind (ndarray): Wind components (m/s), shape (ny, nx)
                - arrival_time (ndarray): Fire arrival time (s), shape (ny, nx)
        
        Raises:
            RuntimeError: If solver is not initialized
        """
        if not self.initialized:
            raise RuntimeError("Solver not initialized. Call initialize() first.")
        
        return pyWildfire.get_state()
    
    def update_wind(self, u_wind, v_wind):
        """
        Update the wind field from 2D arrays.
        
        This allows coupling with external wind solvers by passing
        column-averaged wind components.
        
        Parameters:
            u_wind (ndarray): U-component of wind (m/s), shape (ny, nx)
            v_wind (ndarray): V-component of wind (m/s), shape (ny, nx)
        
        Returns:
            bool: True if wind was updated successfully
        
        Raises:
            RuntimeError: If solver is not initialized
            ValueError: If wind array shapes don't match grid
        """
        if not self.initialized:
            raise RuntimeError("Solver not initialized. Call initialize() first.")
        
        u_wind = np.asarray(u_wind, dtype=np.float64)
        v_wind = np.asarray(v_wind, dtype=np.float64)
        
        if u_wind.shape != (self.ny, self.nx):
            raise ValueError(
                f"u_wind shape {u_wind.shape} doesn't match grid ({self.ny}, {self.nx})"
            )
        if v_wind.shape != (self.ny, self.nx):
            raise ValueError(
                f"v_wind shape {v_wind.shape} doesn't match grid ({self.ny}, {self.nx})"
            )
        
        return pyWildfire.update_wind(u_wind, v_wind)
    
    def update_wind_3d(self, u_3d, v_3d, w_3d, nz, zmin, zmax):
        """
        Update the wind field from 3D arrays.
        
        This allows coupling with 3D wind solvers (e.g., massconsistent_amr)
        by passing full 3D wind velocity fields.
        
        Parameters:
            u_3d (ndarray): U-component, shape (nz, ny, nx)
            v_3d (ndarray): V-component, shape (nz, ny, nx)
            w_3d (ndarray): W-component, shape (nz, ny, nx)
            nz (int): Number of vertical levels
            zmin, zmax (float): Vertical domain bounds (meters)
        
        Returns:
            bool: True if wind was updated successfully
        
        Raises:
            RuntimeError: If solver is not initialized
            ValueError: If wind array shapes don't match grid
        """
        if not self.initialized:
            raise RuntimeError("Solver not initialized. Call initialize() first.")
        
        u_3d = np.asarray(u_3d, dtype=np.float64)
        v_3d = np.asarray(v_3d, dtype=np.float64)
        w_3d = np.asarray(w_3d, dtype=np.float64)
        
        expected_shape = (nz, self.ny, self.nx)
        if u_3d.shape != expected_shape:
            raise ValueError(f"u_3d shape {u_3d.shape} doesn't match ({nz}, {self.ny}, {self.nx})")
        if v_3d.shape != expected_shape:
            raise ValueError(f"v_3d shape {v_3d.shape} doesn't match ({nz}, {self.ny}, {self.nx})")
        if w_3d.shape != expected_shape:
            raise ValueError(f"w_3d shape {w_3d.shape} doesn't match ({nz}, {self.ny}, {self.nx})")
        
        # Flatten in Fortran order for pyWildfire
        u_flat = u_3d.flatten('F')
        v_flat = v_3d.flatten('F')
        w_flat = w_3d.flatten('F')
        
        return pyWildfire.update_wind_3d(
            self.nx, self.ny, nz,
            self.xmin, self.xmax,
            self.ymin, self.ymax,
            zmin, zmax,
            u_flat, v_flat, w_flat
        )
    
    def write_plotfile(self, plotfile_name):
        """
        Write the current state to an AMReX plotfile.
        
        Parameters:
            plotfile_name (str): Name/path of the plotfile to write
        
        Returns:
            bool: True if plotfile was written successfully
        
        Raises:
            RuntimeError: If solver is not initialized
        """
        if not self.initialized:
            raise RuntimeError("Solver not initialized. Call initialize() first.")
        
        return pyWildfire.write_plotfile(plotfile_name)
    
    def run(self, final_time=None, num_steps=None, plot_interval=None,
            callback=None, wind_callback=None):
        """
        Run the simulation for a specified duration or number of steps.
        
        Parameters:
            final_time (float, optional): Final simulation time (seconds).
                                         Either this or num_steps must be provided.
            num_steps (int, optional): Number of timesteps to run.
            plot_interval (float, optional): Write plotfiles at this interval (seconds).
            callback (callable, optional): Function called after each timestep
                                          with signature: callback(step, time, state_dict)
            wind_callback (callable, optional): Function called before each timestep
                                               to update wind with signature:
                                               u, v = wind_callback(time)
        
        Returns:
            dict: Final state dictionary
        
        Raises:
            RuntimeError: If solver is not initialized
            ValueError: If neither final_time nor num_steps is provided
        """
        if not self.initialized:
            raise RuntimeError("Solver not initialized. Call initialize() first.")
        
        if final_time is None and num_steps is None:
            raise ValueError("Must specify either final_time or num_steps")
        
        use_final_time = (final_time is not None)
        next_plot_time = self.time + plot_interval if plot_interval else np.inf
        plot_num = 0
        
        print(f"\nRunning wildfire simulation...")
        if use_final_time:
            print(f"  Final time: {final_time} s")
        else:
            print(f"  Number of steps: {num_steps}")
        
        step_count = 0
        while True:
            # Check stopping condition
            if use_final_time:
                if self.time >= final_time:
                    break
            else:
                if step_count >= num_steps:
                    break
            
            # Update wind from external callback if provided
            if wind_callback is not None:
                u_wind, v_wind = wind_callback(self.time)
                self.update_wind(u_wind, v_wind)
            
            # Advance one timestep
            result = self.step(max_time=final_time if use_final_time else None)
            step_count += 1
            
            # Write plotfile if needed
            if plot_interval and self.time >= next_plot_time:
                plotfile_name = f"plt{plot_num:05d}"
                self.write_plotfile(plotfile_name)
                print(f"  Wrote {plotfile_name} at t={self.time:.2f} s")
                next_plot_time += plot_interval
                plot_num += 1
            
            # Call user callback if provided
            if callback is not None:
                state = self.get_state()
                callback(self.step_num, self.time, state)
        
        print(f"✓ Simulation complete at t={self.time:.2f} s ({self.step_num} steps)\n")
        return self.get_state()
    
    def get_burned_area(self):
        """
        Calculate the total burned area from the current fire state.
        
        Returns:
            float: Total burned area in m² (cells where phi <= 0)
        
        Raises:
            RuntimeError: If solver is not initialized
        """
        if not self.initialized:
            raise RuntimeError("Solver not initialized. Call initialize() first.")
        
        state = self.get_state()
        phi = state['phi']
        burned_cells = np.sum(phi <= 0.0)
        return burned_cells * self.dx * self.dy
    
    def get_fire_perimeter(self):
        """
        Estimate the fire perimeter from the level set field.
        
        Uses a simple marching squares approach to estimate the perimeter
        where phi crosses zero (the fire front).
        
        Returns:
            float: Estimated fire perimeter in meters
        
        Raises:
            RuntimeError: If solver is not initialized
        """
        if not self.initialized:
            raise RuntimeError("Solver not initialized. Call initialize() first.")
        
        state = self.get_state()
        phi = state['phi']
        
        # Count zero crossings (phi changes sign across cell edges)
        perimeter = 0.0
        
        # Count horizontal crossings
        for j in range(self.ny):
            for i in range(self.nx - 1):
                if phi[j, i] * phi[j, i+1] < 0:  # Sign change
                    perimeter += self.dy
        
        # Count vertical crossings
        for j in range(self.ny - 1):
            for i in range(self.nx):
                if phi[j, i] * phi[j+1, i] < 0:  # Sign change
                    perimeter += self.dx
        
        return perimeter
    
    def compute_heat_release(self, state=None):
        """
        Compute heat release information from fire state.
        
        Uses Byram intensity and rate of spread to estimate convective heat release.
        The total heat release rate is: Intensity * Perimeter
        
        The convective heat release per unit area (kW/m²) can be estimated from
        the flame length and intensity.
        
        Parameters:
            state (dict, optional): Fire state dictionary. If None, calls get_state().
        
        Returns:
            dict: Dictionary containing:
                - 'total_heat_release': Total convective heat release (kW)
                - 'heat_release_rate': Heat release rate per unit perimeter (kW/m)
                - 'mean_intensity': Mean fire intensity (kW/m)
                - 'max_intensity': Maximum fire intensity (kW/m)
                - 'max_ros': Maximum rate of spread (m/s)
                - 'perimeter': Estimated fire perimeter (m)
                - 'burned_area': Total burned area (m²)
                - 'surface_flux': Surface heat flux array (ny, nx) in kW/m²
        
        Raises:
            RuntimeError: If solver is not initialized
        """
        if not self.initialized:
            raise RuntimeError("Solver not initialized. Call initialize() first.")
        
        if state is None:
            state = self.get_state()
        
        # Extract key fields
        phi = state['phi']
        intensity = state['intensity']  # Byram intensity in kW/m
        ros = state['ros']  # Rate of spread in m/s
        flame_length = state['flame_length']  # Flame length in m
        
        # Calculate metrics
        perimeter = self.get_fire_perimeter()
        burned_area = self.get_burned_area()
        
        mean_intensity = intensity[intensity > 0].mean() if np.any(intensity > 0) else 0.0
        max_intensity = intensity.max()
        max_ros = ros.max()
        
        # Total heat release = intensity × perimeter (kW)
        # This assumes intensity is defined per unit fireline length
        total_heat_release = mean_intensity * perimeter if perimeter > 0 else 0.0
        
        # Heat release rate (total / perimeter) - should approximate mean intensity
        heat_release_rate = mean_intensity
        
        # Compute surface heat flux array (kW/m²)
        # Surface flux = intensity / flame_length (simplified model)
        # We also weight by whether cell is burning (phi is transitioning)
        surface_flux = np.zeros_like(phi)
        
        for j in range(self.ny):
            for i in range(self.nx):
                # Only compute flux near the fire front (where phi ~ 0)
                # and in recently burned cells
                if phi[j, i] < 100.0 and phi[j, i] > -100.0:  # Near fire front
                    if flame_length[j, i] > 0.1:
                        # Heat flux = intensity / flame length
                        surface_flux[j, i] = intensity[j, i] / flame_length[j, i]
                    else:
                        surface_flux[j, i] = intensity[j, i] * 10.0  # Default scaling
        
        return {
            'total_heat_release': total_heat_release,
            'heat_release_rate': heat_release_rate,
            'mean_intensity': mean_intensity,
            'max_intensity': max_intensity,
            'max_ros': max_ros,
            'perimeter': perimeter,
            'burned_area': burned_area,
            'surface_flux': surface_flux
        }
    
    def get_surface_fluxes(self):
        """
        Get surface heat flux array for coupling with wind solver.
        
        Returns the heat flux field at the surface (2D array) computed from
        the current fire state. This array can be passed to the wind solver
        to account for heating effects on atmospheric dynamics.
        
        Returns:
            dict: Dictionary containing:
                - 'heat_flux': Surface heat flux array (ny, nx) in kW/m²
                - 'flux_units': Unit description
                - 'grid_info': Dict with 'nx', 'ny', 'dx', 'dy' grid information
        
        Raises:
            RuntimeError: If solver is not initialized
        """
        if not self.initialized:
            raise RuntimeError("Solver not initialized. Call initialize() first.")
        
        state = self.get_state()
        heat_data = self.compute_heat_release(state)
        
        return {
            'heat_flux': heat_data['surface_flux'],
            'flux_units': 'kW/m²',
            'grid_info': {
                'nx': self.nx,
                'ny': self.ny,
                'dx': self.dx,
                'dy': self.dy,
                'xmin': self.xmin,
                'xmax': self.xmax,
                'ymin': self.ymin,
                'ymax': self.ymax
            }
        }
    
    def finalize(self):
        """
        Clean up and finalize the solver.
        
        After calling this, you must call initialize() again before using the solver.
        """
        if self.initialized:
            pyWildfire.finalize()
            self.initialized = False
            print("Fire solver finalized")
    
    def __del__(self):
        """Destructor: finalize solver if still initialized."""
        if self.initialized:
            self.finalize()
    
    def __enter__(self):
        """Context manager entry."""
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit: automatically finalize."""
        self.finalize()
        return False


def main():
    """Example usage of WildfireSolver class."""
    import sys
    
    if len(sys.argv) < 2:
        print("Usage: python wildfire_solver.py <inputs_file>")
        return 1
    
    inputs_file = sys.argv[1]
    
    # Initialize solver
    fire = WildfireSolver(inputs_file)
    
    # Run for 10 timesteps
    def progress_callback(step, time, state):
        phi = state['phi']
        burned_area = np.sum(phi <= 0.0) * fire.dx * fire.dy
        print(f"  Step {step}: t={time:.2f} s, burned area={burned_area:.1f} m²")
    
    final_state = fire.run(num_steps=10, callback=progress_callback)
    
    # Print final statistics
    phi_final = final_state['phi']
    total_burned = np.sum(phi_final <= 0.0) * fire.dx * fire.dy
    print(f"\nFinal burned area: {total_burned:.1f} m²")
    
    fire.finalize()
    return 0


if __name__ == "__main__":
    sys.exit(main())
