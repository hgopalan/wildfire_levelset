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
        
        # Performance configuration (placeholders for future API expansion)
        self._cfl = 0.8
        self._max_time = np.inf
        self._nsteps = None
        
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
    
    def get_status(self):
        """
        Get comprehensive solver status information.
        
        Returns:
            dict: Dictionary containing:
                - 'initialized': Whether solver is initialized
                - 'time': Current simulation time (seconds)
                - 'step': Current step number
                - 'nx', 'ny': Grid dimensions
                - 'xmin', 'xmax', 'ymin', 'ymax': Domain bounds (meters)
                - 'dx', 'dy': Cell spacing (meters)
                - 'domain_area': Total domain area (m²)
        
        Raises:
            RuntimeError: If solver is not initialized
        """
        if not self.initialized:
            raise RuntimeError("Solver not initialized. Call initialize() first.")
        
        domain_area = (self.xmax - self.xmin) * (self.ymax - self.ymin)
        
        return {
            'initialized': self.initialized,
            'time': self.time,
            'step': self.step_num,
            'nx': self.nx,
            'ny': self.ny,
            'xmin': self.xmin,
            'xmax': self.xmax,
            'ymin': self.ymin,
            'ymax': self.ymax,
            'dx': self.dx,
            'dy': self.dy,
            'domain_area': domain_area
        }
    
    def get_domain_bounds(self):
        """
        Get domain boundary information.
        
        Returns:
            dict: Dictionary with 'xmin', 'xmax', 'ymin', 'ymax', 'zmin', 'zmax'
                  Note: zmin=0, zmax=0 for 2D fire solver
        
        Raises:
            RuntimeError: If solver is not initialized
        """
        if not self.initialized:
            raise RuntimeError("Solver not initialized. Call initialize() first.")
        
        return {
            'xmin': self.xmin,
            'xmax': self.xmax,
            'ymin': self.ymin,
            'ymax': self.ymax,
            'zmin': 0.0,  # Fire solver is 2D
            'zmax': 0.0
        }
    
    def get_grid_spacing(self):
        """
        Get grid cell spacing.
        
        Returns:
            dict: Dictionary with 'dx', 'dy', 'dz' (dz=0 for 2D)
        
        Raises:
            RuntimeError: If solver is not initialized
        """
        if not self.initialized:
            raise RuntimeError("Solver not initialized. Call initialize() first.")
        
        return {
            'dx': self.dx,
            'dy': self.dy,
            'dz': 0.0  # Fire solver is 2D
        }
    
    def get_grid_dimensions(self):
        """
        Get grid dimensions.
        
        Returns:
            dict: Dictionary with 'nx', 'ny', 'nz' (nz=1 for 2D)
        
        Raises:
            RuntimeError: If solver is not initialized
        """
        if not self.initialized:
            raise RuntimeError("Solver not initialized. Call initialize() first.")
        
        return {
            'nx': self.nx,
            'ny': self.ny,
            'nz': 1  # Fire solver is 2D (single surface)
        }
    
    def get_field(self, field_name):
        """
        Get a specific fire field by name.
        
        Parameters:
            field_name (str): Name of field to retrieve. Options:
                'phi' - Level set field (signed distance)
                'ros' - Rate of spread (m/s)
                'intensity' - Fireline intensity (kW/m)
                'flame_length' - Flame length (m)
                'u_wind' - U wind component (m/s)
                'v_wind' - V wind component (m/s)
                'arrival_time' - Fire arrival time (seconds)
        
        Returns:
            ndarray: 2D array (ny, nx) containing the requested field
        
        Raises:
            RuntimeError: If solver is not initialized
            ValueError: If field_name is not recognized
        """
        if not self.initialized:
            raise RuntimeError("Solver not initialized. Call initialize() first.")
        
        state = self.get_state()
        
        if field_name not in state:
            raise ValueError(
                f"Unknown field '{field_name}'. Available fields: "
                f"{', '.join(state.keys())}"
            )
        
        return state[field_name]
    
    def get_statistics(self):
        """
        Get comprehensive fire simulation statistics.
        
        Returns:
            dict: Dictionary containing:
                - 'time': Current time (seconds)
                - 'burned_area': Total burned area (m²)
                - 'perimeter': Estimated fire perimeter (m)
                - 'max_ros': Maximum rate of spread (m/s)
                - 'mean_ros': Mean rate of spread where fire is active (m/s)
                - 'max_intensity': Maximum fireline intensity (kW/m)
                - 'mean_intensity': Mean fireline intensity (kW/m)
                - 'total_heat_release': Total convective heat released (kW)
                - 'num_burning_cells': Number of cells with active fire
        
        Raises:
            RuntimeError: If solver is not initialized
        """
        if not self.initialized:
            raise RuntimeError("Solver not initialized. Call initialize() first.")
        
        state = self.get_state()
        phi = state['phi']
        ros = state['ros']
        intensity = state['intensity']
        
        # Calculate statistics
        burned_area = self.get_burned_area()
        perimeter = self.get_fire_perimeter()
        heat_data = self.compute_heat_release(state)
        
        # Count burning cells (where ROS > 0)
        burning_mask = ros > 0.0
        num_burning = np.sum(burning_mask)
        
        # Calculate mean values for active fire cells
        mean_ros = ros[burning_mask].mean() if num_burning > 0 else 0.0
        mean_intensity = intensity[burning_mask].mean() if num_burning > 0 else 0.0
        
        return {
            'time': self.time,
            'burned_area': burned_area,
            'perimeter': perimeter,
            'max_ros': heat_data['max_ros'],
            'mean_ros': mean_ros,
            'max_intensity': heat_data['max_intensity'],
            'mean_intensity': mean_intensity,
            'total_heat_release': heat_data['total_heat_release'],
            'num_burning_cells': int(num_burning)
        }
    
    def get_diagnostic_info(self):
        """
        Get diagnostic information about solver state.
        
        Returns:
            dict: Dictionary containing diagnostic data:
                - 'initialized': Solver ready status
                - 'time': Current simulation time
                - 'step': Current step number
                - 'domain_bounds': Domain extent
                - 'grid_spacing': Cell dimensions
                - 'solver_type': Description of solver configuration
        
        Raises:
            RuntimeError: If solver is not initialized
        """
        if not self.initialized:
            raise RuntimeError("Solver not initialized. Call initialize() first.")
        
        return {
            'initialized': self.initialized,
            'time': self.time,
            'step': self.step_num,
            'domain_bounds': self.get_domain_bounds(),
            'grid_spacing': self.get_grid_spacing(),
            'solver_type': 'WildfireSolver (pyWildfire wrapper)',
            'api_version': '2.0'
        }
    
    def set_fuel_model(self, model_number):
        """
        Set fuel model for Rothermel ROS calculations.
        
        Note: Currently requires setting via inputs file. This is a placeholder
        for future direct API support.
        
        Parameters:
            model_number (int): Fuel model number (1-13 for Anderson, 1-40 for Scott-Burgan)
        
        Returns:
            bool: True if fuel model was set (or would be set in future implementation)
        
        Raises:
            ValueError: If model_number is invalid
        """
        if model_number < 1 or model_number > 100:
            raise ValueError(f"Fuel model number must be 1-100, got {model_number}")
        
        # TODO: When fire_solver_api.H exposes fuel model control
        # return pyWildfire.set_fuel_model(model_number)
        
        # For now, provide documentation
        import warnings
        warnings.warn(
            f"set_fuel_model({model_number}) requires rebuilding with fuel model API. "
            "Set fuel model in inputs file instead: fire.fuel_model = {model_number}",
            FutureWarning
        )
        return False
    
    def get_fuel_properties(self):
        """
        Get current fuel model properties for Rothermel calculations.
        
        Returns:
            dict: Dictionary with fuel properties:
                - 'model_number': Fuel model ID
                - 'fuel_load': Oven-dry fuel load (tons/acre)
                - 'fuel_bed_depth': Fuel bed depth (feet)
                - 'moisture_of_extinction': Moisture of extinction (%)
                - 'low_heat_content': Low heat content (BTU/lb)
                - 'description': Fuel model description
        
        Note: Currently returns placeholder values. Requires C++ API expansion.
        
        Raises:
            RuntimeError: If solver is not initialized
        """
        if not self.initialized:
            raise RuntimeError("Solver not initialized. Call initialize() first.")
        
        # TODO: When fire_solver_api.H exposes fuel properties
        # return pyWildfire.get_fuel_properties()
        
        # Placeholder implementation
        return {
            'model_number': 0,  # Unknown
            'fuel_load': 0.0,
            'fuel_bed_depth': 0.0,
            'moisture_of_extinction': 0.0,
            'low_heat_content': 0.0,
            'description': 'Fuel properties currently unavailable through Python API',
            'note': 'Set fuel parameters in inputs file'
        }
    
    def compute_ros(self):
        """
        Compute Rate of Spread field using Rothermel model.
        
        Returns:
            ndarray: ROS field (ny, nx) in m/s
        
        Note: ROS is computed internally during fire.step(). This method
        simply extracts the current ROS field from the fire state.
        
        Raises:
            RuntimeError: If solver is not initialized
        """
        if not self.initialized:
            raise RuntimeError("Solver not initialized. Call initialize() first.")
        
        state = self.get_state()
        return state['ros']
    
    def set_fuel_moisture(self, dead_1hr, dead_10hr, dead_100hr, live_herbaceous):
        """
        Set fuel moisture content for Rothermel ROS calculations.
        
        Note: Currently requires setting via inputs file. This is a placeholder
        for future direct API support.
        
        Parameters:
            dead_1hr (float): 1-hour fuel moisture (%)
            dead_10hr (float): 10-hour fuel moisture (%)
            dead_100hr (float): 100-hour fuel moisture (%)
            live_herbaceous (float): Live herbaceous fuel moisture (%)
        
        Returns:
            bool: True if moisture was set (or would be set in future)
        
        Raises:
            ValueError: If any value is outside valid range (0-300%)
        """
        for val in [dead_1hr, dead_10hr, dead_100hr, live_herbaceous]:
            if val < 0 or val > 300:
                raise ValueError(f"Fuel moisture must be 0-300%, got {val}")
        
        # TODO: When fire_solver_api.H exposes moisture control
        # return pyWildfire.set_fuel_moisture(dead_1hr, dead_10hr, dead_100hr, live_herbaceous)
        
        import warnings
        warnings.warn(
            "set_fuel_moisture() requires API expansion in fire_solver_api.H. "
            "Set moisture in inputs file instead.",
            FutureWarning
        )
        return False
    
    def set_wind_direction(self, direction_degrees):
        """
        Set predominant wind direction for Rothermel calculations.
        
        Note: Wind direction affects elliptical spread pattern in Rothermel model.
        
        Parameters:
            direction_degrees (float): Wind direction in degrees (0=North, 90=East, etc.)
        
        Returns:
            bool: True if direction was set
        
        Raises:
            ValueError: If direction is outside 0-360 range
        """
        if direction_degrees < 0 or direction_degrees > 360:
            raise ValueError(f"Wind direction must be 0-360°, got {direction_degrees}")
        
        # TODO: When fire_solver_api.H exposes wind direction control
        import warnings
        warnings.warn(
            "set_wind_direction() requires API expansion. "
            "Set wind in inputs file or via update_wind() methods.",
            FutureWarning
        )
        return False
    
    def set_ambient_temperature(self, temp_celsius):
        """
        Set ambient air temperature for Rothermel calculations.
        
        Note: Temperature affects fuel moisture calculations.
        
        Parameters:
            temp_celsius (float): Temperature in °C
        
        Returns:
            bool: True if temperature was set
        """
        # TODO: When fire_solver_api.H exposes temperature control
        import warnings
        warnings.warn(
            "set_ambient_temperature() requires API expansion. "
            "Set temperature in inputs file.",
            FutureWarning
        )
        return False
    
    def set_relative_humidity(self, rh_percent):
        """
        Set relative humidity for Rothermel calculations.
        
        Note: RH affects fuel moisture and fire behavior.
        
        Parameters:
            rh_percent (float): Relative humidity (0-100%)
        
        Returns:
            bool: True if RH was set
        
        Raises:
            ValueError: If RH is outside 0-100 range
        """
        if rh_percent < 0 or rh_percent > 100:
            raise ValueError(f"Relative humidity must be 0-100%, got {rh_percent}")
        
        # TODO: When fire_solver_api.H exposes RH control
        import warnings
        warnings.warn(
            "set_relative_humidity() requires API expansion. "
            "Set RH in inputs file.",
            FutureWarning
        )
        return False
    
    def set_ignition(self, x, y, time=0.0, radius=1.0):
        """
        Set ignition point programmatically.
        
        Note: Currently requires configuring via inputs file.
        
        Parameters:
            x (float): X-coordinate of ignition point (meters)
            y (float): Y-coordinate of ignition point (meters)
            time (float, optional): Time of ignition (seconds), default 0.0
            radius (float, optional): Initial ignition radius (meters), default 1.0
        
        Returns:
            bool: True if ignition was set
        """
        # TODO: When fire_solver_api.H exposes ignition control
        import warnings
        warnings.warn(
            "set_ignition() requires API expansion. "
            "Configure ignition in inputs file.",
            FutureWarning
        )
        return False
    
    def get_ignition_state(self):
        """
        Get current ignition configuration.
        
        Returns:
            dict: Dictionary with ignition information:
                - 'configured': Whether ignition is configured
                - 'status': Current ignition status description
        
        Note: Returns placeholder until C++ API exposes ignition state.
        
        Raises:
            RuntimeError: If solver is not initialized
        """
        if not self.initialized:
            raise RuntimeError("Solver not initialized. Call initialize() first.")
        
        # TODO: When fire_solver_api.H exposes ignition state query
        return {
            'configured': True,  # Assumed from inputs file
            'status': 'Configured from inputs file',
            'note': 'Direct ignition state queries require fire_solver_api expansion'
        }
    
    @property
    def current_time(self):
        """Current simulation time (seconds)."""
        return self.time
    
    @property
    def timestep(self):
        """Current timestep number."""
        return self.step_num
    
    def set_cfl(self, cfl_value):
        """
        Set CFL (Courant-Friedrichs-Lewy) criterion for stability.
        
        Note: Currently a placeholder. Actual CFL control requires C++ API expansion.
        
        Parameters:
            cfl_value (float): CFL criterion (typically 0.1-0.9)
        
        Returns:
            bool: True if CFL was set
        
        Raises:
            ValueError: If CFL is outside valid range (0, 1]
        """
        if cfl_value <= 0 or cfl_value > 1.0:
            raise ValueError(f"CFL must be in (0, 1], got {cfl_value}")
        
        self._cfl = cfl_value
        
        # TODO: When fire_solver_api.H exposes CFL control
        # return pyWildfire.set_cfl(cfl_value)
        
        return True
    
    def set_max_time(self, t_max):
        """
        Set maximum simulation time for run() method.
        
        Parameters:
            t_max (float): Maximum simulation time (seconds)
        
        Returns:
            bool: True if max time was set
        
        Raises:
            ValueError: If t_max is negative
        """
        if t_max < 0:
            raise ValueError(f"Max time must be non-negative, got {t_max}")
        
        self._max_time = t_max
        return True
    
    def set_nsteps(self, nsteps):
        """
        Set maximum number of timesteps for run() method.
        
        Parameters:
            nsteps (int): Maximum number of steps
        
        Returns:
            bool: True if nsteps was set
        
        Raises:
            ValueError: If nsteps is not a positive integer
        """
        if nsteps <= 0 or not isinstance(nsteps, int):
            raise ValueError(f"nsteps must be positive integer, got {nsteps}")
        
        self._nsteps = nsteps
        return True
    
    def get_performance_metrics(self):
        """
        Get performance metrics for the solver.
        
        Returns:
            dict: Dictionary containing:
                - 'time_per_step': Average time per timestep (seconds)
                - 'total_time': Total wall-clock time (seconds)
                - 'steps_completed': Number of completed steps
                - 'memory_usage': Estimated memory usage (MB)
                - 'cfl_setting': Current CFL criterion
        
        Note: Returns placeholder values until C++ profiling infrastructure is added.
        
        Raises:
            RuntimeError: If solver is not initialized
        """
        if not self.initialized:
            raise RuntimeError("Solver not initialized. Call initialize() first.")
        
        return {
            'time_per_step': 0.0,  # TODO: Add timing instrumentation
            'total_time': 0.0,
            'steps_completed': self.step_num,
            'memory_usage': 0.0,  # TODO: Add memory tracking
            'cfl_setting': self._cfl,
            'note': 'Performance metrics require C++ instrumentation'
        }
    
    # ========================================================================
    # PHASE 1: Core Configuration & Properties Methods
    # ========================================================================
    
    def get_config(self):
        """
        Get comprehensive solver configuration.
        
        Returns:
            dict: Configuration with keys: nx, ny, xmin, xmax, ymin, ymax, dx, dy, time, step, dt
        """
        if not self.initialized:
            raise RuntimeError("Solver not initialized. Call initialize() first.")
        
        return pyWildfire.get_config()
    
    def get_vertical_domain(self):
        """
        Get vertical domain bounds.
        
        Returns:
            dict: Dict with 'zmin' and 'zmax' in meters
        """
        if not self.initialized:
            raise RuntimeError("Solver not initialized. Call initialize() first.")
        
        return pyWildfire.get_vertical_domain()
    
    def get_rothermel_properties(self):
        """
        Get Rothermel fuel model properties.
        
        Returns:
            dict: Fuel model properties for all size classes and moisture
        """
        if not self.initialized:
            raise RuntimeError("Solver not initialized. Call initialize() first.")
        
        return pyWildfire.get_rothermel_properties()
    
    def get_wind_ros_relationship(self):
        """
        Get wind-ROS interaction parameters.
        
        Returns:
            dict: B-coefficient and other Rothermel wind parameters
        """
        if not self.initialized:
            raise RuntimeError("Solver not initialized. Call initialize() first.")
        
        return pyWildfire.get_wind_ros_relationship()
    
    def get_spread_parameters(self):
        """
        Get spread model parameters.
        
        Returns:
            dict: Parameters for Richards/Balbi model if enabled
        """
        if not self.initialized:
            raise RuntimeError("Solver not initialized. Call initialize() first.")
        
        return pyWildfire.get_spread_parameters()
    
    def update_rothermel_fuel_load(self, dead_load, live_load):
        """
        Update Rothermel fuel loads.
        
        Parameters:
            dead_load (array-like): Dead fuel loads (tons/acre)
            live_load (array-like): Live fuel loads (tons/acre)
        
        Returns:
            bool: True if update succeeded
        """
        if not self.initialized:
            raise RuntimeError("Solver not initialized. Call initialize() first.")
        
        dead_load = np.asarray(dead_load, dtype=np.float64)
        live_load = np.asarray(live_load, dtype=np.float64)
        return pyWildfire.update_rothermel_fuel_load(dead_load, live_load)
    
    def validate_domain_compatibility(self, wind_nx, wind_ny, 
                                      wind_xmin, wind_xmax, 
                                      wind_ymin, wind_ymax):
        """
        Validate compatibility between fire and wind domains.
        
        Parameters:
            wind_nx, wind_ny (int): Wind grid dimensions
            wind_xmin, wind_xmax, wind_ymin, wind_ymax (float): Wind domain bounds
        
        Returns:
            dict: Compatibility flags
        """
        if not self.initialized:
            raise RuntimeError("Solver not initialized. Call initialize() first.")
        
        return pyWildfire.validate_domain_compatibility(
            wind_nx, wind_ny, wind_xmin, wind_xmax, wind_ymin, wind_ymax)
    
    # ========================================================================
    # PHASE 2: Terrain & Spatial Features Methods
    # ========================================================================
    
    def update_terrain(self, elevation, slope, aspect):
        """
        Update terrain data (elevation, slope, aspect).
        
        Parameters:
            elevation (ndarray): Elevation map (m), shape (ny, nx)
            slope (ndarray): Slope (degrees), shape (ny, nx)
            aspect (ndarray): Aspect (degrees), shape (ny, nx)
        
        Returns:
            bool: True if update succeeded
        """
        if not self.initialized:
            raise RuntimeError("Solver not initialized. Call initialize() first.")
        
        elevation = np.asarray(elevation, dtype=np.float64)
        slope = np.asarray(slope, dtype=np.float64)
        aspect = np.asarray(aspect, dtype=np.float64)
        
        if elevation.shape != (self.ny, self.nx):
            raise ValueError(f"elevation shape {elevation.shape} doesn't match ({self.ny}, {self.nx})")
        if slope.shape != (self.ny, self.nx):
            raise ValueError(f"slope shape {slope.shape} doesn't match ({self.ny}, {self.nx})")
        if aspect.shape != (self.ny, self.nx):
            raise ValueError(f"aspect shape {aspect.shape} doesn't match ({self.ny}, {self.nx})")
        
        return pyWildfire.update_terrain(elevation, slope, aspect, self.nx, self.ny)
    
    def get_terrain_info(self):
        """
        Get terrain information fields.
        
        Returns:
            dict: Dict with elevation, slope, aspect arrays
        """
        if not self.initialized:
            raise RuntimeError("Solver not initialized. Call initialize() first.")
        
        return pyWildfire.get_terrain_info()
    
    def get_ros_at_location(self, x, y):
        """
        Query ROS at a specific location.
        
        Parameters:
            x, y (float): Location coordinates (meters)
        
        Returns:
            float: ROS at location (m/s)
        """
        if not self.initialized:
            raise RuntimeError("Solver not initialized. Call initialize() first.")
        
        return pyWildfire.get_ros_at_location(x, y)
    
    def interpolate_field(self, field_name, x, y):
        """
        Interpolate any field to a specific location.
        
        Parameters:
            field_name (str): Field to interpolate ('phi', 'ros', 'intensity', etc.)
            x, y (float): Location coordinates (meters)
        
        Returns:
            float: Field value at location
        """
        if not self.initialized:
            raise RuntimeError("Solver not initialized. Call initialize() first.")
        
        return pyWildfire.interpolate_field(field_name, x, y)
    
    # ========================================================================
    # PHASE 3 & 4: Fire State Fields Methods
    # ========================================================================
    
    def get_ros_x(self):
        """Get ROS x-component field as numpy array."""
        if not self.initialized:
            raise RuntimeError("Solver not initialized. Call initialize() first.")
        
        data = pyWildfire.get_ros_x()
        return np.reshape(data, (self.ny, self.nx), order='F')
    
    def get_ros_y(self):
        """Get ROS y-component field as numpy array."""
        if not self.initialized:
            raise RuntimeError("Solver not initialized. Call initialize() first.")
        
        data = pyWildfire.get_ros_y()
        return np.reshape(data, (self.ny, self.nx), order='F')
    
    def get_ros_wind(self):
        """Get wind-driven ROS component as numpy array."""
        if not self.initialized:
            raise RuntimeError("Solver not initialized. Call initialize() first.")
        
        data = pyWildfire.get_ros_wind()
        return np.reshape(data, (self.ny, self.nx), order='F')
    
    def get_ros_slope(self):
        """Get slope-driven ROS component as numpy array."""
        if not self.initialized:
            raise RuntimeError("Solver not initialized. Call initialize() first.")
        
        data = pyWildfire.get_ros_slope()
        return np.reshape(data, (self.ny, self.nx), order='F')
    
    def get_residence_time(self):
        """Get fuel residence time field as numpy array."""
        if not self.initialized:
            raise RuntimeError("Solver not initialized. Call initialize() first.")
        
        data = pyWildfire.get_residence_time()
        return np.reshape(data, (self.ny, self.nx), order='F')
    
    def get_fuel_consumption(self):
        """Get fuel consumption field as numpy array."""
        if not self.initialized:
            raise RuntimeError("Solver not initialized. Call initialize() first.")
        
        data = pyWildfire.get_fuel_consumption()
        return np.reshape(data, (self.ny, self.nx), order='F')
    
    def get_front_curvature(self):
        """Get fire front curvature field as numpy array."""
        if not self.initialized:
            raise RuntimeError("Solver not initialized. Call initialize() first.")
        
        data = pyWildfire.get_front_curvature()
        return np.reshape(data, (self.ny, self.nx), order='F')
    
    def get_spread_direction(self):
        """Get primary spread direction field as numpy array (radians)."""
        if not self.initialized:
            raise RuntimeError("Solver not initialized. Call initialize() first.")
        
        data = pyWildfire.get_spread_direction()
        return np.reshape(data, (self.ny, self.nx), order='F')
    
    # ========================================================================
    # PHASE 5: Advanced Ignition & Control Methods
    # ========================================================================
    
    def set_ignition_region(self, xmin, xmax, ymin, ymax, time=0.0):
        """
        Set rectangular ignition region.
        
        Parameters:
            xmin, xmax, ymin, ymax (float): Region bounds (meters)
            time (float): Ignition time (seconds)
        
        Returns:
            bool: True if ignition was set
        """
        if not self.initialized:
            raise RuntimeError("Solver not initialized. Call initialize() first.")
        
        return pyWildfire.set_ignition_region(xmin, xmax, ymin, ymax, time)
    
    def set_ignition_from_array(self, phi_init):
        """
        Set custom ignition pattern from level-set array.
        
        Parameters:
            phi_init (ndarray): Initial level-set field, shape (ny, nx)
        
        Returns:
            bool: True if ignition was set
        """
        if not self.initialized:
            raise RuntimeError("Solver not initialized. Call initialize() first.")
        
        phi_init = np.asarray(phi_init, dtype=np.float64)
        if phi_init.shape != (self.ny, self.nx):
            raise ValueError(f"phi_init shape {phi_init.shape} doesn't match ({self.ny}, {self.nx})")
        
        return pyWildfire.set_ignition_from_array(phi_init)
    
    def set_spread_model(self, model_name):
        """
        Set propagation method.
        
        Parameters:
            model_name (str): 'levelset', 'richards', or 'hybrid'
        
        Returns:
            bool: True if model was set
        """
        if not self.initialized:
            raise RuntimeError("Solver not initialized. Call initialize() first.")
        
        valid_models = ['levelset', 'richards', 'hybrid']
        if model_name not in valid_models:
            raise ValueError(f"Invalid model '{model_name}'. Must be one of {valid_models}")
        
        return pyWildfire.set_spread_model(model_name)
    
    def step_with_subcycles(self, target_dt, max_subcycles=10):
        """
        Advance with subcycling control.
        
        Parameters:
            target_dt (float): Target timestep (seconds)
            max_subcycles (int): Maximum number of subcycles
        
        Returns:
            float: Actual timestep used
        """
        if not self.initialized:
            raise RuntimeError("Solver not initialized. Call initialize() first.")
        
        return pyWildfire.step_with_subcycles(target_dt, max_subcycles)
    
    def get_timestep_recommendation(self):
        """
        Get recommended next timestep.
        
        Returns:
            float: Recommended timestep (seconds)
        """
        if not self.initialized:
            raise RuntimeError("Solver not initialized. Call initialize() first.")
        
        return pyWildfire.get_timestep_recommendation()
    
    # ========================================================================
    # PHASE 6: Surface Fluxes & Emissions Methods
    # ========================================================================
    
    def get_all_surface_fluxes(self):
        """
        Get all surface flux components.
        
        Returns:
            dict: Dict with all flux fields as numpy arrays
        """
        if not self.initialized:
            raise RuntimeError("Solver not initialized. Call initialize() first.")
        
        flux_dict = pyWildfire.get_all_surface_fluxes()
        
        # Reshape all arrays from 1D to 2D
        result = {}
        for key, value in flux_dict.items():
            if isinstance(value, list) and len(value) == self.nx * self.ny:
                result[key] = np.reshape(value, (self.ny, self.nx), order='F')
            else:
                result[key] = value
        
        return result
    
    def get_emission_factors(self):
        """
        Get emission factors for species (per unit fuel consumed).
        
        Returns:
            dict: Dict with co2, co, ch4, pm25, nox, so2 factors
        """
        if not self.initialized:
            raise RuntimeError("Solver not initialized. Call initialize() first.")
        
        return pyWildfire.get_emission_factors()
    
    # ========================================================================
    # PHASE 7: Advanced I/O & Checkpointing Methods
    # ========================================================================
    
    def write_checkpoint(self, filename):
        """
        Write checkpoint to file.
        
        Parameters:
            filename (str): Checkpoint file path
        
        Returns:
            bool: True if write succeeded
        """
        if not self.initialized:
            raise RuntimeError("Solver not initialized. Call initialize() first.")
        
        return pyWildfire.write_checkpoint(filename)
    
    def read_checkpoint(self, filename):
        """
        Read checkpoint from file.
        
        Parameters:
            filename (str): Checkpoint file path
        
        Returns:
            bool: True if read succeeded
        """
        if not self.initialized:
            raise RuntimeError("Solver not initialized. Call initialize() first.")
        
        return pyWildfire.read_checkpoint(filename)
    
    def get_checkpoint_data(self):
        """
        Get checkpoint-compatible data dictionary.
        
        Returns:
            dict: Dict with all fields needed to restart
        """
        if not self.initialized:
            raise RuntimeError("Solver not initialized. Call initialize() first.")
        
        return pyWildfire.get_checkpoint_data()
    
    # ========================================================================
    # PHASE 8: Atmosphere Coupling & Diagnostics Methods
    # ========================================================================
    
    def set_fire_atmosphere_feedback_enabled(self, enabled):
        """
        Enable/disable fire-atmosphere feedback.
        
        Parameters:
            enabled (bool): True to enable feedback
        """
        if not self.initialized:
            raise RuntimeError("Solver not initialized. Call initialize() first.")
        
        pyWildfire.set_fire_atmosphere_feedback_enabled(enabled)
    
    def get_buoyancy_driven_winds(self):
        """
        Get induced wind from fire plume.
        
        Returns:
            dict: Dict with u_induced, v_induced, w_induced arrays
        """
        if not self.initialized:
            raise RuntimeError("Solver not initialized. Call initialize() first.")
        
        wind_dict = pyWildfire.get_buoyancy_driven_winds()
        
        # Reshape all arrays from 1D to 2D
        result = {}
        for key, value in wind_dict.items():
            if isinstance(value, list) and len(value) == self.nx * self.ny:
                result[key] = np.reshape(value, (self.ny, self.nx), order='F')
            else:
                result[key] = value
        
        return result
    
    def get_coupling_statistics(self):
        """
        Get fire-atmosphere coupling statistics.
        
        Returns:
            dict: Dict with total_heat_release, max_flame_height, wind_speed_at_fire, etc.
        """
        if not self.initialized:
            raise RuntimeError("Solver not initialized. Call initialize() first.")
        
        return pyWildfire.get_coupling_statistics()
    
    # ========================================================================
    # PHASE 9: GPU & Performance Methods
    # ========================================================================
    
    def set_accelerated_ros_computation(self, enabled):
        """
        Enable/disable GPU acceleration for ROS computation.
        
        Parameters:
            enabled (bool): True to enable GPU acceleration
        """
        if not self.initialized:
            raise RuntimeError("Solver not initialized. Call initialize() first.")
        
        pyWildfire.set_accelerated_ros_computation(enabled)
    
    def profile_ros_calculation(self):
        """
        Profile ROS calculation bottlenecks.
        
        Returns:
            dict: Dict with timing breakdown (milliseconds)
        """
        if not self.initialized:
            raise RuntimeError("Solver not initialized. Call initialize() first.")
        
        return pyWildfire.profile_ros_calculation()
    
    # ========================================================================
    # PHASE 10: Enhanced Diagnostics Methods
    # ========================================================================
    
    def get_wind_at_surface(self):
        """
        Get wind at surface with derived fields.
        
        Returns:
            dict: Dict with u, v, w, wind_speed, wind_direction arrays
        """
        if not self.initialized:
            raise RuntimeError("Solver not initialized. Call initialize() first.")
        
        wind_dict = pyWildfire.get_wind_at_surface()
        
        # Reshape all arrays from 1D to 2D
        result = {}
        for key, value in wind_dict.items():
            if isinstance(value, list) and len(value) == self.nx * self.ny:
                result[key] = np.reshape(value, (self.ny, self.nx), order='F')
            else:
                result[key] = value
        
        return result
    
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
