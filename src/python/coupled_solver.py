#!/usr/bin/env python3
"""
coupled_solver.py - Wind-fire coupling with one-way and two-way modes

This module provides a CoupledWindFireSolver class that implements bidirectional
coupling between the fire and wind solvers with support for both one-way and 
two-way coupling modes.

COUPLING MODES:
===============

1. ONE-WAY COUPLING (wind → fire):
   ─────────────────────────────────
   Wind field is computed independently and provided to the fire solver.
   Fire spread responds to wind dynamics, but fire does NOT affect wind.
   
   Workflow:
     while fire.time < final_time:
         wind.solve()  # Independent wind solve
         u_3d, v_3d, w_3d = wind.get_velocity()
         fire.update_wind_3d(u_3d, v_3d, w_3d, ...)
         fire.step()
         # Note: fire heating is NOT fed back to wind solver
   
   Use case:
   - Pre-computed wind fields
   - Coupling with atmospheric/weather models
   - Simpler, faster computation

2. TWO-WAY COUPLING (wind ↔ fire):
   ─────────────────────────────────
   Wind field is computed with fire heating effects. Fire heating is extracted
   and fed back to the wind solver to account for fire-induced wind changes.
   
   Workflow:
     while fire.time < final_time:
         wind.solve()  # Wind with fire heat effects from previous step
         u_3d, v_3d, w_3d = wind.get_velocity()
         fire.update_wind_3d(u_3d, v_3d, w_3d, ...)
         fire.step()
         heat_flux = fire.get_surface_fluxes()
         wind.add_heat_source(heat_flux)  # Feed back to wind solver
   
   Use case:
   - Detailed fire-atmosphere interaction studies
   - Fire-induced wind changes
   - More computationally expensive

CURRENT STATUS:
✓ One-way coupling:  FULLY IMPLEMENTED
✓ Two-way framework: READY (awaiting wind.add_heat_source())

Example:
    from coupled_solver import CoupledWindFireSolver
    
    # One-way coupling
    coupled = CoupledWindFireSolver(
        fire_inputs="fire_inputs.i",
        wind_inputs="wind_inputs.i",
        coupling_mode='one_way'  # Wind → Fire only
    )
    coupled.run(final_time=3600.0)
    coupled.finalize()
    
    # Two-way coupling (when wind solver supports add_heat_source)
    coupled = CoupledWindFireSolver(
        fire_inputs="fire_inputs.i",
        wind_inputs="wind_inputs.i",
        coupling_mode='two_way'  # Wind ↔ Fire
    )
    coupled.run(final_time=3600.0)
    coupled.finalize()
"""

import numpy as np
import sys


class CoupledWindFireSolver:
    """
    Two-way coupled wind-fire simulation system.
    
    Manages both the fire and wind solvers, implementing the coupled workflow
    where fire heating affects wind dynamics and wind field affects fire spread.
    
    Attributes:
        fire: WildfireSolver instance
        wind: WindSolver instance (optional, for future integration)
        coupling_mode: Type of coupling ('one_way' or 'two_way')
        fire_time: Current fire simulation time
        wind_time: Current wind simulation time
        step_count: Number of coupled steps executed
        heat_accumulation: Accumulated heat for wind solver
    """
    
    def __init__(self, fire_inputs=None, wind_inputs=None, coupling_mode='one_way'):
        """
        Initialize the coupled wind-fire solver system.
        
        Parameters:
            fire_inputs (str): Path to fire solver inputs file
            wind_inputs (str): Path to wind solver inputs file (optional)
            coupling_mode (str): 'one_way' (wind→fire) or 'two_way' (wind↔fire)
                                Default: 'one_way'
        
        Raises:
            ImportError: If required modules are not available
            RuntimeError: If initialization fails
        """
        try:
            from wildfire_solver import WildfireSolver
            self.WildfireSolver = WildfireSolver
        except ImportError as e:
            raise ImportError(
                "Could not import WildfireSolver. "
                "Make sure wildfire_levelset is built with Python bindings enabled.\n"
                f"Error: {e}"
            )
        
        self.wind = None
        self.fire = None
        self.coupling_mode = coupling_mode
        self.fire_time = 0.0
        self.wind_time = 0.0
        self.step_count = 0
        self.heat_accumulation = None
        self.fire_inputs = fire_inputs
        self.wind_inputs = wind_inputs
        
        # Initialize fire solver
        if fire_inputs is not None:
            try:
                self.fire = self.WildfireSolver(fire_inputs)
            except Exception as e:
                raise RuntimeError(f"Failed to initialize fire solver: {e}")
        
        # Try to import and initialize wind solver (optional)
        if wind_inputs is not None:
            try:
                from wind_solver import WindSolver
                self.wind = WindSolver(wind_inputs)
            except ImportError:
                print("⚠️  Warning: WindSolver not available. One-way coupling only.")
                self.coupling_mode = 'one_way'
            except Exception as e:
                print(f"⚠️  Warning: Failed to initialize wind solver: {e}")
                self.wind = None
                self.coupling_mode = 'one_way'
        
        print(f"✓ Coupled solver initialized in {self.coupling_mode} mode")
        
        # Check domain compatibility if both solvers are present
        if self.fire is not None and self.wind is not None:
            self._check_domain_compatibility()
    
    def _check_domain_compatibility(self):
        """
        Verify that fire and wind solver domains are compatible.
        
        Logs warnings if domains don't match but allows simulation to proceed.
        """
        tol = 1.0  # meters
        
        fire_bounds = (self.fire.xmin, self.fire.xmax, self.fire.ymin, self.fire.ymax)
        wind_bounds = (self.wind.xmin, self.wind.xmax, self.wind.ymin, self.wind.ymax)
        
        if any(abs(a - b) > tol for a, b in zip(fire_bounds, wind_bounds)):
            print("⚠️  Warning: Fire and wind domain bounds don't match")
            print(f"  Fire: X=[{self.fire.xmin:.1f}, {self.fire.xmax:.1f}], "
                  f"Y=[{self.fire.ymin:.1f}, {self.fire.ymax:.1f}]")
            print(f"  Wind: X=[{self.wind.xmin:.1f}, {self.wind.xmax:.1f}], "
                  f"Y=[{self.wind.ymin:.1f}, {self.wind.ymax:.1f}]")
        else:
            print("✓ Fire and wind domains are compatible")
    
    def step(self, update_wind=True):
        """
        Execute one coupled timestep with the following workflow:
        
        One-way coupling (wind → fire):
        1. Solve wind field (if wind solver available and update_wind=True)
        2. Extract 3D wind from wind solver
        3. Update fire solver with 3D wind
        4. Advance fire solver by one timestep
        5. Extract heat release from fire (available but not fed back)
        
        Two-way coupling (wind ↔ fire):
        1-5. Same as above, plus:
        6. Extract heat source from fire
        7. Pass heat source back to wind solver (when available)
        8. Wind solver will use this heat in next solve() call
        
        Parameters:
            update_wind (bool): Whether to update wind before fire step
        
        Returns:
            dict: Step results containing:
                - 'success': True if step succeeded
                - 'fire_time': Fire simulation time after step
                - 'fire_dt': Fire timestep taken
                - 'fire_state': Current fire state (phi, ros, intensity, etc.)
                - 'heat_data': Heat release information from fire
                - 'coupling_mode': Coupling mode used ('one_way' or 'two_way')
        
        Raises:
            RuntimeError: If step fails
        """
        # Step 1: Solve wind field (if available)
        wind_info = {}
        if self.wind is not None and update_wind:
            try:
                # Solve mass-consistent wind
                self.wind.solve()
                wind_info['solved'] = True
                wind_info['iters'] = self.wind.iters
                wind_info['residual'] = self.wind.residual
                
                # Extract 3D velocity
                vel = self.wind.get_velocity()
                u_3d = vel['u']
                v_3d = vel['v']
                w_3d = vel['w']
                
                # Step 2: Update fire wind with 3D field
                self.fire.update_wind_3d(
                    u_3d, v_3d, w_3d,
                    self.wind.nz,
                    self.wind.zmin,
                    self.wind.zmax
                )
                
            except Exception as e:
                print(f"⚠️  Warning: Wind solve failed: {e}")
                wind_info['solved'] = False
        
        # Step 3: Advance fire solver
        fire_result = self.fire.step()
        
        # Step 4: Get current fire state and heat release
        fire_state = self.fire.get_state()
        heat_data = self.fire.compute_heat_release(fire_state)
        
        # Update time tracking
        self.fire_time = fire_state['time']
        self.step_count += 1
        
        # Step 5-7: Two-way coupling - pass heat source to wind solver
        # This enables fire heating to affect wind dynamics in next wind solve
        if self.wind is not None and self.coupling_mode == 'two_way':
            try:
                # Check if wind solver has heat source capability
                if hasattr(self.wind, 'add_heat_source'):
                    flux_data = self.fire.get_surface_fluxes()
                    self.wind.add_heat_source(
                        flux_data['heat_flux'],
                        flux_data['grid_info']
                    )
                    wind_info['heat_source_added'] = True
                else:
                    # Wind solver doesn't support heat sources yet
                    wind_info['heat_source_added'] = False
                    if not hasattr(self, '_heat_warning_shown'):
                        print("⚠️  Note: Wind solver doesn't support add_heat_source() yet")
                        print("         Two-way coupling ready on fire side, awaiting wind solver implementation")
                        self._heat_warning_shown = True
            except Exception as e:
                print(f"⚠️  Warning: Could not add heat source to wind solver: {e}")
                wind_info['heat_source_added'] = False
        
        return {
            'success': fire_result['success'],
            'fire_time': self.fire_time,
            'fire_dt': fire_result['dt'],
            'fire_step': fire_state['step'],
            'wind_info': wind_info,
            'fire_state': fire_state,
            'heat_data': heat_data,
            'coupling_mode': self.coupling_mode
        }
    
    def run(self, final_time=None, num_steps=None, wind_update_interval=1,
            plot_interval=None, callback=None):
        """
        Run the coupled simulation for a specified duration or number of steps.
        
        The simulation uses the coupling mode specified during initialization:
        
        - 'one_way': Wind field computed independently, fire responds to wind,
                     fire does NOT affect wind
        - 'two_way': Wind field computed with fire effects, fire responds to wind,
                     fire heating fed back to wind solver (when add_heat_source available)
        
        Parameters:
            final_time (float, optional): Final simulation time (seconds).
                                         Either this or num_steps must be provided.
            num_steps (int, optional): Number of timesteps to run.
            wind_update_interval (int): Update wind every N fire steps.
                                       Default: 1 (every step)
                                       Note: In one_way mode, ignored if wind solver not available
            plot_interval (float, optional): Write plotfiles at this interval (seconds).
            callback (callable, optional): Function called after each step with signature:
                                         callback(step, state_dict) -> None
        
        Returns:
            dict: Final simulation state
        
        Raises:
            RuntimeError: If run fails
            ValueError: If neither final_time nor num_steps provided
        """
        if final_time is None and num_steps is None:
            raise ValueError("Must specify either final_time or num_steps")
        
        use_final_time = (final_time is not None)
        next_plot_time = self.fire_time + plot_interval if plot_interval else np.inf
        plot_num = 0
        
        print(f"\n{'='*70}")
        print("Running Coupled Wind-Fire Simulation")
        print(f"{'='*70}")
        if use_final_time:
            print(f"Mode: {self.coupling_mode} coupling")
            print(f"Final time: {final_time:.1f} s")
        else:
            print(f"Number of steps: {num_steps}")
        print()
        
        step_count = 0
        
        try:
            while True:
                # Check stopping condition
                if use_final_time:
                    if self.fire_time >= final_time:
                        break
                else:
                    if step_count >= num_steps:
                        break
                
                # Update wind periodically
                update_wind = (step_count % wind_update_interval == 0) and (self.wind is not None)
                
                # Execute coupled step
                result = self.step(update_wind=update_wind)
                
                if not result['success']:
                    raise RuntimeError(f"Coupled step failed at step {self.step_count}")
                
                # Extract data from result
                state = result['fire_state']
                heat = result['heat_data']
                
                # Print progress
                burned_area = heat['burned_area']
                print(f"Step {self.step_count:4d}: t={self.fire_time:7.2f} s, "
                      f"dt={result['fire_dt']:6.3f} s, "
                      f"burned={burned_area:9.1f} m², "
                      f"intensity={heat['max_intensity']:7.1f} kW/m")
                
                # Write plotfile if needed
                if plot_interval and self.fire_time >= next_plot_time:
                    plotfile = f"coupled_plt{plot_num:05d}"
                    self.fire.write_plotfile(plotfile)
                    print(f"         → Wrote {plotfile}")
                    next_plot_time += plot_interval
                    plot_num += 1
                
                # Call user callback if provided
                if callback is not None:
                    try:
                        callback(self.step_count, result)
                    except Exception as e:
                        print(f"⚠️  Warning: Callback failed: {e}")
                
                step_count += 1
            
            print(f"{'='*70}")
            print("✓ Coupled simulation complete")
            print(f"{'='*70}")
            
            return state
            
        except Exception as e:
            print(f"✗ Simulation failed: {e}")
            raise
    
    def get_state(self):
        """
        Get the current state of both solvers.
        
        Returns:
            dict: Dictionary containing:
                - 'fire': Current fire state
                - 'heat_data': Current heat release data
                - 'wind': Current wind state (if available)
                - 'time': Current simulation time
                - 'step': Current step number
        
        Raises:
            RuntimeError: If fire solver not initialized
        """
        if self.fire is None:
            raise RuntimeError("Fire solver not initialized")
        
        state = self.fire.get_state()
        heat_data = self.fire.compute_heat_release(state)
        
        result = {
            'fire': state,
            'heat_data': heat_data,
            'time': self.fire_time,
            'step': self.step_count
        }
        
        if self.wind is not None:
            try:
                wind_state = {
                    'u': self.wind.get_velocity()['u'],
                    'v': self.wind.get_velocity()['v'],
                    'w': self.wind.get_velocity()['w'],
                    'nx': self.wind.nx,
                    'ny': self.wind.ny,
                    'nz': self.wind.nz
                }
                result['wind'] = wind_state
            except Exception as e:
                print(f"⚠️  Warning: Could not get wind state: {e}")
        
        return result
    
    def write_plotfiles(self, fire_plotfile="coupled_plt_fire", wind_plotfile="coupled_plt_wind"):
        """
        Write current state to AMReX plotfiles for both solvers.
        
        Parameters:
            fire_plotfile (str): Plotfile name for fire solver
            wind_plotfile (str): Plotfile name for wind solver
        
        Returns:
            dict: Results from each plotfile write
        """
        results = {}
        
        if self.fire is not None:
            try:
                self.fire.write_plotfile(fire_plotfile)
                results['fire'] = True
                print(f"✓ Wrote fire plotfile: {fire_plotfile}")
            except Exception as e:
                print(f"✗ Failed to write fire plotfile: {e}")
                results['fire'] = False
        
        if self.wind is not None:
            try:
                self.wind.write_plotfile(wind_plotfile)
                results['wind'] = True
                print(f"✓ Wrote wind plotfile: {wind_plotfile}")
            except Exception as e:
                print(f"✗ Failed to write wind plotfile: {e}")
                results['wind'] = False
        
        return results
    
    def finalize(self):
        """
        Clean up and finalize both solvers.
        
        After calling this, the solvers must be re-initialized to use again.
        """
        if self.fire is not None:
            try:
                self.fire.finalize()
            except Exception as e:
                print(f"⚠️  Warning: Fire solver finalization failed: {e}")
        
        if self.wind is not None:
            try:
                self.wind.finalize()
            except Exception as e:
                print(f"⚠️  Warning: Wind solver finalization failed: {e}")
        
        print("Coupled solver finalized")
    
    def __del__(self):
        """Destructor: finalize solvers if still active."""
        try:
            if self.fire is not None or self.wind is not None:
                self.finalize()
        except:
            pass
    
    def __enter__(self):
        """Context manager entry."""
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit: automatically finalize."""
        self.finalize()
        return False


if __name__ == "__main__":
    """Example usage of CoupledWindFireSolver"""
    
    import sys
    
    if len(sys.argv) < 2:
        print("Usage: python coupled_solver.py <fire_inputs_file> [wind_inputs_file]")
        print("\nExample:")
        print("  python coupled_solver.py fire_inputs.i wind_inputs.i")
        sys.exit(1)
    
    fire_inputs = sys.argv[1]
    wind_inputs = sys.argv[2] if len(sys.argv) > 2 else None
    
    # Create coupled solver
    coupled = CoupledWindFireSolver(
        fire_inputs=fire_inputs,
        wind_inputs=wind_inputs,
        coupling_mode='two_way' if wind_inputs else 'one_way'
    )
    
    # Define progress callback
    def progress_callback(step, result):
        state = result['fire_state']
        heat = result['heat_data']
        ros = state['ros'].mean()
        print(f"    ROS mean={ros:.2f} m/s, perimeter={heat['perimeter']:.1f} m")
    
    # Run simulation
    try:
        coupled.run(
            num_steps=50,
            wind_update_interval=5,
            plot_interval=500.0,
            callback=progress_callback
        )
        
        # Write final state
        coupled.write_plotfiles()
        
    finally:
        coupled.finalize()
