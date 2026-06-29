#!/usr/bin/env python3
"""
levelset_coupling.py - Coupling module for wildfire_levelset fire solver integration

This module provides utilities to couple the massconsistent_amr wind solver with the
wildfire_levelset fire solver, supporting both one-way (wind→fire) and two-way
(wind↔fire) coupling modes.

COUPLING MODES:
===============

1. ONE-WAY COUPLING (wind → fire):
   Wind field is computed independently and provided to the fire solver.
   Fire spread responds to wind dynamics, but fire does NOT affect wind.
   
   Use case:
   - Pre-computed or prescribed wind fields
   - Coupling with atmospheric/weather models
   - Simpler, faster computation

2. TWO-WAY COUPLING (wind ↔ fire):
   Wind field is computed with fire heating effects. Fire heating is extracted
   and fed back to the wind solver for fire-induced wind changes.
   
   Use case:
   - Detailed fire-atmosphere interaction studies
   - Fire-induced wind changes (updrafts, flow deflection)
   - More computationally expensive

Example - One-way coupling:
    from levelset_coupling import CoupledWindFireSimulation
    
    coupled = CoupledWindFireSimulation(
        wind_inputs="wind_inputs.i",
        fire_inputs="fire_inputs.i",
        coupling_mode='one_way'
    )
    coupled.run(final_time=3600.0)
    coupled.finalize()

Example - Two-way coupling:
    coupled = CoupledWindFireSimulation(
        wind_inputs="wind_inputs.i",
        fire_inputs="fire_inputs.i",
        coupling_mode='two_way'
    )
    coupled.run(final_time=3600.0)
    coupled.finalize()
"""

import numpy as np
import os
import sys
from typing import Optional, Dict, Tuple, Callable


class CoupledWindFireSimulation:
    """
    Two-way coupled wind-fire simulation system.
    
    Manages both the fire and wind solvers, implementing the coupled workflow where
    fire heating can affect wind dynamics and wind field affects fire spread.
    
    Supports both one-way (wind→fire) and two-way (wind↔fire) coupling modes.
    
    Attributes:
        wind: WindSolver instance
        fire: WildfireSolver instance
        coupling_mode: Type of coupling ('one_way' or 'two_way')
        fire_time: Current fire simulation time (seconds)
        wind_time: Current wind simulation time (seconds)
        step_count: Number of coupled steps executed
        domain_compatible: Whether fire and wind domains are compatible
    """
    
    def __init__(self, 
                 wind_inputs: str, 
                 fire_inputs: str, 
                 coupling_mode: str = 'one_way'):
        """
        Initialize the coupled wind-fire solver system.
        
        Parameters:
            wind_inputs (str): Path to wind solver inputs file
            fire_inputs (str): Path to fire solver inputs file
            coupling_mode (str): 'one_way' (wind→fire) or 'two_way' (wind↔fire)
                                Default: 'one_way'
        
        Raises:
            ImportError: If required solver modules are not available
            RuntimeError: If initialization fails
        """
        # Import solvers - try external wind_solver first (from massconsistent_amr)
        # then fall back to local wind_solver if available
        try:
            # First try to import from massconsistent_amr
            from wind_solver import WindSolver
            self.WindSolver = WindSolver
            self._wind_source = "massconsistent_amr"
        except ImportError:
            try:
                # Fall back to local wind_solver if it exists
                import sys
                import importlib.util
                spec = importlib.util.find_spec("wind_solver")
                if spec is not None:
                    from wind_solver import WindSolver
                    self.WindSolver = WindSolver
                    self._wind_source = "local"
                else:
                    raise ImportError("No wind solver found")
            except ImportError as e:
                raise ImportError(
                    "Could not import WindSolver. "
                    "Make sure massconsistent_amr is built with Python bindings enabled, "
                    "or a local wind_solver module is available.\n"
                    f"Error: {e}"
                )
        
        try:
            from wildfire_solver import WildfireSolver
            self.WildfireSolver = WildfireSolver
        except ImportError as e:
            raise ImportError(
                "Could not import WildfireSolver. "
                "Make sure wildfire_levelset is built with Python bindings enabled.\n"
                f"Error: {e}"
            )
        
        # Initialize state
        self.wind = None
        self.fire = None
        self.coupling_mode = coupling_mode.lower()
        self.fire_time = 0.0
        self.wind_time = 0.0
        self.step_count = 0
        self.domain_compatible = False
        self.wind_inputs = wind_inputs
        self.fire_inputs = fire_inputs
        
        # Validate coupling mode
        if self.coupling_mode not in ('one_way', 'two_way'):
            raise ValueError(f"Invalid coupling_mode '{coupling_mode}'. Must be 'one_way' or 'two_way'.")
        
        # Initialize wind solver
        try:
            self.wind = self.WindSolver(wind_inputs)
        except Exception as e:
            raise RuntimeError(f"Failed to initialize wind solver: {e}")
        
        # Initialize fire solver
        try:
            self.fire = self.WildfireSolver(fire_inputs)
        except Exception as e:
            raise RuntimeError(f"Failed to initialize fire solver: {e}")
        
        # Check domain compatibility
        self._check_domain_compatibility()
        
        print(f"\n✓ Coupled {self.coupling_mode} wind-fire solver initialized")
        print(f"  Wind: {self.wind.nx} × {self.wind.ny} × {self.wind.nz}")
        print(f"  Fire: {self.fire.nx} × {self.fire.ny}")
        print(f"  Domain compatible: {self.domain_compatible}")
    
    def _check_domain_compatibility(self, tolerance: float = 1.0):
        """
        Verify that fire and wind solver domains are compatible.
        
        For the coupling to work correctly:
        - Horizontal domains should match (within tolerance)
        - Grid spacings (dx, dy) should match (within tolerance)
        
        Parameters:
            tolerance (float): Tolerance for domain matching (meters). Default: 1.0
        """
        # Check bounds
        bounds_match = (
            abs(self.wind.xmin - self.fire.xmin) <= tolerance and
            abs(self.wind.xmax - self.fire.xmax) <= tolerance and
            abs(self.wind.ymin - self.fire.ymin) <= tolerance and
            abs(self.wind.ymax - self.fire.ymax) <= tolerance
        )
        
        # Check grid spacing
        spacing_match = (
            abs(self.wind.dx - self.fire.dx) <= tolerance and
            abs(self.wind.dy - self.fire.dy) <= tolerance
        )
        
        self.domain_compatible = bounds_match and spacing_match
        
        if not bounds_match:
            print("⚠️  Warning: Wind and fire domain bounds don't match")
            print(f"  Wind: X=[{self.wind.xmin:.1f}, {self.wind.xmax:.1f}], "
                  f"Y=[{self.wind.ymin:.1f}, {self.wind.ymax:.1f}]")
            print(f"  Fire: X=[{self.fire.xmin:.1f}, {self.fire.xmax:.1f}], "
                  f"Y=[{self.fire.ymin:.1f}, {self.fire.ymax:.1f}]")
        
        if not spacing_match:
            print("⚠️  Warning: Wind and fire grid spacing don't match")
            print(f"  Wind: dx={self.wind.dx:.2f}, dy={self.wind.dy:.2f}")
            print(f"  Fire: dx={self.fire.dx:.2f}, dy={self.fire.dy:.2f}")
        
        if self.domain_compatible:
            print("✓ Wind and fire domains are compatible")
    
    def step(self, update_wind: bool = True) -> Dict:
        """
        Execute one coupled timestep.
        
        One-way coupling workflow:
        1. Solve wind field
        2. Extract 3D wind from wind solver
        3. Update fire solver with 3D wind
        4. Advance fire solver by one timestep
        
        Two-way coupling workflow:
        1-4. Same as above, plus:
        5. Extract heat source from fire
        6. Pass heat source back to wind solver
        
        Parameters:
            update_wind (bool): Whether to solve wind before fire step. Default: True
        
        Returns:
            dict: Step results containing:
                - 'success': True if step succeeded
                - 'fire_time': Fire simulation time after step
                - 'fire_dt': Fire timestep taken
                - 'fire_step': Current fire step number
                - 'wind_solved': Whether wind was solved in this step
                - 'wind_iters': MLMG iterations (if wind solved)
                - 'wind_residual': MLMG residual (if wind solved)
                - 'heat_source_added': Whether heat source was fed back to wind
                - 'fire_state': Current fire state (phi, ros, intensity, etc.)
        
        Raises:
            RuntimeError: If step fails
        """
        result = {'success': True}
        
        # Step 1-3: Solve wind and update fire wind field
        if update_wind:
            try:
                self.wind.solve()
                result['wind_solved'] = True
                result['wind_iters'] = self.wind.iters
                result['wind_residual'] = self.wind.residual
                
                # Extract 3D velocity
                vel = self.wind.get_velocity()
                u_3d = vel['u']  # shape (nz, ny, nx)
                v_3d = vel['v']
                w_3d = vel['w']
                
                # Update fire with 3D wind
                self.fire.update_wind_3d(
                    u_3d, v_3d, w_3d,
                    self.wind.nz,
                    self.wind.zmin,
                    self.wind.zmax
                )
                
                self.wind_time = self.wind.get_status().get('time', 0.0)
                
            except Exception as e:
                print(f"⚠️  Warning: Wind solve failed: {e}")
                result['wind_solved'] = False
        else:
            result['wind_solved'] = False
        
        # Step 4: Advance fire solver
        try:
            fire_result = self.fire.step()
            result['success'] = fire_result['success']
        except Exception as e:
            raise RuntimeError(f"Fire solver step failed: {e}")
        
        # Get fire state
        fire_state = self.fire.get_state()
        self.fire_time = fire_state['time']
        self.step_count += 1
        
        result['fire_time'] = self.fire_time
        result['fire_dt'] = fire_result.get('dt', 0.0)
        result['fire_step'] = fire_state.get('step', self.step_count)
        result['fire_state'] = fire_state
        result['heat_source_added'] = False
        
        # Step 5-6: Two-way coupling - extract and pass heat source to wind
        if self.coupling_mode == 'two_way' and update_wind:
            try:
                # Extract heat source from fire
                heat_data = self.fire.get_surface_fluxes()
                
                if heat_data is not None:
                    heat_flux = heat_data.get('heat_flux')
                    
                    if heat_flux is not None:
                        # Pass heat source to wind solver
                        grid_info = {
                            'xmin': self.wind.xmin,
                            'xmax': self.wind.xmax,
                            'ymin': self.wind.ymin,
                            'ymax': self.wind.ymax,
                            'dx': self.wind.dx,
                            'dy': self.wind.dy
                        }
                        
                        self.wind.add_heat_source(heat_flux, grid_info)
                        result['heat_source_added'] = True
            
            except AttributeError:
                # Fire solver doesn't have get_surface_fluxes method
                pass
            except Exception as e:
                print(f"⚠️  Warning: Could not extract/add heat source: {e}")
        
        return result
    
    def run(self,
            final_time: Optional[float] = None,
            num_steps: Optional[int] = None,
            wind_update_interval: int = 1,
            plot_interval: Optional[float] = None,
            callback: Optional[Callable] = None) -> Dict:
        """
        Run the coupled simulation.
        
        The simulation uses the coupling mode specified during initialization:
        
        - 'one_way': Wind computed independently, fire responds to wind,
                     fire does NOT affect wind
        - 'two_way': Wind computed with fire effects, fire responds to wind,
                     fire heating fed back to wind solver
        
        Parameters:
            final_time (float, optional): Final simulation time (seconds).
                                         Either this or num_steps must be provided.
            num_steps (int, optional): Number of timesteps to run.
            wind_update_interval (int): Update wind every N fire steps.
                                       Default: 1 (every step)
            plot_interval (float, optional): Write plotfiles at this interval (seconds).
            callback (callable, optional): Function called after each step with signature:
                                          callback(step, result_dict) -> None
        
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
        print(f"Coupling mode: {self.coupling_mode}")
        if use_final_time:
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
                update_wind = (step_count % wind_update_interval == 0)
                
                # Execute coupled step
                result = self.step(update_wind=update_wind)
                
                if not result['success']:
                    raise RuntimeError(f"Coupled step failed at step {self.step_count}")
                
                # Extract data from result
                state = result['fire_state']
                
                # Write plotfile if needed
                if plot_interval is not None and self.fire_time >= next_plot_time:
                    plotfile = f"plt_{plot_num:05d}"
                    self.fire.write_plotfile(plotfile)
                    print(f"  Wrote plotfile: {plotfile} at t={self.fire_time:.1f}s")
                    next_plot_time += plot_interval
                    plot_num += 1
                
                # Print progress
                if step_count % 10 == 0:
                    wind_status = "✓" if result['wind_solved'] else "✗"
                    heat_status = "✓" if result.get('heat_source_added', False) else "✗"
                    
                    if self.coupling_mode == 'one_way':
                        print(f"Step {step_count:4d}: t={self.fire_time:8.1f}s  "
                              f"Wind{wind_status}  "
                              f"Burned={np.sum(state['phi'] <= 0):.0f} cells")
                    else:
                        print(f"Step {step_count:4d}: t={self.fire_time:8.1f}s  "
                              f"Wind{wind_status}  "
                              f"Heat{heat_status}  "
                              f"Burned={np.sum(state['phi'] <= 0):.0f} cells")
                
                # Call user callback if provided
                if callback is not None:
                    callback(step_count, result)
                
                step_count += 1
        
        except Exception as e:
            print(f"✗ Simulation failed: {e}")
            raise
        
        print(f"\n✓ Simulation completed at t={self.fire_time:.1f}s")
        print(f"  Total steps: {self.step_count}")
        
        return {
            'final_time': self.fire_time,
            'final_step': self.step_count,
            'fire_state': state
        }
    
    def finalize(self):
        """
        Clean up and finalize both solvers.
        """
        if self.fire is not None:
            try:
                self.fire.finalize()
            except Exception as e:
                print(f"⚠️  Warning: Error finalizing fire solver: {e}")
        
        if self.wind is not None:
            try:
                self.wind.finalize()
            except Exception as e:
                print(f"⚠️  Warning: Error finalizing wind solver: {e}")
        
        print("✓ Simulation finalized")
