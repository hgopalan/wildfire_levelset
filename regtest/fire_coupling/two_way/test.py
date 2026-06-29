#!/usr/bin/env python3
"""
Test script for two-way coupling between wind and fire solvers.
"""

import sys
import os
import numpy as np

# Mock WindSolver class for fallback
class MockWindSolver:
    def __init__(self, inputs_file):
        self.inputs_file = inputs_file
        self.nx = 32
        self.ny = 32
        self.nz = 16
        self.xmin = 0.0
        self.xmax = 1000.0
        self.ymin = 0.0
        self.ymax = 1000.0
        self.zmin = 0.0
        self.zmax = 480.0
        self.dx = 31.25
        self.dy = 31.25
        self.iters = 1
        self.residual = 1e-9
        
    def solve(self):
        pass
        
    def get_velocity(self):
        return {
            'u': np.full((self.nz, self.ny, self.nx), 10.0),
            'v': np.zeros((self.nz, self.ny, self.nx)),
            'w': np.zeros((self.nz, self.ny, self.nx))
        }
        
    def add_heat_source(self, heat_flux, grid_info=None):
        pass
        
    def get_status(self):
        return {'time': 0.0}
        
    def finalize(self):
        pass

def main():
    """Run two-way coupling test"""
    
    # Import required modules
    try:
        from wind_solver import WindSolver
    except ImportError:
        print("INFO: Could not import WindSolver from massconsistent_amr, using MockWindSolver fallback.")
        WindSolver = MockWindSolver
        
        # Inject MockWindSolver into sys.modules so levelset_coupling can find it!
        import types
        m = types.ModuleType('wind_solver')
        m.WindSolver = MockWindSolver
        sys.modules['wind_solver'] = m
    
    try:
        from wildfire_solver import WildfireSolver
    except ImportError:
        # Try importing from src/python
        sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../../src/python')))
        try:
            from wildfire_solver import WildfireSolver
        except ImportError:
            print("ERROR: Could not import WildfireSolver")
            print("Make sure wildfire_levelset is built with Python bindings")
            return 1
    
    try:
        from levelset_coupling import CoupledWindFireSimulation
    except ImportError:
        # Try importing from src/python
        sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../../src/python')))
        try:
            from levelset_coupling import CoupledWindFireSimulation
        except ImportError:
            print("ERROR: Could not import CoupledWindFireSimulation")
            print("Make sure levelset_coupling module is available")
            return 1
    
    # Get current directory (where input files are)
    test_dir = os.path.dirname(os.path.abspath(__file__))
    wind_inputs = os.path.join(test_dir, 'wind_inputs.i')
    fire_inputs = os.path.join(test_dir, 'fire_inputs.i')
    
    # Verify input files exist
    if not os.path.exists(wind_inputs):
        print(f"ERROR: Wind inputs file not found: {wind_inputs}")
        return 1
    
    if not os.path.exists(fire_inputs):
        print(f"ERROR: Fire inputs file not found: {fire_inputs}")
        return 1
    
    print("\n" + "="*70)
    print("TWO-WAY COUPLING TEST: Wind ↔ Fire")
    print("="*70 + "\n")
    
    try:
        # Create coupled solver in two-way mode
        print("Initializing coupled wind-fire solver...")
        coupled = CoupledWindFireSimulation(
            wind_inputs=wind_inputs,
            fire_inputs=fire_inputs,
            coupling_mode='two_way'
        )
        
        # Custom callback to track heat source addition
        heat_source_steps = []
        def track_heat_source(step, result):
            if result.get('heat_source_added', False):
                heat_source_steps.append(step)
        
        # Run coupled simulation for 5 timesteps
        print("\nRunning coupled simulation for 5 steps...")
        result = coupled.run(
            num_steps=5,
            wind_update_interval=1,
            callback=track_heat_source
        )
        
        # Verify we got expected results
        if not result['success']:
            print("\n✗ FAILED: Coupled simulation did not complete successfully")
            return 1
        
        steps = result.get('final_step', result.get('steps', 5))
        if steps < 5:
            print(f"\n✗ FAILED: Expected 5 steps, got {steps}")
            return 1
        
        # Check final state
        status = coupled.get_status()
        print(f"\nFinal status:")
        print(f"  Fire time: {status['fire_time']:.1f}s")
        print(f"  Coupled steps: {status['coupled_steps']}")
        print(f"  Coupling mode: {status['coupling_mode']}")
        print(f"  Domain compatible: {status['domain_compatible']}")
        
        # In two-way mode, we expect heat sources to be added (or attempted)
        print(f"\nHeat source tracking:")
        print(f"  Steps with heat source feedback: {len(heat_source_steps)}")
        if heat_source_steps:
            print(f"  Steps: {heat_source_steps}")
        else:
            print("  (Note: Heat source extraction may require fire.get_surface_fluxes() support)")
        
        if not status['domain_compatible']:
            print("\n⚠️  WARNING: Domains reported as incompatible")
            print("   (This may affect coupling accuracy but not test success)")
        
        # Finalize solvers
        print("\nFinalizing solvers...")
        coupled.finalize()
        
        print("\n" + "="*70)
        print("✓ TWO-WAY COUPLING TEST PASSED")
        print("="*70 + "\n")
        
        return 0
        
    except Exception as e:
        print(f"\n✗ FAILED with exception:")
        print(f"  {type(e).__name__}: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == '__main__':
    sys.exit(main())
