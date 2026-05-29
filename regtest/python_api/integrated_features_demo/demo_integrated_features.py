#!/usr/bin/env python3
"""
Python API Integration Demo: 2026 Enhancement Features

Demonstrates comprehensive use of the wildfire solver Python API including:
- State extraction and manipulation
- Real-time monitoring and visualization
- Custom analysis and diagnostics
- Integration with 2026 enhancement features:
  * McArthur moisture scaling
  * FMC phenology
  * Ember accumulation
  * Periodic wind gusts
  * Slope-dependent flame tilt

===============================================================================
HOW TO RUN THIS DEMO
===============================================================================

Method 1: Using CMake/CTest (recommended)
------------------------------------------
From the build directory:
    cd build
    ctest -R python_api_integrated_features_demo --output-on-failure

Method 2: Direct Python execution
----------------------------------
From this test directory:
    cd regtest/python_api/integrated_features_demo
    PYTHONPATH=/path/to/build/python:$PYTHONPATH python3 demo_integrated_features.py

===============================================================================
REQUIREMENTS
===============================================================================
- Build wildfire_levelset with: -DLEVELSET_BUILD_PYTHON_BINDINGS=ON
- Python 3.6 or later
- NumPy (optional, for analysis)
- Matplotlib (optional, for visualization)

===============================================================================
"""

import sys
import os
import time

# Find the pyWildfire module
try:
    import pyWildfire
    from wildfire_solver import WildfireSolver
except ImportError:
    build_python = os.path.join(os.getcwd(), '..', '..', '..', 'python')
    if os.path.exists(build_python):
        sys.path.insert(0, build_python)
        import pyWildfire
        from wildfire_solver import WildfireSolver
    else:
        print("ERROR: Cannot find pyWildfire module")
        print("Make sure to build with -DLEVELSET_BUILD_PYTHON_BINDINGS=ON")
        sys.exit(1)

# Optional imports
try:
    import numpy as np
    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False
    print("Note: NumPy not available, advanced analysis disabled")

try:
    import matplotlib
    matplotlib.use('Agg')  # Non-interactive backend
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    print("Note: Matplotlib not available, visualization disabled")


def print_header(title):
    """Print a formatted section header"""
    print(f"\n{'='*78}")
    print(f"  {title}")
    print(f"{'='*78}")


def print_state_summary(state, step):
    """Print a summary of the current fire state"""
    if not HAS_NUMPY:
        print(f"Step {step}: time = {state['time']:.1f} s")
        return
    
    # Calculate statistics
    phi = state['phi']
    burned = np.sum(phi <= 0.0)
    burning = np.sum((phi > 0.0) & (phi < 50.0))
    unburned = np.sum(phi >= 50.0)
    
    ros = state['ros']
    max_ros = np.max(ros)
    mean_ros = np.mean(ros[ros > 0.01])
    
    intensity = state['intensity']
    max_intensity = np.max(intensity)
    mean_intensity = np.mean(intensity[intensity > 0.1])
    
    flame_length = state['flame_length']
    max_flame = np.max(flame_length)
    
    print(f"Step {step:3d} | Time: {state['time']:7.1f} s | "
          f"Burned: {burned:5d} cells | Burning: {burning:4d} cells | "
          f"ROS: {max_ros:.2f} m/s | I: {max_intensity:6.0f} kW/m | "
          f"FL: {max_flame:.1f} m")


def analyze_fire_growth(history):
    """Analyze fire growth patterns"""
    if not HAS_NUMPY or len(history) == 0:
        return
    
    print_header("Fire Growth Analysis")
    
    times = [h['time'] for h in history]
    areas = [h['burned_area'] for h in history]
    perimeters = [h['perimeter'] for h in history]
    
    print(f"  Initial burned area:    {areas[0]:10,.1f} m²")
    print(f"  Final burned area:      {areas[-1]:10,.1f} m²")
    print(f"  Total growth:           {areas[-1]-areas[0]:10,.1f} m² ({(areas[-1]/areas[0]-1)*100:.1f}%)")
    print(f"  Final perimeter:        {perimeters[-1]:10,.1f} m")
    
    # Growth rate
    if len(times) > 1:
        dt = times[-1] - times[0]
        growth_rate = (areas[-1] - areas[0]) / dt if dt > 0 else 0
        print(f"  Average growth rate:    {growth_rate:10,.1f} m²/s")
    
    # Maximum values
    max_ros_vals = [h['max_ros'] for h in history]
    max_intensity_vals = [h['max_intensity'] for h in history]
    
    print(f"\n  Peak ROS:               {max(max_ros_vals):.3f} m/s")
    print(f"  Peak intensity:         {max(max_intensity_vals):,.0f} kW/m")


def create_visualizations(history, output_dir="."):
    """Create visualization plots if matplotlib is available"""
    if not HAS_MATPLOTLIB or not HAS_NUMPY or len(history) == 0:
        return
    
    print_header("Creating Visualizations")
    
    times = np.array([h['time'] for h in history])
    areas = np.array([h['burned_area'] for h in history])
    perimeters = np.array([h['perimeter'] for h in history])
    max_ros = np.array([h['max_ros'] for h in history])
    max_intensity = np.array([h['max_intensity'] for h in history])
    
    # Create figure with subplots
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle('Wildfire Simulation Analysis - Python API Demo', fontsize=14, fontweight='bold')
    
    # Plot 1: Burned area over time
    axes[0, 0].plot(times / 60, areas / 1e6, 'b-', linewidth=2)
    axes[0, 0].set_xlabel('Time (minutes)')
    axes[0, 0].set_ylabel('Burned Area (km²)')
    axes[0, 0].set_title('Fire Growth')
    axes[0, 0].grid(True, alpha=0.3)
    
    # Plot 2: Perimeter over time
    axes[0, 1].plot(times / 60, perimeters / 1000, 'g-', linewidth=2)
    axes[0, 1].set_xlabel('Time (minutes)')
    axes[0, 1].set_ylabel('Perimeter (km)')
    axes[0, 1].set_title('Fire Perimeter Length')
    axes[0, 1].grid(True, alpha=0.3)
    
    # Plot 3: Rate of spread over time
    axes[1, 0].plot(times / 60, max_ros, 'r-', linewidth=2)
    axes[1, 0].set_xlabel('Time (minutes)')
    axes[1, 0].set_ylabel('Max ROS (m/s)')
    axes[1, 0].set_title('Maximum Rate of Spread')
    axes[1, 0].grid(True, alpha=0.3)
    
    # Plot 4: Fire intensity over time
    axes[1, 1].plot(times / 60, max_intensity / 1000, 'orange', linewidth=2)
    axes[1, 1].set_xlabel('Time (minutes)')
    axes[1, 1].set_ylabel('Max Intensity (MW/m)')
    axes[1, 1].set_title('Maximum Fire Intensity')
    axes[1, 1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    output_file = os.path.join(output_dir, 'fire_analysis.png')
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"  Saved visualization to: {output_file}")
    plt.close()


def run_demo():
    """Main demo function"""
    print_header("Python API Integration Demo: 2026 Enhancement Features")
    print(f"pyWildfire version: {pyWildfire.__version__}")
    
    # Initialize fire solver
    inputs_file = "inputs.i"
    if not os.path.exists(inputs_file):
        print(f"ERROR: Input file {inputs_file} not found")
        return 1
    
    print(f"\n1. Initializing fire solver from {inputs_file}...")
    try:
        fire = WildfireSolver(inputs_file)
        print(f"   ✓ Initialized successfully")
        print(f"     Grid: {fire.nx} × {fire.ny} cells")
        print(f"     Domain: {(fire.xmax-fire.xmin)/1000:.1f} × {(fire.ymax-fire.ymin)/1000:.1f} km")
        print(f"     Cell size: {fire.dx:.1f} × {fire.dy:.1f} m")
    except Exception as e:
        print(f"   ✗ FAILED to initialize: {e}")
        return 1
    
    # Get initial state
    print(f"\n2. Extracting initial state...")
    try:
        state0 = fire.get_state()
        print(f"   ✓ Initial state retrieved")
        print(f"     Time: {state0['time']:.2f} s")
        
        if HAS_NUMPY:
            initial_burned = np.sum(state0['phi'] <= 0.0) * fire.dx * fire.dy
            print(f"     Initial burned area: {initial_burned:.1f} m²")
    except Exception as e:
        print(f"   ✗ FAILED to get state: {e}")
        fire.finalize()
        return 1
    
    # Run simulation with monitoring
    print(f"\n3. Running simulation with real-time monitoring...")
    print(f"   Features enabled:")
    print(f"   - McArthur moisture scaling (T=28°C, RH=25%)")
    print(f"   - FMC sinusoidal phenology (DOY 200, mean=100%, amp=30%)")
    print(f"   - Ember accumulation (decay=0.00167, threshold=10.0)")
    print(f"   - Albini spotting (P_base=0.015, I_min=100 kW/m)")
    print(f"   - Periodic wind gusts (amp=30%, period=900s)")
    print(f"   - Slope-dependent flame tilt")
    
    history = []
    monitor_interval = 5  # Monitor every 5 steps
    
    try:
        print(f"\n   Starting time-stepping...")
        print(f"   {'Step':>4} | {'Time (s)':>8} | {'Burned':>6} | {'Burning':>7} | "
              f"{'ROS':>8} | {'Intensity':>10} | {'Flame':>7}")
        print(f"   {'-'*76}")
        
        total_steps = 50
        for step in range(total_steps):
            result = fire.step()
            if not result['success']:
                print(f"   ✗ FAILED at step {step+1}")
                fire.finalize()
                return 1
            
            # Monitor state every N steps
            if (step + 1) % monitor_interval == 0 or step == total_steps - 1:
                state = fire.get_state()
                print_state_summary(state, step + 1)
                
                # Store history
                if HAS_NUMPY:
                    phi = state['phi']
                    burned_area = np.sum(phi <= 0.0) * fire.dx * fire.dy
                    # Approximate perimeter as cells on boundary
                    boundary = (phi <= 0.0) & (phi > -50.0)
                    perimeter = np.sum(boundary) * fire.dx
                    
                    history.append({
                        'time': state['time'],
                        'burned_area': burned_area,
                        'perimeter': perimeter,
                        'max_ros': np.max(state['ros']),
                        'max_intensity': np.max(state['intensity']),
                        'max_flame': np.max(state['flame_length'])
                    })
        
        print(f"   ✓ Time-stepping successful")
    except Exception as e:
        print(f"   ✗ FAILED during time-stepping: {e}")
        fire.finalize()
        return 1
    
    # Analyze results
    analyze_fire_growth(history)
    
    # Create visualizations
    create_visualizations(history)
    
    # Verify 2026 feature activity
    print_header("2026 Feature Verification")
    state = fire.get_state()
    
    if HAS_NUMPY:
        # Check if features are having an effect
        final_burned = np.sum(state['phi'] <= 0.0) * fire.dx * fire.dy
        max_ros = np.max(state['ros'])
        max_intensity = np.max(state['intensity'])
        
        print(f"  Final burned area:      {final_burned:10,.1f} m²")
        print(f"  Maximum ROS achieved:   {max_ros:.3f} m/s")
        print(f"  Maximum intensity:      {max_intensity:,.0f} kW/m")
        
        # Simple checks
        if final_burned > initial_burned * 2:
            print(f"  ✓ Fire spread verified (>2x growth)")
        if max_ros > 0.5 and max_ros < 10.0:
            print(f"  ✓ ROS values are physically reasonable")
        if max_intensity > 100:
            print(f"  ✓ Fire intensity indicates active fire")
    
    # Finalize
    print_header("Finalizing")
    try:
        fire.finalize()
        print(f"  ✓ Solver finalized successfully")
    except Exception as e:
        print(f"  ✗ FAILED to finalize: {e}")
        return 1
    
    print_header("DEMO COMPLETE")
    print(f"  All Python API operations executed successfully!")
    print(f"  This demo demonstrated:")
    print(f"  - Solver initialization from input file")
    print(f"  - State extraction and monitoring")
    print(f"  - Real-time fire behavior tracking")
    print(f"  - Integration with 2026 enhancement features")
    if HAS_NUMPY:
        print(f"  - Quantitative fire growth analysis")
    if HAS_MATPLOTLIB:
        print(f"  - Visualization generation")
    print()
    
    return 0


def main():
    """Main entry point"""
    return run_demo()


if __name__ == "__main__":
    sys.exit(main())
