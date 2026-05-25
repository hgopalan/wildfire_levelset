#!/usr/bin/env python3
"""
example_integration.py - Example of integrating massconsistent_amr with wildfire_levelset

This script demonstrates how to:
1. Generate 3D wind data (simulating massconsistent_amr output)
2. Load it into wildfire_levelset format using pyWildfire
3. Visualize the column-averaged wind field

In the future, step 1 would be replaced with actual pyAMReX MultiFab from
massconsistent_amr.

Prerequisites:
    pip install numpy matplotlib

Usage:
    PYTHONPATH=build/python python3 src/python/example_integration.py
"""

import sys
import numpy as np
import matplotlib.pyplot as plt

try:
    import pyWildfire
except ImportError as e:
    print(f"Error: Could not import pyWildfire")
    print(f"  {e}")
    print("\nBuild with: cmake -S . -B build -DLEVELSET_BUILD_PYTHON_BINDINGS=ON")
    print("Then set: PYTHONPATH=build/python")
    sys.exit(1)


def generate_synthetic_wind_field(nx, ny, nz, domain):
    """
    Generate a synthetic 3D wind field simulating massconsistent_amr output.
    
    Creates a logarithmic wind profile with spatial variation to simulate
    terrain effects.
    
    Parameters:
    -----------
    nx, ny, nz : int
        Grid dimensions
    domain : dict
        Domain bounds with keys: xmin, xmax, ymin, ymax, zmin, zmax
    
    Returns:
    --------
    u_3d, v_3d, w_3d : ndarray
        3D wind velocity components, shape (nz, ny, nx)
    """
    print("\nGenerating synthetic 3D wind field...")
    print(f"  Grid: {nx} × {ny} × {nz}")
    
    # Create coordinate arrays
    x = np.linspace(domain['xmin'], domain['xmax'], nx)
    y = np.linspace(domain['ymin'], domain['ymax'], ny)
    z = np.linspace(domain['zmin'], domain['zmax'], nz)
    
    # Initialize 3D arrays
    u_3d = np.zeros((nz, ny, nx))
    v_3d = np.zeros((nz, ny, nx))
    w_3d = np.zeros((nz, ny, nx))
    
    # Reference wind at 10m height
    u_ref = 8.0  # m/s westerly
    v_ref = 2.0  # m/s southerly
    z_ref = 10.0  # reference height
    z0 = 0.1  # roughness length
    
    # Create log-law profile with spatial variation
    for k in range(nz):
        z_height = z[k] if z[k] > z0 else z0
        
        # Log-law scaling
        height_factor = np.log(z_height / z0) / np.log(z_ref / z0)
        
        for j in range(ny):
            for i in range(nx):
                # Add spatial variation (simulating terrain effects)
                # Create a gentle gradient and some perturbations
                x_norm = (x[i] - domain['xmin']) / (domain['xmax'] - domain['xmin'])
                y_norm = (y[j] - domain['ymin']) / (domain['ymax'] - domain['ymin'])
                
                # Spatial modulation (represents terrain channeling, etc.)
                spatial_mod = 1.0 + 0.2 * np.sin(3 * np.pi * x_norm) * np.cos(2 * np.pi * y_norm)
                
                u_3d[k, j, i] = u_ref * height_factor * spatial_mod
                v_3d[k, j, i] = v_ref * height_factor * (1.0 + 0.1 * np.cos(2 * np.pi * x_norm))
                w_3d[k, j, i] = 0.1 * np.sin(np.pi * z_height / domain['zmax'])  # Small vertical motion
    
    print(f"  u range: [{u_3d.min():.2f}, {u_3d.max():.2f}] m/s")
    print(f"  v range: [{v_3d.min():.2f}, {v_3d.max():.2f}] m/s")
    print(f"  w range: [{w_3d.min():.2f}, {w_3d.max():.2f}] m/s")
    
    return u_3d, v_3d, w_3d


def load_wind_data(u_3d, v_3d, w_3d, nx, ny, nz, domain):
    """
    Load 3D wind data into wildfire_levelset format using pyWildfire.
    
    Parameters:
    -----------
    u_3d, v_3d, w_3d : ndarray
        3D wind velocity components
    nx, ny, nz : int
        Grid dimensions
    domain : dict
        Domain bounds
    
    Returns:
    --------
    dict
        Result from pyWildfire.load_wind_from_arrays
    """
    print("\nLoading wind data through pyWildfire...")
    
    # Flatten arrays in Fortran (column-major) order
    # This matches the AMReX plotfile format
    u_flat = u_3d.flatten('F')
    v_flat = v_3d.flatten('F')
    w_flat = w_3d.flatten('F')
    
    result = pyWildfire.load_wind_from_arrays(
        nx, ny, nz,
        domain['xmin'], domain['xmax'],
        domain['ymin'], domain['ymax'],
        domain['zmin'], domain['zmax'],
        u_flat, v_flat, w_flat
    )
    
    print(f"  ✓ Loaded {result['n_points']} 3D points")
    print(f"  ✓ Computed {result['nx_2d']} × {result['ny_2d']} 2D wind field")
    
    return result


def visualize_wind_field(result, domain):
    """
    Visualize the column-averaged 2D wind field.
    
    Parameters:
    -----------
    result : dict
        Result from load_wind_from_arrays
    domain : dict
        Domain bounds
    """
    print("\nGenerating visualization...")
    
    u2d = result['u2d']
    v2d = result['v2d']
    nx_2d = result['nx_2d']
    ny_2d = result['ny_2d']
    
    # Create coordinate arrays for plotting
    x = np.linspace(domain['xmin'], domain['xmax'], nx_2d)
    y = np.linspace(domain['ymin'], domain['ymax'], ny_2d)
    X, Y = np.meshgrid(x, y)
    
    # Compute wind speed and direction
    speed = np.sqrt(u2d**2 + v2d**2)
    
    # Create figure
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # Left plot: Wind speed contours with vectors
    ax1 = axes[0]
    contour = ax1.contourf(X, Y, speed, levels=15, cmap='viridis')
    
    # Subsample for quiver plot (every 2nd point)
    skip = 2
    ax1.quiver(X[::skip, ::skip], Y[::skip, ::skip],
               u2d[::skip, ::skip], v2d[::skip, ::skip],
               color='white', alpha=0.7, scale=100)
    
    plt.colorbar(contour, ax=ax1, label='Wind Speed (m/s)')
    ax1.set_xlabel('UTM X (m)')
    ax1.set_ylabel('UTM Y (m)')
    ax1.set_title('Column-Averaged Wind Speed & Direction')
    ax1.grid(True, alpha=0.3)
    
    # Right plot: U and V components
    ax2 = axes[1]
    
    # Plot both components as lines across middle of domain
    mid_y = ny_2d // 2
    ax2.plot(x, u2d[mid_y, :], 'b-', label='u-component', linewidth=2)
    ax2.plot(x, v2d[mid_y, :], 'r-', label='v-component', linewidth=2)
    ax2.plot(x, speed[mid_y, :], 'k--', label='Speed', linewidth=1.5)
    
    ax2.set_xlabel('UTM X (m)')
    ax2.set_ylabel('Wind Velocity (m/s)')
    ax2.set_title(f'Wind Components at Y = {y[mid_y]:.0f} m')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save figure
    output_file = 'wind_field_visualization.png'
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"  ✓ Saved visualization to {output_file}")
    
    plt.show()
    
    # Print statistics
    print("\nWind Field Statistics:")
    print(f"  Mean speed: {speed.mean():.2f} m/s")
    print(f"  Min speed: {speed.min():.2f} m/s")
    print(f"  Max speed: {speed.max():.2f} m/s")
    print(f"  Mean u-component: {u2d.mean():.2f} m/s")
    print(f"  Mean v-component: {v2d.mean():.2f} m/s")


def main():
    """Main integration example"""
    print("=" * 70)
    print("pyWildfire Integration Example")
    print("Simulating massconsistent_amr → wildfire_levelset workflow")
    print("=" * 70)
    print(f"\npyWildfire version: {pyWildfire.__version__}")
    
    # Define domain (UTM Zone 11N, Southern California region)
    domain = {
        'xmin': 329000.0,
        'xmax': 331000.0,
        'ymin': 3774000.0,
        'ymax': 3776000.0,
        'zmin': 0.0,
        'zmax': 100.0  # 100 meters above ground level
    }
    
    # Grid resolution
    nx, ny, nz = 32, 32, 8
    
    # Step 1: Generate synthetic 3D wind (simulates massconsistent_amr output)
    u_3d, v_3d, w_3d = generate_synthetic_wind_field(nx, ny, nz, domain)
    
    # Step 2: Load into wildfire_levelset format
    result = load_wind_data(u_3d, v_3d, w_3d, nx, ny, nz, domain)
    
    # Step 3: Visualize
    try:
        visualize_wind_field(result, domain)
    except ImportError:
        print("\nNote: matplotlib not available, skipping visualization")
        print("Install with: pip install matplotlib")
    
    print("\n" + "=" * 70)
    print("Integration example completed successfully!")
    print("=" * 70)
    print("\nNext steps:")
    print("  1. Replace synthetic wind generation with actual massconsistent_amr")
    print("  2. Use pyAMReX to extract MultiFab data directly")
    print("  3. Pass result to wildfire_levelset simulation")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
