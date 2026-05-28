.. _python_api:

Python API for Coupled Simulations
===================================

The Python API (``pyWildfire``) provides complete programmatic control of the wildfire solver from Python, enabling coupled wind-fire simulations with external wind solvers.

Overview
--------

The Python bindings expose the fire solver's core functionality through two interfaces:

1. **Low-level API** (``pyWildfire`` module): Direct C++ function bindings
2. **High-level API** (``WildfireSolver`` class): Object-oriented Python wrapper

Key capabilities:

* Initialize fire solver from inputs file
* Time-step the fire simulation
* Extract all fire fields as NumPy arrays (phi, ROS, intensity, flame length, etc.)
* Update wind fields from 2D or 3D arrays
* Write AMReX plotfiles
* Zero-copy data transfer between C++ and Python

Building with Python Bindings
------------------------------

Enable Python bindings during CMake configuration:

.. code-block:: bash

   cmake -S . -B build \
     -DLEVELSET_DIM_2D=ON \
     -DLEVELSET_BUILD_PYTHON_BINDINGS=ON
   cmake --build build -j

Set the Python path:

.. code-block:: bash

   export PYTHONPATH=$PWD/build/python:$PYTHONPATH

Requirements:

* Python 3.6 or later
* NumPy
* pybind11 (included as submodule)

Quick Start
-----------

Basic fire simulation:

.. code-block:: python

   from wildfire_solver import WildfireSolver
   
   # Initialize
   fire = WildfireSolver("inputs.i")
   
   # Run simulation
   for i in range(100):
       fire.step()
       state = fire.get_state()
       burned_area = (state['phi'] <= 0).sum() * fire.dx * fire.dy
       print(f"t={state['time']:.1f}s, burned={burned_area:.0f}m²")
   
   # Finalize
   fire.finalize()

High-Level API (WildfireSolver)
--------------------------------

The ``WildfireSolver`` class provides a Pythonic interface to the fire solver.

Initialization
~~~~~~~~~~~~~~

.. code-block:: python

   from wildfire_solver import WildfireSolver
   
   fire = WildfireSolver("inputs.i")

The constructor initializes AMReX and the fire solver from an inputs file.

Properties
~~~~~~~~~~

Access solver properties:

.. code-block:: python

   nx, ny = fire.nx, fire.ny              # Grid dimensions
   xmin, xmax = fire.xmin, fire.xmax      # Domain bounds (x)
   ymin, ymax = fire.ymin, fire.ymax      # Domain bounds (y)
   dx, dy = fire.dx, fire.dy              # Cell spacing
   t = fire.time                           # Current simulation time

Time-Stepping
~~~~~~~~~~~~~

Advance the simulation by one timestep:

.. code-block:: python

   result = fire.step()
   # result = {'success': True, 'dt': 0.123, 'time': 1.234}

State Extraction
~~~~~~~~~~~~~~~~

Extract all fire fields as NumPy arrays:

.. code-block:: python

   state = fire.get_state()
   
   # Available fields (all shape (ny, nx)):
   phi = state['phi']                   # Level set (m)
   ros = state['ros']                   # Rate of spread (m/s)
   intensity = state['intensity']       # Fire line intensity (kW/m)
   flame_length = state['flame_length'] # Flame length (m)
   u_wind = state['u_wind']             # Wind u-component (m/s)
   v_wind = state['v_wind']             # Wind v-component (m/s)
   time = state['time']                 # Simulation time (s)

Wind Updates
~~~~~~~~~~~~

Update wind from 2D arrays:

.. code-block:: python

   import numpy as np
   
   # Create new wind field (shape: (ny, nx))
   u_new = np.full((fire.ny, fire.nx), 10.0)  # 10 m/s easterly
   v_new = np.zeros((fire.ny, fire.nx))
   
   fire.update_wind(u_new, v_new)

Update wind from 3D arrays (for coupled simulations):

.. code-block:: python

   # 3D wind field from external solver
   nz = 10
   zmin, zmax = 0.0, 100.0
   
   u_3d = np.zeros((nz, fire.ny, fire.nx))  # Shape: (nz, ny, nx)
   v_3d = np.zeros((nz, fire.ny, fire.nx))
   w_3d = np.zeros((nz, fire.ny, fire.nx))
   
   # Fill with wind data...
   
   fire.update_wind_3d(u_3d, v_3d, w_3d, nz, zmin, zmax)

Plotfile Writing
~~~~~~~~~~~~~~~~

Write AMReX plotfile:

.. code-block:: python

   fire.write_plotfile()

Finalization
~~~~~~~~~~~~

Clean up and finalize:

.. code-block:: python

   fire.finalize()

Or use as a context manager:

.. code-block:: python

   with WildfireSolver("inputs.i") as fire:
       while fire.time < final_time:
           fire.step()

Low-Level API (pyWildfire)
---------------------------

The ``pyWildfire`` module provides direct access to C++ functions.

.. code-block:: python

   import pyWildfire
   
   # Initialize
   result = pyWildfire.initialize("inputs.i")
   # result = {'success': True, 'nx': 64, 'ny': 64, ...}
   
   # Advance
   step_result = pyWildfire.advance()
   # step_result = {'success': True, 'dt': 0.123, 'time': 1.234}
   
   # Get state
   state = pyWildfire.get_state()
   
   # Update wind (2D)
   pyWildfire.update_wind(u_2d, v_2d)
   
   # Update wind (3D)
   pyWildfire.update_wind_3d(u_3d, v_3d, w_3d, nz, zmin, zmax)
   
   # Write plotfile
   pyWildfire.write_plotfile()
   
   # Finalize
   pyWildfire.finalize()
   
   # Check initialization status
   is_init = pyWildfire.is_initialized()

Coupled Wind-Fire Simulations
------------------------------

Integration with massconsistent_amr
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The primary use case for the Python API is coupling with external wind solvers like
`massconsistent_amr <https://github.com/hgopalan/massconsistent_amr>`_.

.. code-block:: python

   from wildfire_solver import WildfireSolver
   from pyWindSolver import WindSolver  # massconsistent_amr module
   
   # Initialize both solvers
   fire = WildfireSolver("fire_inputs.i")
   wind = WindSolver("wind_inputs.txt")
   
   # Coupled time loop
   final_time = 3600.0  # 1 hour
   
   while fire.time < final_time:
       # 1. Solve wind field
       wind.solve(fire.time)
       u_3d, v_3d, w_3d = wind.get_velocity_arrays()
       
       # 2. Update fire wind
       fire.update_wind_3d(u_3d, v_3d, w_3d, wind.nz, wind.zmin, wind.zmax)
       
       # 3. Advance fire
       fire.step()
       
       # 4. Optional: Extract state for analysis
       state = fire.get_state()
       burned = (state['phi'] <= 0).sum() * fire.dx * fire.dy
       print(f"t={fire.time:.1f}s, burned={burned:.0f}m²")
   
   # Finalize
   fire.finalize()
   wind.finalize()

One-Way vs Two-Way Coupling
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**One-way coupling** (wind → fire):

* Wind solver runs independently
* Fire receives wind but doesn't affect it
* Simpler, faster, suitable for most applications
* Current implementation

**Two-way coupling** (wind ↔ fire):

* Fire heat release affects wind (buoyancy, updrafts)
* More physically realistic
* Computationally expensive
* Future enhancement

Example two-way coupling (conceptual):

.. code-block:: python

   while fire.time < final_time:
       # Wind → Fire
       wind.solve(fire.time)
       u_3d, v_3d, w_3d = wind.get_velocity_arrays()
       fire.update_wind_3d(u_3d, v_3d, w_3d, wind.nz, wind.zmin, wind.zmax)
       
       # Advance fire
       fire.step()
       
       # Fire → Wind (future feature)
       state = fire.get_state()
       heat_release = compute_heat_release(state)
       wind.add_heat_source(heat_release)  # Affects next wind solve

Replacing the Wind Solver
--------------------------

The Python API is designed to work with any wind solver that can provide 3D velocity fields.

Wind Solver Requirements
~~~~~~~~~~~~~~~~~~~~~~~~~

Your wind solver must provide:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Requirement
     - Description
   * - **Array shape**
     - ``(nz, ny, nx)`` where ``nz`` = vertical levels, ``ny`` = latitude, ``nx`` = longitude
   * - **Array order**
     - Fortran (column-major) order
   * - **Units**
     - Velocities in m/s
   * - **Coordinate system**
     - Same as fire solver (typically UTM)
   * - **Domain overlap**
     - Wind domain must cover fire domain
   * - **Vertical extent**
     - ``zmin`` to ``zmax`` in meters above ground level

Implementation Patterns
~~~~~~~~~~~~~~~~~~~~~~~

**Option 1: AMReX-based solver with Python bindings**

If your wind solver is AMReX-based, follow the massconsistent_amr pattern:

1. Add pybind11 bindings to expose solver functions
2. Implement ``get_velocity_arrays()`` returning NumPy arrays
3. Match the interface shown in the coupled simulation example

**Option 2: External executable (WindNinja, etc.)**

Wrap external solvers with Python:

.. code-block:: python

   import subprocess
   import numpy as np
   
   def run_windninja(fire_domain, time):
       """Run WindNinja and load results"""
       # Call WindNinja
       subprocess.run([
           'WindNinja_cli',
           '--domain', f'{fire_domain.xmin},{fire_domain.xmax},...',
           '--time', str(time),
           '--output', 'wind_output.asc'
       ])
       
       # Read ASCII grid output
       u_3d, v_3d, w_3d = read_windninja_output('wind_output.asc')
       
       return u_3d, v_3d, w_3d
   
   # Use in simulation
   fire = WildfireSolver("inputs.i")
   while fire.time < final_time:
       u_3d, v_3d, w_3d = run_windninja(fire, fire.time)
       fire.update_wind_3d(u_3d, v_3d, w_3d, nz, zmin, zmax)
       fire.step()

**Option 3: WRF or other NetCDF-based models**

Extract wind from model output files:

.. code-block:: python

   from netCDF4 import Dataset
   import wrf  # wrf-python package
   
   def extract_wrf_wind(wrfout_file, time_idx, fire_grid):
       """Extract and interpolate WRF wind to fire grid"""
       ncfile = Dataset(wrfout_file)
       
       # Extract wind variables
       u = wrf.getvar(ncfile, 'ua', timeidx=time_idx)
       v = wrf.getvar(ncfile, 'va', timeidx=time_idx)
       w = wrf.getvar(ncfile, 'wa', timeidx=time_idx)
       
       # Interpolate to fire grid (use scipy.interpolate or similar)
       u_interp = interpolate_to_grid(u, fire_grid)
       v_interp = interpolate_to_grid(v, fire_grid)
       w_interp = interpolate_to_grid(w, fire_grid)
       
       return u_interp, v_interp, w_interp

**Option 4: Custom Python wind solver**

Implement your own wind solver in Python:

.. code-block:: python

   def solve_wind_field(fire, time):
       """Simple log-law wind profile"""
       nz = 10
       zmin, zmax = 0.0, 100.0
       
       u_ref, v_ref = 5.0, 1.0  # Reference wind at 10m
       z_ref, z0 = 10.0, 0.1     # Reference height, roughness
       
       u_3d = np.zeros((nz, fire.ny, fire.nx))
       v_3d = np.zeros((nz, fire.ny, fire.nx))
       w_3d = np.zeros((nz, fire.ny, fire.nx))
       
       for k in range(nz):
           z = zmin + (k + 0.5) * (zmax - zmin) / nz
           factor = np.log(z / z0) / np.log(z_ref / z0)
           u_3d[k, :, :] = u_ref * factor
           v_3d[k, :, :] = v_ref * factor
       
       return u_3d, v_3d, w_3d, nz, zmin, zmax

Example: Creating a Wind Solver Interface
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Template for wrapping any wind solver:

.. code-block:: python

   class CustomWindSolver:
       """Interface template for wind solvers"""
       
       def __init__(self, inputs):
           self.nz = 10
           self.zmin = 0.0
           self.zmax = 100.0
           # Initialize your wind solver
       
       def solve(self, time):
           """Solve wind field at given time"""
           # Run your wind solver
           pass
       
       def get_velocity_arrays(self):
           """Return wind as NumPy arrays"""
           # Must return (u_3d, v_3d, w_3d) with shape (nz, ny, nx)
           # in Fortran order
           return u_3d, v_3d, w_3d
       
       def finalize(self):
           """Clean up"""
           pass

Use with fire solver:

.. code-block:: python

   fire = WildfireSolver("fire_inputs.i")
   wind = CustomWindSolver("wind_inputs.txt")
   
   while fire.time < final_time:
       wind.solve(fire.time)
       u_3d, v_3d, w_3d = wind.get_velocity_arrays()
       fire.update_wind_3d(u_3d, v_3d, w_3d, wind.nz, wind.zmin, wind.zmax)
       fire.step()

Current Limitations
-------------------

Known limitations of the Python API:

1. **Column-Averaging**
   
   The fire solver is 2D (horizontal). 3D wind fields are column-averaged before use.
   Planned: Height-dependent wind influence on fire spread.

2. **One-Way Coupling Only**
   
   Fire heat release does not currently affect wind solver.
   Planned: Two-way coupling with buoyancy feedback.

3. **Domain Matching**
   
   Wind and fire domains must overlap; the user must ensure spatial consistency.
   Planned: Automatic domain intersection and interpolation.

4. **Time Synchronization**
   
   User must manually synchronize wind and fire solver times.
   Planned: Built-in time coordination.

5. **MPI Support**
   
   MPI-parallel simulations not fully tested with Python bindings.
   Planned: Complete MPI support in Python API.

6. **GPU Data Transfer**
   
   GPU builds work but data transfer is not optimized.
   Planned: GPU-aware bindings with zero-copy when possible.

7. **AMR Not Exposed**
   
   Adaptive mesh refinement is not accessible via Python.
   Planned: AMR control from Python.

Performance Considerations
--------------------------

Data Transfer Overhead
~~~~~~~~~~~~~~~~~~~~~~

Copying data between C++ and Python has overhead. Best practices:

* **Minimize state extractions**: Only call ``get_state()`` when needed
* **Use appropriate update frequency**: Don't update wind every timestep if not necessary
* **Batch operations**: Group multiple timesteps between wind updates

.. code-block:: python

   # Good: Update wind periodically
   wind_update_interval = 10  # Update every 10 fire steps
   
   for i in range(1000):
       if i % wind_update_interval == 0:
           u_3d, v_3d, w_3d = wind_solver.solve(fire.time)
           fire.update_wind_3d(u_3d, v_3d, w_3d, nz, zmin, zmax)
       fire.step()

Sub-Cycling
~~~~~~~~~~~

Wind and fire solvers may have different optimal timesteps:

.. code-block:: python

   # Sub-cycling example
   dt_wind = 1.0   # Wind solver timestep (s)
   dt_fire = 0.1   # Fire solver timestep (s)
   
   t = 0.0
   while t < final_time:
       # Update wind
       wind.solve(t)
       u_3d, v_3d, w_3d = wind.get_velocity_arrays()
       fire.update_wind_3d(u_3d, v_3d, w_3d, nz, zmin, zmax)
       
       # Advance fire multiple steps
       for _ in range(int(dt_wind / dt_fire)):
           fire.step()
       
       t += dt_wind

Memory Management
~~~~~~~~~~~~~~~~~

For large simulations, be aware of memory usage:

.. code-block:: python

   # Clear state arrays when done
   state = fire.get_state()
   # Use state...
   del state  # Free memory
   
   # Or extract only needed fields
   phi = fire.get_state()['phi']
   # Don't hold reference to full state dict

Applications
------------

The Python API enables various advanced workflows:

Ensemble Simulations
~~~~~~~~~~~~~~~~~~~~

Run multiple fire scenarios with different wind conditions:

.. code-block:: python

   import numpy as np
   from wildfire_solver import WildfireSolver
   
   wind_speeds = [3, 5, 7, 10]  # m/s
   results = []
   
   for u_wind in wind_speeds:
       fire = WildfireSolver("inputs.i")
       
       # Set constant wind
       u_2d = np.full((fire.ny, fire.nx), u_wind)
       v_2d = np.zeros((fire.ny, fire.nx))
       fire.update_wind(u_2d, v_2d)
       
       # Run simulation
       while fire.time < 3600.0:
           fire.step()
       
       # Record results
       state = fire.get_state()
       burned = (state['phi'] <= 0).sum() * fire.dx * fire.dy
       results.append({'wind': u_wind, 'burned': burned})
       
       fire.finalize()
   
   print("Ensemble results:", results)

Machine Learning Training
~~~~~~~~~~~~~~~~~~~~~~~~~

Generate training data for ML models:

.. code-block:: python

   import h5py
   from wildfire_solver import WildfireSolver
   
   # Generate training dataset
   with h5py.File('training_data.h5', 'w') as f:
       fire = WildfireSolver("inputs.i")
       
       for i in range(100):
           fire.step()
           state = fire.get_state()
           
           # Save snapshot
           grp = f.create_group(f'step_{i:04d}')
           grp['phi'] = state['phi']
           grp['ros'] = state['ros']
           grp['wind_u'] = state['u_wind']
           grp['wind_v'] = state['v_wind']
           grp['time'] = state['time']
       
       fire.finalize()

Custom Analysis
~~~~~~~~~~~~~~~

Implement custom fire behavior analysis:

.. code-block:: python

   from wildfire_solver import WildfireSolver
   import matplotlib.pyplot as plt
   
   fire = WildfireSolver("inputs.i")
   
   times = []
   areas = []
   max_ros = []
   
   while fire.time < final_time:
       fire.step()
       state = fire.get_state()
       
       # Track metrics
       times.append(state['time'])
       areas.append((state['phi'] <= 0).sum() * fire.dx * fire.dy)
       max_ros.append(state['ros'].max())
   
   fire.finalize()
   
   # Plot results
   fig, (ax1, ax2) = plt.subplots(2, 1)
   ax1.plot(times, areas)
   ax1.set_ylabel('Burned Area (m²)')
   ax2.plot(times, max_ros)
   ax2.set_ylabel('Max ROS (m/s)')
   ax2.set_xlabel('Time (s)')
   plt.savefig('fire_metrics.png')

Testing and Validation
-----------------------

Regression tests for the Python API are in ``regtest/python_api/``:

* ``basic_fire_solver/`` - Basic API functionality
* ``coupled_wind_fire/`` - Coupled simulation with synthetic wind

Run tests:

.. code-block:: bash

   cd build
   ctest -L python_api --output-on-failure

See Also
--------

* :ref:`building` - Build instructions
* :ref:`usage` - Input parameters
* :ref:`tools` - Python analysis tools
* `PYTHON_API_IMPLEMENTATION.md <https://github.com/hgopalan/wildfire_levelset/blob/main/PYTHON_API_IMPLEMENTATION.md>`_ - Implementation details
* `massconsistent_amr <https://github.com/hgopalan/massconsistent_amr>`_ - Wind solver with Python bindings
