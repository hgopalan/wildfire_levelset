Overview
========

The Wildfire Level-Set Solver is a computational tool for simulating wildfire spread using level-set methods and the FARSITE fire growth model. The code is built on top of AMReX, a software framework for block-structured adaptive mesh refinement.

Summary Flow
------------

The wildfire simulation follows these key steps for each timestep:

1. **Setup inputs** (landscape, fuel, weather, wind)

   Parse configuration parameters for terrain, fuel properties, weather conditions, and wind field.

2. **Compute surface ROS via Rothermel/Level Set**

   Calculate rate of spread using Rothermel fire spread equations with terrain and wind corrections. When a landscape file is provided, per-cell fuel models are used from a pre-built lookup table.

3. **Generate elliptical wavelets per vertex**

   Create Huygens wavelets with elliptical shapes at each fire front vertex based on local spread rates.

4. **Merge to new perimeter**

   Combine wavelets to form new fire perimeter using level set or FARSITE methods.

5. **Apply crown/spotting sub-models**

   Evaluate crown fire initiation criteria and firebrand spotting (stochastic and/or physics-based Albini) to generate new ignition points.

6. **Simulate post-frontal burnout**

   Compute bulk fuel consumption fraction for areas behind the fire front.

7. **Update states, record outputs, step time**

   Save plotfiles with fire state (including fireline intensity and flame length diagnostics), advance simulation time, update data structures.

Prerequisites
-------------

* C++17 compiler
* CMake (3.20+)
* Git

The project supports both 2D and 3D configurations. The default is 3D, but you can build for 2D using CMake options. Time-dependent wind fields are only available in 2D builds.
