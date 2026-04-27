Overview
========

The Wildfire Level-Set Solver is a computational tool for simulating wildfire spread using level-set methods and the FARSITE fire growth model. The code is built on top of AMReX, a software framework for block-structured adaptive mesh refinement.

Summary Flow
------------

The wildfire simulation follows these key steps for each timestep:

1. **Setup inputs** (landscape, fuel, weather, wind)
   
   Parse configuration parameters for terrain, fuel properties, weather conditions, and wind field.

2. **Compute surface ROS via Rothermel/Level Set**
   
   Calculate rate of spread using Rothermel fire spread equations with terrain and wind corrections.

3. **Generate elliptical wavelets per vertex**
   
   Create Huygens wavelets with elliptical shapes at each fire front vertex based on local spread rates.

4. **Merge to new perimeter**
   
   Combine wavelets to form new fire perimeter using level set or FARSITE methods.

5. **Apply crown/spotting sub-models**
   
   Evaluate crown fire initiation criteria and firebrand spotting probability to generate new ignition points.

6. **Simulate post-frontal burnout**
   
   Compute bulk fuel consumption fraction for areas behind the fire front.

7. **Update states, record outputs, step time**
   
   Save plotfiles with fire state, advance simulation time, update data structures.

Fire Spread Models
------------------

Anderson L/W Ratio (FARSITE Model)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The Length-to-Width ratio calculation is based on Anderson (1983):

.. math::

   L/W = 0.936 \exp(0.2566 U) + 0.461 \exp(-0.1548 U) - 0.397

where :math:`U` is the wind speed. This accounts for elliptical fire shape under wind influence and is used when ``use_farsite_model=true``.

Rothermel Terrain Effects
^^^^^^^^^^^^^^^^^^^^^^^^^^

Slope correction factor based on Rothermel (1972) accounts for uphill/downhill fire spread acceleration/deceleration:

.. math::

   \phi_s = 5.275 \beta^{-0.3} \tan^2(\theta)

where :math:`\beta` is the packing ratio and :math:`\theta` is the slope angle. This is used when ``use_terrain_effects=true`` and ``use_farsite_model=false``.

FARSITE Combined Wind and Terrain Model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Vectorial combination of wind and slope effects accounts for alignment between wind direction and slope aspect with enhanced rate of spread calculations. Used when both ``use_terrain_effects=true`` and ``use_farsite_model=true``.

Prerequisites
-------------

* C++17 compiler
* CMake (3.20+)
* Git

The project supports both 2D and 3D configurations. The default is 3D, but you can build for 2D using CMake options.
