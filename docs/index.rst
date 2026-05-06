.. Wildfire Level-Set Solver documentation master file

Wildfire Level-Set Solver Documentation
========================================

Welcome to the documentation for the Wildfire Level-Set Solver, an AMReX-based C++ level-set solver for wildfire-front style advection based on the community fire model.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   overview
   mathematical_models
   code_structure
   building
   usage
   tools
   wind_solver
   regtests
   comparison
   api_reference

Overview
========

This repository contains a small AMReX-based C++ level-set solver for wildfire-front style advection based on the community fire model. The code implements advanced fire spread models including:

* Rothermel (1972) fire spread equations
* FARSITE elliptical expansion model (Richards 1990)
* Anderson (1983) length-to-width ratio
* Van Wagner (1977) crown fire initiation
* Terrain slope effects
* Firebrand spotting models

Key Features
------------

* **Multiple Fire Spread Models**: Level-set advection (WENO5-Z + RK3) and FARSITE elliptical expansion
* **Multi-Class Rothermel**: Full per-size-class fuel moisture and loading path (1-hr, 10-hr, 100-hr dead; live herbaceous; live woody)
* **Fuel Database**: Anderson 13 (FBFM13) and Scott & Burgan 40 (FBFM40) standard fuel models
* **Terrain Effects**: Slope-based spread rate modifications; terrain and FARSITE landscape file support
* **Crown Fire**: Van Wagner (1977) crown fire initiation criteria
* **Spotting — Stochastic**: Probability-based firebrand generation with lognormal/exponential distance distributions
* **Spotting — Physics-Based**: Albini (1983) thermal-plume lofting with 2-D trajectory integration
* **Bulk Fuel Consumption**: Post-frontal burnout simulation
* **Fire Behavior Diagnostics**: Byram fireline intensity and flame length computed at every time step
* **Time-Dependent Wind**: Sequential wind field snapshots with temporal interpolation (2D)
* **Multiple Ignition Types**: Sphere, box, ellipse, embedded boundary implicit function, CSV fire points
* **2D/3D Support**: Flexible spatial dimensions
* **Embedded Boundaries**: Complex geometry support with AMReX EB

Quick Start
-----------

Clone with submodules::

    git clone --recurse-submodules https://github.com/hgopalan/wildfire_levelset.git
    cd wildfire_levelset

Build the project::

    cmake -S . -B build
    cmake --build build -j

Run a test case::

    cd regtest/basic_levelset
    ../../build/levelset inputs.i

Indices and tables
==================

* :ref:`genindex`
* :ref:`search`
