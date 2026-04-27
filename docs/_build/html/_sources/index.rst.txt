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

* **Multiple Fire Spread Models**: Level-set advection and FARSITE elliptical expansion
* **Terrain Effects**: Slope-based spread rate modifications using Rothermel corrections
* **Crown Fire**: Van Wagner crown fire initiation criteria
* **Spotting**: Firebrand spotting probability models
* **Bulk Fuel Consumption**: Post-frontal burnout simulation
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
    ../../build/levelset inputs

Indices and tables
==================

* :ref:`genindex`
* :ref:`search`
