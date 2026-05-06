Comparison with Other Wildfire Simulation Tools
===============================================

.. note::

   **Note on naming**: WRF-Fire and WRF-SFIRE refer to the same coupled
   fire–atmosphere system built on the WRF mesoscale model. The official
   package name is WRF-SFIRE; "WRF-Fire" is a common shorthand used in the
   community.

.. warning::

   **Disclaimer**: The comparisons below are based on publicly available
   documentation and peer-reviewed literature as of 2025. Capabilities of
   third-party tools may have evolved or may differ across versions. Consult
   the official references linked in the `Tool Documentation References`_
   section for authoritative and up-to-date information. This table is
   intended for high-level orientation only and should not be taken as an
   exhaustive or definitive characterisation of any tool.

Capability Summary
------------------

.. list-table::
   :header-rows: 1
   :widths: 20 20 10 10 10 10 10 10 10

   * - Capability
     - wildfire_levelset
     - FARSITE
     - WRF-SFIRE
     - FlamMap
     - BehavePlus
     - FSPRO
     - QUIC-Fire
     - FIRETEC
   * - **Fire spread model**
     - Rothermel (1972) or Balbi (2009) ROS; Cheney & Gould (1995) grassland;
       elliptical directional spread (Richards 1990); Eulerian level-set
     - Rothermel ROS; explicit Huygens wavelet propagation
     - Rothermel ROS (level-set)
     - Rothermel ROS; minimum travel time (MTT) or Huygens wavelet
     - Rothermel (1972) point-based ROS; no spatial propagation
     - Ensemble of FARSITE simulations; probabilistic spread
     - Coupled semi-empirical fire model with QUIC-URB 3-D wind solver
     - High-fidelity LES with coupled combustion physics; no empirical ROS
   * - **Wind adjustment**
     - Optional WAF (20-ft → midflame) and MEWS cap; 7 wind-terrain feedback
       models (Rothermel 1983 canyon, Viegas & Neto 1994, Pimont et al. 2009,
       WindNinja ridge/canyon)
     - WAF applied internally; MEWS cap
     - Wind derived from WRF atmospheric model; WAF handled in coupling
     - WindNinja or constant/gridded input; WAF optional
     - User-specified WAF option
     - Applied through each FARSITE instance
     - 3-D mass-consistent wind field computed by QUIC-URB
     - Wind computed by Navier-Stokes LES solver
   * - **Terrain representation**
     - Per-cell elevation, slope, aspect, and fuel from ``.lcp`` or XYZ; 2-D;
       3-D terrain and canopy layer can be added in future
     - Full 2-D landscape (``.lcp``)
     - Full 3-D terrain from WRF grid
     - Full 2-D landscape (``.lcp``)
     - Point-based (user-specified slope/aspect); no raster landscape
     - Full 2-D landscape (``.lcp``)
     - Full 3-D terrain and vegetation canopy structure
     - Full 3-D terrain and vegetation canopy
   * - **Fire-atmosphere coupling**
     - None (prescribed wind fields)
     - None (prescribed gridded wind)
     - Two-way (fire ↔ WRF atmosphere, heat/momentum flux)
     - None (prescribed wind)
     - None
     - None (prescribed weather)
     - One-way/quasi-coupled (QUIC-URB wind; no fire-driven feedback)
     - Two-way (fully coupled LES fire–atmosphere)
   * - **Atmospheric model**
     - None
     - None
     - Full WRF mesoscale atmosphere (NWP)
     - None
     - None
     - None
     - Fast-response QUIC-URB 3-D mass-consistent wind model
     - Full 3-D Navier-Stokes LES atmosphere
   * - **Spatial resolution**
     - User-defined (AMReX grid); adaptive mesh refinement (AMR)
     - Landscape raster resolution
     - WRF grid resolution (typically 100 m–1 km)
     - Landscape raster resolution
     - Point-based (no spatial grid)
     - Landscape raster resolution (via FARSITE)
     - User-defined 3-D grid (typically 1–10 m)
     - User-defined 3-D grid (typically 1–10 m)
   * - **Fuel models**
     - Anderson 13 + Scott & Burgan 40
     - Anderson 13 + Scott & Burgan 40
     - Anderson 13
     - Anderson 13 + Scott & Burgan 40
     - Anderson 13 + Scott & Burgan 40
     - Anderson 13 + Scott & Burgan 40
     - Custom 3-D bulk fuel density and arrangement
     - Continuous 3-D fuel density distribution
   * - **Fuel moisture**
     - Per size class (1-hr, 10-hr, 100-hr dead; live herbaceous and live
       woody); time-varying ``.fmd`` schedule
     - Dead/live moisture per size class
     - Dead/live moisture (prescribed)
     - Dead/live per size class; fuel moisture conditioning option
     - Dead/live per size class
     - Dead/live per size class (from historical RAWS)
     - User-specified bulk moisture per cell
     - User-specified moisture content per cell
   * - **Crown fire**
     - Van Wagner (1977) initiation
     - Van Wagner (1977) initiation
     - Van Wagner (1977) initiation
     - Scott & Reinhardt (2001) crown fire assessment
     - Van Wagner (1977)
     - Via FARSITE (Van Wagner)
     - Physics-based ignition via heat flux
     - Governed by physical combustion model
   * - **Firebrand spotting**
     - Stochastic + Albini (1983) 2-D trajectory model
     - Albini empirical spotting
     - Not included by default
     - Not standard
     - Albini empirical (point calculation)
     - Via FARSITE (Albini)
     - Not standard
     - Physics-based firebrand transport (select versions)
   * - **Bulk fuel burnout**
     - Residence-time model
     - Per-cell burnout tracking
     - Post-frontal fuel consumption
     - Not tracked separately
     - Not applicable
     - Via FARSITE
     - Physics-based fuel consumption
     - Physical combustion with full burnout
   * - **Multi-ignition**
     - CSV ignition file
     - Interactive ignition map
     - WPS/WRF fire restart
     - Ignition map or polygon
     - Single point or simple polygon
     - Multiple ignition points
     - User-defined ignition locations
     - User-defined ignition locations
   * - **GPU acceleration**
     - Yes (AMReX CUDA/HIP/SYCL kernels)
     - No
     - No
     - No
     - No
     - No
     - Partial (QUIC-URB GPU support in select versions)
     - No (primarily MPI CPU)
   * - **MPI parallelism**
     - Yes (AMReX domain decomposition)
     - No
     - Yes (WRF MPI)
     - No
     - No
     - No (ensemble via independent FARSITE instances)
     - Yes (QUIC MPI/OpenMP)
     - Yes (HPC clusters, MPI)
   * - **Embedded boundaries / obstacles**
     - Yes (AMReX EB for buildings, fuel breaks)
     - No
     - No
     - No
     - No
     - No
     - Yes (QUIC-URB handles buildings and obstacles)
     - Yes (complex 3-D geometry)

Tool Documentation References
------------------------------

For authoritative and up-to-date information on each tool, refer to the
official sources:

* **wildfire_levelset**: https://hgopalan.github.io/wildfire_levelset/
* **FARSITE**: https://www.firelab.org/project/farsite
* **WRF-SFIRE**: https://github.com/openwfm/WRF-SFIRE
* **FlamMap**: https://www.firelab.org/project/flammap
* **BehavePlus**: https://www.firelab.org/project/behaveplusfiremodeling
* **FSPRO**: https://www.firelab.org/project/fspro
* **QUIC-Fire**: https://www.lanl.gov/projects/quic-fire/
* **FIRETEC**: https://www.lanl.gov/org/padwp/adcles/fluid-dynamics-solid-mechanics/index.php
  (see also Linn et al. 2002)

Key Differences from FARSITE
-----------------------------

* **Eulerian level-set vs. explicit Huygens wavelets**: This solver implements
  the same elliptical directional spread (Richards 1990) as FARSITE, but
  embeds it in an Eulerian level-set framework. The fire perimeter is tracked
  as the zero contour of a smooth signed-distance function rather than a chain
  of individually advected vertex wavelets. The Eulerian approach handles
  complex topology changes (merging fronts, islands) automatically without
  explicit connectivity management.

* **Andrews (2018) wind adjustments are optional**: FARSITE applies the Wind
  Adjustment Factor (WAF) and maximum effective wind speed (MEWS) cap
  internally. This solver exposes both via ``rothermel.use_waf`` and
  ``rothermel.use_wind_limit`` so users can enable them when wind input is at
  20-ft (NWP or WRF) height rather than midflame height.

* **Full 2-D landscape terrain**: Both solvers support per-cell elevation,
  slope, aspect, and fuel model from FARSITE ``.lcp`` landscape files
  (LANDFIRE data). Canopy fuel layers can be added to this solver in future
  without architectural changes.

* **Per size-class fuel moisture**: Dead fuel moisture is specified
  independently for 1-hr, 10-hr, and 100-hr size classes; live herbaceous and
  live woody moisture are specified separately. This mirrors FARSITE's moisture
  inputs and drives the multi-class Rothermel (1972) reaction intensity
  calculation.

* **Embedded Boundary**: AMReX EB allows irregular obstacles (buildings, fuel
  breaks) on the Cartesian grid without remeshing.

* **GPU and MPI**: The AMReX backend supports CUDA/HIP/SYCL and MPI for large
  domains; FARSITE is serial and CPU-only.

Key Differences from WRF-Fire (WRF-SFIRE)
------------------------------------------

* **No atmospheric coupling**: This solver uses prescribed wind fields
  (constant, CSV, or WRF output); WRF-SFIRE fully couples the fire with the
  WRF atmospheric model, including fire-induced wind, heat flux, and smoke
  transport.

* **Wind adjustment control**: When using WRF output wind via
  ``wrf_wind_reader.py``, enable ``rothermel.use_waf = 1`` to convert the NWP
  wind to midflame height. WRF-SFIRE handles this internally as part of the
  coupled framework.

* **Albini spotting vs. none**: WRF-SFIRE does not include firebrand spotting
  by default.

* **Simpler setup**: Running this solver requires only CMake and an AMReX
  build; WRF-SFIRE requires a full WRF stack (NetCDF, MPI, WPS, WRF
  pre-processing).

* **GPU-native kernels**: WRF-SFIRE is primarily CPU-MPI; this solver's
  hotspot loops use AMReX GPU kernels.
