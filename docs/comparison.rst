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
   :widths: 22 26 10 10 10 10 12

   * - Capability
     - **Wildfire-AMR** (this solver)
     - FARSITE
     - WRF-SFIRE
     - FlamMap
     - BehavePlus
     - QUIC-Fire / FIRETEC
   * - **Surface spread models**
     - Rothermel (1972); Balbi (2009); Cheney–Gould (1995) grassland;
       Cruz–Alexander–Wakimoto (2005) crown; FBP O1a/O1b/S1–S3 (Forestry
       Canada 1992); Lautenberger (2013) physics-based
     - Rothermel (1972)
     - Rothermel (1972)
     - Rothermel (1972)
     - Rothermel (1972)
     - Semi-empirical (QUIC); physics LES (FIRETEC)
   * - **Propagation method**
     - Eulerian level-set (WENO5-Z/RK3); FARSITE Huygens ellipse (Richards
       1990); Minimum Travel Time (Finney 2002)
     - Huygens wavelet (explicit vertex propagation)
     - Eulerian level-set
     - MTT or Huygens wavelet
     - Point-based; no spatial propagation
     - Eulerian CFD grid
   * - **Crown fire**
     - Van Wagner (1977) initiation + passive CF blend;
       Rothermel (1991) active ROS (3.34 × R\ :sub:`surface`);
       Cruz et al. (2005) active ROS;
       Scott & Reinhardt (2001) full TI/CI (bisection)
     - Van Wagner (1977); active crown via R\ :sub:`c` = 3.34 R\ :sub:`s`
     - Van Wagner (1977)
     - Scott & Reinhardt (2001) full assessment
     - Van Wagner (1977)
     - Physics-based (heat flux / combustion)
   * - **Wind adjustment**
     - WAF (Andrews 2018, 20-ft → midflame); MEWS cap; 7 wind-terrain
       feedback models (canyon, Viegas, Pimont, WindNinja ridge/canyon)
     - WAF + MEWS cap (internal)
     - WRF-derived; WAF in coupling layer
     - WindNinja or gridded; WAF optional
     - User-specified WAF
     - 3-D mass-consistent (QUIC-URB) or LES
   * - **Fuel models**
     - FBFM13 (Anderson 13); FBFM40 (Scott–Burgan 40); custom inline;
       FBP O1a/O1b grass + S1–S3 slash; per-cell from LANDFIRE LCP
     - FBFM13 + FBFM40
     - FBFM13
     - FBFM13 + FBFM40
     - FBFM13 + FBFM40
     - Custom 3-D bulk density
   * - **Fuel moisture**
     - 1-hr/10-hr/100-hr dead + live herb/woody per cell; FMD schedule;
       Nelson (2000) diurnal EMC; precipitation wetting;
       per-cell from FARSITE .fms; solar shade adjustment
     - Dead/live per size class; FMD schedule
     - Dead/live (prescribed)
     - Dead/live per class; conditioning option
     - Dead/live per class
     - Bulk moisture per cell
   * - **Flame diagnostics**
     - Byram (1959) fireline intensity & flame length;
       Thomas (1963) inverse (I\ :sub:`B` from observed L\ :sub:`f`);
       scorch height; tree mortality; crown activity class
     - Fireline intensity; flame length
     - Fireline intensity
     - Full fire behavior outputs
     - Full fire behavior outputs
     - Physics-derived heat release
   * - **Firebrand spotting**
     - Albini (1983) 2-D trajectory; Albini (1979) torching-tree;
       stochastic lognormal/exponential distance model
     - Albini empirical
     - Not included
     - Not standard
     - Albini (point)
     - Physics-based (select FIRETEC versions)
   * - **Terrain**
     - Per-cell elevation, slope, aspect from LCP or XYZ terrain file;
       2-D domain
     - Full 2-D LCP landscape
     - Full 3-D WRF terrain
     - Full 2-D LCP landscape
     - Point-based (user slope/aspect)
     - Full 3-D terrain + canopy
   * - **Fire–atmosphere coupling**
     - None (prescribed wind; heat-flux plume correction optional)
     - None
     - Two-way (fire ↔ WRF; heat/momentum flux)
     - None
     - None
     - One-way (QUIC-URB) or two-way LES (FIRETEC)
   * - **Fuel burnout / depletion**
     - Residence-time burnout; per-class exponential decay;
       per-cell residual fuel tracking
     - Per-cell burnout
     - Post-frontal consumption
     - Not tracked
     - Not applicable
     - Physics-based combustion
   * - **Ignition types**
     - Point CSV; sphere/box/ellipse; closed polygon; polyline (line fire);
       dynamic fire-points polling; Gaussian sigma blending
     - Interactive ignition map; polygon
     - WPS/WRF restart
     - Ignition map or polygon
     - Single point / simple polygon
     - User-defined locations
   * - **Barrier / suppression**
     - Polyline firebreaks (CSV); multiple barrier files;
       extinguish-on-contact logic
     - Dozer/hand lines; aerial retardant
     - Not included
     - Not included
     - Not applicable
     - Not included
   * - **Spatial resolution**
     - User-defined; AMReX AMR-ready; uniform 2-D default
     - Landscape raster resolution
     - WRF grid (100 m – 1 km)
     - Landscape raster resolution
     - Point-based
     - 1–10 m 3-D grid
   * - **GPU acceleration**
     - ✓ AMReX CUDA / HIP / SYCL kernels
     - ✗
     - ✗
     - ✗
     - ✗
     - Partial (QUIC-URB select versions)
   * - **MPI parallelism**
     - ✓ AMReX domain decomposition
     - ✗
     - ✓ WRF MPI
     - ✗
     - ✗
     - ✓ MPI/OpenMP
   * - **Embedded boundaries**
     - ✓ AMReX EB (buildings, fuel breaks, complex geometry)
     - ✗
     - ✗
     - ✗
     - ✗
     - ✓ (QUIC-URB / FIRETEC 3-D geometry)
   * - **Open source**
     - ✓ (MIT licence)
     - ✗ (USFS proprietary binary)
     - ✓ (WRF open source)
     - ✗ (USFS proprietary binary)
     - ✓ (open source)
     - Partial (QUIC-Fire: research licence; FIRETEC: research only)

Tool Documentation References
------------------------------

For authoritative and up-to-date information on each tool, refer to the
official sources:

* **Wildfire-AMR**: https://hgopalan.github.io/wildfire_levelset/
* **FARSITE**: https://www.firelab.org/project/farsite
* **WRF-SFIRE**: https://github.com/openwfm/WRF-SFIRE
* **FlamMap**: https://www.firelab.org/project/flammap
* **BehavePlus**: https://www.firelab.org/project/behaveplusfiremodeling
* **QUIC-Fire**: https://www.lanl.gov/projects/quic-fire/
* **FIRETEC**: https://www.lanl.gov/org/padwp/adcles/fluid-dynamics-solid-mechanics/index.php

Key Differences from FARSITE
-----------------------------

* **Eulerian level-set vs. explicit Huygens wavelets**: Wildfire-AMR embeds
  the same Richards (1990) elliptical directional spread as FARSITE inside an
  Eulerian level-set (WENO5-Z/RK3). The fire perimeter is the zero contour of
  a signed-distance function; merging fronts and islands are handled
  automatically without explicit connectivity management.

* **Extended spread model library**: In addition to Rothermel (1972),
  Wildfire-AMR includes Balbi (2009) physics-based, Cheney–Gould (1995)
  Australian grassland, Cruz et al. (2005) crown fire, Canadian FBP
  O1a/O1b/S1–S3 (Forestry Canada 1992), and Lautenberger (2013)
  physics-based. FARSITE ships only with Rothermel (1972).

* **Richer crown fire pipeline**: Wildfire-AMR supports Van Wagner (1977)
  passive blending via the crowning fraction CF = (I\ :sub:`B`/I\ :sub:`o`)\ :sup:`2/3`,
  Rothermel (1991) active crown ROS (3.34 × R\ :sub:`surface`), Cruz et al.
  (2005) wind-dependent crown ROS, and full Scott & Reinhardt (2001) TI/CI via
  bisection. FARSITE uses the Rothermel (1991) multiplier internally.

* **Wind adjustments are optional**: FARSITE applies WAF and MEWS internally.
  Wildfire-AMR exposes both via ``rothermel.use_waf`` and
  ``rothermel.use_wind_limit`` so users can match FARSITE behaviour or supply
  midflame-height wind directly.

* **GPU and MPI**: AMReX CUDA/HIP/SYCL kernels and MPI domain decomposition.
  FARSITE is serial and CPU-only.

* **Embedded Boundary**: AMReX EB allows buildings and fuel breaks on the
  Cartesian grid without remeshing. FARSITE has no EB support.

Key Differences from WRF-Fire (WRF-SFIRE)
------------------------------------------

* **No atmospheric coupling**: Wildfire-AMR uses prescribed wind fields
  (constant, CSV, or WRF output); WRF-SFIRE fully couples fire with WRF,
  including fire-induced wind, heat flux, and smoke transport.

* **Wind adjustment**: When driving with WRF output, enable
  ``rothermel.use_waf = 1`` to convert NWP wind to midflame height.
  WRF-SFIRE handles this inside the coupled framework.

* **Richer fire behaviour models**: WRF-SFIRE uses Rothermel (1972) only.
  Wildfire-AMR supports seven additional spread models and a richer crown fire
  pipeline (see above).

* **Simpler setup**: Wildfire-AMR requires only CMake and AMReX; WRF-SFIRE
  requires a full WRF stack (NetCDF, MPI, WPS, WRF pre-processing).

* **GPU-native kernels**: WRF-SFIRE is CPU-MPI; Wildfire-AMR uses AMReX GPU
  kernels throughout.

Key Differences from FlamMap
-----------------------------

* **Time-dependent propagation**: FlamMap computes static fire behaviour maps
  (no time stepping); Wildfire-AMR evolves the fire front dynamically via
  level-set, FARSITE Huygens, or MTT.

* **Crown fire depth**: FlamMap provides the full Scott & Reinhardt (2001)
  crown fire assessment. Wildfire-AMR now matches this with bisection-based
  TI/CI plus Van Wagner (1977) passive blending and Rothermel (1991) active
  crown ROS.

* **GPU / open source**: FlamMap is a closed-source Windows binary; Wildfire-AMR
  is MIT-licensed and GPU-accelerated.
