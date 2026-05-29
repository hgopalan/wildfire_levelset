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

Tools are grouped into three columns to keep the table compact:

* **Group A — Operational fire behavior tools**: FARSITE · FlamMap · BehavePlus
* **Group B — Physics-based LES/CFD tools**: QUIC-Fire (QUIC-URB) · FIRETEC
* **Group C — Coupled fire–atmosphere NWP**: WRF-Fire (WRF-SFIRE)

.. list-table::
   :header-rows: 1
   :widths: 22 28 20 15 15

   * - Capability
     - **Wildfire-AMR** (this solver)
     - **Group A** — FARSITE / FlamMap / BehavePlus
     - **Group B** — QUIC-Fire / FIRETEC
     - **Group C** — WRF-Fire (WRF-SFIRE)
   * - **Surface spread model**
     - Rothermel (1972), Balbi (2009), Cheney–Gould (1995),
       Cruz (2005), FBP (1992), Lautenberger (2013),
       plus Viegas eruptive option
     - Rothermel (1972) in all three
     - Semi-empirical QUIC; physics LES FIRETEC
     - Rothermel (1972)
   * - **Propagation method**
     - Eulerian level-set (WENO5-Z / RK3);
       FARSITE Huygens ellipse; MTT
     - FARSITE: Huygens wavelet;
       FlamMap: MTT or Huygens;
       BehavePlus: point (no spatial propagation)
     - Eulerian CFD (QUIC-URB / HIGRAD-FIRETEC)
     - Eulerian level-set on WRF grid
   * - **Crown fire**
     - Van Wagner (1977), Rothermel (1991), Cruz (2005);
       Scott–Reinhardt (2001) TI/CI; CFSA *(new)*
     - Van Wagner (1977); FARSITE/FlamMap active
       via R\ :sub:`c`\ = 3.34 R\ :sub:`s`;
       BehavePlus: Van Wagner point calc.
     - Physics-based combustion
     - Van Wagner (1977)
   * - **Wind adjustment**
     - WAF (Andrews 2018); MEWS cap; multi-layer canopy *(new)*;
       8 wind-terrain models
     - FARSITE/FlamMap: WAF + MEWS (internal);
       BehavePlus: user-specified WAF
     - 3-D mass-consistent (QUIC-URB) or LES
     - WRF-derived; WAF in coupling layer
   * - **Fuel models**
     - FBFM13, FBFM40, FBP grass/slash, Lautenberger,
       per-cell LCP; two-fuel blending *(new)*
     - FARSITE/FlamMap: FBFM13 + FBFM40;
       BehavePlus: FBFM13 + FBFM40
     - Custom 3-D bulk density per cell
     - FBFM13
   * - **Fuel moisture**
     - All size classes; FMD/Nelson diurnal; precipitation;
       FWI system; time-lag DE *(new)*; duff/smolder *(new)*
     - FARSITE/FlamMap: dead/live + FMD + conditioning;
       BehavePlus: dead/live per class
     - Bulk moisture per cell
     - Dead/live (prescribed)
   * - **Non-burnable masking**
     - ✓ Codes 91–99 / NB1–NB9 → ROS = 0
     - ✓ (all three)
     - N/A (3-D grid)
     - Partial (fuel mask)
   * - **Flame diagnostics**
     - Byram intensity + flame length + tilt + convective number;
       packing ratio (β/β_opt); flame depth; scorch height; tree mortality;
       TI/CI; NFDRS ERC; McArthur FFDI; fire acceleration *(new)*
     - FARSITE: intensity + flame length;
       FlamMap: full outputs;
       BehavePlus: full outputs
     - Physics heat release
     - Intensity + flame length
   * - **Firebrand spotting**
     - Albini (1983) 2-D trajectory + torching;
       stochastic distance model
     - FARSITE: Albini empirical;
       FlamMap: not standard;
       BehavePlus: Albini (point)
     - Select FIRETEC versions
     - Not included
   * - **Terrain & landscape**
     - Per-cell elev./slope/aspect/fuel from LCP
       or XYZ terrain file; 2-D domain
     - All three: full 2-D LCP landscape;
       BehavePlus: user slope/aspect
     - Full 3-D terrain + canopy
     - Full 3-D WRF terrain
   * - **Weather input**
     - Single .wtr; FMD schedule;
       **multi-station IDW** spatial interpolation
     - FARSITE: per-station .wtr;
       FlamMap: gridded or single station;
       BehavePlus: single point
     - Prescribed per-cell
     - WRF atmospheric profiles
   * - **Fire–atmosphere coupling**
     - None (prescribed wind; heat-flux plume
       correction optional)
     - None in all three
     - One-way QUIC-URB; two-way LES FIRETEC
     - Two-way (fire ↔ WRF)
   * - **Barrier / suppression**
     - Polyline firebreaks (CSV);
       aerial retardant (ROS + spotting suppressed)
     - FARSITE: dozer/hand lines + retardant;
       FlamMap/BehavePlus: not included
     - Not included
     - Not included
   * - **GPU acceleration**
     - ✓ AMReX CUDA / HIP / SYCL
     - ✗ (all three are serial / CPU-only)
     - Partial (QUIC-URB select versions)
     - ✗
   * - **MPI parallelism**
     - ✓ AMReX domain decomposition
     - ✗ (all three serial)
     - ✓ MPI/OpenMP
     - ✓ WRF MPI
   * - **Embedded boundaries**
     - ✓ AMReX EB (buildings, fuel breaks)
     - ✗ (all three)
     - ✓ QUIC-URB / FIRETEC 3-D geometry
     - ✗
   * - **Open source**
     - ✓ MIT licence
     - FARSITE/FlamMap: proprietary USFS binary;
       BehavePlus: open source
     - Research licence (QUIC-Fire);
       restricted (FIRETEC)
     - ✓ WRF open source

Tool Documentation References
------------------------------

For authoritative and up-to-date information on each tool, refer to the
official sources:

* **Wildfire-AMR**: https://hgopalan.github.io/wildfire_levelset/
* **FARSITE**: https://www.firelab.org/project/farsite
* **FlamMap**: https://www.firelab.org/project/flammap
* **BehavePlus**: https://www.firelab.org/project/behaveplusfiremodeling
* **QUIC-Fire**: https://www.lanl.gov/projects/quic-fire/
* **FIRETEC**: https://www.lanl.gov/org/padwp/adcles/fluid-dynamics-solid-mechanics/index.php
* **WRF-SFIRE**: https://github.com/openwfm/WRF-SFIRE

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

* **Non-burnable cell masking**: Fuel model codes 91–99 and NB1–NB9 (water,
  rock, urban, bare ground) are explicitly zeroed in the ROS kernel so fire
  cannot creep through sparse-fuel numerical noise into non-burnable areas.

* **Multiple weather stations**: ``multi_wtr_file`` loads per-station .wtr
  files and produces spatially-varying wind and T/RH via IDW interpolation,
  matching FARSITE's multi-station weather capability.

* **Retardant spotting suppression**: ``retardant_file`` now suppresses both
  ROS and spotting probability inside active drop zones, consistent with
  FARSITE's aerial retardant model.

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
  crown fire assessment. Wildfire-AMR matches this with bisection-based
  TI/CI plus Van Wagner (1977) passive blending and Cruz et al. (2005) active
  crown ROS.

* **GPU / open source**: FlamMap is a closed-source Windows binary; Wildfire-AMR
  is MIT-licensed and GPU-accelerated.

