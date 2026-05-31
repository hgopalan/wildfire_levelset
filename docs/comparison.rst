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
     - Rothermel (1972) plus six alternative models
       (Balbi, Cheney–Gould, Cruz, FBP, Lautenberger, Viegas)
     - Rothermel (1972) in all three
     - Semi-empirical QUIC; physics LES FIRETEC
     - Rothermel (1972)
   * - **Propagation method**
     - Eulerian level-set (WENO5-Z/RK3) • FARSITE Huygens • MTT
     - FARSITE: Huygens wavelet • FlamMap: MTT/Huygens •
       BehavePlus: point only
     - Eulerian CFD (QUIC-URB / FIRETEC LES)
     - Eulerian level-set on WRF grid
   * - **Crown fire**
     - Van Wagner (1977) • Rothermel (1991) • Cruz (2005) •
       Scott–Reinhardt TI/CI • CFSA
     - Van Wagner (1977) • FARSITE/FlamMap: R\ :sub:`c`\ = 3.34 R\ :sub:`s` •
       BehavePlus: point calculation
     - Physics-based combustion
     - Van Wagner (1977)
   * - **Wind adjustment**
     - WAF (Andrews/BehavePlus) • MEWS cap • Multi-layer canopy •
       Wind-terrain models
     - FARSITE/FlamMap: WAF + MEWS • BehavePlus: user WAF
     - 3-D mass-consistent (QUIC-URB) or LES
     - WRF-derived with WAF coupling
   * - **Fuel models**
     - FBFM13 • FBFM40 • FBP • Lautenberger •
       Per-cell LCP • Two-fuel blending
     - FARSITE/FlamMap/BehavePlus: FBFM13 + FBFM40
     - Custom 3-D bulk density fields
     - FBFM13
   * - **Fuel moisture**
     - All size classes • FMD/Nelson diurnal • Precipitation •
       FWI system • Time-lag DE • Duff/smolder
     - FARSITE/FlamMap: dead/live + FMD + conditioning •
       BehavePlus: per size class
     - Bulk moisture fields
     - Dead/live (prescribed)
   * - **Non-burnable masking**
     - ✓ Codes 91–99 / NB1–NB9 → ROS = 0
     - ✓ (all three)
     - N/A (3-D grid)
     - Partial (fuel mask)
   * - **Flame diagnostics**
     - Byram intensity • Flame length/tilt • Convective number •
       Packing ratio • Scorch height • TI/CI • NFDRS ERC • McArthur FFDI
     - FARSITE: intensity + flame length •
       FlamMap/BehavePlus: full outputs
     - Physics-based heat release
     - Intensity + flame length
   * - **Firebrand spotting**
     - Albini (1983) trajectory • Ember cascade •
       Probabilistic ignition
     - FARSITE: Albini empirical •
       FlamMap: optional • BehavePlus: point only
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
     - Single/multi-station .wtr • FMD schedule • IDW interpolation
     - FARSITE: per-station .wtr • FlamMap: gridded/single •
       BehavePlus: single point
     - Prescribed per cell
     - WRF atmospheric profiles
   * - **Fire–atmosphere coupling**
     - None (prescribed wind with optional plume correction)
     - None
     - One-way QUIC-URB • Two-way FIRETEC LES
     - Two-way (fire ↔ WRF)
   * - **Barrier / suppression**
     - Polyline firebreaks • Aerial retardant
     - FARSITE: dozer/hand lines + retardant •
       FlamMap/BehavePlus: none
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

* **Propagation method**: Eulerian level-set (WENO5-Z/RK3) with automatic merging vs. FARSITE's explicit Huygens wavelets

* **Extended spread models**: Six alternative models beyond Rothermel (1972)

* **Non-burnable masking**: Explicit zeroing of fuel codes 91–99 and NB1–NB9

* **Multi-station weather**: IDW interpolation for spatially-varying wind/T/RH

* **Optional wind adjustments**: User control of WAF and MEWS vs. FARSITE's internal application

* **GPU and MPI**: AMReX CUDA/HIP/SYCL kernels with domain decomposition vs. serial CPU-only

* **Embedded Boundary**: AMReX EB support for buildings and fuel breaks

Key Differences from WRF-Fire (WRF-SFIRE)
------------------------------------------

* **Atmospheric coupling**: Prescribed wind fields vs. WRF-SFIRE's fully coupled fire–atmosphere system

* **Fire behavior models**: Seven spread models and extended crown fire pipeline vs. Rothermel (1972) only

* **Setup complexity**: CMake + AMReX vs. full WRF stack (NetCDF, MPI, WPS)

* **GPU support**: Native AMReX GPU kernels vs. CPU-MPI only

Key Differences from FlamMap
-----------------------------

* **Propagation**: Dynamic time-stepping (level-set/Huygens/MTT) vs. static fire behavior maps

* **Crown fire**: Scott & Reinhardt TI/CI with Van Wagner and Cruz models vs. FlamMap's static assessment

* **Licensing**: MIT open source with GPU support vs. proprietary Windows binary

