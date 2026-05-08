FARSITE Topographic Horizon Shading
=====================================

.. contents::
   :depth: 2
   :local:

Overview
--------

``solar_radiation.H`` already computes per-cell terrain shading from the local
surface normal (incidence-angle test) and canopy shading from cover fraction.
Neither method detects shadows cast by terrain **tens or hundreds of cells
away** — e.g. a deep canyon cell that lies in the shadow of a flanking ridge
for much of the morning.

FARSITE addresses this with a **full topographic horizon scan**: for each cell,
march 8 compass directions to find the maximum elevation angle to any visible
terrain point, then shade the cell whenever the solar elevation is below that
angle.  This is now available in ``wildfire_levelset`` via the
``solar_radiation.use_topographic_horizon`` input parameter.

Algorithm
---------

**Precomputation (one-time)**

The function ``compute_topographic_horizon_angles`` in ``src/solar_radiation.H``
executes once after terrain setup:

1. The elevation ``MultiFab`` is copied to pinned (host-accessible) memory.
2. A flat ``std::vector<Real>`` covering the full domain is built from the
   per-tile contributions and summed across all MPI ranks with
   ``ParallelDescriptor::ReduceRealSum``.
3. For each locally-owned cell ``(i, j)`` and each of the 8 compass directions
   ``d``, the function ray-marches outward:

   .. math::

      \theta_\text{hz}^{(d)}(i,j) =
      \max_{s=1,\ldots,S}
      \arctan\!\left(
        \frac{z(i + s\,\delta i_d,\; j + s\,\delta j_d) - z(i,j)}{s \cdot \Delta_d}
      \right)

   where :math:`\Delta_d` is the physical step distance in direction
   :math:`d` (cell spacing for axis-aligned directions,
   :math:`\sqrt{\Delta x^2 + \Delta y^2}` for diagonals).

4. The result is stored in an 8-component ``MultiFab`` (one angle per direction,
   in radians).

**Per-timestep shading check (GPU)**

Inside ``compute_shade_fraction_mf``, the GPU kernel interpolates the horizon
angle for the current solar azimuth between the two bounding compass directions:

.. math::

   d_\text{lo} &= \lfloor \varphi_s / 45°\rfloor \bmod 8 \\
   d_\text{hi} &= (d_\text{lo} + 1) \bmod 8 \\
   t &= (\varphi_s / 45° - d_\text{lo}) \\
   \theta_\text{hz} &= \theta_\text{hz}^{(d_\text{lo})} \cdot (1 - t)
                     + \theta_\text{hz}^{(d_\text{hi})} \cdot t

A cell is marked as fully shaded (:math:`f = 1`) when

.. math::

   \theta_s < \theta_\text{hz}

where :math:`\theta_s` is the solar elevation angle.

Input Parameters
----------------

All parameters are under the ``solar_radiation.`` prefix.

.. list-table::
   :header-rows: 1
   :widths: 38 12 50

   * - Parameter
     - Default
     - Description
   * - ``use_topographic_horizon``
     - ``0``
     - Set to ``1`` to enable the FARSITE 8-direction horizon scan.
       **Disabled by default** because the scan involves an MPI global
       elevation gather and an O(N²) CPU sweep.
   * - ``horizon_scan_max_dist_m``
     - ``0.0``
     - Maximum ray-march distance [m]. ``0`` scans the full domain.
       Reduce this to bound cost on large domains (e.g. ``800.0`` for a
       1 km × 1 km domain at 15 m resolution).

Example::

    solar_radiation.enable                  = 1
    solar_radiation.use_topographic_horizon = 1
    solar_radiation.horizon_scan_max_dist_m = 800.0

Performance Notes
-----------------

.. warning::

   The topographic horizon scan is **expensive** and should be disabled
   when not needed.

+---------------------------+------------------------------------------+
| Cost source               | Mitigation                               |
+===========================+==========================================+
| MPI ``ReduceRealSum`` on  | Proportional to domain cells.  Large     |
| global elevation array    | domains (> 10⁶ cells) will incur         |
|                           | significant communication overhead.      |
+---------------------------+------------------------------------------+
| O(N³) CPU sweep           | Use ``horizon_scan_max_dist_m`` to limit  |
| (N² cells × N steps/ray)  | the scan radius to the relevant distance. |
+---------------------------+------------------------------------------+
| One-time precomputation   | Cost paid once at startup; negligible    |
| (not per timestep)        | relative to a long simulation.           |
+---------------------------+------------------------------------------+

Set ``solar_radiation.use_topographic_horizon = 0`` (the default) to skip
the scan entirely and use only local surface-normal and canopy shading.

Interaction with Other Shading Sources
---------------------------------------

The total shade fraction applied to fuel temperature is:

.. code-block:: none

    f = max(f_local_terrain,  f_canopy,  f_topographic_horizon,  cloud_cover)

where:

- **f_local_terrain**: 1 when the local surface normal faces away from the sun
  (cos_i ≤ 0), 0 otherwise.
- **f_canopy**: canopy_cover [%] / 100 (requires ``use_canopy_shading = 1``
  and spatial crown data from an LCP file).
- **f_topographic_horizon**: 1 when the solar elevation < interpolated ridge
  horizon angle (requires ``use_topographic_horizon = 1``).
- **cloud_cover**: domain-uniform overcast fraction (``solar_radiation.cloud_cover``).

Shade-Adjusted EMC
------------------

When ``diurnal_moisture.enable = 1``, the shade fraction modulates the solar
heating increment applied to fuel temperature:

.. math::

   T_\text{fuel}(i,j) = T_\text{air} + \Delta T_\text{solar} \cdot (1 - f)

The equilibrium moisture content (EMC) is then recomputed with
:math:`T_\text{fuel}` using the Simard (1968) equations, and the result is
written to components 0–2 of the spatial-moisture ``MultiFab``
(:math:`M_{d1}`, :math:`M_{d10}`, :math:`M_{d100}`).

Cells in deep canyon shadow receive *less* solar heating and therefore retain
*higher* fuel moisture throughout the day, reducing their fire-weather hazard.

Regression Test
---------------

``regtest/terrain/solar_horizon_shading/`` — a 1 km × 1 km synthetic
ridge-canyon-ridge terrain:

- W/E ridges at ~200 m, canyon floor at ~50 m.
- Nominal E/W horizon angle from the floor ≈ 16.7°.
- Start time 06:30 AM PDT, July 1, 2024 (lat 34.1° N) → solar elevation ≈ 10°.
- The horizon scan correctly shadows the canyon floor; local surface-normal
  shading would leave it unshaded.

Run with::

    cd build
    ctest -R solar_horizon_shading --output-on-failure

Python Analysis Tool
--------------------

``tools/topographic_horizon_analysis.py`` is a pure-Python implementation of
the same algorithm.  It is useful for:

- Pre-computing and visualising horizon maps before a full C++ run.
- Diagnosing whether canyon cells will be shaded at a given time of day.
- Understanding the algorithm without reading C++ code.

Usage::

    # Synthetic canyon (no terrain file needed)
    python tools/topographic_horizon_analysis.py --synthetic \
        --solar-elevation 10.0 --solar-azimuth 90.0 --plot

    # Real terrain file
    python tools/topographic_horizon_analysis.py terrain.csv \
        --nx 64 --ny 64 --max-scan-dist 800 \
        --solar-elevation 10.0 --solar-azimuth 90.0

References
----------

- Finney, M.A. (2004). FARSITE: Fire Area Simulator—Model Evaluation,
  Modifications and Testing. USDA FS Research Paper RMRS-RP-4 (Revised).
- Nelson, R.M. Jr. (2000). Prediction of diurnal change in 10-h fuel stick
  moisture content. *Canadian Journal of Forest Research*, 30(7), 1071–1087.
- Simard, A.J. (1968). The moisture content of forest fuels.
  Forestry Branch Information Report FF-X-14, Dept. of Forestry, Canada.
