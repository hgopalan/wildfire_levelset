New FARSITE-Parity Features (2025)
===================================

This page documents the easy- and medium-difficulty FARSITE features added in the 2025
update.  All new C++ headers live in ``src/`` and are included automatically by ``main.cpp``.

.. contents::
   :depth: 2
   :local:

Scott & Reinhardt (2001) Torching / Crowning Index Diagnostics
---------------------------------------------------------------

**Header**: ``src/fire_ecology_model.H``

The ``ECOLOGY_NCOMP`` constant was increased from 4 to 6.  Two new per-cell
diagnostic components are now written to every plotfile:

* **Component 4** — ``torching_ratio`` = :math:`I_B / I_o`

  Ratio of Byram fireline intensity to the Van Wagner (1977) crown initiation threshold.
  A value ≥ 1.0 indicates that passive or active crown fire conditions are met.

  .. math::

     I_o \ [\text{kW/m}] = 0.010 \times \text{CBH} \times (460 + 25.9 \times \text{FMC})

* **Component 5** — ``crowning_ratio`` = :math:`R / R^\prime_{SA}`

  Ratio of surface or crown ROS to the Van Wagner (1977) critical active-crown ROS.
  A value ≥ 1.0 indicates active crown fire conditions.

  .. math::

     R^\prime_{SA} \ [\text{m/s}] = \frac{3.0}{\text{CBD} \times 60}

These are dimensionless proxies for the Scott & Reinhardt (2001) Torching Index (TI)
and Crowning Index (CI), computed without the expensive nonlinear inversion of the
Rothermel wind function.

**Parameters**: No new input parameters required.  Uses ``crown.CBH``, ``crown.CBD``,
and ``crown.FMC`` from the existing crown initiation parameters.

**Test**: ``regtest/diagnostics/scott_reinhardt_indices/``

FMC Seasonal / Phenological Schedule
--------------------------------------

**Header**: ``src/fmc_schedule.H``

The foliar moisture content (FMC) used by the Van Wagner (1977) crown fire initiation
threshold is now updated each timestep from a seasonal schedule.  Two modes:

1. **CSV file** (``fmc_schedule.file``): two-column CSV of (day-of-year, fmc_pct).
2. **Built-in parametric curve** (``fmc_schedule.use_farsite_curve = 1``): piecewise linear
   phenology with configurable spring green-up and autumn curing breakpoints.

Default breakpoints match Southern California chaparral (Rothermel 1983):

+----------------+------------+-------+
| Phase          | DOY range  | FMC   |
+================+============+=======+
| Winter dormant | 1–89       | 85 %  |
+----------------+------------+-------+
| Green-up       | 90–149     | rises |
+----------------+------------+-------+
| Summer peak    | 150–239    | 140 % |
+----------------+------------+-------+
| Curing         | 240–300    | falls |
+----------------+------------+-------+
| Post-cure      | 301–365    | 85 %  |
+----------------+------------+-------+

**Input parameters** (prefix ``fmc_schedule.``):

.. list-table::
   :header-rows: 1

   * - Parameter
     - Default
     - Description
   * - ``enable``
     - 0
     - 1 to enable FMC seasonal schedule
   * - ``file``
     - ""
     - Path to two-column CSV (doy, fmc_pct)
   * - ``use_farsite_curve``
     - 0
     - 1 to use built-in parametric curve
   * - ``start_doy``
     - 1
     - Day-of-year at simulation t = 0
   * - ``spring_start``
     - 90
     - DOY when green-up begins
   * - ``summer_peak``
     - 150
     - DOY when FMC reaches maximum
   * - ``fall_start``
     - 240
     - DOY when curing begins
   * - ``fall_end``
     - 300
     - DOY when curing ends
   * - ``fmc_min``
     - 85.0
     - Dormant / cured FMC [%]
   * - ``fmc_max``
     - 140.0
     - Peak green FMC [%]

**Test**: ``regtest/crown_fire/fmc_seasonal/``

Precipitation-Driven Dead Fuel Moisture Wetting
------------------------------------------------

**Header**: ``src/precipitation_moisture.H``

When ``diurnal_moisture.enable = 1`` and a rain rate is specified, dead fuel
moisture is driven by an exponential wetting / drying model:

* **Wetting** (rain > threshold): :math:`M_\text{new} = M_\text{sat}(1 - e^{-\Delta t / \tau_\text{wet}}) + M_\text{old}\,e^{-\Delta t / \tau_\text{wet}}`
* **Drying** (no rain): :math:`M_\text{new} = M_\text{EMC} + (M_\text{old} - M_\text{EMC})\,e^{-\Delta t / \tau_\text{dry}}`

Per size-class time constants (Rothermel 1983):

+----------+----------+----------+
| Class    | τ_wet    | τ_dry    |
+==========+==========+==========+
| 1-hr     | 1 hr     | 4 hr     |
+----------+----------+----------+
| 10-hr    | 2 hr     | 24 hr    |
+----------+----------+----------+
| 100-hr   | 4 hr     | 10 days  |
+----------+----------+----------+

**Input parameters** (prefix ``diurnal_moisture.``):

.. list-table::
   :header-rows: 1

   * - Parameter
     - Default
     - Description
   * - ``precip_rain_rate_mm_hr``
     - 0.0
     - Constant rain rate [mm/hr]
   * - ``precip_schedule_file``
     - ""
     - CSV of (time_s, rain_mm_hr) for time-varying rain
   * - ``precip_threshold_mm_hr``
     - 0.25
     - Minimum rain rate to trigger wetting [mm/hr]
   * - ``M_sat``
     - 1.20
     - Saturation moisture content [fraction]

**Test**: ``regtest/moisture/precip_wetting/``

Polygon and Polyline Fire Ignition
-----------------------------------

**Header**: ``src/polygon_ignition.H``

Two new ``source_type`` values allow fires to be initialised from GIS-style vertex files:

* ``source_type = polygon``: Closed ignition polygon.  Uses the 2-D winding-number
  algorithm to determine interior cells (GPU-safe ``AMREX_GPU_HOST_DEVICE`` function).
  Interior receives ``phi = −distance_to_boundary``; exterior receives the positive
  signed distance.

* ``source_type = polyline``: Line-fire ignition.  Cells within ``polyline_width``
  of any segment receive ``phi = −polyline_width``; other cells receive the distance
  to the nearest segment.

Both modes read vertices from a CSV file (X, Y[, Z] columns; ``#`` comment lines allowed).
The GPU kernel uses ``amrex::Gpu::DeviceVector`` to store vertices on the device and
``amrex::ParallelFor`` for GPU-safe rasterisation.

**Input parameters**:

.. list-table::
   :header-rows: 1

   * - Parameter
     - Default
     - Description
   * - ``source_type``
     - ``sphere``
     - ``polygon`` or ``polyline``
   * - ``fire_polygon_file``
     - ""
     - Path to CSV vertex file (required for polygon/polyline)
   * - ``polyline_width``
     - 10.0
     - Half-width [m] of line-ignition zone
   * - ``fire_polygon_z_level``
     - 0.5
     - Z-centre of ignition layer in 3-D runs [m]

**Tests**: ``regtest/ignition/polygon_ignition/``, ``regtest/ignition/polyline_ignition/``

Per-Cell Live Canopy Moisture from FMS File
--------------------------------------------

**Header**: ``src/fms_moisture.H``

A FARSITE fuel moisture scenario (``.fms``) file can be supplied via ``fms_file`` to
provide per-fuel-model dead and live moisture values.  These are stored in a spatially-
varying ``MultiFab`` (``FMS_MOISTURE_NCOMP = 5`` components) and override the global
``rothermel.M_d1``, …, ``rothermel.M_lw`` on a per-cell basis.

**Requires**: A landscape file (``rothermel.landscape_file``) to supply per-cell fuel
model codes.

**Input parameter**: ``fms_file = path/to/fire.fms``

The GPU kernel uses device-side lookup arrays populated from the FMS map before the
kernel launch — no host functions or ``std::map`` are used inside ``ParallelFor``.

Compact Wind Direction Schedule
---------------------------------

**Header**: ``src/wind_dir_schedule.H``

A compact three-column CSV schedule (``time_s, speed_ms, dir_deg``) updates the uniform
wind field at every timestep without requiring a full spatial wind grid per snapshot.

Meteorological convention: ``dir_deg`` is the direction *from* which the wind blows
(270° = westerly wind blowing eastward).  The schedule is linearly interpolated with
circular direction averaging to avoid 0°/360° wrap-around artefacts.

**Input parameter**: ``wind_dir_schedule_file = path/to/wind_schedule.csv``

When this file is set it overrides ``u_x``/``u_y`` and the ``use_time_dependent_wind``
file-grid path.  Compatible with turbulent wind perturbation (perturbations are added
on top of the scheduled base wind).

**Test**: ``regtest/wind/wind_dir_schedule/``

GPU Safety Notes
-----------------

All new kernels follow the AMReX GPU-safe programming model:

* ``AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE`` decorators on helper functions called
  from ``ParallelFor`` lambdas (``segment_dist_sq``, ``polygon_winding_number``).
* ``amrex::Gpu::DeviceVector`` for per-cell or per-vertex data passed to kernels.
* ``amrex::Gpu::copy(hostToDevice, ...)`` for transferring host ``std::vector`` data.
* No ``std::vector``, ``std::map``, ``std::string``, or host-only library calls inside
  ``ParallelFor`` lambdas.
* Host-only functions (schedule loading, file parsing, Print/Abort) are called outside
  kernels and store results in device-accessible arrays before kernel launch.
* All kernel captures are by value (scalars or device pointers).

The new headers compile cleanly for CUDA, HIP, and SYCL backends.
