New FARSITE-Parity Features (2025)
===================================

This page documents the easy- and medium-difficulty FARSITE features added in the 2025
update.  All new C++ headers live in ``src/`` and are included automatically by ``main.cpp``.

.. contents::
   :depth: 2
   :local:

WAF Formula Selection (Andrews vs. BehavePlus)
-----------------------------------------------

**Header**: ``src/andrews_model.H``

The Wind Adjustment Factor (WAF) formula used when ``rothermel.use_waf = 1`` is now
selectable via ``rothermel.waf_formula``.  Two fuel-model–aware formulas are provided:

**"andrews"** (default) — logarithmic Albini & Baughman (1979)

  *Open / shrub fuel beds* (unsheltered):

  .. math::

     \text{WAF} = \frac{1.83}{\ln\!\left(\dfrac{20 + 0.36\,h}{0.13\,h}\right)}

  *Forest / closed-canopy fuel beds* (sheltered, FARSITE-style):

  .. math::

     \text{WAF}_\text{sheltered} = \frac{0.555}{\sqrt{f_c\,h_c}\,
       \ln\!\left(\dfrac{20 + 0.36\,h_c}{0.13\,h_c}\right)}

  The sheltered formula is applied automatically when canopy cover fraction
  :math:`f_c \geq 0.5` and canopy height :math:`h_c >` fuel depth :math:`h`.

**"behaviorplus"** — BehavePlus linear / exponential

  *Open / shrub fuel beds* (unsheltered):

  .. math::

     \text{WAF} = 0.36 + 0.004\,h_\text{in}   \quad [h_\text{in} = 12 h_\text{ft}]

  This linear approximation is used by BehavePlus for grass, shrub, and open
  timber-litter fuel models.

  *Forest / closed-canopy fuel beds* (sheltered):

  .. math::

     \text{WAF}_\text{canopy} = (0.36 + 0.004\,h_{c,\text{in}}) \times
       \exp(-\alpha_c\,f_c)

  The exponential Beer–Lambert term accounts for the wind speed attenuation
  through the canopy column.  The attenuation strength is controlled by
  ``rothermel.waf_canopy_alpha`` (:math:`\alpha_c`, default 1.5).

**Input parameters** (prefix ``rothermel.``):

.. list-table::
   :header-rows: 1

   * - Parameter
     - Default
     - Description
   * - ``use_waf``
     - 0
     - 1 to enable Wind Adjustment Factor
   * - ``waf_formula``
     - ``"andrews"``
     - Formula: ``"andrews"`` or ``"behaviorplus"``
   * - ``waf_canopy_alpha``
     - 1.5
     - Canopy attenuation coefficient α_c for BehavePlus exponential canopy WAF

**Example — BehavePlus open/shrub WAF**::

    rothermel.use_waf       = 1
    rothermel.waf_formula   = behaviorplus

**Example — BehavePlus with custom canopy attenuation**::

    rothermel.use_waf          = 1
    rothermel.waf_formula      = behaviorplus
    rothermel.waf_canopy_alpha = 2.0   # stronger sheltering

**Tests**: ``regtest/wind/waf_andrews/``, ``regtest/wind/waf_behaviorplus/``



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

Non-Burnable Cell Masking
--------------------------

**Header**: ``src/compute_rothermel_R.H``

Fuel model codes that correspond to non-burnable landscape types are now
explicitly zeroed in the Rothermel ROS kernel.  This prevents fire from
numerically creeping across water bodies, bare rock, or urban areas due to
floating-point noise in zero-fuel cells.

Affected codes:

* **FBFM40 / FBFM13** landscape integer codes: 91 (Urban/Developed), 92
  (Snow/Ice), 93 (Agriculture), 98 (Open Water), 99 (Bare Ground).

When ``rothermel.landscape_file`` is provided and the per-cell fuel code
maps to one of these values, the kernel returns ``R = 0`` for that cell
without entering the wind-slope-fuel computation.

**Parameters**: None.  Automatic when ``rothermel.landscape_file`` is used.

**Test**: ``regtest/misc/nonburnable_mask/`` (setup: ``create_landscape.py``).

ROS Floor (FARSITE Stall Threshold)
-------------------------------------

**Header**: ``src/constants.H``, ``src/compute_rothermel_R.H``

A FARSITE-style stall threshold is applied: cells where the computed ROS is
below ``WildfireConst::ROS_MIN_MS`` (1 × 10⁻⁴ m/s) are forced to zero.
This prevents numerical drift in extremely low-moisture conditions and
matches FARSITE's internal practice of not propagating the fire front when
spread is negligible.

**Parameters**: None (constant in ``constants.H``).

Per-Fuel-Model Burnout (Residence) Time
-----------------------------------------

**Header**: ``src/burnout_time.H``

``compute_burnout_time`` now accepts an optional per-cell fuel model MultiFab
and a device-side ``FuelResidenceTime`` table populated by
``build_fuel_tau_table``.  The residence time for each fuel code is derived
from the Rothermel (1983) formula:

.. math::

   \tau_i \ [\text{s}] = \frac{\alpha \rho_p}{\sigma_i}

where :math:`\alpha` = ``WildfireConst::BURNOUT_ALPHA`` = 3600 s·ft⁻¹/(lb/ft³),
:math:`\rho_p` = 32 lb/ft³ (particle density), and :math:`\sigma_i` = 1739 ft⁻¹
(1-hr fine fuel representative SAV).  When no landscape file is active the
global ``farsite.tau_residence`` scalar is used (backward compatible).

**Parameters**: Automatic when ``rothermel.landscape_file`` is present.

Live Fuel Moisture Conditioning
---------------------------------

**Header**: ``src/main.cpp`` (conditioning block)

During the pre-simulation conditioning period (``conditioning.n_days > 0``),
live fuel moistures (M_lh, M_lw) are now ramped linearly from their initial
values toward equilibrium dead-fuel targets over the conditioning steps.
This matches FARSITE's behaviour of gradually updating live fuel moisture
during the spin-up phase rather than holding it constant.

The ramp is:

.. code-block:: text

   M_lh_target = M_d100 × 1.5   (live herbaceous ≈ 150% of 100-hr dead)
   M_lw_target = M_d100 × 2.0   (live woody     ≈ 200% of 100-hr dead)

At the end of conditioning, updated M_lh and M_lw are logged alongside the
dead fuel values.

**Parameters**: Uses ``conditioning.n_days`` (existing parameter).

Spotting Suppression Inside Retardant Drop Zones
--------------------------------------------------

**Header**: ``src/retardant_drop.H``

A new function ``apply_retardant_to_spotting_probability`` is called after
``apply_retardant_to_ros`` in the time loop.  Inside active drop zones the
spotting probability (component 0 of ``spotting_data``) is scaled by the
same factor as ROS:

.. math::

   f_{\text{ret}} = 1 - \varepsilon \cdot e^{-(t - t_{\text{drop}}) / \tau_{\text{decay}}}

This ensures that cells covered by retardant do not generate long-range spot
fires, consistent with FARSITE's suppression model.

**Parameters**: Existing ``retardant_file`` (no new parameters).

Spatial Fuel Moisture in Plotfiles
------------------------------------

**Header**: ``src/main.cpp``

The ``spatial_moisture_mf`` MultiFab (5 components: 1-hr, 10-hr, 100-hr dead
plus live herb / woody) is now written to every AMReX plotfile as five named
variables:

.. list-table::
   :header-rows: 1
   :widths: 22 50

   * - Variable name
     - Description
   * - ``moisture_d1``
     - 1-hr dead fuel moisture [fraction]
   * - ``moisture_d10``
     - 10-hr dead fuel moisture [fraction]
   * - ``moisture_d100``
     - 100-hr dead fuel moisture [fraction]
   * - ``moisture_lh``
     - Live herbaceous fuel moisture [fraction]
   * - ``moisture_lw``
     - Live woody fuel moisture [fraction]

When the Nelson (2000) diurnal EMC or FARSITE .fms spatial scenario file is
active, these fields reflect the spatially-varying per-cell values; otherwise
they contain the global Rothermel moisture values uniformly.

**Parameters**: None (automatic).

**Test**: ``regtest/moisture/spatial_moisture_output/``.

Multiple Weather Stations with Spatial IDW Interpolation
----------------------------------------------------------

**Header**: ``src/multi_wtr_weather.H``

A new ``multi_wtr_file`` input enables loading multiple FARSITE .wtr weather
station files and producing spatially-varying wind, temperature, and relative
humidity via inverse-distance-weighting (IDW).

Station list CSV format (``multi_wtr_file``)::

   # station_id, x_m, y_m, wtr_file
   1, 330000.0, 3775000.0, station1.wtr
   2, 335000.0, 3775000.0, station2.wtr
   3, 332500.0, 3780000.0, station3.wtr

At each timestep, the U and V wind components from each station are IDW-
interpolated to every grid cell, and the domain-mean T/RH from all stations
are used to update the global diurnal moisture model.

IDW formula:

.. math::

   V_{\text{cell}} = \frac{\sum_i w_i V_i}{\sum_i w_i},
   \quad w_i = d_i^{-p}

where :math:`d_i` is the distance from the cell to station :math:`i`,
:math:`p` is the IDW power exponent (``multi_wtr_idw_power``, default 2.0),
and :math:`V_i` is the wind component at station :math:`i` at the current time.

**Parameters**:

.. list-table::
   :header-rows: 1
   :widths: 28 50

   * - Parameter
     - Description
   * - ``multi_wtr_file``
     - Path to station list CSV (default: ``""`` = disabled)
   * - ``multi_wtr_idw_power``
     - IDW exponent :math:`p` (default: 2.0)

When ``multi_wtr_file`` is set, the diurnal moisture model is automatically
enabled and the domain-mean T/RH tracks the station-averaged IDW centroid.

Burn-Period Controls (Diurnal Active-Spread Window)
----------------------------------------------------

**Header**: ``src/parse_inputs.H`` / ``src/main.cpp``

In operational fire modelling, fire spread is often restricted to a daily
*burn period* when weather conditions (low humidity, warm temperatures) are
most conducive.  This mirrors the FARSITE / FSPro burn-period concept used in
operational forecasting.

When ``burn_period.enable = 1``, the computed rate-of-spread field (``R_mf``)
is zeroed to zero outside the specified local clock window
[``start_hour``, ``end_hour``), preventing any level-set, FARSITE, or MTT
advance during the inactive hours.  All other processes (moisture evolution,
spotting diagnostics, ecology metrics) continue normally.

The current local clock hour is computed as::

    clock_hour = (sim_start_hour + time_s / 3600) mod 24

where ``sim_start_hour`` defaults to ``solar_radiation.sim_start_hour`` when
solar radiation is enabled, or to ``burn_period.sim_start_hour`` otherwise.

Midnight-crossing windows (e.g. ``start_hour = 22``, ``end_hour = 6``)
are handled correctly: fire is active when ``clock_hour >= 22`` **or**
``clock_hour < 6``.

**Parameters**:

.. list-table::
   :header-rows: 1
   :widths: 30 50

   * - Parameter
     - Description
   * - ``burn_period.enable``
     - 1 to enable burn-period gating (default: 0)
   * - ``burn_period.start_hour``
     - Local hour (decimal) when fire becomes active (default: 10.0 = 10:00 AM)
   * - ``burn_period.end_hour``
     - Local hour (decimal) when fire becomes inactive (default: 20.0 = 8:00 PM)
   * - ``burn_period.sim_start_hour``
     - Local clock hour at simulation t=0 (default: inherited from
       ``solar_radiation.sim_start_hour`` when solar is enabled, else 0.0)

**Example** – active spread 10:00 AM to 8:00 PM::

    burn_period.enable     = 1
    burn_period.start_hour = 10.0
    burn_period.end_hour   = 20.0
    burn_period.sim_start_hour = 8.0   # simulation starts at 8 AM

**Example** – overnight window (active 9 PM to 7 AM)::

    burn_period.enable     = 1
    burn_period.start_hour = 21.0
    burn_period.end_hour   = 7.0
    burn_period.sim_start_hour = 0.0

