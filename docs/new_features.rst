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



Smoke Plume-Rise Model (Briggs 1965)
-------------------------------------

**Header**: ``src/smoke_plume_rise.H``

Computes the effective smoke plume-rise height Δh [m] above each fire-front
cell using the Briggs (1965/1969) buoyancy-dominated formula.  The result is
written to every plotfile as the variable ``plume_rise_m`` and the domain
maximum is printed in per-step log output.

**Theory**

1. Buoyancy flux per unit fire-line length:

   .. math::

      F_B = \frac{g \, I_B}{\pi \, \rho_a \, c_p \, T_a}  \quad [\text{m}^4/\text{s}^3]

   where :math:`I_B` is the Byram fireline intensity [W/m] (= 1000 × kW/m),
   :math:`g = 9.81\,\text{m/s}^2`, :math:`\rho_a` is air density [kg/m³],
   :math:`c_p` is specific heat of air [J/(kg·K)], and :math:`T_a` is ambient
   temperature [K].

2. Distance to final rise :math:`x_f = 3.5\,x^*`:

   .. math::

      x^* = \begin{cases}
        49\,F_B^{5/8}   & F_B < 55\,\text{m}^4/\text{s}^3 \\
        120\,F_B^{2/5}  & F_B \geq 55\,\text{m}^4/\text{s}^3
      \end{cases}

3. Final plume rise:

   .. math::

      \Delta h = 1.6 \, F_B^{1/3} \, x_f^{2/3} / u  \quad [\text{m}]

   where :math:`u` is wind speed [m/s] (clamped to 0.5 m/s to avoid
   division by zero).

**Input parameters** (prefix ``smoke_plume.``):

.. list-table::
   :header-rows: 1

   * - Parameter
     - Default
     - Description
   * - ``enable``
     - 0
     - 1 to compute plume rise (written to plotfile as ``plume_rise_m``)
   * - ``T_a``
     - 303.15
     - Ambient temperature [K] (default = 30 °C)
   * - ``rho_a``
     - 1.20
     - Air density [kg/m³]
   * - ``Cp_a``
     - 1005.0
     - Specific heat of air [J/(kg·K)]

**Example**::

    smoke_plume.enable = 1
    smoke_plume.T_a    = 308.0   # 35 °C hot day
    smoke_plume.rho_a  = 1.18    # slightly lower density at altitude
    smoke_plume.Cp_a   = 1005.0

**References**: Briggs (1965) JAPCA 15:433; Briggs (1969) USAEC TID-25075.


KML Perimeter Export
--------------------

**Function**: ``write_fire_perimeter_kml()`` in ``src/write_xy_data.H``

Writes the fire perimeter as a KML (Keyhole Markup Language) document that
can be opened directly in Google Earth or any GIS tool supporting KML.
Perimeter coordinates are converted from simulation UTM [m] to WGS-84 lon/lat
using the standard inverse Transverse Mercator formula (Redfearn 1948) when a
UTM zone number is specified.

A ``.kml`` file is written at every plotfile step and for every isochrone,
alongside the existing ``.csv`` and ``.geojson`` perimeter files.

**Input parameters**:

.. list-table::
   :header-rows: 1

   * - Parameter
     - Default
     - Description
   * - ``write_perimeter_kml``
     - 0
     - 1 to write ``perimeter_NNNN.kml`` at each plotfile step
   * - ``kml_utm_zone``
     - 0
     - UTM zone number 1–60; set 0 to write raw UTM coords (not WGS-84)
   * - ``kml_utm_northern``
     - 1
     - 1 for northern hemisphere, 0 for southern

**Example** – Southern California (UTM Zone 11N)::

    write_perimeter_csv     = 1
    write_perimeter_geojson = 1
    write_perimeter_kml     = 1
    kml_utm_zone            = 11
    kml_utm_northern        = 1


Simulation Date/Time Display in Logs and Report
------------------------------------------------

When a simulation start date is configured, each per-step log line includes
the calendar date and time alongside the simulation time in seconds:

.. code-block:: text

    Time:3600.0 (2024-08-15 09:00) with timestep:45.2

The simulation start date is taken from either:

* ``sim_datetime.year/month/day`` (explicit standalone fields), or
* ``solar_radiation.year/month/day`` when ``solar_radiation.enable = 1``
  (the two sets are automatically linked).

The start hour comes from ``solar_radiation.sim_start_hour`` (decimal hours).
If no start date is set, the plain "Time:T s" format is used as before.

The HTML fire report (``fire_report_file``) also shows start and end
calendar dates in the header.

**Input parameters** (prefix ``sim_datetime.``):

.. list-table::
   :header-rows: 1

   * - Parameter
     - Default
     - Description
   * - ``sim_datetime.year``
     - 0
     - Calendar year at t=0 (0 = inherit from ``solar_radiation``)
   * - ``sim_datetime.month``
     - 0
     - Calendar month at t=0 (1–12)
   * - ``sim_datetime.day``
     - 0
     - Calendar day at t=0 (1–31)

**Example**::

    sim_datetime.year  = 2024
    sim_datetime.month = 8
    sim_datetime.day   = 15
    # Start hour comes from solar_radiation.sim_start_hour (default 8.0)


Post-fire Fuel Adjustment for Re-entry Spots
---------------------------------------------

When firebrand spotting (basic or Albini) is active alongside per-cell fuel
depletion tracking (``fuel_depletion.enable = 1``), this option reduces the
catching probability for spot fires that land in previously-burned areas where
the residual fuel load is depleted.

The adjustment works as follows:

* For each spot-fire landing point, the local ``residual_fuel_mf`` value
  :math:`f_r \in [0,1]` is looked up (1 = fully loaded, 0 = completely burned).
* If :math:`f_r < \text{spotting\_fuel\_threshold}`, the spot is suppressed
  (no ignition regardless of P_catch).
* Otherwise, the effective catching probability is scaled:
  :math:`P_\text{eff} = P_\text{catch} \times f_r`.

This prevents re-burned spots from growing as fast as fresh-fuel ignitions
and is physically consistent with the exponential burnout model used to
compute ``residual_fuel_mf``.

**Input parameters** (prefix ``fuel_depletion.``):

.. list-table::
   :header-rows: 1

   * - Parameter
     - Default
     - Description
   * - ``adjust_spotting_reentry``
     - 0
     - 1 to scale P_catch by residual fuel for firebrand spots (requires ``enable=1``)
   * - ``spotting_fuel_threshold``
     - 0.05
     - Minimum residual fuel fraction below which spots cannot ignite [0–1]

**Example**::

    fuel_depletion.enable                  = 1
    fuel_depletion.tau_burnout             = 1800.0
    fuel_depletion.adjust_spotting_reentry = 1
    fuel_depletion.spotting_fuel_threshold = 0.05

    albini_spotting.enable = 1
    albini_spotting.P_catch = 0.8

FARSITE Vectorial Slope/Wind Combination
-----------------------------------------

**Headers**: ``src/compute_rothermel_R.H``, ``src/advection.H``

The standard Rothermel (1972) scalar ROS formula adds wind and slope factors
directly:

.. math::

   R = R_0 \,(1 + \varphi_w + \varphi_s)

This overestimates spread when wind and slope point in different directions
(e.g. a cross-slope wind or a downslope wind that opposes the slope factor).

FARSITE (Finney 1998) avoids this by treating the two factors as *vectors*
before combining them into a scalar magnitude:

.. math::

   \boldsymbol{\varphi}_w &= \varphi_w \,\hat{u}  \quad (\text{upwind unit vector})\\
   \boldsymbol{\varphi}_s &= \varphi_s \,\hat{n}  \quad (\text{upslope unit vector})\\
   \varphi_\text{combined} &= \left|\boldsymbol{\varphi}_w + \boldsymbol{\varphi}_s\right|\\
   R &= R_0 \,(1 + \varphi_\text{combined})

This is now available as an option with ``rothermel.use_slope_wind_vectors = 1``.

**Comparison of slope/wind interaction modes** (all mutually exclusive):

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Mode
     - Formula
   * - Default scalar additive (mode C)
     - :math:`R = R_0(1 + \varphi_w + \varphi_s)`
   * - Scalar cross-term (mode B)
     - :math:`R = R_0(1 + \varphi_w + \varphi_s + k\,\varphi_w\varphi_s)`
   * - FARSITE vectorial (mode A)
     - :math:`R = R_0(1 + |\varphi_w\hat{u} + \varphi_s\hat{n}|)`

Modes A (vectorial) and B (cross-term) require ``rothermel.terrain_file`` (or a landscape file) to
supply per-cell slope vectors.  When terrain is not provided the vectorial
mode falls back to scalar additive.  Modes B and C (cross-term) and A
(vectorial) cannot be active simultaneously; setting both raises an error at
startup.

**Input parameters**:

.. list-table::
   :header-rows: 1

   * - Parameter
     - Default
     - Description
   * - ``rothermel.use_slope_wind_vectors``
     - 0
     - 1 to enable FARSITE-style vectorial slope/wind combination
   * - ``rothermel.use_slope_wind_cross``
     - 0
     - 1 to enable :math:`k\,\varphi_w\varphi_s` cross-term (existing option)
   * - ``rothermel.k_slope_wind_cross``
     - 1.0
     - Cross-term coupling coefficient :math:`k` (only used when ``use_slope_wind_cross = 1``)

**Example — vectorial mode (FARSITE-parity)**::

    rothermel.terrain_file          = terrain.xyz
    rothermel.use_slope_wind_vectors = 1

**Example — cross-term with custom coupling**::

    rothermel.terrain_file         = terrain.xyz
    rothermel.use_slope_wind_cross = 1
    rothermel.k_slope_wind_cross   = 0.5   # partial coupling

**References**: Finney (1998) FARSITE RMRS-RP-4; Rothermel (1983) GTR INT-143;
Anderson (1982) Res. Note INT-328.

Burn-Period (Daytime Burning Window) Gate
------------------------------------------

**Headers**: ``src/parse_inputs.H``, ``src/main.cpp``

The rate-of-spread field ``R_mf`` is zeroed outside the configured local clock
window ``[burn_period.start_hour, burn_period.end_hour)``, pausing all spread
paths (level-set advection, FARSITE ellipse, MTT) during inactive hours.
Moisture evolution, spotting diagnostics, and all other sub-models continue
normally.  This mirrors the FARSITE / FSPro burn-period concept used in
operational forecasting.

The gate is applied as a single ``R_mf.setVal(0.0)`` call immediately after
the retardant-suppression step, so no changes are required to any of the
individual spread-path implementations.

**Input parameters** (prefix ``burn_period.``):

.. list-table::
   :header-rows: 1

   * - Parameter
     - Default
     - Description
   * - ``enable``
     - 0
     - 1 to activate the burn-period gate
   * - ``start_hour``
     - 10.0
     - Local decimal hour when spread becomes active (inclusive)
   * - ``end_hour``
     - 20.0
     - Local decimal hour when spread becomes inactive (exclusive)
   * - ``sim_start_hour``
     - 0.0
     - Local clock hour at simulation :math:`t = 0` (inherits from
       ``solar_radiation.sim_start_hour`` when solar radiation is enabled)

Midnight-crossing windows (e.g. ``start_hour = 22``, ``end_hour = 6``) are
handled correctly.

**Example**::

    burn_period.enable         = 1
    burn_period.start_hour     = 10.0
    burn_period.end_hour       = 20.0
    burn_period.sim_start_hour = 8.0   # simulation starts at 08:00

Post-frontal Fuel Consumption Raster Output
--------------------------------------------

**Header**: ``src/main.cpp``; GIS export: ``tools/plotfile_to_geotiff.py``

Two complementary rasters quantify post-frontal fuel state:

* **``fuel_consumption``** — instantaneous bulk fuel consumption fraction
  :math:`f_c \in [0, 1]` computed during the FARSITE spread step when
  ``farsite.use_bulk_fuel_consumption = 1``.  Based on fire intensity and
  residence time (Albini 1976 / Rothermel 1983 per-class exponential burnout).

* **``residual_fuel``** — fraction of fuel remaining behind the fire front.
  Updated every timestep for burned cells via exponential decay:

  .. math::

     f_r(t) = \exp\!\left(-\frac{t - t_\text{arrive}}{\tau_\text{burnout}}\right)

  Requires ``fuel_depletion.enable = 1`` and ``fuel_depletion.tau_burnout``
  [s].  :math:`f_r = 1` (unburned), :math:`f_r \to 0` (fully consumed).

Both variables are written to every AMReX plotfile and are exported as
GeoTIFF by ``plotfile_to_geotiff.py`` (now included in the default
``FIRE_VARS`` export set alongside ``phi``, ``fireline_intensity``, etc.).

**Input parameters**:

.. list-table::
   :header-rows: 1

   * - Parameter
     - Default
     - Description
   * - ``farsite.use_bulk_fuel_consumption``
     - 0
     - 1 to compute ``fuel_consumption`` during FARSITE spread
   * - ``fuel_depletion.enable``
     - 0
     - 1 to track exponential post-frontal burnout (``residual_fuel``)
   * - ``fuel_depletion.tau_burnout``
     - 3600.0
     - Characteristic burnout time :math:`\tau` [s]

**Example**::

    farsite.use_bulk_fuel_consumption = 1
    fuel_depletion.enable             = 1
    fuel_depletion.tau_burnout        = 1800.0

**GIS export**::

    python3 tools/plotfile_to_geotiff.py plt0200 \
        -v residual_fuel fuel_consumption arrival_time \
        --utm-origin 500000 3700000 --epsg 32611 --outdir gis_out
