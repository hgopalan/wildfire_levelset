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
interpolated to every grid cell. Temperature and relative humidity are also
spatially interpolated using IDW when ``diurnal_moisture.enable = 1``, creating
spatially-varying fuel moisture content across the domain. When the solar
radiation model is enabled (``solar_radiation.enable = 1``), the shade-adjusted
EMC computation uses the per-cell interpolated T and RH values instead of the
diurnal sinusoidal model.

IDW formula:

.. math::

   V_{\text{cell}} = \frac{\sum_i w_i V_i}{\sum_i w_i},
   \quad w_i = d_i^{-p}

where :math:`d_i` is the distance from the cell to station :math:`i`,
:math:`p` is the IDW power exponent (``multi_wtr_idw_power``, default 2.0),
and :math:`V_i` is the value (wind component, temperature, or RH) at station
:math:`i` at the current time.

**Spatial Interpolation Details**:

The following fields are spatially interpolated via IDW:

- **Wind U and V components**: Interpolated in Cartesian components then
  converted to speed and direction.
- **Temperature [°C]**: Spatially interpolated as a scalar field.
- **Relative Humidity [%]**: Spatially interpolated as a scalar field.

When ``diurnal_moisture.enable = 1`` is set, the spatially-interpolated T and RH
values are used to compute per-cell equilibrium moisture content (EMC) via the
Nelson (2000) model, creating spatially-varying fuel moisture across the domain.

When combined with ``solar_radiation.enable = 1``, the shade-adjusted EMC
computation uses the per-cell T/RH MultiFabs instead of the global diurnal
sinusoidal model, enabling realistic spatial moisture gradients from weather
station data.

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

When ``multi_wtr_file`` is set:

- The diurnal moisture model is automatically enabled.
- Spatial T and RH interpolation is performed each timestep via
  ``apply_multi_wtr_TRH_to_spatial()``.
- The domain-mean T/RH (simple average of all stations) is also computed as a
  fallback for any global moisture parameters that may still use scalar values.

**Functions**:

- ``load_multi_wtr_stations()`` — Parses station list and loads .wtr files
- ``apply_multi_wtr_to_vel()`` — IDW interpolation of wind to velocity MultiFab
- ``apply_multi_wtr_TRH_to_spatial()`` — IDW interpolation of T/RH to MultiFabs
- ``get_domain_TRH_at_time()`` — Domain-mean T/RH (fallback for scalar use)

**Test**: ``regtest/moisture/multi_wtr_spatial/``

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

Fire Size and Growth Rate Statistics
-------------------------------------

**Header**: ``src/write_xy_data.H``; **Tool**: ``tools/fire_size_summary.py``

The fire statistics CSV output (``fire_stats_file``) has been extended with
additional growth metrics and fire ellipse geometry estimates:

**New CSV columns**:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Column
     - Description
   * - ``growth_rate_ha_min``
     - Instantaneous fire growth rate [ha/min]; (current_area - previous_area) / time_delta
   * - ``major_axis_m``
     - Fire ellipse semi-major axis [m]; estimated from perimeter and area using Ramanujan's ellipse approximation
   * - ``minor_axis_m``
     - Fire ellipse semi-minor axis [m]; estimated from perimeter and area

These metrics characterize fire expansion speed and shape evolution over the
simulation duration.  The ellipse axes are computed from the perimeter and
burned area by solving the simultaneous equations:

.. math::

   A = \pi a b, \quad P \approx \pi(a + b)

where :math:`a, b` are the semi-major and semi-minor axes respectively.

**Processing tool updates**:

The ``fire_size_summary.py`` tool now:

* Displays growth rate and ellipse axes in ASCII summary table
* Includes them in percentile statistics (10th, 50th, 90th percentiles)
* Generates additional plots: ``fire_growth_rate.{fmt}`` and ``fire_ellipse_axes.{fmt}``

**Usage**::

    python3 tools/fire_size_summary.py fire_stats.csv --plot

This generates updated plots including fire growth rate dynamics and ellipse
geometry evolution.  Growth rate can be negative during controlled burns or when
fire front contracts (e.g., due to retardant or suppression).

**CSV format example**::

    step,time_s,...,growth_rate_ha_min,major_axis_m,minor_axis_m
    0,0.0,...,0.0,0.0,0.0
    20,200.0,...,0.5,125.3,85.2
    40,400.0,...,0.45,180.6,120.1

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


Scott & Reinhardt (2001) Crown Fire Surface Area (CFSA)
--------------------------------------------------------

**Header**: ``src/scott_reinhardt_cfsa.H``

The Crown Fire Surface Area (CFSA) model computes the effective burning surface
area per unit ground area during active crown fire, accounting for 3-D canopy
structure. This provides improved crown fire heat release and emissions estimates.

**Model equation:**

.. math::

   \mathrm{CFSA} = \min(\mathrm{CBD} \times (\mathrm{CH} - \mathrm{CBH}) \times k_{sa},\, \mathrm{CFSA}_{\max})

where:

* :math:`\mathrm{CBD}` = canopy bulk density [kg/m³]
* :math:`\mathrm{CH}` = canopy height [m]
* :math:`\mathrm{CBH}` = canopy base height [m]
* :math:`k_{sa}` = surface area coefficient [m²/kg] (typical: 4-8)
* :math:`\mathrm{CFSA}_{\max}` = maximum CFSA [dimensionless] (typical: 2-4)

**Functions:**

* ``compute_crown_fire_surface_area(CBD, CH, CBH, k_sa, max_cfsa)`` — returns CFSA
* ``compute_crown_heat_release_with_cfsa(cfsa, w_crown, h_crown, tau_crown)`` — heat release [kW/m²]
* ``estimate_crown_fuel_loading_from_cbd(CBD, canopy_depth)`` — crown fuel loading [kg/m²]

**Input parameters** (prefix ``crown.``):

.. list-table::
   :header-rows: 1

   * - Parameter
     - Default
     - Description
   * - ``cfsa_k_sa``
     - 6.0
     - Surface area coefficient [m²/kg]
   * - ``cfsa_max``
     - 3.0
     - Maximum CFSA [dimensionless]

**Applications**: Crown fire heat release calculations, emissions modeling,
radiation heat flux, plume rise and smoke production.

**Test**: ``regtest/crown_fire/scott_reinhardt_cfsa/``

**Reference**: Scott, J.H. & Reinhardt, E.D. (2001). *Assessing Crown Fire Potential
by Linking Models of Surface and Crown Fire Behavior.* USDA Forest Service
Research Paper RMRS-RP-29.


Multi-Layer Canopy Wind Profile
--------------------------------

**Header**: ``src/multi_layer_canopy_wind.H``

Implements vertical wind profile through multi-layer canopy following exponential
attenuation (Massman 1997, Inoue 1963). Provides realistic wind speed variation
with height for improved surface vs. crown fire modeling.

**Within-canopy exponential profile** (z < h):

.. math::

   u(z) = u_h \times \exp\!\left(-\alpha \left(1 - \frac{z}{h}\right)\right)

**Above-canopy logarithmic profile** (z > h):

.. math::

   u(z) = u_h \times \frac{\ln\!\left(\frac{z-d}{z_0}\right)}{\ln\!\left(\frac{h-d}{z_0}\right)}

where :math:`\alpha` is the attenuation coefficient related to LAI via
:math:`\alpha \approx 2.5 \times \mathrm{LAI}^{0.5}`.

**Functions:**

* ``compute_canopy_wind_attenuation_coefficient(LAI)`` — returns α from LAI
* ``compute_wind_speed_at_height_within_canopy(u_h, z, h, alpha)`` — wind at height z
* ``compute_wind_speed_at_height_above_canopy(u_h, z, h, d_ratio, z0_ratio)`` — wind above canopy
* ``compute_multi_layer_wind_profile(...)`` — computes wind at N vertical layers
* ``compute_effective_wind_speed_for_layer(...)`` — wind for surface/crown fire layers

**Input parameters** (prefix ``canopy_wind.``):

.. list-table::
   :header-rows: 1

   * - Parameter
     - Default
     - Description
   * - ``enable``
     - 0
     - 1 to activate multi-layer canopy wind
   * - ``n_layers``
     - 5
     - Number of vertical layers
   * - ``alpha``
     - 2.5
     - Attenuation coefficient [dimensionless]
   * - ``LAI``
     - 4.0
     - Leaf area index [m²/m²]
   * - ``d_ratio``
     - 0.7
     - Displacement height ratio d/h
   * - ``z0_ratio``
     - 0.1
     - Roughness length ratio z₀/h
   * - ``z_ref``
     - 10.0
     - Reference height for input wind [m]

**Test**: ``regtest/wind/multi_layer_canopy/``

**References**: Massman, W.J. (1997). *An analytical one-dimensional model of
momentum transfer by vegetation of arbitrary structure.* Boundary-Layer
Meteorology, 83(3), 407-421. Inoue, E. (1963). *On the turbulent structure of
airflow within crop canopies.* Journal of the Meteorological Society of Japan, 41, 317-326.


Fine Fuel Moisture Time-Lag Differential Equations
---------------------------------------------------

**Header**: ``src/fuel_moisture_timelag_de.H``

Implements physically-based fuel moisture dynamics using time-lag differential
equations (Nelson 2000, Viney 1991). Provides continuous moisture evolution with
hysteresis and temperature correction.

**Governing equation:**

.. math::

   \frac{dM}{dt} = \frac{M_e - M}{\tau_{\text{eff}}} + P(t)

where :math:`M` is fuel moisture [fraction], :math:`M_e` is equilibrium moisture
content (EMC), :math:`\tau_{\text{eff}}` is the effective time-lag [hours], and
:math:`P(t)` is precipitation wetting [1/hour].

**Equilibrium moisture with hysteresis** (Nelson 2000):

* **Adsorption** (wetting): :math:`M_e = 0.03229 + 0.2810 H + 0.4093 H^2 - 1.3560 H^3 + 1.6596 H^4`
* **Desorption** (drying): :math:`M_e = 0.05800 + 0.1985 H + 0.6250 H^2 - 1.1830 H^3 + 1.0570 H^4`

where :math:`H = \mathrm{RH}/100`.

**Temperature correction:**

.. math::

   \tau_{\text{eff}} = \tau \times \exp\!\left(-0.015 (T - 20)\right)

Higher temperature accelerates drying by reducing effective time lag.

**Functions:**

* ``compute_emc_adsorption(RH)`` — EMC for wetting [fraction]
* ``compute_emc_desorption(RH)`` — EMC for drying [fraction]
* ``compute_emc_with_hysteresis(RH, M_current)`` — EMC with hysteresis
* ``compute_temperature_correction_factor(T)`` — temperature factor
* ``compute_precipitation_wetting_rate(rain_rate)`` — wetting term P(t)
* ``update_fuel_moisture_timelag_de(M, RH, T, rain, dt, tau)`` — forward Euler update
* ``update_multi_class_fuel_moisture(...)`` — update 1hr/10hr/100hr/1000hr simultaneously

**Input parameters** (prefix ``fuel_moisture_de.``):

.. list-table::
   :header-rows: 1

   * - Parameter
     - Default
     - Description
   * - ``enable``
     - 0
     - 1 to activate time-lag DE model
   * - ``method``
     - ``"timelag_de"``
     - Solver method
   * - ``T_ambient``
     - 25.0
     - Ambient temperature [°C]
   * - ``RH``
     - 40.0
     - Relative humidity [%]
   * - ``tau_1hr``
     - 1.0
     - 1-hour time lag [hours]
   * - ``tau_10hr``
     - 10.0
     - 10-hour time lag [hours]
   * - ``tau_100hr``
     - 100.0
     - 100-hour time lag [hours]
   * - ``tau_1000hr``
     - 1000.0
     - 1000-hour time lag [hours]
   * - ``use_hysteresis``
     - 1
     - Enable adsorption/desorption hysteresis
   * - ``temp_correction``
     - 1
     - Enable temperature-dependent drying
   * - ``M_1hr_init``
     - 0.12
     - Initial 1-hr moisture [fraction]
   * - ``M_10hr_init``
     - 0.15
     - Initial 10-hr moisture [fraction]
   * - ``M_100hr_init``
     - 0.18
     - Initial 100-hr moisture [fraction]
   * - ``rain_rate``
     - 0.0
     - Precipitation rate [mm/h]

**Test**: ``regtest/moisture/fuel_moisture_timelag_de/``

**References**:

* Nelson, R.M. (2000). *Prediction of diurnal change in 10-h fuel stick moisture content.* Canadian Journal of Forest Research, 30, 1071-1087.
* Viney, N.R. (1991). *A review of fine fuel moisture modelling.* International Journal of Wildland Fire, 1(4), 215-234.


Integrated Fire Behavior Features
==================================

This section describes advanced fire behavior features integrated into the AMReX-based level-set wildfire framework with GPU support.

Fuel Continuity Factor
----------------------

**Header**: ``src/rothermel_model.H``

Accounts for patchy or discontinuous fuel beds by applying a multiplicative factor (0-1) to computed rate of spread (ROS).

**Input Parameters:**

.. list-table::
   :header-rows: 1

   * - Parameter
     - Default
     - Description
   * - ``rothermel.fuel_continuity``
     - 1.0
     - Fuel continuity factor: 1.0 = continuous, 0.5 = 50% coverage, 0.0 = gaps

**References**: Finney, M.A. (2006). *An overview of FlamMap fire modeling capabilities.* USDA Forest Service, Rocky Mountain Research Station, RMRS-P-41, 213-220.


NFDRS Fire Danger Class
-----------------------

**Header**: ``src/fire_intensity_class.H``

Classifies fireline intensity into 5 operational fire danger categories (Low/Moderate/High/Very High/Extreme) following National Fire Danger Rating System (NFDRS) conventions.

**Implementation:**

* Automatically computed from fireline intensity
* Output in plotfile as ``nfdrs_danger_class``

**Danger Classes:**

* 1 = Low (< 100 kW/m)
* 2 = Moderate (100-500 kW/m)
* 3 = High (500-2000 kW/m)
* 4 = Very High (2000-5000 kW/m)
* 5 = Extreme (> 5000 kW/m)

**References**: Deeming, J.E., et al. (1977). *National Fire Danger Rating System - 1978.* USDA GTR INT-39.


Crown Fraction Burned (CFB) Diagnostic
--------------------------------------

**Header**: ``src/crown_initiation.H``

Diagnostic metric distinguishing passive vs active crown fire based on fire intensity.

**Input Parameters:**

.. list-table::
   :header-rows: 1

   * - Parameter
     - Default
     - Description
   * - ``crown_fraction.enable``
     - 0
     - Set to 1 to enable CFB diagnostic

**Output Variables:**

* ``crown_fraction_burned`` - CFB ratio (0-1)

**References**: Scott, J.H. & Reinhardt, E.D. (2001). *Assessing canopy fire behavior in a wildland fire information system.* USDA RMRS-RP-29.


Effective Wind Speed
--------------------

**Header**: ``src/effective_wind_speed.H``

Combines ambient wind and slope effects into a single scalar effective wind speed via vector addition.

**Formula:**

.. math::

   U_{\text{eff}} = \sqrt{U_{\text{wind}}^2 + U_{\text{slope\_equiv}}^2}

**Input Parameters:**

.. list-table::
   :header-rows: 1

   * - Parameter
     - Default
     - Description
   * - ``effective_wind.enable``
     - 0
     - Set to 1 to enable effective wind speed calculation

**Output Variables:**

* ``effective_wind_speed`` - Combined wind speed [m/s]

**Requirements:** Terrain slopes must be available (landscape file)

**References**: Rothermel, R.C. (1983). *How to predict the spread and intensity of forest and range fires.* GTR INT-143.


Thomas Flame Length Model
-------------------------

Alternative flame length formula :math:`L = 0.0266 \times I^{0.667}` (vs default Byram: :math:`L = 0.0775 \times I^{0.46}`)

**Input Parameters:**

.. list-table::
   :header-rows: 1

   * - Parameter
     - Default
     - Description
   * - ``flame_length_model.model``
     - ``"byram"``
     - Flame length model: ``"byram"`` or ``"thomas"``

**Selection:**

* ``"byram"`` - Default Byram (1959) formula
* ``"thomas"`` - Thomas (1963) formula

**References**: Thomas, P.H. (1963). *Size of flames from natural fires.* Ninth Symposium (International) on Combustion, 9(1):844-859.


Fuel Boundary Smoothing
-----------------------

**Header**: ``src/fuel_boundary_smoothing.H``

Distance-weighted ROS blending at fuel boundaries to eliminate unrealistic discontinuities.

**Input Parameters:**

.. list-table::
   :header-rows: 1

   * - Parameter
     - Default
     - Description
   * - ``fuel_boundary.enable``
     - 0
     - Set to 1 to enable fuel boundary smoothing
   * - ``fuel_boundary.transition_cells``
     - 2.0
     - Number of cells for transition zone

**Requirements:** Landscape file with fuel model data must be present

**References**: Finney, M.A. (1998, 2006). *FARSITE: Fire Area Simulator.* USDA RMRS-RP-4, RMRS-P-41.


CSIRO Grassfire Acceleration
----------------------------

**Header**: ``src/fire_acceleration.H``

Models non-equilibrium fire growth during initial spread phase using exponential acceleration function.

**Input Parameters:**

.. list-table::
   :header-rows: 1

   * - Parameter
     - Default
     - Description
   * - ``grassfire_accel.enable``
     - 0
     - Set to 1 to enable grassfire acceleration
   * - ``grassfire_accel.t_accel``
     - 600.0
     - Acceleration time constant [seconds]

**Formula:**

.. math::

   \text{acceleration\_factor} &= 1 - \exp(-dt / t_{\text{accel}}) \\
   R_{\text{accelerated}} &= R \times (1 + \text{acceleration\_factor})

**References**: Cheney, N.P. & Gould, J.S. (1995). *Fire growth in grassland fuels.* International Journal of Wildland Fire, 5(4):237-247.


Burnout Time Separation
-----------------------

Splits total residence time into flaming and smoldering phase durations by fuel type.

**Input Parameters:**

.. list-table::
   :header-rows: 1

   * - Parameter
     - Default
     - Description
   * - ``burnout_separation.enable``
     - 1
     - Enable burnout time separation
   * - ``burnout_separation.flaming_fraction_fine``
     - 0.70
     - Fine fuels: 70% flaming
   * - ``burnout_separation.flaming_fraction_medium``
     - 0.40
     - Medium: 40% flaming
   * - ``burnout_separation.flaming_fraction_heavy``
     - 0.20
     - Heavy: 20% flaming
   * - ``burnout_separation.flaming_fraction_duff``
     - 0.10
     - Duff: 10% flaming

**Output Variables:**

* ``burnout_flaming_time`` - Flaming phase duration [seconds]
* ``burnout_smoldering_time`` - Smoldering phase duration [seconds]

**References**:

* Anderson, H.E. (1969). *Sustained burning of piled coarse forest fuels.* USDA Research Paper INT-69.
* Frandsen, W.H. (1997). *Ignition probability of organic soils.* Canadian Journal of Forest Research, 27(9):1471-1477.


Simard Moisture Model
---------------------

**Header**: ``src/simard_moisture.H``

Exponential time-lag moisture update based on equilibrium moisture approach with size-dependent time constants.

**Input Parameters:**

.. list-table::
   :header-rows: 1

   * - Parameter
     - Default
     - Description
   * - ``simard_moisture.enable``
     - 0
     - Enable Simard moisture model
   * - ``simard_moisture.tau_1hr``
     - 1.0
     - 1-hour fuel lag [hours]
   * - ``simard_moisture.tau_10hr``
     - 10.0
     - 10-hour fuel lag [hours]
   * - ``simard_moisture.tau_100hr``
     - 100.0
     - 100-hour fuel lag [hours]

**Formula:**

.. math::

   M(t+dt) = M_{\text{eq}} + (M(t) - M_{\text{eq}}) \times \exp(-dt/\tau)

**References**: Simard, A.J. (1968). *The moisture content of forest fuels - I.* Forestry Branch Info. Rep. FF-X-14.


Post-Frontal Smoldering
-----------------------

**Header**: ``src/duff_moisture_smoldering.H``

Tracks residual combustion after flame front passage with exponential decay for smoke/air quality applications.

**Input Parameters:**

.. list-table::
   :header-rows: 1

   * - Parameter
     - Default
     - Description
   * - ``post_frontal.enable``
     - 0
     - Enable post-frontal smoldering
   * - ``post_frontal.tau_fine``
     - 1800.0
     - Fine fuels decay [seconds] = 30 min
   * - ``post_frontal.tau_medium``
     - 3600.0
     - Medium fuels decay = 1 hour
   * - ``post_frontal.tau_heavy``
     - 7200.0
     - Heavy fuels decay = 2 hours
   * - ``post_frontal.tau_duff``
     - 21600.0
     - Duff decay = 6 hours

**Output Variables:**

* ``time_since_burn`` - Elapsed time since cell burned [seconds]
* ``residual_heat_release`` - Residual smoldering intensity [kW/m²]

**References**: Frandsen, W.H. (1997). *Ignition probability of organic soils.* Canadian Journal of Forest Research, 27(9):1471-1477. Urbanski, S.P. (2014). *Wildland fire emissions, carbon, and climate: Emission factors.* Forest Ecology and Management, 317:1-8.


Integration Architecture
------------------------

**Files Modified/Created:**

1. **parse_inputs.H/cpp** - Added parameter structures for all features
2. **multifab_setup.H** - Added MultiFab fields for diagnostic output:

   * ``crown_fraction_burned_mf``
   * ``effective_wind_speed_mf``
   * ``burnout_phases_mf``
   * ``residual_heat_release_mf``
   * ``time_since_burn_mf``

3. **plot_results.H** - Added new diagnostic variables to plotfile output
4. **main.cpp** - Integrated feature calls into main simulation loop
5. **wildfire_includes.H** - Added includes for feature header files

**New Header Files:**

* ``effective_wind_speed.H`` - Effective wind speed computation
* ``fire_intensity_class.H`` - NFDRS danger classification
* ``crown_initiation.H`` - CFB computation
* ``fuel_boundary_smoothing.H`` - Boundary smoothing
* ``fire_acceleration.H`` - CSIRO acceleration
* ``simard_moisture.H`` - Simard moisture model
* ``duff_moisture_smoldering.H`` - Post-frontal tracking

**Performance Notes:**

* **GPU Acceleration:** All features use AMREX_GPU_HOST_DEVICE macros for GPU compatibility
* **Computational Cost:** Features are only computed when explicitly enabled
* **Memory Overhead:** Additional MultiFabs for diagnostics add ~50 MB per feature for 512³ domain
* **Physics-Based:** All formulas based on peer-reviewed wildfire literature

Easy-to-Implement Operational Fire Danger Indices (2026)
=========================================================

This section documents additional fire danger and fire behavior indices added in 2026
based on operational wildfire management tools. These features require minimal code changes
and provide high operational value for fire danger rating and decision support.

Keetch-Byram Drought Index (KBDI)
-----------------------------------

**Header**: ``src/kbdi_index.H``

The Keetch-Byram Drought Index is a widely-used measure of soil moisture deficit
and fire potential, particularly in the southeastern United States. The KBDI ranges
from 0 (no moisture deficit) to 800 (extreme drought).

**Model Formulation**

The KBDI represents the net effect of evapotranspiration and precipitation on soil
moisture depletion. Daily increment (Keetch & Byram 1968, Janis et al. 2002):

.. math::

   dQ = (800 - Q) \times DF \times ET \times 10^{-3}

where:

.. math::

   ET = 0.968 \times \exp(0.0875 \times T_{\max} + 1.5552) - 8.3

   DF = \frac{1}{1 + 10.88 \times \exp(-0.001736 \times R_A)}

* :math:`Q` = current KBDI value (0-800)
* :math:`T_{\max}` = daily maximum temperature [°C]
* :math:`R_A` = mean annual precipitation [mm] (default: 1000 mm)
* :math:`ET` = evapotranspiration factor (temperature-dependent)
* :math:`DF` = drought factor (precipitation-dependent)

Precipitation reduces KBDI:

.. math::

   Q_{\text{new}} = \max(0, Q_{\text{old}} - (P - 5))  \quad \text{if } P > 5 \text{ mm}

where :math:`P` is daily precipitation [mm].

**KBDI Interpretation**

.. list-table::
   :header-rows: 1

   * - KBDI Range
     - Soil Condition
     - Fire Potential
   * - 0-200
     - Moist soil and duff layers
     - Low
   * - 200-400
     - Soil and duff drying
     - Increasing
   * - 400-600
     - Significant moisture depletion
     - High
   * - 600-800
     - Severe drought
     - Extreme

**Drought Factor Conversion**

The KBDI can be converted to the McArthur Drought Factor (DF) used in the
Australian Forest Fire Danger Index (Noble et al. 1980):

.. math::

   DF = 10.5 \times [1 - \exp(-Q/800)] \times [1 - 0.25 \times \exp(-Q/200)]

This provides a 0-10 drought factor scale compatible with ``mcarthur_ffdi.H``.

**API Functions**

Scalar KBDI update (daily):

.. code-block:: cpp

   KBDIState kbdi_state;
   kbdi_state.Q = 100.0;  // Initial KBDI
   update_kbdi_scalar(kbdi_state, T_max, precip_mm, annual_precip);

Spatially-varying KBDI (MultiFab):

.. code-block:: cpp

   update_kbdi_field(kbdi_mf, temperature_mf, precip_mf, annual_precip);

Drought factor conversion:

.. code-block:: cpp

   Real drought_factor = compute_drought_factor_from_kbdi(kbdi_state.Q);

**Input Parameters**

KBDI should be updated once per day (not every timestep). Typical usage in a
daily weather update loop:

.. code-block:: text

   kbdi.enable = 1                # Enable KBDI calculation
   kbdi.initial_value = 100.0     # Initial KBDI (0-800)
   kbdi.annual_precip_mm = 1000.0 # Mean annual precipitation [mm]

**References:**

* Keetch, J.J. & Byram, G.M. (1968). "A drought index for forest fire control."
  USDA Forest Service Research Paper SE-38.
* Janis, M.J., et al. (2002). "Monitoring the Keetch-Byram Drought Index."
  International Journal of Wildland Fire, 11:217-222.
* Noble, I.R., Bary, G.A.V., & Gill, A.M. (1980). "McArthur's fire-danger
  meters expressed as equations." Australian Journal of Ecology, 5(2):201-203.

Haines Index (Lower Atmospheric Severity Index)
------------------------------------------------

**Header**: ``src/haines_index.H``

The Haines Index indicates the potential for large plume-dominated wildfires
and erratic fire behavior based on atmospheric stability and moisture content.
Values range from 2 (very low potential) to 6 (high potential for extreme behavior).

**Model Formulation**

The Haines Index combines two components:

.. math::

   HI = A + B

where :math:`A` = stability term (1-3) and :math:`B` = moisture term (1-3).

**Elevation Variants**

Three variants exist for different terrain elevations:

.. list-table::
   :header-rows: 1
   :widths: 15 20 35 30

   * - Variant
     - Elevation
     - Stability (A)
     - Moisture (B)
   * - LOW
     - < 3000 ft
     - :math:`\Delta T` = T₉₅₀ - T₈₅₀
     - T-Td at 850 mb
   * - MID
     - 3000-7000 ft
     - :math:`\Delta T` = T₈₅₀ - T₇₀₀
     - T-Td at 850 mb
   * - HIGH
     - > 7000 ft
     - :math:`\Delta T` = T₇₀₀ - T₅₀₀
     - T-Td at 700 mb

**Stability Term (A)**

.. list-table::
   :header-rows: 1

   * - Variant
     - A = 1 (Stable)
     - A = 2 (Moderate)
     - A = 3 (Unstable)
   * - LOW
     - ΔT < 4°C
     - 4°C ≤ ΔT < 8°C
     - ΔT ≥ 8°C
   * - MID
     - ΔT < 6°C
     - 6°C ≤ ΔT < 11°C
     - ΔT ≥ 11°C
   * - HIGH
     - ΔT < 18°C
     - 18°C ≤ ΔT < 22°C
     - ΔT ≥ 22°C

**Moisture Term (B)**

.. list-table::
   :header-rows: 1

   * - Variant
     - B = 1 (Moist)
     - B = 2 (Moderate)
     - B = 3 (Dry)
   * - LOW
     - T-Td < 6°C
     - 6°C ≤ T-Td < 10°C
     - T-Td ≥ 10°C
   * - MID
     - T-Td < 6°C
     - 6°C ≤ T-Td < 13°C
     - T-Td ≥ 13°C
   * - HIGH
     - T-Td < 15°C
     - 15°C ≤ T-Td < 21°C
     - T-Td ≥ 21°C

**Haines Index Interpretation**

.. list-table::
   :header-rows: 1

   * - HI Value
     - Fire Potential
   * - 2-3
     - Very low
   * - 4
     - Low
   * - 5
     - Moderate
   * - 6
     - High (extreme fire behavior likely)

**API Functions**

From upper-air sounding data:

.. code-block:: cpp

   int HI = compute_haines_index(T_lower, T_upper, Td_lower, HainesVariant::MID);

From surface observations (approximation):

.. code-block:: cpp

   int HI = compute_haines_index_surface(T_surface, RH_surface, elevation_m);

Spatially-varying Haines Index:

.. code-block:: cpp

   compute_haines_field(haines_mf, T_lower_mf, T_upper_mf, Td_lower_mf, 
                        HainesVariant::MID);

**Input Parameters**

When upper-air data is available:

.. code-block:: text

   haines.enable = 1
   haines.variant = MID        # LOW, MID, or HIGH
   haines.T_lower = 15.0       # Temperature at lower level [°C]
   haines.T_upper = 5.0        # Temperature at upper level [°C]
   haines.Td_lower = 8.0       # Dewpoint at lower level [°C]

When only surface data is available:

.. code-block:: text

   haines.enable = 1
   haines.use_surface_approximation = 1
   # Uses existing temperature and humidity fields

**References:**

* Haines, D.A. (1988). "A lower atmosphere severity index for wildland fires."
  National Weather Digest, 13:23-27.
* Werth, P.A. & Ochoa, R. (1993). "The evaluation of Idaho wildfire growth
  using the Haines Index." Weather and Forecasting, 8(2):223-234.
* Mills, G.A. & McCaw, L. (2010). "Atmospheric stability environments and
  fire weather in Australia." CAWCR Technical Report No. 20.

Backing Fire Rate of Spread
----------------------------

**Header**: ``src/backing_fire_ros.H``

The backing fire is the portion of the fire perimeter spreading against the wind
direction. It exhibits significantly lower rate of spread (ROS) than the heading
fire, typically 10-40% depending on wind speed and fuel conditions.

**Empirical Backing Fire Ratios**

The backing fire ROS ratio :math:`c/a` (backing coefficient / heading coefficient)
decreases with increasing wind speed:

.. list-table::
   :header-rows: 1

   * - Wind Speed
     - Backing Ratio (c/a)
     - Typical Range
   * - Light (< 5 mph)
     - 0.30-0.40
     - 30-40% of heading ROS
   * - Moderate (5-15 mph)
     - 0.15-0.25
     - 15-25% of heading ROS
   * - High (> 15 mph)
     - 0.10-0.15
     - 10-15% of heading ROS

**Model Formulation**

In the Richards (1990) ellipse model, fire spread is computed as:

Heading direction (cos θ > 0):

.. math::

   R = (a \cos\theta + b \sin\theta) \times R_{\text{head}}

Backing direction (cos θ < 0):

.. math::

   R = (b \sin\theta - c \cos\theta) \times R_{\text{head}}

where:

* :math:`a` = head fire coefficient (maximum spread, downwind)
* :math:`b` = flank fire coefficient (crosswind spread)
* :math:`c` = backing fire coefficient (minimum spread, upwind)
* :math:`\theta` = angle between fire normal and wind direction

**Backing Coefficient Adjustment**

The backing coefficient :math:`c` can be computed from wind-dependent empirical
ratios:

.. math::

   c = \frac{c}{a} \times a

where :math:`c/a` is the backing fire ratio from field observations.

**API Functions**

Compute backing fire ratio from wind speed:

.. code-block:: cpp

   Real backing_ratio = compute_backing_fire_ratio(wind_speed_mph);

Compute backing ROS directly:

.. code-block:: cpp

   Real R_backing = compute_backing_ros_from_heading(R_heading, wind_speed_mph);

Adjust ellipse coefficient:

.. code-block:: cpp

   Real coeff_c = adjust_backing_coefficient(coeff_a, wind_speed_mph);

**Integration with FARSITE Ellipse Model**

The backing fire formulation is already integrated into ``farsite_ellipse.H``
through the Richards (1990) ellipse coefficients. The ``backing_fire_ros.H``
header provides documentation and helper functions for ensuring coefficients
match empirical observations.

**References:**

* Rothermel, R.C. (1972). "A mathematical model for predicting fire spread
  in wildland fuels." USDA Forest Service Research Paper INT-115.
* Anderson, H.E. (1983). "Predicting wind-driven wildland fire size and shape."
  USDA Forest Service Research Paper INT-305.
* Richards, G.D. (1990). "An elliptical growth model of forest fire fronts."
  International Journal of Numerical Methods in Engineering, 30(6):1163-1179.
* Finney, M.A. (1998). "FARSITE: Fire Area Simulator." USDA Forest Service
  Research Paper RMRS-RP-4.

Fire Spread Direction Output
-----------------------------

**Header**: ``src/plot_results.H``

Fire spread direction vectors (``spread_dir_x``, ``spread_dir_y``) are written
to every plotfile for visualization and validation. These components represent
the direction of maximum fire spread based on the ROS gradient field.

**Model Formulation**

The spread direction is computed from the level-set advection velocity:

.. math::

   \mathbf{d} = \frac{\nabla \phi}{|\nabla \phi|}

where :math:`\phi` is the level-set function. The spread direction points
toward the unburned region (increasing :math:`\phi`).

**Plotfile Variables**

.. list-table::
   :header-rows: 1

   * - Variable Name
     - Description
     - Units
   * - ``spread_dir_x``
     - X-component of spread direction
     - Dimensionless (unit vector)
   * - ``spread_dir_y``
     - Y-component of spread direction
     - Dimensionless (unit vector)

**Visualization**

Spread direction arrows can be visualized in ParaView or VisIt using the
vector field ``(spread_dir_x, spread_dir_y)``. This is valuable for:

* Validating elliptical fire spread models
* Identifying wind-terrain interaction effects
* Debugging fire front propagation
* Operational fire behavior briefings

**Note:** Fire spread direction is automatically computed and written to
plotfiles. No additional input parameters are required.

Summary of 2026 Easy-Implementation Features
---------------------------------------------

Four new operational fire danger indices have been added:

1. **Keetch-Byram Drought Index (KBDI)** - Soil moisture deficit (0-800)
2. **Haines Index** - Atmospheric stability for plume-dominated fires (2-6)
3. **Backing Fire ROS** - Empirical backing fire rate reductions (10-40%)
4. **Fire Spread Direction** - Spread direction vectors for visualization

**Key Characteristics:**

* **Minimal Code Changes:** Each feature implemented as single header file (~100-250 lines)
* **High Operational Value:** Used by fire management agencies worldwide
* **GPU Compatible:** All functions use ``AMREX_GPU_HOST_DEVICE`` macros
* **Literature-Based:** Formulas from peer-reviewed wildfire research
* **Easy Integration:** Drop-in additions to existing simulation workflow

**Implementation Pattern:**

All features follow the established codebase conventions:

* Single-header implementation (``src/*.H``)
* GPU-compatible device functions
* No external dependencies beyond AMReX
* Minimal computational overhead
* Well-documented with references


Phase 1 & Phase 2 Operational Enhancements (2026)
==================================================

Five additional operational fire behavior and fire danger features have been 
implemented to enhance multi-day simulations and operational fire management support.

Grass Curing Model
------------------

**Header**: ``src/grass_curing_model.H``

Grass curing percentage models affect fuel availability and moisture in grassland
and savanna fuel types. Curing describes the fraction of herbaceous fuel that is
dead (cured) versus live (green).

**Model Options:**

1. **Seasonal Model** (Luke & McArthur 1978):

   .. math::

      C(t) = C_{\text{mean}} + C_{\text{amp}} \times \cos\left(\frac{2\pi(d - d_0)}{365}\right)

   where :math:`d` is day of year, :math:`d_0` is peak green day (default: day 15).

2. **Moisture-Dependent Model** (Cheney & Sullivan 2008):

   .. math::

      C = 100 \times \left(1 - \exp\left(-k \times \frac{\text{KBDI}}{200}\right)\right)

   Curing increases exponentially with drought severity.

3. **Growing Degree Day Model**:

   .. math::

      C = 100 \times \min\left(1, \frac{\text{GDD}_{\text{accum}}}{\text{GDD}_{\text{threshold}}}\right)

   Curing increases linearly with accumulated heat units.

**Functions:**

* ``compute_seasonal_curing()`` - Sinusoidal annual cycle
* ``compute_moisture_dependent_curing()`` - KBDI-driven curing
* ``compute_gdd_curing()`` - Phenology-based curing
* ``apply_curing_to_fuel_load()`` - Partition herbaceous fuel between dead/live

**Interpretation:**

.. list-table::
   :header-rows: 1

   * - Curing %
     - Description
     - Fire Danger
   * - 0-30%
     - Predominantly green grass
     - Low (limited fire spread)
   * - 30-70%
     - Mixed green/cured
     - Moderate
   * - 70-90%
     - Mostly cured
     - High (rapid grass fire spread)
   * - 90-100%
     - Fully cured/dry
     - Very High to Extreme


Diurnal Weather Cycles
-----------------------

**Header**: ``src/diurnal_weather.H``

Implements daily (diurnal) variation of temperature, relative humidity, and wind
speed for realistic multi-day fire simulations.

**Temperature:**

.. math::

   T(t) = T_{\text{mean}} + T_{\text{amp}} \times \sin\left(\frac{2\pi(t - t_{\text{min}})}{24}\right)

where :math:`t_{\text{min}}` is time of minimum temperature (default: 6:00 AM).

**Relative Humidity (anti-phase with temperature):**

.. math::

   \text{RH}(t) = \text{RH}_{\text{mean}} - \text{RH}_{\text{amp}} \times \sin\left(\frac{2\pi(t - t_{\text{max}})}{24}\right)

**Wind Speed:**

.. math::

   U(t) = U_{\text{mean}} \times \left[1 + \alpha \times \sin\left(\frac{2\pi(t - 6)}{24}\right)\right]

where :math:`\alpha` is amplitude factor (default: 0.4 for ±40% variation).

**Functions:**

* ``compute_diurnal_temperature()`` - Temperature at given hour
* ``compute_diurnal_rh()`` - Relative humidity at given hour
* ``compute_diurnal_wind_speed()`` - Wind speed scaling factor
* ``apply_diurnal_weather()`` - Update all fields

**Typical Diurnal Pattern:**

.. list-table::
   :header-rows: 1

   * - Time
     - Temperature
     - RH
     - Wind Speed
     - Fire Activity
   * - 06:00 (sunrise)
     - Minimum
     - Maximum
     - Minimum
     - Low spread, potential burnout
   * - 15:00 (afternoon)
     - Maximum
     - Minimum
     - Maximum
     - Peak fire activity
   * - 21:00 (evening)
     - Decreasing
     - Increasing
     - Decreasing
     - Moderating fire behavior


Elevation Temperature Lapse Rate
---------------------------------

**Header**: ``src/elevation_lapse_rate.H``

Applies elevation-based temperature and pressure corrections for mountainous terrain.

**Temperature Lapse:**

.. math::

   T(z) = T_{\text{ref}} - \Gamma \times (z - z_{\text{ref}})

where :math:`\Gamma` is environmental lapse rate (default: 0.0065 °C/m = 6.5 °C/km,
standard atmosphere).

**Relative Humidity Lapse:**

.. math::

   \text{RH}(z) = \text{RH}_{\text{ref}} \times \exp(k \times \Delta z)

where :math:`k \approx 0.0001` m⁻¹ (empirical coefficient).

**Barometric Pressure:**

.. math::

   P(z) = P_0 \times \left(1 - \frac{L \times z}{T_0}\right)^{\frac{gM}{RL}}

Standard atmosphere barometric formula.

**Functions:**

* ``apply_temperature_lapse()`` - Temperature correction for elevation
* ``apply_rh_lapse()`` - RH correction for elevation
* ``compute_barometric_pressure()`` - Pressure at elevation
* ``apply_elevation_lapse_fields()`` - Update temperature and RH fields
* ``compute_air_density_correction()`` - Density ratio for combustion

**Effects:**

For every 1000 m elevation gain:
* Temperature decreases ~6.5 °C
* RH increases ~10-15%
* Air density decreases ~12%
* Fuel moisture typically increases


NFDRS Spread Component (SC)
----------------------------

**Header**: ``src/nfdrs_spread_component.H``

Implements the U.S. National Fire Danger Rating System Spread Component, a
numerical rating of fire spread potential.

**Formulation (Deeming et al. 1977):**

.. math::

   \text{SC} = a \times (1 + b \times U) \times (1 + c \times S) \times \exp(-d \times \text{FM}/100)

where:
* :math:`U` - wind speed [mph]
* :math:`S` - slope [%]
* :math:`FM` - 1-hour dead fuel moisture [%]
* :math:`a,b,c,d` - fuel model coefficients

**Functions:**

* ``compute_nfdrs_spread_component()`` - SC from weather and fuel
* ``compute_nfdrs_fuel_model_factor()`` - Fuel-specific adjustment
* ``compute_nfdrs_sc_field()`` - SC for entire domain
* ``compute_nfdrs_burning_index()`` - BI = (SC × ERC) / 10

**Interpretation:**

.. list-table::
   :header-rows: 1

   * - SC Value
     - Spread Potential
     - Management Implication
   * - 0-10
     - Low
     - Limited fire spread
   * - 10-20
     - Moderate
     - Normal fire behavior
   * - 20-40
     - High
     - Aggressive initial attack needed
   * - 40-60
     - Very High
     - Difficult control, rapid spread
   * - 60+
     - Extreme
     - Potential for large fire growth

**Relationship to Other NFDRS Components:**

* **Burning Index (BI)** = SC × ERC / 10
* **Fire Load Index (FLI)** = BI × manning level / 100
* Used in dispatch planning and resource allocation


Chandler Burning Index (CBI)
-----------------------------

**Header**: ``src/chandler_burning_index.H``

Chandler Burning Index is a fire danger rating used by the U.S. National Weather
Service for fire weather forecasts and Red Flag Warning criteria.

**Formulation (Chandler et al. 1983):**

.. math::

   \text{CBI} = \left[(110 - 1.373 \times \text{RH}) - 0.54 \times (10.20 - T_F)\right] \times (1 + 0.0345 \times U_{\text{mph}})

where:
* :math:`T_F` - temperature [°F]
* :math:`\text{RH}` - relative humidity [%]
* :math:`U_{\text{mph}}` - wind speed [mph]

**Functions:**

* ``compute_chandler_burning_index()`` - CBI from T, RH, wind
* ``get_cbi_danger_class()`` - Fire danger classification (0-4)
* ``compute_fosberg_fire_weather_index()`` - Related FFWI index
* ``compute_cbi_field()`` - CBI for entire domain
* ``compute_cbi_scalar()`` - Single-point CBI calculation

**Fire Weather Thresholds:**

.. list-table::
   :header-rows: 1

   * - CBI Value
     - Danger Class
     - Fire Weather Action
   * - 0-50
     - Low
     - Routine operations
   * - 50-75
     - Moderate
     - Increased awareness
   * - 75-90
     - High
     - Fire Weather Watch
   * - 90-97.5
     - Very High
     - Red Flag Warning
   * - 97.5+
     - Extreme
     - Extreme Red Flag Warning

**Operational Use:**

* Fire weather forecasts and watches/warnings
* Public fire danger communication
* Burn ban decision support
* Often paired with KBDI for comprehensive fire danger assessment


Summary of All 2025 FARSITE-Parity Features
============================================

This summary table lists all 64 features and enhancements added in the 2025 update,
organized by component type with implementation difficulty and source file references.

.. list-table:: Complete Feature Reference
   :header-rows: 1
   :widths: 5 35 25 15 40

   * - # 
     - Feature / Function
     - Source File
     - Difficulty
     - Description
   * - 1
     - WAF Andrews Formula
     - ``andrews_model.H``
     - Easy
     - Logarithmic Albini & Baughman wind adjustment
   * - 2
     - WAF BehavePlus Formula
     - ``andrews_model.H``
     - Easy
     - Linear/exponential BehavePlus formula
   * - 3
     - Torching Ratio (I_B/I_o)
     - ``fire_ecology_model.H``
     - Easy
     - Crown initiation diagnostic field
   * - 4
     - Crowning Index (CI)
     - ``fire_ecology_model.H``
     - Easy
     - Crown fire initiation threshold
   * - 5
     - FMC CSV Schedule
     - ``fmc_schedule.H``
     - Easy
     - Day-of-year foliar moisture content
   * - 6
     - FMC FARSITE Curve
     - ``fmc_schedule.H``
     - Easy
     - Built-in phenological FMC schedule
   * - 7
     - Dead Fuel Wetting
     - ``diurnal_moisture.H``
     - Medium
     - Exponential precipitation wetting
   * - 8
     - 1hr Fuel Time Constant
     - ``diurnal_moisture.H``
     - Easy
     - Time-lag constant for 1-hour fuel
   * - 9
     - 10hr Fuel Time Constant
     - ``diurnal_moisture.H``
     - Easy
     - Time-lag constant for 10-hour fuel
   * - 10
     - 100hr Fuel Time Constant
     - ``diurnal_moisture.H``
     - Easy
     - Time-lag constant for 100-hour fuel
   * - 11
     - Polygon Fire Ignition
     - ``ignition.H``
     - Medium
     - Closed ignition polygon from CSV vertices
   * - 12
     - Polyline Fire Ignition
     - ``ignition.H``
     - Medium
     - Line-fire ignition from CSV polyline
   * - 13
     - Per-Fuel FMS Moisture
     - ``fms_moisture.H``
     - Medium
     - FARSITE FMS file dead/live moisture
   * - 14
     - Spatial FMS Override
     - ``fms_moisture.H``
     - Medium
     - Per-cell fuel moisture MultiFab
   * - 15
     - Wind Speed Schedule
     - ``wind_schedule.H``
     - Easy
     - 3-column CSV time-wind schedule
   * - 16
     - Wind Direction Schedule
     - ``wind_schedule.H``
     - Easy
     - Direction from CSV schedule
   * - 17
     - DeviceVector GPU Storage
     - ``Multiple headers``
     - Medium
     - amrex::Gpu::DeviceVector for device data
   * - 18
     - GPU Copy hostToDevice
     - ``Multiple headers``
     - Medium
     - Host-to-device transfer primitives
   * - 19
     - Non-Burnable Cell Masking
     - ``rothermel.H``
     - Easy
     - Fuel codes 91-99 zero ROS masking
   * - 20
     - Non-Burnable NB Codes
     - ``rothermel.H``
     - Easy
     - NB1-NB9 fuel codes masking
   * - 21
     - ROS Minimum Threshold
     - ``rothermel.H``
     - Easy
     - 1×10⁻⁴ m/s FARSITE-style stall
   * - 22
     - Fuel Residence Time
     - ``rothermel.H``
     - Medium
     - Rothermel burnout residence time
   * - 23
     - Live Fuel Conditioning Ramp
     - ``conditioning.H``
     - Easy
     - Linear ramp to equilibrium during spin-up
   * - 24
     - Retardant Spotting Suppression
     - ``spotting.H``
     - Medium
     - Aerial retardant suppresses spotting
   * - 25
     - Spatial Fuel Moisture Plotfile
     - ``plot_results.H``
     - Easy
     - Write 1hr/10hr/100hr/live moisture fields
   * - 26
     - Crown Fraction Burned (CFB)
     - ``crown_dynamics.H``
     - Hard
     - Scott & Reinhardt crown fire dynamics
   * - 27
     - Active Crown ROS Scaling
     - ``crown_dynamics.H``
     - Hard
     - Cruz et al. active crown rate
   * - 28
     - Passive Crown ROS Blending
     - ``crown_dynamics.H``
     - Hard
     - Van Wagner passive crown blending
   * - 29
     - Crown CFSA Wind Profile
     - ``crown_dynamics.H``
     - Hard
     - Multi-layer wind for CFSA
   * - 30
     - FBP Model Initialization
     - ``fbp_model.H``
     - Hard
     - Canadian Forest Fire Weather Index
   * - 31
     - FBP Fuel Types (C, D, M)
     - ``fbp_model.H``
     - Hard
     - Conifer, Deciduous, Mixed fuel types
   * - 32
     - FBP Fine Fuel Moisture Code
     - ``fbp_model.H``
     - Hard
     - FFMC fine fuel moisture response
   * - 33
     - FBP Duff Moisture Code
     - ``fbp_model.H``
     - Hard
     - DMC organic layer moisture
   * - 34
     - FBP Drought Code
     - ``fbp_model.H``
     - Hard
     - DC long-term drought index
   * - 35
     - FBP Initial Spread Index
     - ``fbp_model.H``
     - Hard
     - ISI fire intensity component
   * - 36
     - FBP Build-up Index
     - ``fbp_model.H``
     - Hard
     - BUI fuel consumption index
   * - 37
     - FBP Fire Weather Index
     - ``fbp_model.H``
     - Hard
     - FWI overall danger rating
   * - 38
     - FBP Rate of Spread
     - ``fbp_model.H``
     - Hard
     - Canadian FBP ROS calculation
   * - 39
     - Seasonal Curing Curve
     - ``grass_curing_model.H``
     - Medium
     - Luke & McArthur sinusoidal annual cycle
   * - 40
     - Moisture-Dependent Curing
     - ``grass_curing_model.H``
     - Medium
     - Cheney & Sullivan KBDI-driven curing
   * - 41
     - Growing Degree Day (GDD)
     - ``grass_curing_model.H``
     - Medium
     - Phenology-based curing
   * - 42
     - Curing to Fuel Load Partition
     - ``grass_curing_model.H``
     - Medium
     - Dead/live herbaceous fuel split
   * - 43
     - Diurnal Temperature Cycle
     - ``diurnal_weather.H``
     - Medium
     - Daily temperature variation
   * - 44
     - Diurnal Relative Humidity
     - ``diurnal_weather.H``
     - Medium
     - Daily RH variation
   * - 45
     - Diurnal Wind Speed Factor
     - ``diurnal_weather.H``
     - Medium
     - Daily wind speed scaling
   * - 46
     - Diurnal Weather Application
     - ``diurnal_weather.H``
     - Medium
     - Update all fields for hour
   * - 47
     - Temperature Lapse Rate
     - ``elevation_lapse.H``
     - Easy
     - Elevation-based temperature correction
   * - 48
     - RH Lapse Rate
     - ``elevation_lapse.H``
     - Easy
     - Elevation-based RH correction
   * - 49
     - Barometric Pressure
     - ``elevation_lapse.H``
     - Easy
     - ICAO standard atmosphere
   * - 50
     - Air Density Correction
     - ``elevation_lapse.H``
     - Easy
     - Density ratio for combustion
   * - 51
     - NFDRS Spread Component (SC)
     - ``nfdrs_indices.H``
     - Medium
     - US National Fire Danger rating
   * - 52
     - NFDRS Fuel Model Factor
     - ``nfdrs_indices.H``
     - Easy
     - Fuel-specific NFDRS adjustment
   * - 53
     - NFDRS Burning Index
     - ``nfdrs_indices.H``
     - Medium
     - BI = (SC × ERC) / 10
   * - 54
     - NFDRS SC Field Output
     - ``nfdrs_indices.H``
     - Easy
     - SC for entire domain
   * - 55
     - Chandler Burning Index (CBI)
     - ``fire_danger_indices.H``
     - Medium
     - CBI from T, RH, wind
   * - 56
     - CBI Danger Classification
     - ``fire_danger_indices.H``
     - Easy
     - Fire danger class 0-4
   * - 57
     - Fosberg Fire Weather Index
     - ``fire_danger_indices.H``
     - Medium
     - Related FFWI index
   * - 58
     - CBI Field Output
     - ``fire_danger_indices.H``
     - Easy
     - CBI for entire domain
   * - 59
     - CBI Scalar Calculation
     - ``fire_danger_indices.H``
     - Easy
     - Single-point CBI
   * - 60
     - compute_crown_fire_surface_area()
     - ``crown_fire_model.H``
     - Medium
     - CFSA from CBD and canopy height
   * - 61
     - compute_emc_with_hysteresis()
     - ``fuel_moisture.H``
     - Medium
     - Equilibrium moisture content hysteresis
   * - 62
     - update_multi_class_fuel_moisture()
     - ``fuel_moisture.H``
     - Medium
     - Simultaneous 4-class moisture update
   * - 63
     - Multi-layer Wind Profile
     - ``wind_model.H``
     - Hard
     - N-layer exponential wind profile
   * - 64
     - Solar Radiation Calculation
     - ``solar_radiation.H``
     - Hard
     - Radiation preheating diagnostics

**Feature Statistics:**

* **Total Features:** 64
* **Easy Implementation:** 25 features
* **Medium Implementation:** 30 features
* **Hard Implementation:** 9 features

**Primary Application Areas:**

* **Wind & Fuel Models:** WAF formulas, multi-layer profiles, fuel conditioning (14 features)
* **Fire Behavior:** Crown fire dynamics, FBP model, spotting suppression (19 features)
* **Fire Danger Indices:** CBI, NFDRS, diurnal weather, curing (26 features)
* **GPU/Technical:** Device vectors, compute kernels, plotfile output (5 features)

All features are GPU-compatible, header-only implementation, based on peer-reviewed 
literature, and backward compatible (optional).
