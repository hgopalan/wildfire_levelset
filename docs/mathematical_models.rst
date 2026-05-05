Mathematical Models
===================

This section describes the mathematical models implemented in the wildfire level-set solver.

Fire Spread Models
------------------

Three primary fire spread model families are implemented, each with optional add-on modules.
Select the active model with ``fire_spread_model`` and the propagation method with
``propagation_method``.

**Rothermel (1972) Surface Fire Spread Model** (``fire_spread_model = rothermel``)

  An empirical model based on fuel properties, reaction intensity, and propagating flux ratio.
  Full details: `Rothermel (1972) Fire Spread Model`_.

  Optional add-ons:

  * *Andrews (2018) wind adjustments* — Wind Adjustment Factor (WAF) and Maximum Effective
    Wind Speed (MEWS) cap. See `Andrews (2018) Wind Adjustments for Rothermel`_.
  * *Viegas (2004) eruptive fire diagnostics* — exponential slope enhancement factor and
    eruptive-regime flag for steep terrain. See `Viegas (2004) Eruptive Fire Diagnostics`_.
  * *Wind-terrain feedback models* — seven options coupling terrain or fire-induced winds into
    the effective wind seen by the spread model (canyon wind, buoyancy upslope, Pimont
    exponential correction, WindNinja ridge/canyon speed-up).
    See `Wind-Terrain Feedback Models (Options 1–6)`_ and
    `WindNinja Ridge/Canyon Terrain Speed-up (Option 7)`_.
  * *Heat-flux-driven wind corrections* — fire-plume buoyancy and horizontal inflow based on a
    spatially-varying or uniform heat-release field. See `Heat Flux MultiFab and Fire-Induced Wind`_.

**Cheney & Gould (1995/1998) Grassland Fire Spread Model** (``fire_spread_model = cheney_gould``)

  An empirical ROS formula calibrated on Australian grassland fires; activates with
  ``fire_spread_model = cheney_gould``.  No add-ons; terrain slope is intentionally omitted.
  Full details: `Cheney & Gould (1995 / 1998) Grassland Fire Spread Model`_.

**Balbi (2009) Physical Fire Spread Model** (``fire_spread_model = balbi``)

  A radiation-driven model that derives ROS from an energy balance between the tilted flame and
  the unburned fuel ahead; fully replaces Rothermel.
  Full details: `Balbi (2009) Physical Fire Spread Model`_.

  Optional add-ons:

  * *Viegas (2004) eruptive fire diagnostics* — uses the Balbi amplitude coefficient as the
    ROS baseline instead of Rothermel :math:`R_0`.  See `Viegas (2004) Eruptive Fire Diagnostics`_.
  * *All six wind-terrain feedback options* — fully compatible; wind and buoyancy velocity
    augmented by heat flux when ``heat_flux.*`` is active.
    See `Wind-Terrain Feedback Models (Options 1–6)`_ and `Heat Flux MultiFab and Fire-Induced Wind`_.

**Fire propagation methods** (``propagation_method``)

  * *Level-set method* (``propagation_method = levelset``) — fire perimeter tracked as the zero
    contour of a signed-distance function :math:`\phi`, advanced with WENO5-Z spatial
    discretization and RK3 time integration.  See `Level-Set Method`_.
  * *FARSITE elliptical expansion* (``propagation_method = farsite``) — Richards (1990) elliptical
    model embedded in the Eulerian level-set framework; implements the Huygens wavelet principle.
    See `FARSITE Elliptical Expansion Model (Richards 1990)`_.

Rothermel (1972) Fire Spread Model
-----------------------------------

The Rothermel fire spread model is the foundation for computing the rate of spread (ROS) of wildfire. The model computes the no-wind, no-slope rate of spread and then applies wind and slope correction factors.

Basic Rate of Spread
^^^^^^^^^^^^^^^^^^^^^

The no-wind, no-slope rate of spread :math:`R_0` (ft/min) is given by:

.. math::

   R_0 = \frac{I_R \xi}{\rho_b \epsilon_h Q_{ig}}

where:

* :math:`I_R` is the reaction intensity (BTU/ft²/min)
* :math:`\xi` is the propagating flux ratio
* :math:`\rho_b` is the oven-dry bulk density (lb/ft³)
* :math:`\epsilon_h` is the effective heating number
* :math:`Q_{ig}` is the heat of preignition (BTU/lb)

Reaction Intensity
^^^^^^^^^^^^^^^^^^

The reaction intensity :math:`I_R` is computed as:

.. math::

   I_R = \Gamma' w_n h \eta_M \eta_s

where:

* :math:`\Gamma'` is the optimum reaction velocity (1/min)
* :math:`w_n` is the net fuel loading (lb/ft²)
* :math:`h` is the heat content of fuel (BTU/lb)
* :math:`\eta_M` is the moisture damping coefficient
* :math:`\eta_s` is the mineral damping coefficient

The optimum reaction velocity is:

.. math::

   \Gamma' = \Gamma_{max} \left(\frac{\beta}{\beta_{op}}\right)^A \exp\left[A\left(1 - \frac{\beta}{\beta_{op}}\right)\right]

where:

.. math::

   \Gamma_{max} = \frac{\sigma^{1.5}}{495 + 0.0594\sigma^{1.5}}

   \beta_{op} = 3.348 \sigma^{-0.8189}

   A = 133 \sigma^{-0.7913}

* :math:`\sigma` is the surface-area-to-volume ratio (ft²/ft³)
* :math:`\beta` is the packing ratio (dimensionless)
* :math:`\beta_{op}` is the optimum packing ratio

Packing Ratio
^^^^^^^^^^^^^

The packing ratio is defined as:

.. math::

   \beta = \frac{\rho_b}{\rho_p}

where:

* :math:`\rho_b = w_0 / \delta` is the oven-dry bulk density (lb/ft³)
* :math:`\rho_p` is the oven-dry particle density (lb/ft³)
* :math:`w_0` is the oven-dry fuel loading (lb/ft²)
* :math:`\delta` is the fuel bed depth (ft)

Moisture Damping
^^^^^^^^^^^^^^^^

**Single-class path (aggregate)**

When no per-class loads are specified (the default), the moisture damping
coefficient :math:`\eta_M` is computed from a single composite moisture content
:math:`M_f` and extinction moisture :math:`M_x`:

.. math::

   \eta_M = \max\left(1 - 2.59r_m + 5.11r_m^2 - 3.52r_m^3, 0\right)

where :math:`r_m = \min(M_f / M_x, 1)`, :math:`M_f` is the fuel moisture content, and :math:`M_x` is the moisture of extinction.

**Multi-class path (per size class)**

When per-class fuel loads are provided (via a fuel model database entry or
explicit ``rothermel.w_d1``, etc. inputs), the full Rothermel (1972)
multi-class formulation is used.  Each size class has its own oven-dry fuel
load :math:`w_i`, surface-area-to-volume ratio :math:`\sigma_i`, and moisture
content :math:`M_i`.

The fixed SAV values for coarser dead-fuel classes are:
:math:`\sigma_{d10} = 109\ \text{ft}^{-1}` and :math:`\sigma_{d100} = 30\ \text{ft}^{-1}`.

The per-class weighting factor is:

.. math::

   A_i = w_i \, \sigma_i \, \exp\!\left(-\frac{138}{\sigma_i}\right)

**Dead-fuel composite moisture ratio** (Rothermel 1972, Eq. 86):

.. math::

   r_{m,\text{dead}} = \frac{\sum_{i \in \text{dead}} A_i M_i}{M_x \sum_{i \in \text{dead}} A_i}

**Live-fuel moisture of extinction** (Rothermel 1972, Eq. 88 approximation):

.. math::

   M_{x,\text{live}} = \max\!\left(
       2.9 \,\frac{A_\text{dead}}{A_\text{live}}
       \left(1 - \frac{\eta_{M,\text{dead}}}{0.3}\right) - 0.226,\; 0.30
   \right)

**Live-fuel composite moisture ratio**:

.. math::

   r_{m,\text{live}} = \frac{\sum_{i \in \text{live}} A_i M_i}{M_{x,\text{live}} \sum_{i \in \text{live}} A_i}

Both dead and live moisture damping coefficients use the cubic polynomial (Eq. 29).

**Reaction intensity** sums dead and live contributions:

.. math::

   I_R = \Gamma' \eta_s \left(
       w_{n,\text{dead}} \, h \, \eta_{M,\text{dead}}
     + w_{n,\text{live}} \, h \, \eta_{M,\text{live}}
   \right)

where :math:`w_{n,\text{dead}} = \sum_{i \in \text{dead}} w_i (1-S_T)` and
:math:`w_{n,\text{live}} = \sum_{i \in \text{live}} w_i (1-S_T)`.

The characteristic SAV :math:`\sigma_c` (load-weighted mean across all
classes) is used in place of the single aggregate :math:`\sigma` for
:math:`\Gamma'`, :math:`\xi`, :math:`\epsilon_h`, :math:`Q_{ig}`, and the
wind-factor coefficients :math:`C`, :math:`B`, :math:`E`.

Mineral Damping
^^^^^^^^^^^^^^^

The mineral damping coefficient :math:`\eta_s` is:

.. math::

   \eta_s = 0.174 S_e^{-0.19}

where :math:`S_e` is the effective mineral content.

Propagating Flux Ratio
^^^^^^^^^^^^^^^^^^^^^^^

The propagating flux ratio :math:`\xi` is:

.. math::

   \xi = \frac{\exp[(0.792 + 0.681\sqrt{\sigma})(\beta + 0.1)]}{192 + 0.2595\sigma}

Heat of Preignition
^^^^^^^^^^^^^^^^^^^

The heat of preignition :math:`Q_{ig}` (BTU/lb) depends on fuel moisture:

.. math::

   Q_{ig} = 250 + 1116 M_f

Wind Factor
^^^^^^^^^^^

The wind factor :math:`\phi_w` multiplies the base rate of spread to account for wind effects:

.. math::

   \phi_w = C \left(\frac{\beta}{\beta_{op}}\right)^{-E} U^B

where:

* :math:`U` is the midflame wind speed (ft/min)
* :math:`C = 7.47 \exp(-0.133 \sigma^{0.55})`
* :math:`B = 0.02526 \sigma^{0.54}`
* :math:`E = 0.715 \exp(-3.59 \times 10^{-4} \sigma)`

The rate of spread with wind is:

.. math::

   R_w = R_0 (1 + \phi_w)

Slope Factor
^^^^^^^^^^^^

The slope factor :math:`\phi_s` accounts for terrain slope:

.. math::

   \phi_s = 5.275 \beta^{-0.3} \tan^2(\theta)

where :math:`\theta` is the slope angle. The slope can be decomposed into x and y components:

.. math::

   \tan^2(\theta) = \left(\frac{\partial z}{\partial x}\right)^2 + \left(\frac{\partial z}{\partial y}\right)^2

Combined Wind and Slope
^^^^^^^^^^^^^^^^^^^^^^^^

When both wind and slope effects are present, the total rate of spread is:

.. math::

   R = R_0 (1 + \phi_w + \phi_s)

FARSITE Elliptical Expansion Model (Richards 1990)
---------------------------------------------------

The FARSITE model uses an elliptical expansion pattern to represent the anisotropic spread of fire under the influence of wind and terrain.

Richards' Coefficients
^^^^^^^^^^^^^^^^^^^^^^^

The elliptical fire spread is characterized by three coefficients:

* :math:`a` - head fire coefficient (maximum spread, downwind)
* :math:`b` - flank fire coefficient (crosswind spread)
* :math:`c` - backing fire coefficient (minimum spread, upwind)

These coefficients relate the directional rate of spread to the base rate:

.. math::

   R(\theta) = R_0 f(a, b, c, \theta)

where :math:`\theta` is the angle between the fire normal and the wind direction.

Anderson Length-to-Width Ratio
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The Anderson (1983) model relates the ellipse shape to wind speed:

.. math::

   L/W = 0.936 \exp(0.2566 U) + 0.461 \exp(-0.1548 U) - 0.397

where :math:`L/W` is the length-to-width ratio and :math:`U` is the wind speed (mph). This ratio is then converted to the Richards' coefficients for the elliptical spread calculation.

Van Wagner (1977) Crown Fire Initiation
----------------------------------------

Crown fire initiation occurs when sufficient heat is available to raise the canopy fuel to ignition temperature.

Crown Fire Initiation Criterion
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The critical surface fire intensity :math:`I_0` (kW/m) required for crown fire initiation is:

.. math::

   I_0 = \frac{CBH \times (460 + 25.9 M_c)}{18 h}

where:

* :math:`CBH` is the canopy base height (m)
* :math:`M_c` is the canopy moisture content (%)
* :math:`h` is the crown fuel consumption depth (m)

If the actual surface fire intensity :math:`I > I_0`, then crown fire initiation occurs.

Active Crown Fire Criterion
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For active crown fire spread, the critical wind speed :math:`U_{active}` is:

.. math::

   U_{active} = \frac{3}{\sqrt{CBD}}

where :math:`CBD` is the canopy bulk density (kg/m³).

Firebrand Spotting
------------------

Firebrand spotting generates new ignition points ahead of the main fire front.
Two independent spotting models are available; they can be enabled separately or
together.

Probability-Based Spotting (``spotting.*``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A stochastic model driven by wind speed, fire intensity, and fuel moisture.
For each fire-front cell a spotting probability is computed and a Bernoulli draw
decides whether a firebrand is generated.  If so, the landing distance is
sampled from a log-normal or exponential distribution and the spot is placed
downwind with optional lateral dispersion.

Albini (1983) Spotting with 2-D Trajectory (``albini_spotting.*``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A physics-based spotting model that couples Albini's thermal-plume lofting
formula with a 2-D horizontal particle trajectory integrated through the
simulation wind field.

**Stage 1 – Lofting height (Albini 1983)**

Byram's fire line intensity is computed from the Rothermel rate-of-spread:

.. math::

   I_B \;[\text{kW/m}] = \frac{R \;[\text{ft/min}]\; w_0 \;[\text{lb/ft}^2]\; h \;[\text{BTU/lb}]}{60} \times 3.459

The maximum height reached by a firebrand above the fire front is:

.. math::

   H_z \;[\text{m}] = 12.2 \; I_B^{1/3}

where :math:`I_B` is in kW/m.  Cells with :math:`I_B < I_{B,\min}` do not
generate firebrands.

**Stage 2 – 2-D horizontal trajectory**

Each lofted firebrand descends at a constant terminal velocity
:math:`v_t` (m/s).  The flight time is:

.. math::

   t_f = \frac{H_z}{v_t}

During the flight the horizontal position is integrated with forward Euler over
:math:`n_{\text{traj}}` sub-steps using bilinear interpolation of the
cell-centred velocity field :math:`(u, v)`:

.. math::

   x(t + \Delta t) &= x(t) + u\!\bigl(x(t),\,y(t)\bigr)\,\Delta t \\
   y(t + \Delta t) &= y(t) + v\!\bigl(x(t),\,y(t)\bigr)\,\Delta t

where :math:`\Delta t = t_f / n_{\text{traj}}`.

The landing position is the final :math:`(x, y)` after the full trajectory.
A circular ignition zone of radius ``albini_spotting.spot_radius`` is then
imposed by setting :math:`\phi` to a negative signed-distance value in that
neighbourhood.

**Diagnostic output fields**

Each plot file contains four Albini-specific scalar fields:

* ``albini_Hz`` – lofting height :math:`H_z` at fire-front source cells.
* ``albini_count`` – number of firebrands launched from each source cell at the
  last spotting step.
* ``albini_dist`` – maximum landing distance from each source cell.
* ``albini_active`` – flag (1) at cells that received a spot ignition.

Bulk Fuel Consumption
----------------------

The bulk fuel consumption fraction :math:`f_c` represents the fraction of available fuel consumed behind the fire front:

.. math::

   f_c = f_c(I_R, \tau_{res})

where :math:`\tau_{res}` is the residence time. The consumption affects the fire intensity and heat release rate.

The consumption fraction is bounded by configurable minimum and maximum values
(``farsite.f_consumed_min`` and ``farsite.f_consumed_max``):

.. math::

   f_c = f_{c,\min} + (f_{c,\max} - f_{c,\min}) \left(1 - \exp\left(-\frac{t}{\tau_{res}}\right)\right)

Byram Fire Behavior Diagnostics
---------------------------------

Fireline Intensity
^^^^^^^^^^^^^^^^^^

Byram's (1959) fire line intensity :math:`I_B` quantifies the rate of heat release per unit length of fire front:

.. math::

   I_B \;[\text{kW/m}] = H \;[\text{kJ/kg}] \times w_a \;[\text{kg/m}^2] \times R \;[\text{m/s}]

where :math:`H` is the low heat of combustion, :math:`w_a` is the available (net) fuel load
per unit area, and :math:`R` is the rate of spread. The available fuel load
accounts for the mineral content correction:

.. math::

   w_a = w_0 (1 - S_T) \times 4.8824 \;[\text{kg/m}^2]

(where :math:`w_0` is in lb/ft² and 4.8824 is the lb/ft² to kg/m² conversion).

Flame Length
^^^^^^^^^^^^^

Byram's (1959) empirical relationship between fireline intensity and flame length:

.. math::

   L_f \;[\text{m}] = 0.0775 \times I_B^{0.46}

These fields (``fireline_intensity`` and ``flame_length``) are computed at every time
step and written to each plotfile.

Cheney & Gould (1995 / 1998) Grassland Fire Spread Model
----------------------------------------------------------

The Cheney–Gould model is a purely empirical rate-of-spread (ROS) formula calibrated against a large number of experimental grassland fires conducted in Australia. It is activated by setting ``fire_spread_model = cheney_gould`` and generally outperforms the Rothermel model in open grassland fuels.

Head-Fire Rate of Spread
^^^^^^^^^^^^^^^^^^^^^^^^^

The head-fire rate of spread :math:`R` (km/h) is a piecewise-linear function of
the 10-m open wind speed :math:`U_{10}` (km/h), a moisture correction factor
:math:`f_{MC}`, and a curing factor :math:`CF`:

.. math::

   R = \begin{cases}
       (0.165 + 0.534\,U_{10}) \times f_{MC} \times CF & U_{10} \le 5\ \text{km/h} \\
       (-0.020 + 0.640\,U_{10}) \times f_{MC} \times CF & U_{10} > 5\ \text{km/h}
   \end{cases}

The result is clamped to :math:`R \ge 0` — in the high-wind regime the
intercept :math:`-0.020` makes the base ROS slightly negative for extremely
low wind speeds just above 5 km/h, which is physically meaningless — and
then converted to m/s for the simulation:

.. math::

   R\ [\text{m/s}] = \frac{R\ [\text{km/h}]}{3.6}

The 10-m wind speed :math:`U_{10}` is obtained by converting the simulation wind
field (m/s) via :math:`U_{10} = |\mathbf{u}| \times 3.6`.

Moisture Correction Factor
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. math::

   f_{MC} = \exp(-0.108 \times MC)

where :math:`MC` is the dead fine fuel moisture content [%]. At :math:`MC = 0`
(bone-dry fuel) :math:`f_{MC} = 1`; at :math:`MC \approx 22\%` the ROS is
halved relative to dry conditions.

Curing Factor
^^^^^^^^^^^^^

:math:`CF \in [0,\,1]` represents the degree of curing of the grass:
:math:`CF = 1` is fully cured (dry standing grass), while :math:`CF = 0` is
completely green (no spread). The factor is clamped to the :math:`[0,1]`
interval at runtime.

Note on Terrain Slope
^^^^^^^^^^^^^^^^^^^^^

Terrain slope is **not** accounted for in the original Cheney–Gould empirical
formulation. The model was calibrated on flat or gently sloping Australian
grasslands and the slope correction is intentionally omitted in this
implementation. For slope effects, use the Rothermel or Balbi models instead.

Viegas (2004) Eruptive Fire Diagnostics
-----------------------------------------

The Viegas (2004) model is an **optional parallel diagnostic** that characterises
eruptive (blow-up) fire behaviour on steep terrain.  Unlike the primary spread
models it does **not** alter the fire front propagation; instead it writes five
diagnostic fields to every plotfile so that potential under-prediction by
Rothermel's quadratic slope factor can be identified.  Enable with
``viegas.enable = 1``.

Physical Basis
^^^^^^^^^^^^^^

On steep slopes the Rothermel quadratic slope factor :math:`\phi_s` significantly
under-predicts the rate of spread compared to field and laboratory observations of
eruptive fires (Viegas 2004).  Viegas introduces an *exponential* slope enhancement
factor that captures the runaway acceleration observed at critical slope angles.

Exponential Slope Enhancement Factor
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. math::

   \Phi_{s,V} = \exp\!\bigl(a_V \,\tan\varphi\bigr)

where:

* :math:`a_V` is the Viegas slope coefficient (default 1.83, calibrated to
  laboratory fuels; dimensionless)
* :math:`\varphi` is the terrain slope angle; :math:`\tan\varphi =
  \sqrt{(\partial z/\partial x)^2 + (\partial z/\partial y)^2}`

Viegas Rate of Spread
^^^^^^^^^^^^^^^^^^^^^

The Viegas ROS uses the Rothermel no-wind, no-slope rate :math:`R_0` and the
Rothermel wind factor :math:`\phi_w` as a baseline and replaces only the slope
factor:

.. math::

   R_V = R_0 \,(1 + \phi_w) \,\Phi_{s,V}
       = R_0 \,(1 + \phi_w) \,\exp\!\bigl(a_V \,\tan\varphi\bigr)

Eruptive Regime Flag
^^^^^^^^^^^^^^^^^^^^

A cell is flagged as being in the eruptive regime when the terrain slope exceeds
a critical value :math:`\tan\varphi_c` (Viegas 2004, Section 3):

.. math::

   \text{flag} =
   \begin{cases}
     1 & \tan\varphi > \tan\varphi_c \\
     0 & \text{otherwise}
   \end{cases}

The default critical slope is :math:`\tan\varphi_c = 0.4` (approximately
22°).  Flagged cells indicate locations where the Rothermel model is expected
to under-predict the rate of spread.

ROS Excess Ratio
^^^^^^^^^^^^^^^^

The signed excess ratio quantifies the relative difference between the two models:

.. math::

   \varepsilon = \frac{R_V - R_{\text{Rothermel}}}{R_{\text{Rothermel}}}

Positive values indicate that Viegas predicts a higher spread rate than
Rothermel; the larger :math:`\varepsilon`, the greater the potential
under-prediction by Rothermel on steep terrain.

Flame-Tilt Angle (Hazard Assessment)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The combined wind–slope flame-tilt angle :math:`\alpha` is computed for hazard
assessment purposes only; it is **not** fed back into the rate-of-spread
calculation:

.. math::

   \tan\alpha = \frac{U}{v_b} + \tan\varphi

where the buoyancy velocity scale :math:`v_b` (m/s) is:

.. math::

   v_b = \sqrt{\frac{g \,\delta_m \,(T_f - T_a)}{T_a}}

and:

* :math:`U` – wind speed magnitude (m/s)
* :math:`g = 9.81` m/s² – gravitational acceleration
* :math:`\delta_m` – fuel bed depth (m), converted from the Rothermel database value in ft
* :math:`T_f` – mean flame temperature (K); default 1000 K
* :math:`T_a` – ambient temperature (K); default 300 K

Merging Strategy
^^^^^^^^^^^^^^^^

The following table summarises which couplings between the Viegas and Rothermel
models are permitted in this implementation:

.. list-table::
   :header-rows: 1
   :widths: 55 10 35

   * - Coupling
     - Status
     - Rationale
   * - Rothermel :math:`R_0` used as Viegas baseline
     - ✅ compatible
     - Shared no-wind, no-slope ROS
   * - Same slope/wind/fuel inputs as Rothermel
     - ✅ compatible
     - Consistent environmental state per cell
   * - Viegas flame-tilt for hazard assessment only
     - ✅ compatible
     - Read-only diagnostic output
   * - Viegas induced wind fed into Rothermel
     - ❌ excluded
     - Would invalidate Rothermel calibration
   * - Viegas eruptive acceleration fed into Rothermel ROS
     - ❌ excluded
     - Incompatible ROS formulations
   * - Viegas critical slope replacing Rothermel :math:`\phi_s`
     - ❌ excluded
     - Different slope-factor definitions
   * - Viegas dynamic ROS inside Rothermel's static formula
     - ❌ excluded
     - Conceptually incompatible

Diagnostic Output Fields
^^^^^^^^^^^^^^^^^^^^^^^^^

Five fields are written to every plotfile when ``viegas.enable = 1``:

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   * - Field
     - Units
     - Description
   * - ``viegas_ROS``
     - m/s
     - :math:`R_0 (1+\phi_w)\exp(a_V \tan\varphi)`
   * - ``viegas_eruptive_flag``
     - –
     - 1 where :math:`\tan\varphi > \tan\varphi_c`
   * - ``viegas_ROS_excess``
     - –
     - :math:`(R_V - R_{\text{Rothermel}})/R_{\text{Rothermel}}`
   * - ``viegas_flame_tilt``
     - rad
     - :math:`\arctan(U/v_b + \tan\varphi)` (hazard only)
   * - ``viegas_slope_factor``
     - –
     - :math:`\exp(a_V \tan\varphi)`

Input Parameters
^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   * - Parameter
     - Default
     - Description
   * - ``viegas.enable``
     - 0
     - Enable (1) or disable (0) the Viegas diagnostics
   * - ``viegas.a_V``
     - 1.83
     - Exponential slope coefficient (dimensionless)
   * - ``viegas.tan_phi_c``
     - 0.4
     - Critical slope :math:`\tan\varphi_c` for eruptive-regime flag (≈ 22°)
   * - ``viegas.T_a``
     - 300.0 K
     - Ambient temperature for buoyancy velocity
   * - ``viegas.T_f``
     - 1000.0 K
     - Mean flame temperature for buoyancy velocity

Reference
^^^^^^^^^

Viegas, D.X. (2004). "Slope and wind effects on fire propagation."
*International Journal of Wildland Fire*, 13(2), 143–156.
https://doi.org/10.1071/WF03046

Viegas-Balbi Coupled Diagnostic
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When ``fire_spread_model = balbi``, the Viegas diagnostic uses the Balbi
amplitude coefficient :math:`A` and buoyancy velocity :math:`v_b` instead of
the Rothermel :math:`R_0` and wind factor :math:`\phi_w`.

The Viegas-Balbi rate of spread is:

.. math::

   R_V^{\text{Balbi}} = A \,\bigl(1 + \sin\alpha_w - \cos\alpha_w\bigr) \,
                        \Phi_{s,V}

where:

* :math:`\alpha_w` is the wind-only flame tilt angle:
  :math:`\tan\alpha_w = U / v_b` (no terrain slope in the tilt angle — the
  slope enters only through :math:`\Phi_{s,V}`)
* :math:`A` is the Balbi amplitude coefficient [m/s] pre-computed from fuel
  properties (radiation term, moisture correction, and fuel geometry)
* :math:`\Phi_{s,V} = \exp(a_V \tan\varphi)` is the same Viegas exponential
  slope enhancement factor used in the Rothermel baseline

This formulation is physically consistent: Balbi's radiation-driven ROS is
modulated by the same exponential slope term that Viegas calibrated from
laboratory and field eruptive fire experiments.  The ROS excess ratio becomes

.. math::

   \varepsilon = \frac{R_V^{\text{Balbi}} - R_{\text{Balbi}}}{R_{\text{Balbi}}}

Balbi (2009) Physical Fire Spread Model
----------------------------------------

The Balbi model provides an alternative rate-of-spread calculation rooted in
radiation physics rather than empirical curve fits.  It is activated by setting
``fire_spread_model = balbi`` and automatically replaces the Rothermel model in
both the level-set advection path and the pre-computed :math:`R` field used by
FARSITE.  The model is fully compatible with all six Viegas wind-terrain options
and FARSITE landscape file inputs.

Physical Basis
^^^^^^^^^^^^^^^

The rate of spread :math:`R` [m/s] is derived from the balance between radiant
heat flux from a tilted flame and the energy required to ignite the unburned fuel
ahead:

.. math::

   R = A_{\text{coeff}} \,(1 + \sin\alpha - \cos\alpha)

where:

* :math:`A_{\text{coeff}} = \chi\,\sigma_m\,\delta_m / (2\,\tau_0\,B^*)` [m/s] — fuel radiative amplitude
* :math:`\chi = r_{00}\,\sigma_m / (1 + r_{00}\,\sigma_m)` — forward radiation fraction
* :math:`B^* = (C_{pf}\,(T_i - T_a) + M_f\,\Lambda) / h` — dimensionless ignition energy
* :math:`\tan\alpha = U / v_b + \tan\theta` — flame tilt from wind and slope
* :math:`v_b = \sqrt{g\,\delta_m\,(T_f - T_a) / T_a}` — buoyancy velocity scale [m/s]

Inputs from the LCP / Fuel Database
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The Balbi model shares fuel parameters with the Rothermel model:

.. list-table::
   :header-rows: 1
   :widths: 25 20 55

   * - Parameter
     - Source
     - Description
   * - :math:`\sigma` [ft⁻¹→m⁻¹]
     - Fuel database / LCP
     - Surface-area-to-volume ratio
   * - :math:`\delta` [ft→m]
     - Fuel database / LCP
     - Fuel bed depth
   * - :math:`w_0` [lb/ft²→kg/m²]
     - Fuel database / LCP
     - Oven-dry fuel load
   * - :math:`\rho_p` [lb/ft³→kg/m³]
     - Fuel database / LCP
     - Particle density
   * - :math:`h` [BTU/lb→J/kg]
     - Fuel database / LCP
     - Heat of combustion
   * - :math:`M_f` [−]
     - ``rothermel.M_f``
     - Fuel moisture fraction
   * - slope, aspect
     - LCP / terrain file
     - Terrain slope components
   * - Wind :math:`U` [m/s]
     - ``u_x``, ``u_y``, wind file
     - Wind speed

Extra Inputs via Parmparse (``balbi.*``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   * - Parameter
     - Default
     - Description
   * - ``balbi.T_a``
     - 300.0 K
     - Ambient temperature
   * - ``balbi.T_f``
     - 1000.0 K
     - Mean flame temperature
   * - ``balbi.T_i``
     - 600.0 K
     - Ignition temperature
   * - ``balbi.delta_H``
     - 2.26e6 J/kg
     - Latent heat of water vaporisation
   * - ``balbi.C_pf``
     - 1800.0 J/(kg·K)
     - Specific heat of dry fuel
   * - ``balbi.r_00``
     - 2.5e-4 m
     - Radiation length scale
   * - ``balbi.tau_0``
     - 75591.0 s/m
     - Residence-time coefficient

Auto-Generated Balbi Fuel Parameter Table
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When ``fire_spread_model = balbi``, the solver automatically prints a table of
pre-computed Balbi parameters for every fuel model in the active database
(FBFM13 or FBFM40) at startup.  Columns are: fuel code, short name,
:math:`\sigma_m`, :math:`\delta_m`, :math:`\chi`, :math:`B^*`, :math:`v_b`,
:math:`A_{\text{coeff}}`, and the predicted ROS at 5 m/s wind on flat ground.

Reference
^^^^^^^^^^

Balbi, J.-H., Rossi, J.-L., Marcelli, T., and Santoni, P.-A. (2009).
"A physical model for wildland fires."
*Combustion and Flame*, 156(12), 2217–2230.
https://doi.org/10.1016/j.combustflame.2009.07.010

Andrews (2018) Wind Adjustments for Rothermel
----------------------------------------------

Andrews (2018) documents two critical wind adjustments that improve the physical
accuracy of the Rothermel surface fire spread model when wind input is from NWP
or WRF (20-ft / 10-m height).  Both are optional flags that can be enabled
independently.

Wind Adjustment Factor (WAF)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The Rothermel model requires wind speed at *midflame height*, but field
measurements and WRF/NWP output are typically at 20 ft (6.1 m) above open
terrain.  The WAF converts 20-ft open wind to midflame height:

.. math::

   \text{WAF} = \frac{1.83}{\ln\!\left(\dfrac{20 + 0.36\,h}{0.13\,h}\right)}

where :math:`h` is the fuel bed depth [ft].  Typical values: 0.42 for FM4
(:math:`h = 2` ft), 0.35 for FM9 (:math:`h = 3` ft).

Enable via ``rothermel.use_waf = 1``.  When a landscape file is active, WAF is
computed per cell using each fuel model's depth.

Maximum Effective Wind Speed (MEWS)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Rothermel's empirical wind factor :math:`\phi_w` grows without bound, producing
unrealistically large ROS at high wind speeds outside the model's calibration
range.  The MEWS cap limits the effective wind factor by inverting Rothermel's
wind-factor formula at a maximum wind factor value equal to 90% of the reaction
intensity (matching units: :math:`I_R` in BTU/ft²/min, :math:`\phi_w`
dimensionless via the original Rothermel calibration):

.. math::

   \phi_{w,\max} &= 0.9\,I_R \quad [I_R \text{ in BTU/ft}^2\text{/min}] \\
   U_{\max}      &= \left(\frac{\phi_{w,\max}}{C\,\beta_{\text{ratio}}^{-E}}\right)^{1/B}
                    \quad [\text{ft/min}]

Wind speed used in the ROS formula is capped at :math:`U_{\max}` before
evaluating the wind factor.  Enable via ``rothermel.use_wind_limit = 1``.

Reference
^^^^^^^^^^

Andrews, P.L. (2018). *The Rothermel Surface Fire Spread Model and Associated
Developments: A Comprehensive Explanation.* Gen. Tech. Rep. RMRS-GTR-371.
USDA Forest Service.
https://doi.org/10.2737/RMRS-GTR-371

Weise & Biging (1996) Fire Whirl Model
----------------------------------------

The Weise & Biging (1996) fire whirl model is an **optional diagnostic sub-model**
that computes fire whirl characteristics from existing fire behavior outputs.  It
does **not** modify the primary fire spread model; instead it runs alongside any
spread model to predict whirl formation when enabled via ``weise_biging.enable = 1``.

Physical Basis
^^^^^^^^^^^^^^^

Fire whirls are columnar vortices generated by fire-induced updrafts interacting
with ambient wind and terrain.  The model characterises whirl geometry and
kinematics from Byram fireline intensity and flame length.

Flame Tilt Angle
^^^^^^^^^^^^^^^^^

The Weise & Biging (1996) empirical flame tilt angle :math:`\theta`:

.. math::

   \tan\theta = 0.382\,Fr^{0.432} + \tan\beta

where :math:`Fr = U^2 / (g\,L_f)` is the modified Froude number, :math:`U` is
wind speed [m/s], :math:`g = 9.81` m/s², :math:`L_f` is Byram flame length [m],
and :math:`\beta` is the terrain slope angle.

Vertical Flame Height
^^^^^^^^^^^^^^^^^^^^^^^

.. math::

   H_f = L_f \cos\theta

Fire Whirl Geometry (Rankine-Vortex Scaling)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. math::

   H_w &= H_f \\
   r_w &= \max(c_r\,H_w,\; 0.01\ \text{m}) \\
   \Gamma &= U\,H_w \\
   \Omega &= \frac{\Gamma}{2\pi\,r_w^2} \\
   v_\theta &= \Omega\,r_w

Configuration via Parmparse (``weise_biging.*``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   * - Parameter
     - Default
     - Description
   * - ``weise_biging.enable``
     - 0
     - Enable (1) or disable (0) the fire whirl model
   * - ``weise_biging.c_r``
     - 0.1
     - Whirl core radius-to-height ratio (dimensionless)
   * - ``weise_biging.I_B_min``
     - 1.0 kW/m
     - Minimum Byram fireline intensity threshold

Diagnostic Output Fields
^^^^^^^^^^^^^^^^^^^^^^^^^^

Six fields are written to every plotfile when enabled:

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   * - Field
     - Units
     - Description
   * - ``weise_flame_height``
     - m
     - Vertical component of tilted flame
   * - ``weise_flame_tilt``
     - rad
     - Flame tilt angle from vertical
   * - ``weise_whirl_height``
     - m
     - Fire whirl column height
   * - ``weise_whirl_radius``
     - m
     - Fire whirl core radius
   * - ``weise_angular_velocity``
     - rad/s
     - Whirl angular velocity
   * - ``weise_max_tang_vel``
     - m/s
     - Maximum tangential velocity at core edge

Reference
^^^^^^^^^^

Weise, D.R. and Biging, G.S. (1996). "Effects of wind velocity and slope on
flame properties." *Canadian Journal of Forest Research*, 26(10), 1849–1858.
https://doi.org/10.1139/x26-210

Wind-Terrain Feedback Models (Options 1–6)
-------------------------------------------

The ``wind_terrain.model`` ParmParse key selects how terrain-induced or
fire-feedback winds modify the effective wind seen by the fire spread model.
**Option 1 (``none``) is the default and preserves all existing behaviour.**
Options 2–6 progressively couple Viegas or terrain-channel wind physics into
the actual fire propagation.  Option 7 (``windninja_ridge_canyon``) is described
in the following section.

Options Summary
^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 10 25 65

   * - Option
     - ``wind_terrain.model``
     - Effect on fire spread model
   * - 1
     - ``none``
     - No modification — default behaviour
   * - 2
     - ``viegas_ros``
     - :math:`R_{\text{final}} = \max(R_\text{Rothermel}, R_\text{Viegas})` in eruptive cells (:math:`\tan\varphi > \tan\varphi_c`)
   * - 3
     - ``viegas_wind``
     - Add buoyancy-induced upslope wind :math:`U_\text{ind} = v_b\tan\varphi` in eruptive cells only
   * - 4
     - ``canyon_wind``
     - Scale ambient wind by :math:`(1 + k_\text{canyon}\tan\varphi)` everywhere (Rothermel 1983)
   * - 5
     - ``viegas_neto``
     - Add upslope buoyancy wind :math:`U_\text{ind} = v_b\tan\varphi` at every cell (no threshold)
   * - 6
     - ``pimont``
     - Scale ambient wind by :math:`\exp(k_\text{pimont}\tan\varphi)` everywhere (Pimont et al. 2009)

Selecting ``viegas_ros``, ``viegas_wind``, or ``viegas_neto`` automatically enables
the Viegas (2004) diagnostic model (``viegas.enable = 1``).  Terrain slopes must be
provided (via ``rothermel.terrain_file`` or ``rothermel.landscape_file``) for any
slope-dependent option to have an effect.

Physical Equations
^^^^^^^^^^^^^^^^^^^

**Option 2 — Viegas ROS override** (requires terrain slopes):

.. math::

   R_{\text{final}}(i,j) =
   \begin{cases}
     \max\!\bigl(R_\text{Rothermel}(i,j),\; R_\text{Viegas}(i,j)\bigr)
       & \tan\varphi(i,j) > \tan\varphi_c \\
     R_\text{Rothermel}(i,j)
       & \text{otherwise}
   \end{cases}

**Options 3 and 5 — Buoyancy-induced upslope wind**:

.. math::

   v_b                  &= \sqrt{g\,\delta_m\,(T_f - T_a) / T_a} \\
   U_\text{ind}         &= v_b\,\tan\varphi \\
   \Delta\mathbf{U}     &= U_\text{ind}\,\hat{\mathbf{n}}_s \\
   \mathbf{U}_\text{eff}&= \mathbf{U}_\text{ambient} + \Delta\mathbf{U}

where :math:`\hat{\mathbf{n}}_s = \nabla z / |\nabla z|` is the unit upslope
vector, :math:`\delta_m` is fuel bed depth [m], and :math:`T_a`, :math:`T_f`
are taken from ``viegas.T_a`` / ``viegas.T_f``.  Option 3 applies
:math:`\Delta\mathbf{U}` only where :math:`\tan\varphi > \tan\varphi_c`;
Option 5 applies it everywhere.

**Option 4 — Canyon wind** (Rothermel 1983):

.. math::

   \mathbf{U}_\text{eff} = \mathbf{U}_\text{ambient}
   \,\max\!\bigl(1,\; 1 + k_\text{canyon}\tan\varphi\bigr)

**Option 6 — Pimont exponential correction** (Pimont et al. 2009):

.. math::

   \mathbf{U}_\text{eff} = \mathbf{U}_\text{ambient}
   \,\exp\!\bigl(k_\text{pimont}\tan\varphi\bigr)

Configuration via Parmparse (``wind_terrain.*``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   * - Parameter
     - Default
     - Description
   * - ``wind_terrain.model``
     - ``none``
     - Model: ``none``, ``viegas_ros``, ``viegas_wind``, ``canyon_wind``, ``viegas_neto``, ``pimont``, or ``windninja_ridge_canyon``
   * - ``wind_terrain.k_canyon``
     - 1.0
     - Terrain channeling coefficient for ``canyon_wind`` (Option 4)
   * - ``wind_terrain.k_pimont``
     - 0.5
     - Exponential slope coefficient for ``pimont`` (Option 6)

Viegas parameters (``viegas.a_V``, ``viegas.tan_phi_c``, ``viegas.T_a``,
``viegas.T_f``) are shared between the Viegas diagnostic model and the
buoyancy-wind options (2, 3, 5); see
`Viegas (2004) Eruptive Fire Diagnostics`_ above.

References
^^^^^^^^^^^

Rothermel, R.C. (1983). *How to Predict the Spread and Intensity of Forest and
Range Fires.* USDA Forest Service General Technical Report INT-143.

Viegas, D.X. & Neto, L.P.S. (1994). "Wind tunnel study of the convective air
flow of slope fires." Annual Report, Project COMOESTAS.

Pimont, F., Dupuy, J.-L., Linn, R.R. & Dupont, S. (2009). "Validation of
FIRETEC wind-flows over a canopy and fuel-break." *International Journal of
Wildland Fire*, 18(7), 775–790.
https://doi.org/10.1071/WF07130

WindNinja Ridge/Canyon Terrain Speed-up (Option 7)
---------------------------------------------------

The WindNinja empirical ridge/canyon model (Option 7,
``wind_terrain.model = windninja_ridge_canyon``) accounts for terrain-driven
wind acceleration observed in ridge and canyon topography (Forthofer 2007).

Wind-Slope Alignment
^^^^^^^^^^^^^^^^^^^^^

The alignment between the ambient wind and the upslope direction is:

.. math::

   a = \frac{\mathbf{U} \cdot \hat{\mathbf{n}}_s}{|\mathbf{U}|}

where :math:`\hat{\mathbf{n}}_s = \nabla z / |\nabla z|` is the unit upslope
vector and :math:`\mathbf{U}` is the horizontal wind velocity.  The alignment
satisfies :math:`a \in [-1,\,1]`:

* :math:`a > 0`: wind has an upslope component → **ridge** acceleration
* :math:`a < 0`: wind has a downslope component → **canyon** channeling
* :math:`a = 0`: wind is perpendicular to slope, no modification

Ridge Speed-up
^^^^^^^^^^^^^^^

When the wind climbs a slope, terrain convergence amplifies the wind speed:

.. math::

   f_{\text{ridge}} = 1 + k_{\text{ridge}} \,\tan\varphi \, a \quad (a > 0)

Canyon Channeling
^^^^^^^^^^^^^^^^^

When the wind flows down a slope (drainage flow channeling):

.. math::

   f_{\text{canyon}} = 1 + k_{\text{canyon,WN}} \,\tan\varphi \, |a| \quad (a < 0)

In both cases the effective wind is :math:`\mathbf{U}_{\text{eff}} = f \,\mathbf{U}`
with :math:`f \ge 1` (the model never suppresses wind speed).

Parameters:

* :math:`k_{\text{ridge}}` – ridge speed-up coefficient (default 1.0)
* :math:`k_{\text{canyon,WN}}` – canyon channeling coefficient (default 0.5)

Reference: Forthofer, J.M. (2007). *Modeling Wind in Complex Terrain for Use
in Fire Spread Prediction.* Colorado State University MS thesis.

Heat Flux MultiFab and Fire-Induced Wind
-----------------------------------------

A spatially-varying (or uniform) heat flux field :math:`Q` [W/m²] represents
the fire heat release rate at each cell.  It drives two WindNinja-style wind
corrections that are applied after any terrain-based modification.

Upward Convective Velocity
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The fire plume drives a vertical velocity proportional to the buoyancy of the
hot gas column above the fuel bed:

.. math::

   w_{\uparrow} = k_{\uparrow} \sqrt{\frac{g \, Q \, H_{\text{ref}}}{\rho_{\text{air}} \, C_p \, T_a}}

where:

* :math:`g = 9.81` m/s² – gravitational acceleration
* :math:`Q` – heat flux [W/m²] at the cell
* :math:`H_{\text{ref}}` – reference height [m] (``heat_flux.ref_height``, default 10 m)
* :math:`\rho_{\text{air}}` – air density [kg/m³] (``heat_flux.rho_air``, default 1.2)
* :math:`C_p` – specific heat of air [J/(kg·K)] (``heat_flux.Cp_air``, default 1005)
* :math:`T_a` – ambient temperature [K] (``heat_flux.T_a``, default 300)
* :math:`k_{\uparrow}` – upward velocity coefficient (``heat_flux.k_upward``, default 1.0)

In a 3-D build, :math:`w_{\uparrow}` is added directly to the vertical wind
component.  In a 2-D build the upward motion is projected onto the horizontal
plane as an inflow opposite to the ambient wind direction.

Induced Horizontal Inflow
^^^^^^^^^^^^^^^^^^^^^^^^^^

The fire-plume entrainment also creates a horizontal inflow toward the fire
perimeter:

.. math::

   U_{\text{ind}} = k_{\text{ind}} \, w_{\uparrow}

directed anti-parallel to the level-set gradient :math:`\nabla\phi`
(i.e., toward decreasing :math:`\phi`, which is the interior of the fire).
This term is only applied outside the fire perimeter (:math:`\phi \ge 0`).

Balbi Buoyancy Velocity Augmentation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For the Balbi spread model, the fire heat flux additionally augments the
fuel-derived buoyancy velocity:

.. math::

   v_{b,Q} = k_{\uparrow} \sqrt{\frac{g \, Q \, H_{\text{ref}}}{\rho_{\text{air}} \, C_p \, T_a}}

.. math::

   v_{b,\text{eff}} = \sqrt{v_{b,\text{fuel}}^2 + v_{b,Q}^2}

The combined :math:`v_{b,\text{eff}}` is then used in the Balbi tilt-angle
calculation.  A larger buoyancy velocity makes the flame more vertical (smaller
tilt angle), reducing radiant heat transfer to unburned fuel and thereby reducing
ROS.  This is physically consistent with observations of strongly buoyant plumes
above high-intensity fires.

Input Parameters
^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 30 10 60

   * - Parameter
     - Default
     - Description
   * - ``heat_flux.value``
     - 0.0
     - Uniform heat flux [W/m²]. Set > 0 to activate.
   * - ``heat_flux.file``
     - ""
     - ASCII file (X Y Q columns) for spatially-varying heat flux (2D only).
   * - ``heat_flux.rho_air``
     - 1.2
     - Air density [kg/m³]
   * - ``heat_flux.Cp_air``
     - 1005.0
     - Specific heat of air [J/(kg·K)]
   * - ``heat_flux.T_a``
     - 300.0
     - Ambient temperature [K]
   * - ``heat_flux.ref_height``
     - 10.0
     - Reference height for buoyancy velocity [m]
   * - ``heat_flux.k_upward``
     - 1.0
     - Upward velocity coefficient
   * - ``heat_flux.k_induced``
     - 0.5
     - Induced horizontal inflow coefficient
   * - ``heat_flux.enable_upward``
     - 0
     - 1 to enable upward velocity term
   * - ``heat_flux.enable_induced``
     - 0
     - 1 to enable induced horizontal inflow term

Level-Set Method
----------------

The level-set function :math:`\phi(x,y,t)` represents the fire front as the zero level set :math:`\phi = 0`. The evolution equation is:

.. math::

   \frac{\partial \phi}{\partial t} + V|\nabla\phi| = 0

where :math:`V` is the local rate of spread computed from the fire models above.

Numerical Scheme
^^^^^^^^^^^^^^^^

The level-set equation is discretized using:

* Upwind finite differences for spatial derivatives
* Explicit time stepping (e.g., forward Euler or RK2)
* Reinitialization to maintain signed distance property

The gradient :math:`|\nabla\phi|` is computed using:

.. math::

   |\nabla\phi| = \sqrt{\left(\frac{\partial\phi}{\partial x}\right)^2 + \left(\frac{\partial\phi}{\partial y}\right)^2}

Post-Processing Diagnostics
-----------------------------

In addition to fireline intensity and flame length, the solver always writes four fire-ecology
fields and three emissions fields to every plotfile.  The equations below describe how each
field is computed.

Scorch Height (Van Wagner 1973)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Van Wagner's (1973) empirical relationship between fireline intensity and the height to which
foliage is killed by convective heat:

.. math::

   H_s \;[\text{m}] = 0.1483 \; I_B^{2/3}

where :math:`I_B` [kW/m] is the Byram fireline intensity.  Cells where :math:`I_B = 0` receive
:math:`H_s = 0`.  Output field: ``scorch_height``.

Probability of Ignition (Anderson 1970)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The probability that a firebrand landing on fine fuel will produce sustained ignition (Anderson
1970, Rothermel 1983):

.. math::

   P_i = \min\!\left(1,\; \max\!\left(0,\; 4.8 \times 10^{-5} \, T_{\text{fuel,F}}^{1.4}
         \exp\!\left(-0.07 \, M_{\%}\right)\right)\right)

where:

* :math:`T_{\text{fuel,F}} = T_{a,\text{F}} + \Delta T_{\text{solar}}` — fine-fuel temperature
  in °F (ambient temperature converted from °C plus a solar-heating increment; default :math:`+25` °F)
* :math:`M_{\%}` — 1-hr dead fuel moisture content [%]

Output field: ``prob_ignition``.

Tree Mortality (Ryan & Reinhardt 1988)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A logistic mortality model based on the fraction of the live crown that is scorched.

**Crown Scorch Ratio (CSR)**:

.. math::

   \text{CSR} = \text{clip}\!\left(\frac{H_s - H_{\text{CBH}}}{H_{\text{tree}} - H_{\text{CBH}}},\;0,\;1\right)

where :math:`H_{\text{CBH}}` is the canopy base height [m] (``crown.CBH``) and
:math:`H_{\text{tree}}` is the tree height [m] (``fire_ecology.tree_height``).

**Mortality probability** (only when :math:`H_s > H_{\text{CBH}}`):

.. math::

   P_{\text{kill}} = \frac{1}{1 + \exp\!\bigl(-4 \,(\text{CSR} - 0.6)\bigr)}

The logistic curve gives ≈50 % mortality when 60 % of the live crown is scorched.
Output field: ``tree_mortality``.

Crown Activity Class (Van Wagner 1977)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The crown activity class categorises each cell as surface, passive crown, or active crown fire:

.. list-table::
   :header-rows: 1
   :widths: 10 25 65

   * - Class
     - Label
     - Criterion
   * - 0
     - Surface fire
     - :math:`I_B < I_o`
   * - 1
     - Passive crown
     - :math:`I_B \ge I_o` **and** :math:`R < R'_{SA}`
   * - 2
     - Active crown
     - :math:`I_B \ge I_o` **and** :math:`R \ge R'_{SA}`

**Crown initiation intensity** (Van Wagner 1977):

.. math::

   I_o \;[\text{kW/m}] = 0.010 \; H_{\text{CBH}} \,(460 + 25.9 \, F_{\text{MC}})

**Critical active-crown ROS** (Van Wagner 1977):

.. math::

   R'_{SA} \;[\text{m/s}] = \frac{3.0}{C_{\text{BD}} \times 60}

where :math:`C_{\text{BD}}` [kg/m³] is the canopy bulk density (``crown.CBD``) and
:math:`F_{\text{MC}}` [%] is the foliar moisture content (``crown.FMC``).
Output field: ``crown_activity``.

Fire Emissions (Seiler & Crutzen 1980)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Per-cell cumulative emissions for three species, following the WRF-Fire convention:

.. math::

   E_{\text{species}} \;[\text{kg/m}^2] = W_{\text{fuel}} \;[\text{kg/m}^2]
   \times f_c \;[-] \times \text{EF}_{\text{species}} \;[\text{kg/kg}]

where:

* :math:`W_{\text{fuel}}` — oven-dry fuel load [kg/m²] (converted from lb/ft²)
* :math:`f_c` — fuel consumption fraction from the bulk model, or ``emissions.default_consumed_frac``
  (default 0.9) for burned cells without the bulk model active
* :math:`\text{EF}` — emission factors [kg species per kg dry fuel consumed]:

  - CO₂: 1.570 kg/kg (Seiler & Crutzen 1980)
  - CO:  0.102 kg/kg
  - PM₂.₅: 0.0162 kg/kg (Andreae & Merlet 2001)

Output fields: ``co2_emissions``, ``co_emissions``, ``pm25_emissions``.

References
^^^^^^^^^^^

* Van Wagner, C.E. (1973). "Height of crown scorch in forest fires." *Canadian Journal of
  Forest Research*, 3(3), 373–378.
* Anderson, H.E. (1970). "Forest fuel ignitibility." *Fire Technology*, 6(4), 312–319.
* Ryan, K.C. & Reinhardt, E.D. (1988). "Predicting postfire mortality of seven western
  conifers." *Canadian Journal of Forest Research*, 18(10), 1291–1297.
* Seiler, W. & Crutzen, P.J. (1980). "Estimates of gross and net fluxes of carbon between the
  biosphere and the atmosphere from biomass burning." *Climatic Change*, 2(3), 207–247.
* Andreae, M.O. & Merlet, P. (2001). "Emission of trace gases and aerosols from biomass
  burning." *Global Biogeochemical Cycles*, 15(4), 955–966.
