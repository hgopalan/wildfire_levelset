Mathematical Models
===================

This section describes the mathematical models implemented in the wildfire level-set solver.

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
