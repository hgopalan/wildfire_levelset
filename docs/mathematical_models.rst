Mathematical Models
===================

This section describes the mathematical models implemented in the wildfire level-set solver.


Fire Spread Models
------------------

Multiple fire spread model families are implemented.  Select the active
model with ``fire_spread_model``:

* ``rothermel`` — Rothermel (1972) empirical surface fire spread (default)
* ``balbi`` — Balbi (2009) physics-based radiation-driven spread
* ``cheney_gould`` — Cheney & Gould (1995/1998) Australian grassland empirical model
* ``fbp_o1a`` / ``fbp_o1b`` — Canadian FBP grass fuel types (matted / standing)
* ``fbp_s1`` / ``fbp_s2`` / ``fbp_s3`` — Canadian FBP slash fuel types
* ``lautenberger`` — Lautenberger (2013) semi-physical Eulerian model
* ``cruz_crown`` — Cruz et al. active crown fire rate of spread

Rothermel (1972) Fire Spread Model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The Rothermel fire spread model is the foundation for computing the rate of spread (ROS) of wildfire. The model computes the no-wind, no-slope rate of spread and then applies wind and slope correction factors.

Basic Rate of Spread
~~~~~~~~~~~~~~~~~~~~~

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
~~~~~~~~~~~~~~~~~~

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
~~~~~~~~~~~~~

The packing ratio is defined as:

.. math::

   \beta = \frac{\rho_b}{\rho_p}

where:

* :math:`\rho_b = w_0 / \delta` is the oven-dry bulk density (lb/ft³)
* :math:`\rho_p` is the oven-dry particle density (lb/ft³)
* :math:`w_0` is the oven-dry fuel loading (lb/ft²)
* :math:`\delta` is the fuel bed depth (ft)

Moisture Damping
~~~~~~~~~~~~~~~~

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
~~~~~~~~~~~~~~~

The mineral damping coefficient :math:`\eta_s` is:

.. math::

   \eta_s = 0.174 S_e^{-0.19}

where :math:`S_e` is the effective mineral content.

Propagating Flux Ratio
~~~~~~~~~~~~~~~~~~~~~~~

The propagating flux ratio :math:`\xi` is:

.. math::

   \xi = \frac{\exp[(0.792 + 0.681\sqrt{\sigma})(\beta + 0.1)]}{192 + 0.2595\sigma}

Heat of Preignition
~~~~~~~~~~~~~~~~~~~

The heat of preignition :math:`Q_{ig}` (BTU/lb) depends on fuel moisture:

.. math::

   Q_{ig} = 250 + 1116 M_f

Wind Factor
~~~~~~~~~~~

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
~~~~~~~~~~~~

The slope factor :math:`\phi_s` accounts for terrain slope:

.. math::

   \phi_s = 5.275 \beta^{-0.3} \tan^2(\theta)

where :math:`\theta` is the slope angle. The slope can be decomposed into x and y components:

.. math::

   \tan^2(\theta) = \left(\frac{\partial z}{\partial x}\right)^2 + \left(\frac{\partial z}{\partial y}\right)^2

Combined Wind and Slope
~~~~~~~~~~~~~~~~~~~~~~~~

When both wind and slope effects are present, the total rate of spread is:

.. math::

   R = R_0 (1 + \phi_w + \phi_s)


Balbi (2009) Physical Fire Spread Model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The Balbi model provides an alternative rate-of-spread calculation rooted in
radiation physics rather than empirical curve fits.  It is activated by setting
``fire_spread_model = balbi`` and automatically replaces the Rothermel model in
both the level-set advection path and the pre-computed :math:`R` field used by
FARSITE.  The model is fully compatible with all six Viegas wind-terrain options
and FARSITE landscape file inputs.

Physical Basis
~~~~~~~~~~~~~~~

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
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When ``fire_spread_model = balbi``, the solver automatically prints a table of
pre-computed Balbi parameters for every fuel model in the active database
(FBFM13 or FBFM40) at startup.  Columns are: fuel code, short name,
:math:`\sigma_m`, :math:`\delta_m`, :math:`\chi`, :math:`B^*`, :math:`v_b`,
:math:`A_{\text{coeff}}`, and the predicted ROS at 5 m/s wind on flat ground.

Reference
~~~~~~~~~~

Balbi, J.-H., Rossi, J.-L., Marcelli, T., and Santoni, P.-A. (2009).
"A physical model for wildland fires."
*Combustion and Flame*, 156(12), 2217–2230.
https://doi.org/10.1016/j.combustflame.2009.07.010



Cheney & Gould (1995 / 1998) Grassland Fire Spread Model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The Cheney–Gould model is a purely empirical rate-of-spread (ROS) formula calibrated against a large number of experimental grassland fires conducted in Australia. It is activated by setting ``fire_spread_model = cheney_gould`` and generally outperforms the Rothermel model in open grassland fuels.

Head-Fire Rate of Spread
~~~~~~~~~~~~~~~~~~~~~~~~~

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
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. math::

   f_{MC} = \exp(-0.108 \times MC)

where :math:`MC` is the dead fine fuel moisture content [%]. At :math:`MC = 0`
(bone-dry fuel) :math:`f_{MC} = 1`; at :math:`MC \approx 22\%` the ROS is
halved relative to dry conditions.

Curing Factor
~~~~~~~~~~~~~

:math:`CF \in [0,\,1]` represents the degree of curing of the grass:
:math:`CF = 1` is fully cured (dry standing grass), while :math:`CF = 0` is
completely green (no spread). The factor is clamped to the :math:`[0,1]`
interval at runtime.

Note on Terrain Slope
~~~~~~~~~~~~~~~~~~~~~

Terrain slope is **not** accounted for in the original Cheney–Gould empirical
formulation. The model was calibrated on flat or gently sloping Australian
grasslands and the slope correction is intentionally omitted in this
implementation. For slope effects, use the Rothermel or Balbi models instead.


Canadian Forest Fire Behaviour Prediction (FBP) System
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The Canadian FBP system provides empirical rate-of-spread equations calibrated for
specific fuel types.  The solver implements five fuel types: O1a (matted grass),
O1b (standing grass), S1 (Jack/Lodgepole Pine slash), S2 (White Spruce/Balsam slash),
and S3 (Coastal Cedar/Hemlock/Douglas Fir slash).

Select with ``fire_spread_model = fbp_o1a`` (or ``fbp_o1b``, ``fbp_s1``, ``fbp_s2``, ``fbp_s3``).

**Grass fuel types (O1a / O1b)**

The fine-fuel moisture factor follows Van Wagner (1985):

.. math::

   f_F = 91.9 \exp(-0.1386\,M_f)
         \left(1 + \frac{M_f^{5.31}}{4.93 \times 10^7}\right)

The wind factor and Initial Spread Index (ISI) are:

.. math::

   f_W = \exp(0.05039\,U_{10}), \qquad \text{ISI} = 0.208 \, f_W \, f_F

where :math:`U_{10}` is the 10-m wind speed (km/h) and :math:`M_f` is the fine fuel moisture content (%).
The curing factor (CF) accounts for the fraction of dead fuel:

.. math::

   \mathrm{CF} = \Bigl(1 - \exp\!\bigl(-2.1\,\mathrm{PC}/100\bigr)\Bigr)^2

Rate of spread [m/s]:

.. math::

   R = C_c \exp(k_c\,\text{ISI}) \times \mathrm{CF}

Coefficients per fuel type:

* **O1a** (matted grass): :math:`C_c = 190/60,\; k_c = 0.031`
* **O1b** (standing grass): :math:`C_c = 250/60,\; k_c = 0.035`

**Slash fuel types (S1 / S2 / S3)**

Same ISI computation; :math:`\mathrm{CF} = 1`. Rate of spread [m/s]:

.. math::

   R = C_s \exp(k_s\,\text{ISI})

Coefficients:

* **S1**: :math:`C_s = 75/60,\; k_s = 0.110`
* **S2**: :math:`C_s = 200/60,\; k_s = 0.062`
* **S3**: :math:`C_s = 320/60,\; k_s = 0.010`

Parameters: ``fbp.fuel_type``, ``fbp.curing`` (curing [%], O1a/O1b only), ``fbp.moisture`` (dead fine fuel moisture [%]).

Reference: Forestry Canada Fire Danger Group (1992). *Development and Structure of
the Canadian Forest Fire Behavior Prediction System.*  Information Report ST-X-3.


Lautenberger (2013) Physics-Based Fire Spread Model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Select with ``fire_spread_model = lautenberger``.

The model derives rate of spread from a semi-physical fuel-energy balance
calibrated against US wildfire data (Lautenberger 2013):

.. math::

   R_0 = A_L \cdot \sigma_m \cdot \exp(-B_L \cdot M_f)

where :math:`\sigma_m` is the fuel SAV ratio [m\ :sup:`-1`] (converted
from ft\ :sup:`-1`), :math:`M_f` is fuel moisture (fraction), :math:`A_L`
is a pre-factor [m²/s], and :math:`B_L` is the moisture sensitivity.

Wind and slope corrections are applied linearly:

.. math::

   \phi_W = C_L \cdot U \qquad \phi_S = D_L \cdot \tan\theta

   R = \max\bigl(R_0 \cdot (1 + \phi_W + \phi_S),\; 0\bigr)

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - Parameter
     - Default
     - Description
   * - ``lautenberger.A_L``
     - ``1.05e-5``
     - Pre-factor :math:`A_L` [m²/s]
   * - ``lautenberger.B_L``
     - ``2.5``
     - Moisture sensitivity :math:`B_L` [-]
   * - ``lautenberger.C_L``
     - ``0.40``
     - Wind correction :math:`C_L` [(m/s)\ :sup:`-1`]
   * - ``lautenberger.D_L``
     - ``0.50``
     - Slope correction :math:`D_L` [-]

Reference: Lautenberger, C. (2013). *Wildland fire modeling with an Eulerian level
set method and automated calibration.* Fire Safety Journal, 62, 289–298.


ROS Floor (FARSITE Stall Threshold)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A FARSITE-style stall threshold is applied across all spread models: cells where
the computed ROS is below ``WildfireConst::ROS_MIN_MS`` (1 × 10⁻⁴ m/s) are
forced to zero.  This prevents numerical drift in extremely low-moisture conditions
and matches FARSITE's internal practice of not propagating the fire front when
spread is negligible.

**Parameters**: None (constant in ``constants.H``).

Non-Burnable Cell Masking
^^^^^^^^^^^^^^^^^^^^^^^^^^

Fuel model codes that correspond to non-burnable landscape types are explicitly
zeroed in the ROS kernel.  This prevents fire from numerically creeping across
water bodies, bare rock, or urban areas due to floating-point noise in zero-fuel
cells.

Affected codes (FBFM40 / FBFM13): 91 (Urban/Developed), 92 (Snow/Ice), 93
(Agriculture), 98 (Open Water), 99 (Bare Ground).

When ``rothermel.landscape_file`` is provided and the per-cell fuel code maps to
one of these values, the kernel returns ``R = 0`` for that cell without entering
the wind-slope-fuel computation.

**Parameters**: None.  Automatic when ``rothermel.landscape_file`` is used.


Propagation Methods
-------------------

Three propagation methods are available and selected with ``propagation_method``:
``levelset`` (WENO5-Z level-set advection), ``farsite`` (FARSITE Huygens-wavelet
expansion), and ``mtt`` (Minimum Travel Time fast-march).

Level-Set Method
^^^^^^^^^^^^^^^^

The level-set function :math:`\phi(x,y,t)` represents the fire front as the zero level set :math:`\phi = 0`. The evolution equation is:

.. math::

   \frac{\partial \phi}{\partial t} + V|\nabla\phi| = 0

where :math:`V` is the local rate of spread computed from the fire models above.

Numerical Scheme
~~~~~~~~~~~~~~~~

The level-set equation is discretized using:

* Upwind finite differences for spatial derivatives
* Explicit time stepping (e.g., forward Euler or RK2)
* Reinitialization to maintain signed distance property

The gradient :math:`|\nabla\phi|` is computed using:

.. math::

   |\nabla\phi| = \sqrt{\left(\frac{\partial\phi}{\partial x}\right)^2 + \left(\frac{\partial\phi}{\partial y}\right)^2}


Terrain-Corrected Gradient
~~~~~~~~~~~~~~~~~~~~~~~~~~

When a terrain file is present, the gradient is computed on the terrain
surface rather than on the flat horizontal plane.  The horizontal grid
spacing :math:`\Delta x` corresponds to a surface arc length of
:math:`\Delta s_x = \Delta x \sqrt{1 + (\partial z/\partial x)^2}`, and
similarly in :math:`y`.  Dividing finite differences by these arc-length
spacings gives the surface gradient components:

.. math::

   \frac{\partial\phi}{\partial s_x} = \frac{\Delta\phi}{\Delta x \sqrt{1 + s_x^2}}, \qquad
   \frac{\partial\phi}{\partial s_y} = \frac{\Delta\phi}{\Delta y \sqrt{1 + s_y^2}}

where :math:`s_x = \partial z / \partial x` and
:math:`s_y = \partial z / \partial y` are the terrain slope components stored
in the two-component slopes ``MultiFab``.

The terrain-corrected gradient magnitude is then:

.. math::

   |\nabla_s\phi| = \sqrt{
       \left(\frac{\partial\phi}{\partial s_x}\right)^2 +
       \left(\frac{\partial\phi}{\partial s_y}\right)^2}

This correction is applied inside ``godunov_norm_grad_phi``
(``src/numerical_schemes.H``) using effective spacings
:math:`\Delta x_\text{eff} = \Delta x \sqrt{1 + s_x^2}` and
:math:`\Delta y_\text{eff} = \Delta y \sqrt{1 + s_y^2}` in both the WENO3
and the first-order fallback stencils.  Reinitialization uses the default
flat-terrain path (slope = 0) because the signed-distance property is defined
in the horizontal plane.

For steep terrain (e.g., slope magnitude 1.2, equivalent to ~50°) the
effective spacing is 57 % larger than the flat value, so the surface gradient
magnitude is correspondingly smaller.  Without this correction the gradient
would be over-estimated on steep flanks, causing the fire to propagate too
slowly upslope.



FARSITE Elliptical Expansion Model (Richards 1990)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The FARSITE model uses an elliptical expansion pattern to represent the anisotropic spread of fire under the influence of wind and terrain.

Richards' Coefficients
~~~~~~~~~~~~~~~~~~~~~~~

The elliptical fire spread is characterized by three coefficients:

* :math:`a` - head fire coefficient (maximum spread, downwind)
* :math:`b` - flank fire coefficient (crosswind spread)
* :math:`c` - backing fire coefficient (minimum spread, upwind)

These coefficients relate the directional rate of spread to the base rate:

.. math::

   R(\theta) = R_0 f(a, b, c, \theta)

where :math:`\theta` is the angle between the fire normal and the wind direction.

Anderson Length-to-Width Ratio
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Anderson (1983) model relates the ellipse shape to wind speed:

.. math::

   L/W = 0.936 \exp(0.2566 U) + 0.461 \exp(-0.1548 U) - 0.397

where :math:`L/W` is the length-to-width ratio and :math:`U` is the wind speed (mph). This ratio is then converted to the Richards' coefficients for the elliptical spread calculation.


Alternative Fire Shape Models
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Four fire shape models are available for the FARSITE propagation path.  All
share the same three user coefficients :math:`a` (head), :math:`b` (flank), and
:math:`c` (backing), and the base rate of spread :math:`R_{head}` from the
active fire spread model (Rothermel, Balbi, etc.).  Select the shape with:

.. code-block:: ini

   farsite.fire_shape_model = richards           # default
   # or: catchpole_demestre | wilson | lemniscate


Richards (1990) – Linear Blend (default)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The default model is a simple linear combination of head, flank, and backing
spread distances:

.. math::

   r(\theta) = \begin{cases}
       (a \cos\theta + b \sin\theta) \, R_{head} \, \Delta t & \cos\theta \ge 0 \\
       (b \sin\theta - c \cos\theta) \, R_{head} \, \Delta t & \cos\theta < 0
   \end{cases}

Enabled by ``farsite.fire_shape_model = richards`` (default).


Catchpole & de Mestre (1986) – True Double-Ellipse
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Each Huygens wavelet is composed of **two half-ellipses** joined at the
crosswind axis.  The head half-ellipse has semi-axes :math:`a` and :math:`b`;
the backing half has semi-axes :math:`c` and :math:`b`.

The polar spread distance from the fire-front origin is:

.. math::

   r(\theta) = \begin{cases}
       \dfrac{a\,b}{\sqrt{b^2\cos^2\theta + a^2\sin^2\theta}} \, R_{head} \, \Delta t
           & \cos\theta \ge 0 \text{ (head half)} \\[8pt]
       \dfrac{c\,b}{\sqrt{b^2\cos^2\theta + c^2\sin^2\theta}} \, R_{head} \, \Delta t
           & \cos\theta < 0 \text{ (backing half)}
   \end{cases}

This model is geometrically exact for the double-ellipse shape and gives:

* :math:`r(0) = a\,R_{head}\,\Delta t` (head fire rate)
* :math:`r(\pi/2) = b\,R_{head}\,\Delta t` (flank rate)
* :math:`r(\pi) = c\,R_{head}\,\Delta t` (backing rate)

Enabled by ``farsite.fire_shape_model = catchpole_demestre``.

**Reference:** Catchpole, E.A. & de Mestre, N.J. (1986). *Australian Forestry*, 49(2), 102–111.


Wilson (1988) – Single Ellipse from Rear Focus
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Wilson's model places the fire origin at the **rear focus** of a single
ellipse.  The ellipse geometry is derived from the head (:math:`a`) and backing
(:math:`c`) coefficients only; the flank spread follows naturally.

Defining:

.. math::

   R_H = a\,R_{head}\,\Delta t, \quad R_B = c\,R_{head}\,\Delta t

the ellipse parameters are:

.. math::

   a_{ell} = \frac{R_H + R_B}{2}, \quad
   e = \frac{R_H - R_B}{R_H + R_B}, \quad
   \ell = a_{ell}(1 - e^2) = \frac{2\,R_H\,R_B}{R_H + R_B}

The conic-section polar equation from the rear focus gives:

.. math::

   r(\theta) = \frac{\ell}{1 - e\cos\theta}

At the principal angles:

* :math:`r(0) = R_H = a\,R_{head}\,\Delta t` (head fire)
* :math:`r(\pi) = R_B = c\,R_{head}\,\Delta t` (backing fire)
* :math:`r(\pi/2) = \ell = \dfrac{2R_H R_B}{R_H + R_B}` (harmonic mean — flank)

The :math:`b` coefficient is not used by this model; the flank spread is
determined solely by :math:`a` and :math:`c`.

Enabled by ``farsite.fire_shape_model = wilson``.

**Reference:** Wilson, A.A.G. (1988). *Canadian Journal of Forest Research*, 18(6), 682–687.


Alexander et al. – Lemniscate (Limaçon)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The lemniscate model represents the fire wavelet as a Limaçon (cardioid
family):

.. math::

   r(\theta) = \left[b + \frac{a - c}{2}\cos\theta\right] R_{head}\,\Delta t

Principal values:

* :math:`r(0) = \left(b + \tfrac{a-c}{2}\right) R_{head}\,\Delta t`
* :math:`r(\pi/2) = b\,R_{head}\,\Delta t` (flank rate equals :math:`b` exactly)
* :math:`r(\pi) = \left(b - \tfrac{a-c}{2}\right) R_{head}\,\Delta t`

When the natural choice :math:`b = (a + c)/2` is used, the formula reduces to:

.. math::

   r(\theta) = \frac{a(1+\cos\theta) + c(1-\cos\theta)}{2}\,R_{head}\,\Delta t

giving :math:`r(0) = a`, :math:`r(\pi) = c`, and a flank rate that is the
arithmetic mean of head and backing.  Values below zero are clamped to zero.

Enabled by ``farsite.fire_shape_model = lemniscate``.

**Reference:** Alexander, M.E., Stocks, B.J., & Lawson, B.D. (1991). *Information Report NOR-X-310*, Canadian Forest Service.

Minimum Travel Time (MTT) Propagation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The Minimum Travel Time method (`propagation_method = mtt`) pre-computes a
scalar **arrival-time field** :math:`T(i,j,k)` — the time at which fire first
reaches cell :math:`(i,j,k)` — using a Dijkstra fast-marching sweep over the
precomputed ROS field.  The level-set field is then updated analytically at
each simulation time :math:`t`:

.. math::

   \phi(i,j,k,t) = T(i,j,k) - t

giving:

* :math:`\phi < 0` — burned (arrival before current time)
* :math:`\phi = 0` — fire front (arriving now)
* :math:`\phi > 0` — unburned (fire not yet arrived)

Arrival Time Computation
~~~~~~~~~~~~~~~~~~~~~~~~~

Starting from all ignition cells (where initial :math:`\phi < 0`,
:math:`T = 0`), the arrival time propagates to neighbouring cells via:

.. math::

   T_{\text{new}} = T_{\text{current}} + \frac{d_s}{R_{\text{interface}}}

where :math:`d_s` is the cell-centre-to-cell-centre distance and
:math:`R_{\text{interface}}` is the **harmonic mean** of the two adjacent
ROS values:

.. math::

   R_{\text{interface}} = \frac{2\,R_a\,R_b}{R_a + R_b}

The harmonic mean is chosen to penalise paths through low-ROS cells and
prevents numerical blow-up when one cell has very small ROS.

A standard min-heap priority queue (Dijkstra) guarantees that the globally
minimum arrival time is always processed next, making the result exact for
any non-negative, spatially-varying ROS field.

Advantages and Limitations
~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Advantages** relative to the level-set and FARSITE methods:

* The arrival-time field is monotonically increasing; no reinitialization is
  needed.
* The phi update is trivially parallel (one subtraction per cell per step).
* The fast-march is run only once, regardless of the number of time steps.

**Limitations**:

* The ROS field is frozen at the start; time-varying wind, moisture, or
  diurnal models are not captured after the initial sweep.
* The method is isotropic (spread depends only on local ROS magnitude, not
  direction). Wind direction effects are included through the directional ROS
  model but fire does not curve around obstacles post-computation.

Select MTT with::

    propagation_method = mtt

It is compatible with all fire spread models (``rothermel``, ``balbi``,
``cheney_gould``, ``cruz_crown``).

**Reference:** Finney, M.A. (2002). "Fire growth using minimum travel time
methods." *Canadian Journal of Forest Research*, 32(8), 1420–1424.
https://doi.org/10.1139/x02-068


Crown Models
------------

Van Wagner (1977) Crown Fire Initiation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Crown fire initiation occurs when sufficient heat is available to raise the canopy fuel to ignition temperature.

Crown Fire Initiation Criterion
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The critical surface fire intensity :math:`I_0` (kW/m) required for crown fire initiation is:

.. math::

   I_0 = \frac{CBH \times (460 + 25.9 M_c)}{18 h}

where:

* :math:`CBH` is the canopy base height (m)
* :math:`M_c` is the canopy moisture content (%)
* :math:`h` is the crown fuel consumption depth (m)

If the actual surface fire intensity :math:`I > I_0`, then crown fire initiation occurs.

Active Crown Fire Criterion
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For active crown fire spread, the critical wind speed :math:`U_{active}` is:

.. math::

   U_{active} = \frac{3}{\sqrt{CBD}}

where :math:`CBD` is the canopy bulk density (kg/m³).


Rothermel (1991) Active Crown Fire ROS
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When ``crown.use_rothermel1991_crown = 1``, active-crown-fire cells use a simple
multiplicative relationship derived by Rothermel (1991):

.. math::

   R_{\text{crown}} = 3.34 \; R_{\text{surface}}

This is applied only in cells where ``crown_activity == 2`` (active crown fire),
and replaces (or is blended with) the surface ROS.

Parameters: ``crown.use_rothermel1991_crown = 1``.

Reference: Rothermel, R.C. (1991). *Predicting Behavior and Size of Crown Fires
in the Northern Rocky Mountains.* USDA Forest Service Research Paper INT-438.


Van Wagner (1977) Passive-Crown Blending
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When ``crown.use_passive_blend = 1`` the transition from surface to active crown
fire is smoothed via Van Wagner's (1977) passive crown-fire factor:

.. math::

   \mathrm{CF} = \left(\frac{I_B}{I_o}\right)^{2/3}

The blended ROS is:

.. math::

   R = (1 - \mathrm{CF})\,R_{\text{surface}} + \mathrm{CF}\,R_{\text{crown}}

where :math:`I_B` is the Byram fireline intensity (kW/m) and :math:`I_o` is the
Van Wagner crown fire initiation threshold (kW/m):

.. math::

   I_o = 0.010\,h_{\text{CBH}}\,(460 + 25.9\,\mathrm{FMC})

Parameters: ``crown.use_passive_blend = 1``, ``crown.CBH``, ``crown.FMC``.

Reference: Van Wagner, C.E. (1977). *Conditions for the start and spread of crown
fire.* Canadian Journal of Forest Research, 7(1), 23–34.


Scott & Reinhardt (2001) Bisection-Based TI and CI
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The full Scott & Reinhardt (2001) Torching Index (TI) and Crowning Index (CI) are
defined as the minimum 20-ft (6.1 m) open-wind speed (km/h) at which passive torching or
active crowning can occur.  These are output as the plotfile fields
``torching_index_kmh`` and ``crowning_index_kmh``.

Enable with ``scott_reinhardt_full.enable = 1``.

**Torching Index** is found by bisection on the wind speed :math:`U` that satisfies:

.. math::

   I_B(U) = I_o

where :math:`I_B(U) = H_l \, w_n \, R(U) / 60` is the Byram fireline intensity
(kW/m) computed from the Rothermel surface ROS :math:`R(U)` [m/min], net fuel
load :math:`w_n` [kg/m²], and heat of combustion :math:`H_l` [kJ/kg].

**Crowning Index** is found by bisection on the wind speed :math:`U` that satisfies:

.. math::

   R(U) = R'_{SA} = \frac{3.0}{\mathrm{CBD}} \; \text{[m/min]}

Reference: Scott, J.H. & Reinhardt, E.D. (2001). *Assessing Crown Fire Potential
by Linking Models of Surface and Crown Fire Behavior.* USDA Forest Service
Research Paper RMRS-RP-29.


Spotting Models
---------------

Firebrand Spotting
^^^^^^^^^^^^^^^^^^

Firebrand spotting generates new ignition points ahead of the main fire front.
Two independent spotting models are available; they can be enabled separately or
together.

Probability-Based Spotting (``spotting.*``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A stochastic model driven by wind speed, fire intensity, and fuel moisture.
For each fire-front cell a spotting probability is computed and a Bernoulli draw
decides whether a firebrand is generated.  If so, the landing distance is
sampled from a log-normal or exponential distribution and the spot is placed
downwind with optional lateral dispersion.

Albini (1983) Spotting with 2-D Trajectory (``albini_spotting.*``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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


Flux-Based Ember Cascade Model (``ember_cascade.*``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The flux-based ember cascade model complements the single-brand Albini trajectory
approach by treating dense firebrand showers as a **continuous landing-flux
field** rather than tracking individual particles.  It is appropriate when fire
intensity is high enough to generate large, simultaneous ember cascades from the
convective plume — a regime where the Monte-Carlo single-brand treatment becomes
computationally or physically inadequate (Sardoy et al. 2007).

**Stage 1 – Plume lofting height (Albini 1983)**

Same as the trajectory model:

.. math::

   H_z \;[\text{m}] = 12.2 \; I_B^{1/3}

where :math:`I_B` [kW/m] is the Byram fireline intensity.  Cells with
:math:`I_B < I_{B,\min}` produce no embers.

**Stage 2 – Ember source flux**

The volumetric firebrand production rate per unit fire-line length is
parameterised as a power-law of intensity:

.. math::

   F_{\text{src}} \;[\text{embers/m/s}] = k_{\text{flux}} \left(\frac{I_B}{I_{B,\text{ref}}}\right)^{\alpha}

where ``k_flux``, ``I_B_ref``, and ``flux_exp`` (:math:`\alpha`) are user
parameters.

**Stage 3 – Wind-advected Gaussian landing kernel**

Flight time and mean transport displacement:

.. math::

   t_f = \frac{H_z}{v_t}, \quad
   \delta x = u \, t_f, \quad
   \delta y = v \, t_f

Turbulent diffusion from the plume broadens the ember cloud:

.. math::

   \sigma = \sigma_{\text{base}} + k_\sigma \, H_z

The landing-flux density at cell :math:`(x, y)` from source :math:`s` is:

.. math::

   \delta F_{\text{land}}(x,y)
     = F_{\text{src},s}
       \cdot \frac{A_{\text{cell}}}{2\pi\sigma_s^2}
       \cdot \exp\!\left(-\frac{(\Delta x_s)^2 + (\Delta y_s)^2}{2\sigma_s^2}\right)

where :math:`\Delta x_s = x - (x_s + \delta x_s)`.  Only cells within
:math:`n_{\sigma,\text{cut}} \times \sigma_s` of the mean landing centre
receive contributions (Gaussian truncation for efficiency).  The total landing
flux is the sum over all active source cells:

.. math::

   F_{\text{land}}(x,y) = \sum_s \delta F_{\text{land},s}(x,y)

**Stage 4 – Stochastic Poisson ignition**

During each check interval of duration :math:`\Delta t_{\text{check}}`, each
cell is treated as an independent Poisson receiver:

.. math::

   \lambda = \frac{F_{\text{land}} \, \Delta t_{\text{check}}}{N_{\text{min}}}

   P_{\text{ign}} = 1 - e^{-\lambda}

A uniform random draw on :math:`[0,1]` is compared against :math:`P_{\text{ign}}`;
cells that pass receive a spot-fire ignition of radius ``spot_radius``.
``N_min_density`` [embers/m²/s] controls the ignition sensitivity.

**3-D wind transport**

When ``ember_cascade.use_3d_wind = 1`` (or ``require_3d_wind = 1``), the
height-averaged horizontal wind from a
`massconsistent_amr <https://github.com/hgopalan/massconsistent_amr>`_ plotfile
is used for :math:`(u, v)` in the transport displacement.  Setting
``require_3d_wind = 1`` causes the simulation to abort with a diagnostic message
if no valid 3-D wind file is provided.

**Diagnostic output fields**

* ``ember_cascade_flux`` – accumulated ember landing-flux density
  :math:`F_{\text{land}}` [embers/m²/s] at each cell.
* ``ember_cascade_ignition`` – flag (1.0) at cells that received a spot-fire
  ignition during the most recent check interval.

**Input parameters** (prefix ``ember_cascade.``):

.. list-table::
   :header-rows: 1
   :widths: 30 10 60

   * - Parameter
     - Default
     - Description
   * - ``enable``
     - 0
     - 1 to activate flux-based ember cascade model
   * - ``I_B_min``
     - 10.0
     - Minimum Byram fireline intensity [kW/m] to emit embers
   * - ``terminal_velocity``
     - 1.0
     - Firebrand terminal descent velocity [m/s]
   * - ``k_flux``
     - 1.0
     - Ember source flux coefficient
   * - ``I_B_ref``
     - 100.0
     - Reference intensity [kW/m] for flux scaling
   * - ``flux_exp``
     - 1.0
     - Intensity exponent :math:`\alpha`
   * - ``sigma_base``
     - 50.0
     - Minimum Gaussian spread radius [m]
   * - ``k_sigma``
     - 0.1
     - Spread growth per metre of plume height [m/m]
   * - ``n_sigma_cutoff``
     - 4.0
     - Gaussian kernel truncation (multiples of :math:`\sigma`)
   * - ``N_min_density``
     - 1.0×10⁻³
     - Ignition threshold [embers/m²/s]
   * - ``spot_radius``
     - 5.0
     - New ignition zone radius [m]
   * - ``check_interval``
     - 5
     - Run cascade check every N timesteps
   * - ``random_seed``
     - 0
     - RNG seed (0 = system time)
   * - ``use_3d_wind``
     - 0
     - 1 to use height-averaged wind from plt file
   * - ``require_3d_wind``
     - 0
     - 1 to abort if no valid 3-D plt wind data provided
   * - ``plt_wind_file``
     - *(empty)*
     - Path to massconsistent_amr plt directory

**References**:
Albini (1983); Sardoy, Consalvi, Porterie & Fernandez-Pello (2007).
*Modeling transport and combustion of firebrands from burning trees.*
Combustion and Flame 150(3):151–169.


Spotting Suppression Inside Retardant Drop Zones
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Inside active retardant drop zones the spotting probability is scaled by the
same suppression factor as ROS:

.. math::

   f_{\text{ret}} = 1 - \varepsilon \cdot e^{-(t - t_{\text{drop}}) / \tau_{\text{decay}}}

This ensures that cells covered by retardant do not generate long-range spot
fires, consistent with FARSITE's suppression model.

**Parameters**: Existing ``retardant_file`` (no new parameters).


Model Enhancements
------------------

These modules augment the primary fire spread models.  They are activated
independently via ParmParse flags and are compatible with any of the fire
spread models unless otherwise noted.

Barrier Polygons / Firebreaks
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The barrier polygon feature (`barrier_files`) allows one or more CSV files of
barrier vertex coordinates to be read at startup.  At every time step, any
grid cell whose centre falls nearest to a barrier vertex and that is currently
burning (:math:`\phi < 0`) is extinguished by setting:

.. math::

   \phi(i,j,k) \leftarrow +\tfrac{1}{2}\min(\Delta x,\, \Delta y)

This models firebreaks, fuel breaks, roads, and other physical barriers
**without** requiring the AMReX Embedded Boundary (EB) framework.

Nearest-Cell Mapping
~~~~~~~~~~~~~~~~~~~~~

For each vertex :math:`(X_v, Y_v)` in a barrier CSV file, the nearest cell
index is found by:

.. math::

   i = i_0 + \left\lfloor \frac{X_v - x_{\text{lo}}}{\Delta x} \right\rfloor, \quad
   j = j_0 + \left\lfloor \frac{Y_v - y_{\text{lo}}}{\Delta y} \right\rfloor

with the result clamped to the domain bounds.  The set of barrier cells is
de-duplicated once at startup and stored on the GPU device for efficiency.

Usage::

    # Two CSV firebreak files (space-separated list)
    barrier_files = road_break.csv ridge_break.csv

CSV format::

    # X Y  (one vertex per line; lines starting with # are ignored)
    330100.0  3775300.0
    330200.0  3775300.0

The barrier is applied after every propagation step regardless of the
propagation method (``levelset``, ``farsite``, or ``mtt``).


Andrews (2018) Wind Adjustments for Rothermel
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Andrews (2018) documents two critical wind adjustments that improve the physical
accuracy of the Rothermel surface fire spread model when wind input is from NWP
or WRF (20-ft / 10-m height).  Both are optional flags that can be enabled
independently.

Wind Adjustment Factor (WAF)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Rothermel model requires wind speed at *midflame height*, but field
measurements and WRF/NWP output are typically at 20 ft (6.1 m) above open
terrain.  The WAF converts 20-ft open wind to midflame height.

Two formulas are available, selected via ``rothermel.waf_formula``:

**Andrews logarithmic formula** (``waf_formula = "andrews"``, default)

Derived from Albini & Baughman (1979).  Used for all fuel types; the sheltered
(closed-canopy) variant is applied automatically when canopy data are available.

*Unsheltered (open / shrub fuel beds)*:

.. math::

   \text{WAF} = \frac{1.83}{\ln\!\left(\dfrac{20 + 0.36\,h}{0.13\,h}\right)}

where :math:`h` is the fuel bed depth [ft].  Typical values: 0.42 for FM4
(:math:`h = 2` ft), 0.35 for FM9 (:math:`h = 3` ft).

*Sheltered (closed-canopy, FARSITE-style)*:

When canopy cover fraction :math:`f_c \geq 0.5` and canopy height
:math:`h_c > h`:

.. math::

   \text{WAF}_\text{sheltered} =
     \frac{0.555}{\sqrt{f_c \, h_c}
     \ln\!\left(\dfrac{20 + 0.36\,h_c}{0.13\,h_c}\right)}

**BehavePlus linear formula** (``waf_formula = "behaviorplus"``)

A simpler linear approximation used in BehavePlus for open and shrub fuel models:

.. math::

   \text{WAF} = 0.36 + 0.004\,h_\text{in}

where :math:`h_\text{in} = 12 \times h_\text{ft}` is the fuel bed depth in
inches.  Typical values for common fuels:

+----------+---------+---------+---------------+
| Fuel     | h [ft]  | h [in]  | WAF (linear)  |
+==========+=========+=========+===============+
| FM1      | 1.0     | 12      | 0.408         |
+----------+---------+---------+---------------+
| FM2/FM3  | 2.5     | 30      | 0.480         |
+----------+---------+---------+---------------+
| FM4      | 6.0     | 72      | 0.648         |
+----------+---------+---------+---------------+
| FM6      | 2.5     | 30      | 0.480         |
+----------+---------+---------+---------------+
| FM9      | 0.2     | 2.4     | 0.370         |
+----------+---------+---------+---------------+

For closed-canopy (sheltered) cells the BehavePlus path uses an exponential
Beer–Lambert-style canopy attenuation:

.. math::

   \text{WAF}_\text{canopy} =
     \text{WAF}_\text{open}(h_c) \times \exp(-\alpha_c \, f_c)

where :math:`\alpha_c` is the canopy attenuation coefficient
(``rothermel.waf_canopy_alpha``, default 1.5) and :math:`f_c` is the canopy
cover fraction.  At :math:`f_c = 1` and :math:`\alpha_c = 1.5` the canopy
reduces the midflame wind to about 22 % of the above-canopy value.

Enable via ``rothermel.use_waf = 1``.  When a landscape file is active and
provides canopy cover and canopy height data, WAF is computed per cell using
each fuel model's depth and the local canopy structure.

Maximum Effective Wind Speed (MEWS)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
~~~~~~~~~~

Andrews, P.L. (2018). *The Rothermel Surface Fire Spread Model and Associated
Developments: A Comprehensive Explanation.* Gen. Tech. Rep. RMRS-GTR-371.
USDA Forest Service.
https://doi.org/10.2737/RMRS-GTR-371

Albini, F.A. & Baughman, R.G. (1979). *Estimating Windspeeds for Predicting
Wildland Fire Behavior.* USDA For. Serv. Res. Pap. INT-221.

Massman, W.J. (1987). A comparative study of some mathematical models of the
mean wind structure and aerodynamic drag of plant canopies.
*Boundary-Layer Meteorology*, 40(1–2), 179–197.



Heat Flux MultiFab and Fire-Induced Wind
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A spatially-varying (or uniform) heat flux field :math:`Q` [W/m²] represents
the fire heat release rate at each cell.  It drives two WindNinja-style wind
corrections that are applied after any terrain-based modification.

Upward Convective Velocity
~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
~~~~~~~~~~~~~~~~~~~~~~~~~~

The fire-plume entrainment also creates a horizontal inflow toward the fire
perimeter:

.. math::

   U_{\text{ind}} = k_{\text{ind}} \, w_{\uparrow}

directed anti-parallel to the level-set gradient :math:`\nabla\phi`
(i.e., toward decreasing :math:`\phi`, which is the interior of the fire).
This term is only applied outside the fire perimeter (:math:`\phi \ge 0`).

Balbi Buoyancy Velocity Augmentation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
~~~~~~~~~~~~~~~~

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


FMC Seasonal / Phenological Schedule
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The foliar moisture content (FMC) used by the Van Wagner (1977) crown fire initiation
threshold is updated each timestep from a seasonal schedule.  Two modes:

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

.. list-table::
   :header-rows: 1

   * - Parameter
     - Default
     - Description
   * - ``fmc_schedule.enable``
     - 0
     - 1 to enable FMC seasonal schedule
   * - ``fmc_schedule.file``
     - ""
     - Path to two-column CSV (doy, fmc_pct)
   * - ``fmc_schedule.use_farsite_curve``
     - 0
     - 1 to use built-in parametric curve
   * - ``fmc_schedule.start_doy``
     - 1
     - Day-of-year at simulation t = 0
   * - ``fmc_schedule.spring_start``
     - 90
     - DOY when green-up begins
   * - ``fmc_schedule.summer_peak``
     - 150
     - DOY when FMC reaches maximum
   * - ``fmc_schedule.fall_start``
     - 240
     - DOY when curing begins
   * - ``fmc_schedule.fall_end``
     - 300
     - DOY when curing ends
   * - ``fmc_schedule.fmc_min``
     - 85.0
     - Dormant / cured FMC [%]
   * - ``fmc_schedule.fmc_max``
     - 140.0
     - Peak green FMC [%]


Precipitation-Driven Dead Fuel Moisture Wetting
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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

.. list-table::
   :header-rows: 1

   * - Parameter
     - Default
     - Description
   * - ``diurnal_moisture.precip_rain_rate_mm_hr``
     - 0.0
     - Constant rain rate [mm/hr]
   * - ``diurnal_moisture.precip_schedule_file``
     - ""
     - CSV of (time_s, rain_mm_hr) for time-varying rain
   * - ``diurnal_moisture.precip_threshold_mm_hr``
     - 0.25
     - Minimum rain rate to trigger wetting [mm/hr]
   * - ``diurnal_moisture.M_sat``
     - 1.20
     - Saturation moisture content [fraction]


Per-Cell Live Canopy Moisture from FMS File
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A FARSITE fuel moisture scenario (``.fms``) file can be supplied via ``fms_file`` to
provide per-fuel-model dead and live moisture values.  These are stored in a spatially-
varying MultiFab and override the global ``rothermel.M_d1``, …, ``rothermel.M_lw`` on
a per-cell basis.

**Requires**: A landscape file (``rothermel.landscape_file``) to supply per-cell fuel
model codes.

**Input parameter**: ``fms_file = path/to/fire.fms``


Per-Fuel-Model Burnout (Residence) Time
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When a landscape file is active, the residence time for each fuel code is derived
from the Rothermel (1983) formula:

.. math::

   \tau_i \ [\text{s}] = \frac{\alpha \rho_p}{\sigma_i}

where :math:`\alpha` = 3600 s·ft⁻¹/(lb/ft³), :math:`\rho_p` = 32 lb/ft³, and
:math:`\sigma_i` = 1739 ft⁻¹.  When no landscape file is active the global
``farsite.tau_residence`` scalar is used (backward compatible).

**Parameters**: Automatic when ``rothermel.landscape_file`` is present.


Live Fuel Moisture Conditioning
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

During the pre-simulation conditioning period (``conditioning.n_days > 0``),
live fuel moistures (M_lh, M_lw) are ramped linearly from their initial
values toward equilibrium dead-fuel targets over the conditioning steps,
matching FARSITE's behaviour:

.. code-block:: text

   M_lh_target = M_d100 × 1.5   (live herbaceous ≈ 150% of 100-hr dead)
   M_lw_target = M_d100 × 2.0   (live woody     ≈ 200% of 100-hr dead)

**Parameters**: Uses ``conditioning.n_days`` (existing parameter).


Burn-Period Controls (Diurnal Active-Spread Window)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When ``burn_period.enable = 1``, the computed rate-of-spread field is zeroed
outside the specified local clock window [``start_hour``, ``end_hour``),
preventing any level-set, FARSITE, or MTT advance during inactive hours.
All other processes (moisture evolution, spotting diagnostics, ecology metrics)
continue normally.

Midnight-crossing windows (e.g. ``start_hour = 22``, ``end_hour = 6``)
are handled correctly.

.. list-table::
   :header-rows: 1

   * - Parameter
     - Description
   * - ``burn_period.enable``
     - 1 to enable burn-period gating (default: 0)
   * - ``burn_period.start_hour``
     - Local hour (decimal) when fire becomes active (default: 10.0)
   * - ``burn_period.end_hour``
     - Local hour (decimal) when fire becomes inactive (default: 20.0)
   * - ``burn_period.sim_start_hour``
     - Local clock hour at simulation t=0 (default: 0.0)


FARSITE Topographic Horizon Shading
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``solar_radiation.H`` computes per-cell terrain shading from the local surface
normal (incidence-angle test) and canopy shading from cover fraction.  Neither
method detects shadows cast by terrain tens or hundreds of cells away.

FARSITE addresses this with a **full topographic horizon scan**: for each cell,
march 8 compass directions to find the maximum elevation angle to any visible
terrain point, then shade the cell whenever the solar elevation is below that
angle.  Enable via ``solar_radiation.use_topographic_horizon = 1``.

**Precomputation (one-time)**

For each locally-owned cell :math:`(i, j)` and each of the 8 compass directions
:math:`d`, the function ray-marches outward:

.. math::

   \theta_\text{hz}^{(d)}(i,j) =
   \max_{s=1,\ldots,S}
   \arctan\!\left(
     \frac{z(i + s\,\delta i_d,\; j + s\,\delta j_d) - z(i,j)}{s \cdot \Delta_d}
   \right)

where :math:`\Delta_d` is the physical step distance in direction :math:`d`
(cell spacing for axis-aligned directions, :math:`\sqrt{\Delta x^2 + \Delta y^2}`
for diagonals).  The result is stored in an 8-component MultiFab.

**Per-timestep shading check (GPU)**

The GPU kernel interpolates the horizon angle for the current solar azimuth
between the two bounding compass directions:

.. math::

   d_\text{lo} &= \lfloor \varphi_s / 45°\rfloor \bmod 8 \\
   d_\text{hi} &= (d_\text{lo} + 1) \bmod 8 \\
   t &= (\varphi_s / 45° - d_\text{lo}) \\
   \theta_\text{hz} &= \theta_\text{hz}^{(d_\text{lo})} \cdot (1 - t)
                     + \theta_\text{hz}^{(d_\text{hi})} \cdot t

A cell is marked as fully shaded (:math:`f = 1`) when
:math:`\theta_s < \theta_\text{hz}`, where :math:`\theta_s` is the solar
elevation angle.

The total shade fraction applied to fuel temperature is:

.. code-block:: none

   f = max(f_local_terrain,  f_canopy,  f_topographic_horizon,  cloud_cover)

.. list-table::
   :header-rows: 1
   :widths: 38 12 50

   * - Parameter
     - Default
     - Description
   * - ``solar_radiation.use_topographic_horizon``
     - ``0``
     - Set to ``1`` to enable the FARSITE 8-direction horizon scan.
   * - ``solar_radiation.horizon_scan_max_dist_m``
     - ``0.0``
     - Maximum ray-march distance [m]. ``0`` scans the full domain.

.. warning::

   The topographic horizon scan is **expensive** (one-time O(N³) CPU sweep plus
   MPI global elevation gather).  Use ``horizon_scan_max_dist_m`` to bound cost
   on large domains.  Disabled by default.


Diagnostic Models
-----------------

The following models can be run alongside any firespread model.  They do not alter fire-front
propagation but write additional diagnostic fields to every plotfile.

Byram Fire Behavior Diagnostics
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Fireline Intensity
~~~~~~~~~~~~~~~~~~

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
~~~~~~~~~~~~~

Byram's (1959) empirical relationship between fireline intensity and flame length:

.. math::

   L_f \;[\text{m}] = 0.0775 \times I_B^{0.46}

These fields (``fireline_intensity`` and ``flame_length``) are computed at every time
step and written to each plotfile.

Viegas (2004) Eruptive Fire Diagnostics
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The Viegas (2004) model is an **optional parallel diagnostic** that characterises
eruptive (blow-up) fire behaviour on steep terrain.  Unlike the primary spread
models it does **not** alter the fire front propagation; instead it writes five
diagnostic fields to every plotfile so that potential under-prediction by
Rothermel's quadratic slope factor can be identified.  Enable with
``viegas.enable = 1``.

Physical Basis
~~~~~~~~~~~~~~

On steep slopes the Rothermel quadratic slope factor :math:`\phi_s` significantly
under-predicts the rate of spread compared to field and laboratory observations of
eruptive fires (Viegas 2004).  Viegas introduces an *exponential* slope enhancement
factor that captures the runaway acceleration observed at critical slope angles.

Exponential Slope Enhancement Factor
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. math::

   \Phi_{s,V} = \exp\!\bigl(a_V \,\tan\varphi\bigr)

where:

* :math:`a_V` is the Viegas slope coefficient (default 1.83, calibrated to
  laboratory fuels; dimensionless)
* :math:`\varphi` is the terrain slope angle; :math:`\tan\varphi =
  \sqrt{(\partial z/\partial x)^2 + (\partial z/\partial y)^2}`

Viegas Rate of Spread
~~~~~~~~~~~~~~~~~~~~~

The Viegas ROS uses the Rothermel no-wind, no-slope rate :math:`R_0` and the
Rothermel wind factor :math:`\phi_w` as a baseline and replaces only the slope
factor:

.. math::

   R_V = R_0 \,(1 + \phi_w) \,\Phi_{s,V}
       = R_0 \,(1 + \phi_w) \,\exp\!\bigl(a_V \,\tan\varphi\bigr)

Eruptive Regime Flag
~~~~~~~~~~~~~~~~~~~~

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
~~~~~~~~~~~~~~~~

The signed excess ratio quantifies the relative difference between the two models:

.. math::

   \varepsilon = \frac{R_V - R_{\text{Rothermel}}}{R_{\text{Rothermel}}}

Positive values indicate that Viegas predicts a higher spread rate than
Rothermel; the larger :math:`\varepsilon`, the greater the potential
under-prediction by Rothermel on steep terrain.

Flame-Tilt Angle (Hazard Assessment)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
~~~~~~~~~~~~~~~~

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
~~~~~~~~~~~~~~~~~~~~~~~~~

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
~~~~~~~~~~~~~~~~

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
~~~~~~~~~

Viegas, D.X. (2004). "Slope and wind effects on fire propagation."
*International Journal of Wildland Fire*, 13(2), 143–156.
https://doi.org/10.1071/WF03046

Viegas-Balbi Coupled Diagnostic
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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


Weise & Biging (1996) Fire Whirl Model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The Weise & Biging (1996) fire whirl model is an **optional diagnostic sub-model**
that computes fire whirl characteristics from existing fire behavior outputs.  It
does **not** modify the primary fire spread model; instead it runs alongside any
spread model to predict whirl formation when enabled via ``weise_biging.enable = 1``.

Physical Basis
~~~~~~~~~~~~~~~

Fire whirls are columnar vortices generated by fire-induced updrafts interacting
with ambient wind and terrain.  The model characterises whirl geometry and
kinematics from Byram fireline intensity and flame length.

Flame Tilt Angle
~~~~~~~~~~~~~~~~~

The Weise & Biging (1996) empirical flame tilt angle :math:`\theta`:

.. math::

   \tan\theta = 0.382\,Fr^{0.432} + \tan\beta

where :math:`Fr = U^2 / (g\,L_f)` is the modified Froude number, :math:`U` is
wind speed [m/s], :math:`g = 9.81` m/s², :math:`L_f` is Byram flame length [m],
and :math:`\beta` is the terrain slope angle.

Vertical Flame Height
~~~~~~~~~~~~~~~~~~~~~~~

.. math::

   H_f = L_f \cos\theta

Fire Whirl Geometry (Rankine-Vortex Scaling)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. math::

   H_w &= H_f \\
   r_w &= \max(c_r\,H_w,\; 0.01\ \text{m}) \\
   \Gamma &= U\,H_w \\
   \Omega &= \frac{\Gamma}{2\pi\,r_w^2} \\
   v_\theta &= \Omega\,r_w

Configuration via Parmparse (``weise_biging.*``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
~~~~~~~~~~~~~~~~~~~~~~~~~~

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
~~~~~~~~~~

Weise, D.R. and Biging, G.S. (1996). "Effects of wind velocity and slope on
flame properties." *Canadian Journal of Forest Research*, 26(10), 1849–1858.
https://doi.org/10.1139/x26-210



Bulk Fuel Consumption
^^^^^^^^^^^^^^^^^^^^^^

The bulk fuel consumption fraction :math:`f_c` represents the fraction of available fuel consumed behind the fire front:

.. math::

   f_c = f_c(I_R, \tau_{res})

where :math:`\tau_{res}` is the residence time. The consumption affects the fire intensity and heat release rate.

The consumption fraction is bounded by configurable minimum and maximum values
(``farsite.f_consumed_min`` and ``farsite.f_consumed_max``):

.. math::

   f_c = f_{c,\min} + (f_{c,\max} - f_{c,\min}) \left(1 - \exp\left(-\frac{t}{\tau_{res}}\right)\right)


Post-Processing Diagnostics
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In addition to fireline intensity and flame length, the solver always writes four fire-ecology
fields and three emissions fields to every plotfile.  The equations below describe how each
field is computed.

Scorch Height (Van Wagner 1973)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Van Wagner's (1973) empirical relationship between fireline intensity and the height to which
foliage is killed by convective heat:

.. math::

   H_s \;[\text{m}] = 0.1483 \; I_B^{2/3}

where :math:`I_B` [kW/m] is the Byram fireline intensity.  Cells where :math:`I_B = 0` receive
:math:`H_s = 0`.  Output field: ``scorch_height``.

Probability of Ignition (Anderson 1970)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
  - PM\ :sub:`2.5`\: 0.0162 kg/kg (Andreae & Merlet 2001)

Output fields: ``co2_emissions``, ``co_emissions``, ``pm25_emissions``.


Always-Present Output Fields
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Every plotfile produced by the solver contains the following fields regardless of which
models are enabled.  Optional-model fields (Albini, Viegas, Weise & Biging, stochastic
spotting) are described in their respective model sections.

**Core simulation fields**

.. list-table::
   :header-rows: 1
   :widths: 22 12 66

   * - Field
     - Units
     - Description
   * - ``phi``
     - m
     - Level-set signed-distance function.  Negative inside the fire perimeter,
       positive outside.  The zero contour :math:`\phi = 0` is the active fire front.
   * - ``velx``
     - m/s
     - Wind velocity component in the *x*-direction at cell centres.
   * - ``vely``
     - m/s
     - Wind velocity component in the *y*-direction at cell centres.

**Terrain and fuel fields**

These fields carry the static landscape inputs and are constant throughout the
simulation.  They are written once per plotfile for easy co-location with fire
behaviour outputs.

.. list-table::
   :header-rows: 1
   :widths: 22 12 66

   * - Field
     - Units
     - Description
   * - ``elevation``
     - m
     - Ground elevation above mean sea level, read from the terrain file or landscape
       file (LCP).  Zero when neither is provided.
   * - ``slope``
     - --
     - Terrain slope magnitude :math:`\sqrt{(\partial z/\partial x)^2 + (\partial z/\partial y)^2}`.
       Dimensionless (rise/run).
   * - ``aspect``
     - rad
     - Terrain aspect (downslope azimuth from north, clockwise).  Used internally for
       anisotropic slope corrections.
   * - ``fuel_model``
     - --
     - Integer fuel model index for each cell.  Values correspond to FBFM13 (1–13)
       or FBFM40 (91–204) identifiers, or the uniform model code when no landscape
       file is provided.

**Crown fire field**

.. list-table::
   :header-rows: 1
   :widths: 22 12 66

   * - Field
     - Units
     - Description
   * - ``crown_fraction``
     - --
     - Fraction of the total rate of spread contributed by crown fire
       (:math:`R_{\text{crown}} / R_{\text{total}}` when crown fire is initiated).
       Zero for pure surface fire, one for pure crown fire.  Computed only when
       the crown fire initiation model is active (``crown.use_crown_initiation = 1``).

**Bulk fuel consumption field**

.. list-table::
   :header-rows: 1
   :widths: 22 12 66

   * - Field
     - Units
     - Description
   * - ``fuel_consumption``
     - --
     - Fraction of oven-dry fuel consumed behind the fire front
       (:math:`f_c \in [0,1]`, see `Bulk Fuel Consumption`_).
       Zero in unburned cells and in cells ahead of the front.

**Stochastic spotting fields**

These four fields are written when the probability-based spotting model is active
(``spotting.enable = 1``; see `Firebrand Spotting`_ for the full model description).

.. list-table::
   :header-rows: 1
   :widths: 22 12 66

   * - Field
     - Units
     - Description
   * - ``spot_prob``
     - --
     - Probability that a firebrand generated at this cell produces a secondary
       ignition.  Computed from wind speed, fuel moisture, and fire intensity
       (Bernoulli spot-probability model).
   * - ``spot_count``
     - --
     - Number of spot ignitions generated from this source cell at the most
       recent spotting step.
   * - ``spot_dist``
     - m
     - Maximum spotting distance from this source cell at the most recent
       spotting step (landing distance of the most distant active firebrand).
   * - ``spot_active``
     - --
     - Flag (1) at cells that received at least one spot ignition during the
       most recent spotting step; 0 elsewhere.

**Spatial Fuel Moisture Fields**

These five fields are written when the Nelson (2000) diurnal EMC or FARSITE
``.fms`` spatial scenario file is active; otherwise they contain the global
Rothermel moisture values uniformly.

.. list-table::
   :header-rows: 1
   :widths: 22 12 66

   * - Field
     - Units
     - Description
   * - ``moisture_d1``
     - --
     - 1-hr dead fuel moisture [fraction]
   * - ``moisture_d10``
     - --
     - 10-hr dead fuel moisture [fraction]
   * - ``moisture_d100``
     - --
     - 100-hr dead fuel moisture [fraction]
   * - ``moisture_lh``
     - --
     - Live herbaceous fuel moisture [fraction]
   * - ``moisture_lw``
     - --
     - Live woody fuel moisture [fraction]

References
~~~~~~~~~~~

* Van Wagner, C.E. (1973). "Height of crown scorch in forest fires." *Canadian Journal of
  Forest Research*, 3(3), 373–378.
* Anderson, H.E. (1970). "Forest fuel ignitibility." *Fire Technology*, 6(4), 312–319.
* Ryan, K.C. & Reinhardt, E.D. (1988). "Predicting postfire mortality of seven western
  conifers." *Canadian Journal of Forest Research*, 18(10), 1291–1297.
* Seiler, W. & Crutzen, P.J. (1980). "Estimates of gross and net fluxes of carbon between the
  biosphere and the atmosphere from biomass burning." *Climatic Change*, 2(3), 207–247.
* Andreae, M.O. & Merlet, P. (2001). "Emission of trace gases and aerosols from biomass
  burning." *Global Biogeochemical Cycles*, 15(4), 955–966.



Wind Models
-----------

Several wind modelling approaches are available, from a simple constant
uniform wind to spatially-varying fields and terrain-following feedback models.
The selection is made through a combination of input parameters as described
below.

Uniform Wind
^^^^^^^^^^^^

The simplest wind specification sets constant velocity components over the
entire domain for the whole simulation:

.. code-block:: ini

   u_x = 2.0   # x-wind component [m/s]
   u_y = 0.5   # y-wind component [m/s]

The default is ``u_x = 0.25``, ``u_y = 0.0``.  This is the
fallback when no wind file or terrain-following solver is used.

Wind File — Steady (Spatially-Varying)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A spatially-varying but time-constant wind field can be supplied as a CSV file:

.. code-block:: ini

   velocity_file = wind_data.csv

The file must contain one sample point per line with columns
**X  Y  U  V** (in metres and m/s respectively).  The solver interpolates
the wind onto the computational grid using inverse-distance weighting (IDW)
at initialisation.  When ``velocity_file`` is set, it overrides
``u_x`` / ``u_y`` / ``u_z``.

CSV format::

    # X [m]  Y [m]  U [m/s]  V [m/s]
    0       0       2.0      0.5
    1000    0       2.1      0.4
    0       1000    1.9      0.6
    1000    1000    2.0      0.5

Wind File — Unsteady (Time-Dependent)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When ``use_time_dependent_wind = 1``, the solver loads a sequence of wind
snapshot files and performs temporal linear interpolation between them.  This
is only available in 2-D builds.

.. code-block:: ini

   velocity_file          = wind_t0.csv      # base snapshot (t = 0)
   use_time_dependent_wind = 1
   wind_time_spacing      = 3600.0           # seconds between snapshots

Snapshot files follow the naming convention:

* ``<velocity_file>`` — time index 0
* ``<base>_1.csv``    — time index 1
* ``<base>_2.csv``    — time index 2

At each simulation step the solver:

1. Determines which two snapshots bracket the current time.
2. Loads both files (IDW spatial interpolation).
3. Performs temporal linear interpolation.

For example, with ``wind_time_spacing = 3600.0`` s:

* At t = 0 s  : uses snapshot 0
* At t = 1800 s : 50 % blend of snapshots 0 and 1
* At t = 3600 s : uses snapshot 1


Turbulent Wind Perturbation (``turb_wind.*``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Adds stochastic velocity perturbations to the ambient wind field at every
timestep.  The unperturbed base wind is preserved in a separate MultiFab;
the perturbed field is always computed as *base + perturbation* so the
perturbation never drifts away from the intended background.

Two models are available.

Ornstein-Uhlenbeck Process (``turb_wind.model = ou_process``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The OU process produces temporally correlated wind gusts whose
auto-correlation decays exponentially in time (Uhlenbeck & Ornstein 1930).

**Discrete-time OU update (exact Euler–Maruyama)**

.. math::

   \alpha &= \exp(-\theta\,\Delta t) \\
   \sigma_{\text{step}} &= \sigma\,\sqrt{1 - \alpha^2} \\
   u'(t + \Delta t) &= \alpha\,u'(t) + \sigma_{\text{step}}\,\eta_u(t) \\
   v'(t + \Delta t) &= \alpha\,v'(t) + \sigma_{\text{step}}\,\eta_v(t)

where:

* :math:`\theta` [s⁻¹] is the reversion rate (``turb_wind.theta``);
  decorrelation time = :math:`1/\theta`.
* :math:`\sigma` [m/s] is the stationary standard deviation of each cell's
  perturbation (``turb_wind.sigma``).
* :math:`\eta_u, \eta_v` are unit-variance noise terms (white when
  ``L_c = 0``, spatially correlated when ``L_c > 0``; see below).

The effective wind used for all ROS computations is

.. math::

   \mathbf{u}_{\text{eff}}(x,y) = \mathbf{u}_{\text{base}}(x,y)
                                 + \bigl(u'(x,y),\; v'(x,y)\bigr).

**Domain-uniform mode (L_c = 0)**

A single pair :math:`(u', v')` is drawn and added to every cell.  This
corresponds to a gust event coherent over the entire domain — appropriate
when the domain is small relative to the turbulent integral scale.

**Spatially correlated mode (L_c > 0)**

Per-cell OU states :math:`u'(x,y)` and :math:`v'(x,y)` are maintained in a
MultiFab.  The noise field :math:`\eta(x,y)` is constructed as:

1. Draw independent :math:`\mathcal{N}(0,1)` white noise :math:`\xi(x,y)` at
   every cell.
2. Apply a separable 1-D Gaussian convolution in :math:`x` then :math:`y`
   with kernel standard deviation

   .. math::

      \sigma_k = \frac{L_c}{\Delta x} \quad\text{[cells]}

   using Neumann (nearest-cell) clamping at domain and subdomain boundaries.
3. Rescale the smoothed field to unit RMS so the per-cell stationary standard
   deviation remains exactly :math:`\sigma`.

The resulting noise field has the approximate 2-D autocovariance

.. math::

   C_\eta(r_x, r_y) \approx
   \exp\!\left(-\frac{r_x^2 + r_y^2}{2 L_c^2}\right)

i.e.\ a Gaussian envelope with e-folding length :math:`L_c` (``turb_wind.L_c``).

Random Fourier Feature Spectral Noise (``turb_wind.model = spectral_noise``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The spectral noise model combines a physically correct spatial power spectrum
with OU temporal evolution, following the Random Fourier Feature (RFF)
methodology (Kraichnan 1970; Rahimi & Recht 2007).

**Initialisation – wavenumber and phase sampling**

At startup :math:`N` wavenumber pairs and phases are drawn once and remain
fixed for the entire simulation:

.. math::

   k_{x,n},\; k_{y,n} &\;\sim\; \mathcal{N}\!\left(0,\; \frac{1}{L_c^2}\right)
   \quad\text{independently,} \\
   \phi^u_n,\; \phi^v_n &\;\sim\; \mathcal{U}[0,\; 2\pi),

where :math:`L_c` is the user-specified spatial correlation length
(``turb_wind.L_c`` [m]) and :math:`N` is the number of modes
(``turb_wind.N_modes``).  Sampling :math:`k_x, k_y` independently from
:math:`\mathcal{N}(0, 1/L_c^2)` corresponds to a 2-D **isotropic Gaussian**
power spectral density :math:`S(k) \propto \exp(-k^2 L_c^2/2)`, whose
Fourier transform is the Gaussian autocorrelation function.

**Temporal update – OU amplitude evolution (CPU, every timestep)**

Each mode carries a pair of scalar OU amplitude coefficients
:math:`(A^u_n, A^v_n)` with unit stationary variance:

.. math::

   \alpha &= \exp(-\theta\,\Delta t), \\
   A^u_n(t+\Delta t) &= \alpha\,A^u_n(t) + \sqrt{1-\alpha^2}\;\xi^u_n(t), \\
   A^v_n(t+\Delta t) &= \alpha\,A^v_n(t) + \sqrt{1-\alpha^2}\;\xi^v_n(t),

where :math:`\xi^u_n, \xi^v_n \sim \mathcal{N}(0,1)` are independent each step.

**Field reconstruction (GPU, every timestep)**

The perturbation field is reconstructed on the GPU as a weighted cosine
superposition:

.. math::

   u'(x,y) &= \sigma\sqrt{\frac{2}{N}}\sum_{n=1}^{N}
              A^u_n\cos\!\bigl(k_{x,n}\,x + k_{y,n}\,y + \phi^u_n\bigr), \\
   v'(x,y) &= \sigma\sqrt{\frac{2}{N}}\sum_{n=1}^{N}
              A^v_n\cos\!\bigl(k_{x,n}\,x + k_{y,n}\,y + \phi^v_n\bigr).

The scale factor :math:`\sigma\sqrt{2/N}` ensures that the stationary
variance of each cell's perturbation equals :math:`\sigma^2` [m²/s²]:

.. math::

   \mathrm{Var}[u'(x,y)]
   = \sigma^2 \cdot \frac{2}{N} \cdot N \cdot 1 \cdot \mathbb{E}[\cos^2] = \sigma^2.

**Properties**

* Temporal decorrelation time = :math:`1/\theta` [s] (same as ``ou_process``).
* Approximate 2-D spatial autocovariance (by the Bochner theorem):

  .. math::

     C(r) \approx \sigma^2 \exp\!\left(-\frac{r^2}{2 L_c^2}\right).

* Energy is distributed across all resolved wavenumbers according to the
  Gaussian PSD — unlike the grid-smoothing approach, no artificial
  truncation of the spectrum occurs.
* GPU-accelerated: the inner-mode loop in the reconstruction kernel runs on
  the device; CPU workload per step is only :math:`O(N)` for the OU updates.

**References**

* Kraichnan, R.H. (1970). Diffusion by a random velocity field.
  *Phys. Fluids* 13(1):22–31.
* Rahimi, A. & Recht, B. (2007). Random features for large-scale kernel
  machines. *NIPS 2007*.

Bounded Direction Random Walk (``turb_wind.model = direction_walk``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

At each timestep a Gaussian angular increment is drawn and accumulated:

.. math::

   \Delta\vartheta &\sim \mathcal{N}(0,\,\sigma_\vartheta^2) \\
   \vartheta(t+\Delta t) &= \operatorname{clamp}
       \bigl(\vartheta(t) + \Delta\vartheta,\;
             {-\vartheta_{\max}},\; {+\vartheta_{\max}}\bigr)

The entire base wind field is then rotated by the cumulative angle
:math:`\vartheta`:

.. math::

   u_{\text{eff}} &= u_{\text{base}} \cos\vartheta - v_{\text{base}} \sin\vartheta \\
   v_{\text{eff}} &= u_{\text{base}} \sin\vartheta + v_{\text{base}} \cos\vartheta

Wind *speed* is preserved exactly; only direction fluctuates.  Parameters:
``turb_wind.sigma_theta`` [rad/step] and ``turb_wind.theta_max`` [rad].

**References**

* Uhlenbeck, G.E. & Ornstein, L.S. (1930). On the theory of the Brownian
  motion. *Phys. Rev.* 36(5):823–841.
* Finney, M.A. et al. (2011). Role of wind, fuel moisture, and terrain in
  controlling fire movement. *Ecosphere* 2(3):art17.
* Cruz, M.G. et al. (2019). Uncertainty in wildfire behaviour research.
  *Current Forestry Reports* 5:155–172.


Wind-Terrain Feedback Models (Options 1–6)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``wind_terrain.model`` ParmParse key selects how terrain-induced or
fire-feedback winds modify the effective wind seen by the fire spread model.
**Option 1 (``none``) is the default and preserves all existing behaviour.**
Options 2–6 progressively couple Viegas or terrain-channel wind physics into
the actual fire propagation.  Option 7 (``windninja_ridge_canyon``) is described
in the following section.

Options Summary
~~~~~~~~~~~~~~~~

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
~~~~~~~~~~~~~~~~~~~

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
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
~~~~~~~~~~~

Rothermel, R.C. (1983). *How to Predict the Spread and Intensity of Forest and
Range Fires.* USDA Forest Service General Technical Report INT-143.

Viegas, D.X. & Neto, L.P.S. (1994). "Wind tunnel study of the convective air
flow of slope fires." Annual Report, Project COMOESTAS.

Pimont, F., Dupuy, J.-L., Linn, R.R. & Dupont, S. (2009). "Validation of
FIRETEC wind-flows over a canopy and fuel-break." *International Journal of
Wildland Fire*, 18(7), 775–790.
https://doi.org/10.1071/WF07130


WindNinja Ridge/Canyon Terrain Speed-up (Option 7)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The WindNinja empirical ridge/canyon model (Option 7,
``wind_terrain.model = windninja_ridge_canyon``) accounts for terrain-driven
wind acceleration observed in ridge and canyon topography (Forthofer 2007).

Wind-Slope Alignment
~~~~~~~~~~~~~~~~~~~~~

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
~~~~~~~~~~~~~~~

When the wind climbs a slope, terrain convergence amplifies the wind speed:

.. math::

   f_{\text{ridge}} = 1 + k_{\text{ridge}} \,\tan\varphi \, a \quad (a > 0)

Canyon Channeling
~~~~~~~~~~~~~~~~~

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

Compact Wind Direction Schedule
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A compact three-column CSV schedule (``wind_dir_schedule_file``) updates the
uniform wind field at every timestep without requiring a full spatial wind grid
per snapshot:

.. code-block:: ini

   wind_dir_schedule_file = wind_schedule.csv

CSV format: ``time_s, speed_ms, dir_deg`` where ``dir_deg`` is the direction
*from* which the wind blows (270° = westerly wind blowing eastward).  The
schedule is linearly interpolated with circular direction averaging to avoid
0°/360° wrap-around artefacts.

When this file is set it overrides ``u_x`` / ``u_y`` and the
``use_time_dependent_wind`` file-grid path.  Compatible with turbulent wind
perturbation (perturbations are added on top of the scheduled base wind).

**Input parameter**: ``wind_dir_schedule_file = path/to/wind_schedule.csv``


Multiple Weather Stations with Spatial IDW Interpolation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A new ``multi_wtr_file`` input enables loading multiple FARSITE ``.wtr`` weather
station files and producing spatially-varying wind, temperature, and relative
humidity via inverse-distance-weighting (IDW).

Station list CSV format (``multi_wtr_file``)::

   # station_id, x_m, y_m, wtr_file
   1, 330000.0, 3775000.0, station1.wtr
   2, 335000.0, 3775000.0, station2.wtr
   3, 332500.0, 3780000.0, station3.wtr

At each timestep the U and V wind components from each station are IDW-
interpolated to every grid cell:

.. math::

   V_{\text{cell}} = \frac{\sum_i w_i V_i}{\sum_i w_i},
   \quad w_i = d_i^{-p}

where :math:`d_i` is the distance from the cell to station :math:`i` and
:math:`p` is the IDW power exponent (default 2.0).

.. list-table::
   :header-rows: 1

   * - Parameter
     - Description
   * - ``multi_wtr_file``
     - Path to station list CSV (default: ``""`` = disabled)
   * - ``multi_wtr_idw_power``
     - IDW exponent :math:`p` (default: 2.0)
