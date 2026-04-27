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

The moisture damping coefficient :math:`\eta_M` is:

.. math::

   \eta_M = \max\left(1 - 2.59r_m + 5.11r_m^2 - 3.52r_m^3, 0\right)

where :math:`r_m = \min(M_f / M_x, 1)`, :math:`M_f` is the fuel moisture content, and :math:`M_x` is the moisture of extinction.

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

Spotting Distance
^^^^^^^^^^^^^^^^^

The maximum spotting distance :math:`d_{max}` depends on:

* Firebrand lofting height
* Wind speed
* Firebrand size and terminal velocity

The spotting probability decreases with distance from the source and is modulated by terrain and fuel moisture.

Bulk Fuel Consumption
----------------------

The bulk fuel consumption fraction :math:`f_c` represents the fraction of available fuel consumed behind the fire front:

.. math::

   f_c = f_c(I_R, \tau_{res})

where :math:`\tau_{res}` is the residence time. The consumption affects the fire intensity and heat release rate.

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
