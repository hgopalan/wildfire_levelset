API Reference
=============

This section provides a reference for the main functions and data structures in the wildfire level-set solver.

Core Functions
--------------

Rothermel Model
^^^^^^^^^^^^^^^

.. code-block:: cpp

   RothermelComputed compute_rothermel_params(const RothermelParams& rp)

Computes all Rothermel (1972) fire spread parameters from fuel properties.

**Parameters:**

* ``rp``: Structure containing fuel properties (see RothermelParams)

**Returns:**

* ``RothermelComputed``: Structure with computed parameters including:
  
  - ``R0``: No-wind, no-slope rate of spread [ft/min]
  - ``I_R``: Reaction intensity [BTU/ft²/min]
  - ``C``, ``B``, ``E``: Wind factor coefficients
  - ``phi_s``: Slope factor
  - ``beta``: Packing ratio

**Equations:**

The function implements the following key equations:

Packing ratio:

.. math::

   \beta = \frac{\rho_b}{\rho_p} = \frac{w_0}{\delta \rho_p}

Optimum reaction velocity:

.. math::

   \Gamma' = \Gamma_{max} \left(\frac{\beta}{\beta_{op}}\right)^A \exp\left[A\left(1 - \frac{\beta}{\beta_{op}}\right)\right]

Reaction intensity:

.. math::

   I_R = \Gamma' w_n h \eta_M \eta_s

No-wind, no-slope ROS:

.. math::

   R_0 = \frac{I_R \xi}{\rho_b \epsilon_h Q_{ig}}

FARSITE Ellipse
^^^^^^^^^^^^^^^

.. code-block:: cpp

   void compute_farsite_spread(
       MultiFab &phi,
       const MultiFab &vel,
       MultiFab &farsite_spread,
       const Geometry &geom,
       Real dt,
       const RothermelParams &rp,
       const InputParameters::FARSITEParams &fp,
       const InputParameters::CrownInitiationParams &cp,
       const MultiFab *slopes = nullptr,
       MultiFab *fuel_consumption = nullptr,
       MultiFab *crown_fire_fraction = nullptr
   )

Computes FARSITE elliptical fire spread using Richards' (1990) model.

**Parameters:**

* ``phi``: Level-set function (modified in-place)
* ``vel``: Velocity field (wind)
* ``farsite_spread``: Output array for spread rates
* ``geom``: AMReX geometry object
* ``dt``: Time step
* ``rp``: Rothermel fuel parameters
* ``fp``: FARSITE model parameters
* ``cp``: Crown fire parameters
* ``slopes``: Optional terrain slope data
* ``fuel_consumption``: Optional output for fuel consumption fraction
* ``crown_fire_fraction``: Optional output for crown fire fraction

**Algorithm:**

1. Find fire front locations where :math:`\phi \approx 0`
2. Compute base ROS using Rothermel model
3. Apply Van Wagner crown fire modification (if enabled)
4. Compute elliptical spread using Richards' coefficients :math:`a, b, c`
5. Update level-set function with new fire front positions

Level-Set Advection
^^^^^^^^^^^^^^^^^^^

.. code-block:: cpp

   void advect_levelset(
       MultiFab &phi,
       const MultiFab &velocity,
       const Geometry &geom,
       Real dt
   )

Advances the level-set function using the advection equation.

**Parameters:**

* ``phi``: Level-set function (modified in-place)
* ``velocity``: Velocity field (fire spread rate)
* ``geom``: AMReX geometry object
* ``dt``: Time step

**Equation:**

.. math::

   \frac{\partial \phi}{\partial t} + V|\nabla\phi| = 0

where :math:`V` is the local fire spread rate.

Terrain Slope
^^^^^^^^^^^^^

.. code-block:: cpp

   void compute_slopes_from_terrain(
       MultiFab& slopes,
       const Geometry& geom,
       const std::string& filename
   )

Reads terrain elevation data and computes slopes for each grid cell.

**Parameters:**

* ``slopes``: Output MultiFab with 2 components (slope_x, slope_y)
* ``geom``: AMReX geometry object
* ``filename``: Path to terrain data file (X Y Z format)

**Algorithm:**

1. Read X Y Z terrain data from file
2. Interpolate elevation to grid cell centers using inverse distance weighting (IDW)
3. Compute slopes using central differences:

.. math::

   s_x = \frac{\partial z}{\partial x} = \frac{z(x+\Delta x) - z(x-\Delta x)}{2\Delta x}

   s_y = \frac{\partial z}{\partial y} = \frac{z(y+\Delta y) - z(y-\Delta y)}{2\Delta y}

Crown Fire Initiation
^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: cpp

   void apply_crown_fire(
       MultiFab &spread_rate,
       const MultiFab &surface_intensity,
       const InputParameters::CrownInitiationParams &cp,
       MultiFab *crown_fraction = nullptr
   )

Applies Van Wagner (1977) crown fire initiation criteria.

**Parameters:**

* ``spread_rate``: Fire spread rate (modified if crown fire occurs)
* ``surface_intensity``: Surface fire intensity [kW/m]
* ``cp``: Crown fire parameters (CBH, CBD, moisture, etc.)
* ``crown_fraction``: Optional output for crown fire active fraction

**Criterion:**

Crown fire initiates when surface intensity exceeds critical value:

.. math::

   I > I_0 = \frac{CBH \times (460 + 25.9 M_c)}{18 h}

Active crown fire occurs when wind speed exceeds:

.. math::

   U > U_{active} = \frac{3}{\sqrt{CBD}}

Firebrand Spotting
^^^^^^^^^^^^^^^^^^

.. code-block:: cpp

   void apply_spotting(
       MultiFab &phi,
       const MultiFab &spread_rate,
       const MultiFab &wind,
       const Geometry &geom,
       const InputParameters::SpottingParams &sp,
       Real time
   )

Generates new ignition points from firebrand spotting.

**Parameters:**

* ``phi``: Level-set function (modified to add spot fires)
* ``spread_rate``: Current fire spread rate
* ``wind``: Wind velocity field
* ``geom``: AMReX geometry object
* ``sp``: Spotting parameters (probability, max distance, etc.)
* ``time``: Current simulation time

**Model:**

Spotting probability decreases with distance :math:`d` from fire front:

.. math::

   P_{spot}(d) = P_0 \exp\left(-\frac{d}{d_{max}}\right)

where :math:`P_0` is the base probability and :math:`d_{max}` is the maximum spotting distance.

Bulk Fuel Consumption
^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: cpp

   void compute_fuel_consumption(
       MultiFab &consumption_fraction,
       const MultiFab &phi,
       const MultiFab &intensity,
       const InputParameters::FuelConsumptionParams &fcp,
       Real dt
   )

Computes the fraction of fuel consumed behind the fire front.

**Parameters:**

* ``consumption_fraction``: Output consumption fraction [0, 1]
* ``phi``: Level-set function
* ``intensity``: Fire intensity [kW/m²]
* ``fcp``: Fuel consumption parameters
* ``dt``: Time step

**Model:**

Consumption fraction depends on residence time :math:`\tau_{res}`:

.. math::

   f_c = 1 - \exp\left(-\frac{t}{\tau_{res}}\right)

where :math:`t` is the time since fire front passage.

Data Structures
---------------

RothermelParams
^^^^^^^^^^^^^^^

Fuel properties for Rothermel fire spread model.

.. code-block:: cpp

   struct RothermelParams {
       Real w0;        // Oven-dry fuel loading [lb/ft²]
       Real sigma;     // Surface-area-to-volume ratio [ft⁻¹]
       Real delta;     // Fuel bed depth [ft]
       Real M_f;       // Fuel moisture content [fraction]
       Real M_x;       // Moisture of extinction [fraction]
       Real h_heat;    // Heat content [BTU/lb]
       Real S_T;       // Total mineral content [fraction]
       Real S_e;       // Effective mineral content [fraction]
       Real rho_p;     // Oven-dry particle density [lb/ft³]
       Real slope_x;   // Terrain slope in x-direction [tan(θ)]
       Real slope_y;   // Terrain slope in y-direction [tan(θ)]
       Real wind_conv; // Wind unit conversion factor
       Real ros_conv;  // ROS unit conversion factor
   };

**Members:**

* ``w0``: Total oven-dry fuel loading per unit area
* ``sigma``: Fuel surface-area-to-volume ratio (determines reaction velocity)
* ``delta``: Depth of the fuel bed
* ``M_f``: Fuel moisture content as a fraction of dry weight
* ``M_x``: Moisture content at which fire will not spread
* ``h_heat``: Heat content (energy per unit mass)
* ``S_T``: Total mineral content (inert material)
* ``S_e``: Effective mineral content (affects reaction intensity)
* ``rho_p``: Density of individual fuel particles
* ``slope_x``, ``slope_y``: Terrain slope components (rise/run)
* ``wind_conv``, ``ros_conv``: Unit conversion factors

RothermelComputed
^^^^^^^^^^^^^^^^^

Computed parameters from Rothermel model.

.. code-block:: cpp

   struct RothermelComputed {
       Real R0;              // No-wind, no-slope rate of spread [ft/min]
       Real I_R;             // Reaction intensity [BTU/ft²/min]
       Real beta_ratio_E;    // Packing-ratio part of wind factor
       Real C;               // Wind factor coefficient C
       Real B;               // Wind factor coefficient B
       Real phi_s;           // Slope factor
       Real beta;            // Packing ratio
       Real wind_conv;       // Wind unit conversion
       Real ros_conv;        // ROS unit conversion
   };

**Members:**

* ``R0``: Base rate of spread without wind or slope effects
* ``I_R``: Reaction intensity (heat release rate per unit area)
* ``beta_ratio_E``: Pre-computed term :math:`(\beta/\beta_{op})^{-E}` for wind factor
* ``C``, ``B``, ``E``: Coefficients for wind factor :math:`\phi_w = C(\beta/\beta_{op})^{-E} U^B`
* ``phi_s``: Slope factor :math:`\phi_s = 5.275\beta^{-0.3}\tan^2(\theta)`
* ``beta``: Packing ratio :math:`\beta = \rho_b/\rho_p`

FARSITEParams
^^^^^^^^^^^^^

Parameters for FARSITE elliptical expansion model.

.. code-block:: cpp

   struct FARSITEParams {
       int enable;              // Enable FARSITE model (1=yes, 0=no)
       int use_anderson_LW;     // Use Anderson L/W ratio (1=yes, 0=no)
       Real coeff_a;            // Richards' head fire coefficient
       Real coeff_b;            // Richards' flank fire coefficient
       Real coeff_c;            // Richards' backing fire coefficient
       Real phi_threshold;      // Fire front detection threshold
   };

**Members:**

* ``enable``: Toggle FARSITE model on/off
* ``use_anderson_LW``: Use Anderson (1983) length-to-width ratio to compute coefficients
* ``coeff_a``: Head fire coefficient (maximum spread, downwind)
* ``coeff_b``: Flank fire coefficient (perpendicular to wind)
* ``coeff_c``: Backing fire coefficient (minimum spread, upwind)
* ``phi_threshold``: Level-set value threshold for identifying fire front

CrownInitiationParams
^^^^^^^^^^^^^^^^^^^^^

Parameters for Van Wagner crown fire model.

.. code-block:: cpp

   struct CrownInitiationParams {
       int enable;                      // Enable crown fire (1=yes, 0=no)
       Real canopy_base_height;         // Canopy base height [m]
       Real canopy_bulk_density;        // Canopy bulk density [kg/m³]
       Real canopy_moisture;            // Canopy moisture content [fraction]
       Real crown_consumption_depth;    // Crown consumption depth [m]
   };

**Members:**

* ``enable``: Toggle crown fire model on/off
* ``canopy_base_height``: Height from ground to bottom of canopy (CBH)
* ``canopy_bulk_density``: Mass of crown fuel per unit canopy volume (CBD)
* ``canopy_moisture``: Moisture content of crown fuel
* ``crown_consumption_depth``: Depth of canopy consumed during crown fire

SpottingParams
^^^^^^^^^^^^^^

Parameters for firebrand spotting model.

.. code-block:: cpp

   struct SpottingParams {
       int enable;                  // Enable spotting (1=yes, 0=no)
       Real probability;            // Base spotting probability
       Real max_distance;           // Maximum spotting distance [m]
       Real wind_speed_threshold;   // Minimum wind speed for spotting [m/s]
   };

**Members:**

* ``enable``: Toggle spotting on/off
* ``probability``: Base probability of firebrand generation per active cell
* ``max_distance``: Maximum distance firebrands can travel
* ``wind_speed_threshold``: Minimum wind speed required for spotting to occur

FuelConsumptionParams
^^^^^^^^^^^^^^^^^^^^^

Parameters for bulk fuel consumption calculation.

.. code-block:: cpp

   struct FuelConsumptionParams {
       int enable;              // Enable fuel consumption (1=yes, 0=no)
       Real residence_time;     // Residence time [seconds]
   };

**Members:**

* ``enable``: Toggle fuel consumption calculation on/off
* ``residence_time``: Time constant for exponential fuel consumption :math:`\tau_{res}`

AMReX Integration
-----------------

MultiFab
^^^^^^^^

The code uses AMReX ``MultiFab`` objects to store multi-dimensional array data:

.. code-block:: cpp

   MultiFab phi(ba, dm, ncomp, nghost);

* ``ba``: BoxArray defining the domain decomposition
* ``dm``: DistributionMapping for parallel execution
* ``ncomp``: Number of components (variables)
* ``nghost``: Number of ghost cells for boundary conditions

Geometry
^^^^^^^^

The ``Geometry`` object defines the physical domain:

.. code-block:: cpp

   Geometry geom(domain, &rb, coord_sys, is_periodic);

* ``domain``: Box defining the index space
* ``rb``: RealBox defining physical coordinates
* ``coord_sys``: Coordinate system (0=Cartesian, 1=cylindrical, 2=spherical)
* ``is_periodic``: Array indicating periodic boundaries

ParallelFor
^^^^^^^^^^^

GPU/CPU parallel loops use ``ParallelFor``:

.. code-block:: cpp

   ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
       // Kernel code executed for each cell (i,j,k)
   });

The ``AMREX_GPU_DEVICE`` decorator enables the same code to run on CPU or GPU.

Units and Conversions
---------------------

The code uses multiple unit systems:

**Rothermel Model (US Customary Units):**

* Distance: feet (ft)
* Fuel loading: pounds per square foot (lb/ft²)
* Fuel depth: feet (ft)
* SAV ratio: inverse feet (ft⁻¹)
* Wind speed: feet per minute (ft/min) or miles per hour (mph)
* Rate of spread: feet per minute (ft/min)
* Heat content: BTU per pound (BTU/lb)
* Intensity: BTU per square foot per minute (BTU/ft²/min)

**SI Units (Used in some sub-models):**

* Distance: meters (m)
* Fuel loading: kilograms per square meter (kg/m²)
* Wind speed: meters per second (m/s)
* Rate of spread: meters per second (m/s)
* Intensity: kilowatts per meter (kW/m)

**Conversion Factors:**

* 1 ft = 0.3048 m
* 1 lb/ft² = 4.8824 kg/m²
* 1 mph = 0.44704 m/s = 26.8224 ft/min
* 1 ft/min = 0.00508 m/s
* 1 BTU/ft²/min = 0.18946 kW/m²

The code handles unit conversions internally through ``wind_conv`` and ``ros_conv`` parameters.
