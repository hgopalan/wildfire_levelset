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
* ``farsite_spread``: Output array for spread displacements
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

   void advect_levelset_weno5z_rk3(
       MultiFab &phi,
       const MultiFab &vel,
       const Geometry &geom,
       Real dt,
       const RothermelParams &rp,
       const MultiFab *terrain_slopes = nullptr
   )

Advances the level-set function using WENO5-Z spatial discretization and
third-order Runge-Kutta (RK3) time integration.

**Parameters:**

* ``phi``: Level-set function (modified in-place)
* ``vel``: Velocity field (wind)
* ``geom``: AMReX geometry object
* ``dt``: Time step
* ``rp``: Rothermel fuel parameters
* ``terrain_slopes``: Optional terrain slope data

**Equation:**

.. math::

   \frac{\partial \phi}{\partial t} + V|\nabla\phi| = 0

where :math:`V` is the local fire spread rate computed from the Rothermel model.

Fire Behavior Diagnostics
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: cpp

   void compute_fire_behavior(
       MultiFab& fireline_intensity_mf,
       MultiFab& flame_length_mf,
       const MultiFab& R_mf,
       const RothermelParams& rp
   )

Computes Byram (1959) fireline intensity and flame length from the Rothermel
rate of spread field.

**Parameters:**

* ``fireline_intensity_mf``: Output fireline intensity field [kW/m]
* ``flame_length_mf``: Output flame length field [m]
* ``R_mf``: Input rate of spread field [m/s]
* ``rp``: Rothermel fuel parameters

**Equations:**

Fireline intensity (Byram 1959):

.. math::

   I_B \;[\text{kW/m}] = H \;[\text{kJ/kg}] \times w_a \;[\text{kg/m}^2] \times R \;[\text{m/s}]

where :math:`H` is the low heat of combustion, :math:`w_a` is the available
(net) fuel load, and :math:`R` is the rate of spread.

Flame length (Byram 1959, SI form):

.. math::

   L_f \;[\text{m}] = 0.0775 \times I_B^{0.46}

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
* ``cp``: Crown fire parameters (CBH, CBD, FMC, etc.)
* ``crown_fraction``: Optional output for crown fire active fraction

**Criterion:**

Crown fire initiates when surface intensity exceeds critical value:

.. math::

   I > I_0 = \frac{CBH \times (460 + 25.9 M_c)}{18 h}

Active crown fire occurs when wind speed exceeds:

.. math::

   U > U_{active} = \frac{3}{\sqrt{CBD}}

Firebrand Spotting (Probability-Based)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: cpp

   void compute_spotting_probability(
       MultiFab &spotting_data,
       const MultiFab &phi,
       const MultiFab &vel,
       const Geometry &geom,
       const RothermelParams &rp,
       const InputParameters::SpottingParams &sp,
       const MultiFab *terrain_slopes = nullptr
   )

   void generate_firebrand_spots(
       MultiFab &phi,
       const MultiFab &spotting_data,
       const MultiFab &vel,
       const Geometry &geom,
       const InputParameters::SpottingParams &sp,
       int step
   )

Generates new ignition points from stochastic firebrand spotting.

**Spotting probability model:**

Spotting probability is weighted by fire intensity and wind speed:

.. math::

   P_{spot}(d) = P_0 \exp\left(-\frac{d}{d_{max}}\right)

where :math:`P_0` is the base probability and :math:`d_{max}` is the maximum
spotting distance. The landing distance is drawn from a lognormal or exponential
distribution.

Albini (1983) Physics-Based Spotting
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: cpp

   void compute_albini_spotting(
       MultiFab& phi,
       MultiFab& albini_data,
       const MultiFab& vel,
       const MultiFab& R_mf,
       const Geometry& geom,
       const RothermelParams& rp,
       const InputParameters::AlbiniSpottingParams& asp,
       int step
   )

Simulates firebrand spotting using Albini's (1983) thermal-plume lofting formula
coupled with a 2-D horizontal particle trajectory through the wind field.

**Parameters:**

* ``phi``: Level-set function (modified to add spot fires)
* ``albini_data``: Diagnostic output (4 components: Hz, count, dist, active flag)
* ``vel``: Wind velocity field
* ``R_mf``: Rothermel rate-of-spread field [m/s]
* ``geom``: AMReX geometry object
* ``rp``: Rothermel fuel parameters
* ``asp``: Albini spotting parameters
* ``step``: Current time step

**Algorithm:**

1. Compute Byram's fire line intensity :math:`I_B` from the Rothermel ROS field
2. Compute Albini lofting height :math:`H_z = 12.2 I_B^{1/3}`
3. Stochastically decide whether a firebrand is launched (intensity-weighted probability)
4. Integrate 2-D trajectory using forward Euler until the particle descends from :math:`H_z` to ground
5. Place a circular spot-fire ignition zone at the landing location

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
       Real w0;          // Oven-dry total fuel loading [lb/ft²]
       Real sigma;       // Surface-area-to-volume ratio [ft⁻¹]
       Real delta;       // Fuel bed depth [ft]
       Real M_f;         // Aggregate fuel moisture content [fraction]
       Real M_x;         // Moisture of extinction [fraction]
       Real h_heat;      // Heat content [BTU/lb]
       Real S_T;         // Total mineral content [fraction]
       Real S_e;         // Effective mineral content [fraction]
       Real rho_p;       // Oven-dry particle density [lb/ft³]
       // Per-class fuel moistures
       Real M_d1;        // 1-hr dead moisture [fraction]
       Real M_d10;       // 10-hr dead moisture [fraction]
       Real M_d100;      // 100-hr dead moisture [fraction]
       Real M_lh;        // Live herbaceous moisture [fraction]
       Real M_lw;        // Live woody moisture [fraction]
       // Per-class fuel loads (from database or user override)
       Real w_d1;        // 1-hr dead load [lb/ft²]
       Real sigma_d1;    // 1-hr dead SAV [ft⁻¹]
       Real w_d10;       // 10-hr dead load [lb/ft²]
       Real w_d100;      // 100-hr dead load [lb/ft²]
       Real w_lh;        // Live herbaceous load [lb/ft²]
       Real sigma_lh;    // Live herbaceous SAV [ft⁻¹]
       Real w_lw;        // Live woody load [lb/ft²]
       Real sigma_lw;    // Live woody SAV [ft⁻¹]
       // Terrain
       Real slope_x;     // Terrain slope in x-direction [tan(θ)]
       Real slope_y;     // Terrain slope in y-direction [tan(θ)]
       std::string terrain_file;    // Path to terrain file (X Y Z)
       std::string landscape_file;  // Path to FARSITE landscape file
       std::string landscape_fuel_type; // "13" (FBFM13) or "40" (FBFM40)
       // Unit conversions
       Real wind_conv;   // ft/min per (simulation velocity unit)
       Real ros_conv;    // simulation length/time per (ft/min)
   };

**Key members:**

* ``w0``: Total oven-dry fuel loading per unit area
* ``sigma``: Fuel surface-area-to-volume ratio (determines reaction velocity)
* ``delta``: Depth of the fuel bed
* ``M_f``: Aggregate fuel moisture content as a fraction of dry weight
* ``M_x``: Moisture content at which fire will not spread
* ``h_heat``: Heat content (energy per unit mass)
* ``S_T``: Total mineral content (inert material)
* ``S_e``: Effective mineral content (affects reaction intensity)
* ``rho_p``: Density of individual fuel particles
* ``M_d1``, ``M_d10``, ``M_d100``: Per-class dead fuel moistures
* ``M_lh``, ``M_lw``: Per-class live fuel moistures
* ``w_d1`` … ``sigma_lw``: Per-class fuel loads and SAVs (populated from fuel model database)
* ``slope_x``, ``slope_y``: Terrain slope components (rise/run)
* ``landscape_fuel_type``: Fuel model system in the landscape file
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
       int enable;                      // Enable FARSITE model (1=yes, 0=no)
       int use_anderson_LW;             // Use Anderson L/W ratio (1=yes, 0=no)
       Real length_to_width_ratio;      // Fixed L/W ratio (default: 3.0)
       Real coeff_a;                    // Richards' head fire coefficient
       Real coeff_b;                    // Richards' flank fire coefficient
       Real coeff_c;                    // Richards' backing fire coefficient
       Real phi_threshold;              // Fire front detection threshold
       // Bulk fuel consumption
       int use_bulk_fuel_consumption;   // Enable bulk fuel consumption (1=yes, 0=no)
       Real tau_residence;              // Residence time [seconds]
       Real f_consumed_max;             // Maximum consumption fraction (0-1)
       Real f_consumed_min;             // Minimum consumption fraction (0-1)
   };

**Members:**

* ``enable``: Toggle FARSITE model on/off
* ``use_anderson_LW``: Use Anderson (1983) length-to-width ratio to compute coefficients
* ``length_to_width_ratio``: Fixed L/W ratio when Anderson formula is not used
* ``coeff_a``: Head fire coefficient (maximum spread, downwind)
* ``coeff_b``: Flank fire coefficient (perpendicular to wind)
* ``coeff_c``: Backing fire coefficient (minimum spread, upwind)
* ``phi_threshold``: Level-set value threshold for identifying fire front
* ``use_bulk_fuel_consumption``: Toggle bulk fuel consumption model
* ``tau_residence``: Residence time constant for exponential fuel consumption
* ``f_consumed_max``, ``f_consumed_min``: Bounds on fuel consumption fraction

CrownInitiationParams
^^^^^^^^^^^^^^^^^^^^^

Parameters for Van Wagner crown fire model.

.. code-block:: cpp

   struct CrownInitiationParams {
       int enable;                      // Enable crown fire (1=yes, 0=no)
       Real CBH;                        // Canopy base height [m]
       Real CBD;                        // Canopy bulk density [kg/m³]
       Real FMC;                        // Foliar moisture content [%]
       Real crown_fraction_weight;      // Crown fire weighting factor (0-2)
       int use_metric_units;            // 1 = metric (m, kW), 0 = imperial
   };

**Members:**

* ``enable``: Toggle crown fire model on/off
* ``CBH``: Height from ground to bottom of canopy (canopy base height)
* ``CBD``: Mass of crown fuel per unit canopy volume (canopy bulk density)
* ``FMC``: Foliar moisture content in percent (50–300 %)
* ``crown_fraction_weight``: Scaling factor for crown fire spread rate contribution
* ``use_metric_units``: Unit system for crown fire intensity computations

SpottingParams
^^^^^^^^^^^^^^

Parameters for probability-based firebrand spotting model.

.. code-block:: cpp

   struct SpottingParams {
       int enable;                      // Enable spotting (1=yes, 0=no)
       Real P_base;                     // Base spotting probability (0-1)
       Real k_wind;                     // Wind speed coefficient
       Real I_critical;                 // Critical intensity threshold
       Real d_mean;                     // Mean spotting distance
       Real d_sigma;                    // Lognormal std deviation
       Real d_lambda;                   // Exponential decay rate
       std::string distance_model;      // "lognormal" or "exponential"
       Real lateral_spread_angle;       // Angular spread [degrees]
       Real spot_radius;                // Spot-fire radius [physical units]
       int random_seed;                 // RNG seed (0 = use clock)
       int check_interval;              // Steps between spotting checks
   };

**Members:**

* ``enable``: Toggle spotting on/off
* ``P_base``: Base probability of firebrand generation per active cell
* ``k_wind``: Wind speed scaling coefficient for probability
* ``I_critical``: Minimum fire intensity for spotting to occur
* ``d_mean``: Mean parameter for spotting distance distribution
* ``d_sigma``, ``d_lambda``: Distribution shape parameters
* ``distance_model``: Statistical distribution for landing distance
* ``lateral_spread_angle``: Angular dispersion perpendicular to wind
* ``spot_radius``: Radius of new ignition zones placed at landing locations
* ``random_seed``: Reproducible random number generation seed
* ``check_interval``: Frequency of spotting evaluation

AlbiniSpottingParams
^^^^^^^^^^^^^^^^^^^^

Parameters for Albini (1983) physics-based firebrand spotting.

.. code-block:: cpp

   struct AlbiniSpottingParams {
       int enable;                      // Enable Albini spotting (1=yes, 0=no)
       Real terminal_velocity;          // Firebrand descent velocity [m/s]
       Real P_base;                     // Max launch probability per front cell
       Real I_B_min;                    // Min Byram intensity to launch [kW/m]
       Real spot_radius;                // Radius of spot-fire zone [m]
       int random_seed;                 // RNG seed (0 = use clock)
       int check_interval;              // Steps between spotting checks
       int n_traj_steps;                // Forward-Euler trajectory sub-steps
   };

**Members:**

* ``enable``: Toggle Albini spotting on/off
* ``terminal_velocity``: Constant descent velocity used to compute flight time
* ``P_base``: Maximum probability of launching a firebrand (intensity-weighted)
* ``I_B_min``: Byram fire line intensity threshold below which no firebrands are generated
* ``spot_radius``: Radius of the circular ignition zone placed at each landing location
* ``random_seed``: Reproducible random number generation seed
* ``check_interval``: Frequency of spotting evaluation
* ``n_traj_steps``: Resolution of the forward-Euler trajectory integration

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
