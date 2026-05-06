.. _wind_solver:

Terrain-Following Mass-Consistent Wind Solver
=============================================

The ``wind_solver`` executable is a stand-alone 3-D wind diagnostic tool that
constructs a terrain-aware, mass-consistent wind field over complex terrain.
It follows the approach pioneered by QUIC-URB (Röckle 1990, Pardyjak & Brown
2001) and the mass-consistent models of Sherman (1978) and Mathiesen (1987).

Overview
--------

Given a terrain file, a reference wind vector at a specified height, and a
surface roughness length, ``wind_solver``:

1. **Reads the terrain** — an arbitrary-density X Y Z point cloud.
2. **Interpolates terrain elevation** onto the computational grid columns
   using inverse-distance weighting (IDW).
3. **Builds a log-law initial wind field** — the von Kármán log-law profile is
   evaluated at the height above local terrain (AGL) for every grid cell.
   Sub-surface cells are set to zero.
4. **Solves the mass-consistent Poisson equation** for the Lagrange multiplier
   λ via AMReX MLMG (``MLABecLaplacian``), enforcing ∇·**u** = 0.
5. **Corrects the wind field** using the Lagrange gradient and writes an
   AMReX plotfile containing all relevant fields.

The solver operates on a **single-level (level-0) Cartesian grid** with
user-specified horizontal (dx, dy) and vertical (dz) spacings.  Vertical
domain extent is 300 m by default, giving a refined near-surface resolution
when dz < dx or dy.

Physical Model
--------------

**Log-law profile**

.. math::

   u(z_\text{agl}) = \frac{u_*}{\kappa}\,\ln\!\left(\frac{z_\text{agl}+z_0}{z_0}\right),
   \qquad u_* = \frac{\kappa\,|\mathbf{U}_\text{ref}|}
                     {\ln\!\left(\dfrac{z_\text{ref}+z_0}{z_0}\right)}

where κ = 0.41 is the von Kármán constant, z₀ is the aerodynamic roughness
length, z_ref is the reference height, and z_agl = z − z_s(x,y) is the height
above the local terrain surface.

**Mass-consistent correction** (Lagrange multiplier method, Sherman 1978)

Minimise

.. math::

   E = \int\!\left[\frac{(u-u_0)^2}{\alpha_h^2}
                  +\frac{(v-v_0)^2}{\alpha_h^2}
                  +\frac{(w-w_0)^2}{\alpha_v^2}\right]\mathrm{d}V

subject to ∇·**u** = 0.  The Euler–Lagrange conditions give

.. math::

   \mathbf{u} = \mathbf{u}_0 - \left(\alpha_h^2\frac{\partial\lambda}{\partial x},\;
                                       \alpha_h^2\frac{\partial\lambda}{\partial y},\;
                                       \alpha_v^2\frac{\partial\lambda}{\partial z}\right)

and substituting into ∇·**u** = 0 yields the anisotropic Poisson equation

.. math::

   -\left(\alpha_h^2\frac{\partial^2\lambda}{\partial x^2}
         +\alpha_h^2\frac{\partial^2\lambda}{\partial y^2}
         +\alpha_v^2\frac{\partial^2\lambda}{\partial z^2}\right)
   = -\nabla\cdot\mathbf{u}_0

solved with AMReX MLMG (``MLABecLaplacian``).

**Boundary conditions for λ**

* x-faces (inflow / outflow): Dirichlet λ = 0
* y-faces (lateral): Neumann ∂λ/∂y = 0
* z-faces (ground, top): Neumann ∂λ/∂z = 0

Building
--------

The wind solver is built automatically when ``LEVELSET_BUILD_WIND_SOLVER=ON``
(the default) and a **3-D build** is selected (``LEVELSET_DIM_2D=OFF``).
Enabling the wind solver also enables AMReX's linear solver library
(``AMReX_LINEAR_SOLVERS``).

.. code-block:: bash

   # Default 3-D build – wind_solver is included automatically
   cmake -S . -B build
   cmake --build build -j

   # Explicitly enable (also the default)
   cmake -S . -B build -DLEVELSET_BUILD_WIND_SOLVER=ON
   cmake --build build -j

   # Disable the wind solver (linear solvers then also disabled for speed)
   cmake -S . -B build -DLEVELSET_BUILD_WIND_SOLVER=OFF
   cmake --build build -j

Usage
-----

.. code-block:: bash

   # From an input file
   ./build/wind_solver inputs.i

   # Inline key=value pairs
   ./build/wind_solver terrain_file=terrain.csv U_ref=8.0 z0=0.05

Input Parameters
----------------

All parameters are read by AMReX ``ParmParse`` and can be placed in an
``inputs.i`` file or passed on the command line as ``key=value`` pairs.

.. list-table::
   :header-rows: 1
   :widths: 20 15 65

   * - Parameter
     - Default
     - Description
   * - ``terrain_file``
     - ``terrain.csv``
     - Path to terrain point-cloud file (X Y Z, whitespace or comma separated; ``#`` comments supported)
   * - ``U_ref``
     - ``10.0``
     - Reference wind x-component [m/s] at height ``z_ref``
   * - ``V_ref``
     - ``0.0``
     - Reference wind y-component [m/s] at height ``z_ref``
   * - ``z_ref``
     - ``10.0``
     - Reference height above local terrain [m]
   * - ``z0``
     - ``0.1``
     - Aerodynamic roughness length [m]
   * - ``dx``
     - ``30.0``
     - Horizontal grid spacing in x [m]
   * - ``dy``
     - ``30.0``
     - Horizontal grid spacing in y [m]
   * - ``dz``
     - ``30.0``
     - Vertical grid spacing [m] (set smaller for near-surface refinement)
   * - ``domain_height``
     - ``300.0``
     - Vertical domain extent above terrain base [m]
   * - ``alpha_h``
     - ``1.0``
     - Horizontal Lagrange anisotropy coefficient (α_h in the Poisson equation)
   * - ``alpha_v``
     - ``1.0``
     - Vertical Lagrange anisotropy coefficient (α_v in the Poisson equation)
   * - ``mlmg_verbose``
     - ``1``
     - MLMG solver verbosity (0 = silent, 4 = maximum)
   * - ``tol_rel``
     - ``1.0e-8``
     - MLMG relative convergence tolerance
   * - ``max_grid_size``
     - ``32``
     - Maximum AMReX box size per dimension
   * - ``plot_file``
     - ``plt_wind``
     - Output plotfile prefix

Terrain File Format
-------------------

The terrain file must contain one data point per line with columns
**X  Y  Z** (in metres, UTM or local coordinates).  Both whitespace and
comma-separated formats are accepted.  Lines beginning with ``#`` are
treated as comments.

.. code-block:: text

   # X [m]  Y [m]  Z [m]
   0.0      0.0    5.2
   30.0     0.0    8.1
   60.0     0.0   12.7
   ...

The horizontal domain extents (x_lo, x_hi, y_lo, y_hi) are derived
automatically from the min/max of the terrain data.  The number of grid cells
is determined from these extents and the requested dx / dy.

Output Plotfile Fields
----------------------

The output AMReX plotfile contains the following cell-centred components:

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Variable
     - Description
   * - ``u``
     - Corrected x-wind [m/s]
   * - ``v``
     - Corrected y-wind [m/s]
   * - ``w``
     - Corrected z-wind [m/s]
   * - ``vel_magnitude``
     - Wind speed |**u**| [m/s]
   * - ``u0``
     - Initial log-law x-wind [m/s]
   * - ``v0``
     - Initial log-law y-wind [m/s]
   * - ``w0``
     - Initial log-law z-wind [m/s]
   * - ``lambda``
     - Lagrange multiplier λ [m²/s]
   * - ``div_before``
     - ∇·**u**₀ before correction [s⁻¹]
   * - ``div_after``
     - ∇·**u** after correction [s⁻¹]
   * - ``terrain_z``
     - Terrain elevation at the column centre [m]

The plotfile can be visualised with VisIt, ParaView (via the ``AMReX``
reader plugin), or converted to GeoTIFF with ``tools/plotfile_to_geotiff.py``.

Typical Workflow
----------------

1. **Prepare terrain** — generate ``terrain.csv`` from a DEM or SRTM data:

   .. code-block:: bash

      python3 tools/srtm_terrain_reader.py \
          --lat-min 40 --lat-max 40.3 --lon-min -105 --lon-max -104.7

2. **Create inputs.i** — choose grid spacing and reference wind:

   .. code-block:: text

      terrain_file  = terrain.csv
      U_ref         = 8.0     # 8 m/s westerly
      V_ref         = 0.0
      z_ref         = 10.0
      z0            = 0.05    # smooth grass
      dx            = 30.0
      dy            = 30.0
      dz            = 10.0    # refined vertical grid
      domain_height = 300.0
      plot_file     = plt_wind

3. **Run the solver**:

   .. code-block:: bash

      ./build/wind_solver inputs.i

4. **Export to GeoTIFF** for GIS:

   .. code-block:: bash

      python3 tools/plotfile_to_geotiff.py plt_wind \
          -v u v w vel_magnitude terrain_z \
          --utm-origin 450000 4400000 --epsg 32613 --outdir gis_out

Vertical Grid Refinement
------------------------

Because dx, dy, and dz are independent parameters, near-surface resolution
can be increased by setting a smaller dz:

.. code-block:: text

   dx = 30.0   # 30 m horizontal
   dy = 30.0
   dz = 5.0    # 5 m vertical near the surface → 60 levels in 300 m

This is especially useful for capturing the log-law profile shape over
complex terrain where the wind gradient is steepest close to the ground.

Anisotropy Coefficients
-----------------------

The α_h and α_v parameters control the relative weight given to horizontal
versus vertical corrections in the mass-consistency step:

* **α_h = α_v = 1** (default): isotropic correction (standard Helmholtz
  projection); the solver adjusts horizontal and vertical winds equally.
* **α_v < α_h** (e.g. α_h = 1, α_v = 0.01): vertical velocity is penalised
  more heavily, so the solver preferentially adjusts horizontal winds —
  the QUIC-URB default which tends to preserve the log-law profile shape.

References
----------

* Sherman, C.A. (1978). A mass-consistent model for wind fields over complex
  terrain.  *Journal of Applied Meteorology*, 17(3), 312–319.
* Mathiesen, M. (1987). Simulation of wind fields in complex terrain.
  *Bound.-Layer Meteorol.*, 38, 213–226.
* Röckle, R. (1990). *Bestimmung der Strömungsverhältnisse im Bereich
  komplexer Bebauungsstrukturen*.  PhD thesis, TH Darmstadt.
* Pardyjak, E.R. & Brown, M.J. (2001). *QUIC-URB v. 1.1: Theory and User's
  Guide*.  Los Alamos National Laboratory, LA-UR-01-4228.
* Forthofer, J.M. (2007). *Modeling Wind in Complex Terrain for Use in
  Fire Spread Prediction*.  MS thesis, Colorado State University.
