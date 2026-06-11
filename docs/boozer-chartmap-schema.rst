Boozer Chartmap Schema
======================

Extended Boozer chartmap files are NetCDF inputs for SIMPLE runs without a
VMEC file at runtime. The file stores geometry and field data on uniform grids,
because the SIMPLE/libneo spline path assumes a constant spacing for each
abscissa.

Radial Grids
------------

The schema uses two radial coordinates:

* ``rho`` is uniform and increasing. It is the chartmap coordinate
  ``rho = sqrt(s)``. Geometry, ``Bmod``, ``B_theta``, and ``B_phi`` use this
  grid.
* ``s`` is uniform and increasing. It spans ``rho(1)**2`` to
  ``rho(n_rho)**2``. ``A_phi`` uses this grid.

``A_phi`` must be a one-dimensional variable over ``s`` and must carry
``radial_abscissa = "s"``. Files that store ``A_phi`` over ``rho`` are invalid.

Dimensions
----------

The logical SIMPLE order is shown below. Fortran writers use this order
directly. Python ``netCDF4`` writers in this repository use reversed dimension
tuples for multidimensional variables and transpose the data before writing.

.. list-table::
   :header-rows: 1

   * - Dimension
     - Meaning
   * - ``rho``
     - Uniform ``sqrt(s)`` grid for geometry, ``Bmod``, ``B_theta``, ``B_phi``.
   * - ``s``
     - Uniform ``s`` grid for ``A_phi``.
   * - ``theta``
     - Endpoint-excluded Boozer poloidal grid.
   * - ``zeta``
     - Endpoint-excluded Boozer toroidal grid.

The file stores no duplicate periodic endpoint planes. Readers append those
planes internally when the spline backend requires a full-period grid.

Variables
---------

.. list-table::
   :header-rows: 1

   * - Variable
     - Logical shape
     - Units
   * - ``rho``
     - ``(rho)``
     - dimensionless
   * - ``s``
     - ``(s)``
     - dimensionless
   * - ``theta``
     - ``(theta)``
     - radians
   * - ``zeta``
     - ``(zeta)``
     - radians
   * - ``x``, ``y``, ``z``
     - ``(rho, theta, zeta)``
     - cm
   * - ``A_phi``
     - ``(s)``
     - G cm^2
   * - ``B_theta``, ``B_phi``
     - ``(rho)``
     - G cm
   * - ``Bmod``
     - ``(rho, theta, zeta)``
     - G
   * - ``num_field_periods``
     - scalar
     - dimensionless

Required Attributes
-------------------

Global attributes:

* ``rho_convention = "rho_tor"``
* ``zeta_convention = "boozer"``
* ``rho_lcfs``
* ``boozer_field = 1``
* ``torflux`` in G cm^2

Variable attributes:

* ``A_phi:radial_abscissa = "s"``
* ``x:units = "cm"``, ``y:units = "cm"``, ``z:units = "cm"``

``rmajor`` is not a schema attribute. The reader derives it from the innermost
geometry surface.
