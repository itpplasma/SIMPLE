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

Runtime input keeps the usual SIMPLE meaning: ``sbeg`` is normalized toroidal
flux ``s`` in both VMEC and chartmap runs. A chartmap run evaluates the start
surface at ``rho = sqrt(sbeg)``.

Angles and Signs
----------------

SIMPLE chartmaps use the same left-handed Boozer orientation as the VMEC path.
The scalar global attribute ``torflux`` is the edge value of ``A_theta``, the
poloidal covariant component of the vector potential:

.. code-block:: text

   A = A_theta grad(theta_B) + A_phi grad(zeta_B)
   A_theta = torflux*s
   A_phi   = -chi

GVEC Boozer coordinates are right-handed for the W7-X/GVEC comparisons, with
``zeta_GVEC = -phi_VMEC``. A GVEC chartmap therefore needs exactly one angle
reversal before SIMPLE reads it:

.. list-table::
   :header-rows: 1

   * - Reversal
     - Resampling
     - Components that change sign
     - Components that keep sign
   * - ``zeta -> -zeta`` (default exporter path)
     - geometry and ``Bmod`` at ``-zeta``
     - ``A_phi``, ``B_phi``
     - ``torflux``/``A_theta``, ``B_theta``
   * - ``theta -> -theta``
     - geometry and ``Bmod`` at ``-theta``
     - ``torflux``/``A_theta``, ``B_theta``
     - ``A_phi``, ``B_phi``

Raw right-handed GVEC chartmaps and double-flipped chartmaps do not match the
SIMPLE/VMEC convention.

Runtime Scaling
---------------

Chartmap files store base-scale CGS quantities. SIMPLE applies the usual
``vmec_B_scale`` and ``vmec_RZ_scale`` settings when it loads a chartmap, so a
chartmap run scales like a VMEC run:

.. list-table::
   :header-rows: 1

   * - Quantity
     - Runtime scale
   * - ``Bmod``
     - ``vmec_B_scale``
   * - ``B_theta``, ``B_phi``
     - ``vmec_B_scale * vmec_RZ_scale``
   * - ``A_phi``, ``torflux``
     - ``vmec_B_scale * vmec_RZ_scale**2``
   * - ``x``, ``y``, ``z``, derived ``rmajor``
     - ``vmec_RZ_scale``

``test_chartmap_scaling`` checks the field object, canonical Boozer splines,
reference-coordinate wrapper, vector potential, covariant field components,
and ``Bmod``.

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
