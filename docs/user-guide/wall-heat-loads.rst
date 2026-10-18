Wall heat loads
===============

SIMPLE gives the position and time at which each lost particle leaves the
confinement region. Heat loads follow by weighting each lost marker with its
share of the alpha power and its remaining energy, and dividing by the area on
which it arrives. Two levels of detail are available.

Losses on the last closed flux surface
--------------------------------------

The default run stops a particle when its guiding centre crosses the last
closed flux surface (LCFS), exit code 1. Its final state in ``times_lost.dat``
and ``results.nc`` (``zend``, ``xend_cart``) is the crossing point. Binning
these points in the poloidal and toroidal angle and dividing by the area of
the bins gives the heat flux on the LCFS:

.. math::

   q_\text{bin} = \frac{P_\alpha}{N A_\text{bin}}
                  \sum_{i \in \text{bin}} \left(\frac{v_i}{v_0}\right)^2 ,

with the alpha power :math:`P_\alpha = P_\text{fusion}/5` for D-T, the number
of markers :math:`N` and the bin area :math:`A_\text{bin}`. Without collisions
:math:`v_i = v_0`. :doc:`/tutorials/wall-heat-loads` (Part A) does this for a
reactor-scale QA configuration. ``DOC/wall-heat-loads.md``
(:doc:`/reference/design-notes`) discusses the slowing-down weighting, the
area computation and typical values.

This estimate is fast and needs no wall geometry. It does not include the
distance between plasma and wall, the finite gyroradius (several centimetres
for 3.5 MeV alphas), or shadowing by wall structures.

Losses on an STL first wall
---------------------------

For losses on the real wall geometry, SIMPLE follows the guiding centres
beyond the LCFS and intersects each integration step with a triangulated wall
surface. This needs

* a build with ``-DSIMPLE_ENABLE_CGAL=ON`` (:doc:`/getting-started/installation`),
* ``coord_input``: a chartmap that extends beyond the LCFS and records its
  position in the attribute ``rho_lcfs`` (:doc:`field-sources`),
* ``field_input``: a field that is valid between LCFS and wall, normally the
  coil field, with ``isw_field_type = 3`` or ``4`` (or ``5`` with
  ``integmode = 0``),
* ``wall_input``: the wall as an STL file, with ``wall_units = 'm'`` or
  ``'mm'``,
* ``output_results_netcdf = .true.`` to write the hit points.

.. code-block:: fortran

   &config
   netcdffile = 'wout.nc'
   field_input = 'coils.simple'
   coord_input = 'chartmap_wall.nc'
   isw_field_type = 3
   wall_input = 'first_wall.stl'
   wall_units = 'mm'
   output_results_netcdf = .true.
   /

SIMPLE queries the wall only after an accepted step has ended outside the
LCFS (``rho > rho_lcfs``), because a straight chord between two points inside
the LCFS can cut a non-convex wall without a physical loss. A particle that
hits the wall ends with exit code 2. ``results.nc`` then contains
``wall_hit = 1``, the intersection point ``wall_hit_cart``, the unit normal
of the hit triangle and the incidence angle. A particle that reaches the outer
chartmap boundary without a wall hit ends according to
``chart_boundary_kind``. The algorithm is described in
:doc:`/model/wall-intersection`; :doc:`/tutorials/wall-heat-loads` (Part B)
runs the complete chain on public data.

The positions are guiding-centre positions. The true impact point of the
particle differs by up to one gyroradius. Full-orbit tracing
(``orbit_model = 7``) does not support walls yet.

Known issues at this revision
-----------------------------

* ``pysimple.engineering.WallHeatMap.from_netcdf`` selects lost particles by
  the ``class_lost`` flag of ``results.nc``. That flag is only set by
  classification runs (``src/simple_main.f90``), so the method reports no
  losses for a standard run. It also indexes ``zend`` as (phase, particle),
  while ``results.nc`` stores (particle, phase). The tutorial script reads
  ``orbit_exit_code`` and ``zend`` directly.
* ``python -m libneo.convert_coils_to_simple`` drops the zero-current points
  that close each coil in the STELLOPT format, so that consecutive coils are
  joined by spurious current segments. Use
  ``test/tests/convert_coils_stellopt_to_simple.py``.
