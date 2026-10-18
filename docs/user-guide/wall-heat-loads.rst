Wall heat loads
===============

SIMPLE records where and when each lost marker leaves the plasma. This page
sets up the two heat-load estimates: on the last closed flux surface (LCFS)
from a standard run, and on a triangulated first wall from an STL file.
Method: :doc:`/model/heat-loads`. Worked example:
:doc:`/tutorials/wall-heat-loads`.

Losses on the last closed flux surface
--------------------------------------

A standard run stops a marker when its guiding centre crosses the LCFS (exit
code 1). Its final state ``zend`` and ``xend_cart`` in ``results.nc``
(``output_results_netcdf = .true.``) is the crossing point. Binning these
points in :math:`(\vartheta, \varphi)`, weighting each by :math:`(v/v_0)^2`
and dividing by the bin area from a chartmap gives the LCFS heat flux; use
``docs/tutorials/scripts/lcfs_heat_flux.py``. The estimate needs no wall
geometry and ignores plasma–wall distance, gyroradius and shadowing.

Limitation: ``pysimple.engineering.WallHeatMap.from_netcdf`` reports no losses
for a standard run (it selects ``class_lost``, set only by classification
runs) and reads ``zend`` with transposed axes; use the tutorial script.

Losses on an STL first wall
---------------------------

SIMPLE follows the guiding centres beyond the LCFS and intersects each
accepted step with the wall surface.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Setting
     - Requirement
   * - build
     - ``-DSIMPLE_ENABLE_CGAL=ON`` (:doc:`/getting-started/installation`)
   * - ``coord_input``
     - chartmap reaching beyond the LCFS, with attribute ``rho_lcfs`` < 1
       (:doc:`field-sources`)
   * - ``field_input``
     - field valid between LCFS and wall, normally the coil field
   * - ``isw_field_type``
     - ``3`` or ``4`` (canonical on the chartmap), or ``5`` with
       ``integmode = 0``
   * - ``wall_input``
     - STL file, closed surface around the LCFS, facet order as in the CAD
       export
   * - ``wall_units``
     - ``'m'`` (default) or ``'mm'``, the STL length unit
   * - ``output_results_netcdf``
     - ``.true.``, to write the hits

The wall is queried only after an accepted step has ended outside the LCFS
(:math:`\rho > \rho_\text{lcfs}`), because a chord between two points inside
the LCFS can cut a non-convex wall. A hit ends the marker with exit code 2 and
writes ``wall_hit = 1``, ``wall_hit_cart``, ``wall_hit_normal_cart``,
``wall_hit_angle_rad`` (0 is normal incidence) and ``wall_hit_triangle``
(zero-based facet index in STL file order, ``-1`` without hit). A marker that
reaches the outer chartmap boundary without a hit ends according to
``chart_boundary_kind`` (:doc:`/model/wall-intersection`). Positions are
guiding-centre positions; full orbits (``orbit_model = 7``) do not support
walls.

With collisions (``swcoll = .true.`` and a ``&profiles`` namelist,
:doc:`running`), hits occur at reduced energy. Choose ``trace_time`` relative
to the slowing-down time of the profiles: much shorter measures prompt losses,
comparable resolves losses at reduced energy. The hit keeps the momentum of
the kicks before the intersecting step, so hit energies have the time
resolution of the microstep.

Heat load per triangle
----------------------

``pysimple.wall_heat_load.TriangleHeatLoad`` deposits each hit on its facet
(formulas and metric definitions in :doc:`/model/heat-loads`). It needs the
original STL with unchanged facet order; the loader checks file name and facet
count.

.. code-block:: python

   from pysimple.wall_heat_load import TriangleHeatLoad

   load = TriangleHeatLoad.from_netcdf("results.nc", "wall.stl",
                                       total_alpha_power_MW=600.0)
   print(load.summary(threshold_W_m2=1e5))  # peak, p95, wetted area, ...
   load.to_vtk("wall_heat.vtk")     # ParaView, cell data per facet
   load.to_netcdf("wall_heat.nc")   # arrays plus vertices(triangle, vertex, xyz)
   load.to_stl("wall_heat.stl")     # geometry plus wall_heat.stl.nc sidecar

Per facet it provides ``power_W``, ``power_density_W_m2``,
``power_density_per_MW`` (W/m² per MW of source alpha power) and
``hit_count``; geometry is in metres. Equal marker weights are power-correct
only for fusion-weighted births (``startmode = 7``, :doc:`/alpha-birth`); for
other start distributions pass ``marker_weights`` (one relative source weight
per marker). No smoothing is applied: check ``hit_count`` and convergence
with marker number.
