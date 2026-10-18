Wall intersection
=================

This page defines the loss boundaries, the segment–triangle test against an
STL wall, and the accuracy of the resulting hit points and times.
Usage: :doc:`/user-guide/wall-heat-loads`.

Loss boundaries
---------------

Without a wall, a marker is lost when its guiding centre crosses the outer
boundary of the integration domain: :math:`s = 1` for VMEC and Boozer
coordinates, or the outer boundary of a chartmap. The crossing is located
inside the step (:doc:`integrators`); the exit code is 1.

For chartmap reference coordinates ``chart_boundary_kind`` sets the meaning
of the outer chartmap boundary (``init_stl_wall_if_enabled`` in
``src/simple_main.f90``):

.. csv-table::
   :header: "Value", "Crossing the outer chartmap boundary"
   :widths: 18 82

   "``auto``", "``lcfs`` if the chartmap attribute ``rho_lcfs`` is 1, else ``domain`` (default)"
   "``lcfs``", "loss at the LCFS, exit code 1"
   "``wall``", "loss at the wall, exit code 2"
   "``domain``", "end of the numerical domain, exit code 101; not a physical loss"

A chartmap with ``rho_lcfs`` :math:`< 1` extends beyond the LCFS; an STL wall
inside that domain is then the physical loss surface.

Segment–triangle test
---------------------

``wall_input`` is read into a CGAL surface mesh with an AABB tree
(``src/wall/stl_wall_cgal.cpp``, exact-predicates/inexact-constructions
kernel), with lengths scaled to metres by ``wall_units`` (``m`` or ``mm``).
A wall requires chartmap reference coordinates.

For every accepted microstep SIMPLE maps the guiding-centre positions at
its start and end, :math:`\mathbf x_0` and :math:`\mathbf x_1`, through the
chartmap to Cartesian coordinates and intersects the straight segment with
the mesh (``stl_wall_first_hit_segment_with_normal``). Of all intersections,
the one closest to :math:`\mathbf x_0` is the hit. ``locate_wall_segment``
then records

* the hit point :math:`\mathbf x_\text{hit}` (``wall_hit_cart``) and the
  triangle index,
* the unit normal :math:`\mathbf n` of the triangle,
* the incidence angle :math:`\alpha`, :math:`\cos\alpha =
  |\hat{\mathbf d}\cdot\mathbf n|` with :math:`\hat{\mathbf d} = (\mathbf x_1
  - \mathbf x_0)/|\mathbf x_1 - \mathbf x_0|` (0 is normal incidence),
* the fraction :math:`f = |\mathbf x_\text{hit} - \mathbf x_0| / |\mathbf x_1
  - \mathbf x_0|`, which interpolates the loss time and :math:`(v/v_0,
  \lambda)` linearly within the microstep,

maps :math:`\mathbf x_\text{hit}` back to reference coordinates (by linear
interpolation in :math:`u` if the chartmap inversion fails) and ends the
orbit with exit code 2.

Query gate
----------

The wall is non-convex: a chord between two points inside the LCFS can cut
it, e.g. on the inboard side of a bean-shaped section. SIMPLE therefore
queries the wall only for microsteps ending at :math:`\rho >
\rho_\text{LCFS}`. If a microstep also crosses the outer chartmap boundary,
the segment up to that crossing is tested and an earlier wall hit wins.

Accuracy
--------

* The guiding centre, not the particle, is intersected: impact points are
  off by up to one gyroradius (a few cm for 3.5 MeV alphas in a reactor),
  and :math:`\alpha` is the angle of the guiding-centre chord, not of the
  particle velocity.
* The chord deviates from an orbit of curvature :math:`k` by about
  :math:`k\ell^2/8` for chord length :math:`\ell`; increase ``npoiper2``
  until hit statistics converge, especially at grazing incidence.

Albert et al. (2026) apply the method to stellarator reactors and compare
with BEAMS3D.

References
----------

* C. G. Albert, G. Grassler, J. Lion, J.-F. Lobsien, J. Schatzlmayr,
  F. Warmer, Symplectic methods for optimization of wall heat loads by fusion
  alpha particles in stellarator reactors, Nucl. Mater. Energy 48, 102200
  (2026), https://doi.org/10.1016/j.nme.2026.102200
