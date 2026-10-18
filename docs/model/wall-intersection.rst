Wall intersection
=================

Loss boundaries
---------------

By default a particle is lost when its guiding centre crosses the last closed
flux surface (LCFS): :math:`s = 1` for VMEC and Boozer coordinates, or the
outer boundary of a chartmap whose ``rho_lcfs`` is 1. The crossing is located
within the step as described in :doc:`integrators` (exit code 1).

When the reference coordinates extend beyond the LCFS (a chartmap with
``rho_lcfs`` < 1), the outer chartmap boundary has the meaning set by
``chart_boundary_kind``: a loss at the LCFS (code 1), a loss at the wall
(code 2), or the end of the numerical domain (code 101, the default ``auto``
for such chartmaps). A triangulated wall gives the physical loss surface
inside this domain.

Segment-triangle intersection
-----------------------------

The wall is read from an STL file into a CGAL surface mesh with an axis-aligned
bounding-box (AABB) tree for fast queries (``src/wall/stl_wall_cgal.cpp``,
exact predicates and inexact constructions kernel). Lengths are converted to
metres with ``wall_units``.

During the integration of a particle, SIMPLE keeps the Cartesian position of
the guiding centre at the start and the end of every accepted microstep,
computed from the reference coordinates (chartmap). The straight segment
between them is intersected with all triangles of the wall. If the segment
intersects the wall, the intersection closest to the start point is the hit.
SIMPLE then records

* the hit point :math:`\mathbf{x}_\text{hit}` (``wall_hit_cart``),
* the unit normal :math:`\mathbf{n}` of the hit triangle,
* the incidence angle :math:`\alpha` from
  :math:`\cos\alpha = |\hat{\mathbf{d}} \cdot \mathbf{n}|`, with
  :math:`\hat{\mathbf{d}}` the direction of the segment (0 is normal
  incidence),
* the loss time, interpolated linearly in the fraction
  :math:`|\mathbf{x}_\text{hit} - \mathbf{x}_0| / |\mathbf{x}_1 - \mathbf{x}_0|`
  of the segment,

and ends the orbit with exit code 2. The final state is mapped back to
reference coordinates (``src/simple_main.f90``, ``locate_wall_segment``).

Query gate
----------

A wall is a non-convex surface that encloses the plasma. A straight chord
between two points inside the LCFS can intersect such a surface, for example
on the inboard side of a bean-shaped cross-section, although the orbit itself
does not reach the wall. SIMPLE therefore queries the wall only after an
accepted microstep has ended outside the LCFS, :math:`\rho > \rho_\text{LCFS}`.
If a microstep also reaches the outer boundary of the chartmap, the partial
segment up to the boundary is tested, and an earlier wall intersection takes
precedence.

Accuracy
--------

* The guiding centre, not the particle, is intersected with the wall. The
  particle impact point differs by up to a gyroradius, which smooths the
  heat-load pattern on the scale of a few centimetres for 3.5 MeV alphas in a
  reactor.
* The chord of a microstep approximates the curved orbit. For a smooth orbit
  with local curvature :math:`\kappa` and segment length :math:`\ell`, the
  maximum distance from the arc to its chord is approximately
  :math:`\kappa\ell^2/8`. Increase ``npoiper2`` to reduce the step length;
  check convergence of wall hits, especially near tangential incidence.
* The incidence angle is the angle of the guiding-centre chord with the wall,
  not of the particle velocity.

The method and its application to wall heat loads in stellarator reactors,
including a comparison with BEAMS3D, are published in Albert et al. (2026).

Reference
---------

* C. G. Albert, G. Grassler, J. Lion, J.-F. Lobsien, J. Schatzlmayr,
  F. Warmer, Symplectic methods for optimization of wall heat loads by fusion
  alpha particles in stellarator reactors, Nucl. Mater. Energy 48, 102200
  (2026), https://doi.org/10.1016/j.nme.2026.102200
