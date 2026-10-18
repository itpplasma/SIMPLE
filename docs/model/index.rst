Model and methods
=================

These pages state the equations SIMPLE solves, as the code solves them, with
the implementing routine for each.

.. toctree::
   :maxdepth: 2

   guiding-centre
   coordinates
   integrators
   collisions
   classification
   wall-intersection
   heat-loads
   extensions

.. figure:: /_static/figures/pipeline.svg
   :alt: Data flow of SIMPLE: an equilibrium or coil field is splined,
         transformed to canonical coordinates, particles are sampled on a flux
         surface, traced in parallel with a symplectic integrator, classified,
         intersected with the wall, and written as loss statistics.
   :width: 100%

   Phases of a ``simple.x`` run.

Phases of a run
---------------

1. **Field.** Read a VMEC equilibrium, a chartmap, a coil file or a SPECTRE
   equilibrium; construct coordinates in which the guiding-centre Lagrangian
   is canonical and spline :math:`A_\vartheta, A_\varphi, h_\vartheta,
   h_\varphi, B` (:doc:`coordinates`).
2. **Sampling.** Draw start positions in reference coordinates with speed
   :math:`v_0` and pitch :math:`\lambda` (:doc:`/user-guide/running`).
3. **Tracing.** Advance each marker in its own OpenMP iteration with the
   selected integrator (:doc:`integrators`), optionally followed by a Monte
   Carlo collision step per microstep (:doc:`collisions`), until it crosses
   the loss boundary or the wall (:doc:`wall-intersection`) or reaches
   ``trace_time``.
4. **Classification** (optional): label regular and chaotic orbits from
   their Poincaré footprints (:doc:`classification`).
5. **Output** (:doc:`/user-guide/outputs`).

Assumptions
-----------

* Test particles in a static magnetic field; no electric field.
* Guiding-centre motion to first order in :math:`\rho_L/L`; the guiding
  centre, not the particle, defines losses and wall hits. Full-orbit tracing
  is an alternative model (:doc:`extensions`).
* Non-relativistic (``rmu = 1e8`` in ``src/params.f90``).
* Collisions are a test-particle operator with a Maxwellian background of
  electrons and two ion species, using alpha-particle mass and charge
  (:doc:`collisions`).
* A marker is lost when it leaves the last closed flux surface or the
  chartmap domain, or hits an STL wall.
