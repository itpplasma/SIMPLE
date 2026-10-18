Model and methods
=================

This part documents what SIMPLE computes and how. It is a technical reference
for the physical model, the coordinates, the numerical integrators, the orbit
classification and the wall intersection. The derivations are published in
the papers listed in :doc:`/reference/publications`; the pages below state the
equations as implemented and point to the source files.

.. toctree::
   :maxdepth: 2

   guiding-centre
   coordinates
   integrators
   classification
   wall-intersection
   extensions

.. figure:: /_static/figures/pipeline.svg
   :alt: Data flow of SIMPLE: an equilibrium or coil field is splined,
         transformed to canonical coordinates, particles are sampled on a flux
         surface, traced in parallel with a symplectic integrator, classified,
         intersected with the wall, and written as loss statistics.
   :width: 100%

   Pipeline of a SIMPLE run.

Overview of a run
-----------------

1. **Field setup.** The magnetic field is read from a VMEC equilibrium, a
   Boozer chartmap, a coil file or a SPECTRE equilibrium. SIMPLE constructs
   coordinates in which the guiding-centre Lagrangian is canonical
   (:doc:`coordinates`) and splines the field components needed by the
   integrator.
2. **Sampling.** Particles start on a flux surface or in a volume with an
   isotropic pitch distribution and the birth energy.
3. **Tracing.** Each guiding-centre orbit is integrated with a symplectic
   method in canonical coordinates (:doc:`integrators`) until it leaves the
   plasma, hits the wall, or reaches the tracing time. Particles are
   distributed over OpenMP threads (or GPU threads).
4. **Classification** (optional). Regular and chaotic orbits are identified
   from a short piece of trajectory (:doc:`classification`).
5. **Output.** Loss times, exit positions, exit codes and confined fractions
   are written per particle and per time.

Scope and assumptions
---------------------

* Test particles: the energetic particles do not modify the field. The field
  is static.
* Guiding-centre motion, to first order in the ratio of gyroradius to the
  scale length of the field (:doc:`guiding-centre`). A full-orbit option and
  hybrid full-orbit treatment for SPECTRE interfaces are described in
  :doc:`extensions`.
* No electric field in the standard model.
* Optional Coulomb collisions with a Maxwellian background as a Monte Carlo
  operator.
* The default loss boundary is the last closed flux surface; a triangulated
  wall can be used instead (:doc:`wall-intersection`).
