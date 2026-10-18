Extensions of the model
=======================

Full-orbit tracing
------------------

``orbit_model = 7`` replaces the guiding-centre model by the full Lorentz
motion of the particle in Cartesian position and velocity
(``src/orbit_fo_boris.f90``). It uses the Boris drift-rotate-drift scheme,
which is explicit, rotates the velocity exactly for a constant field over a
step, and conserves phase-space volume, but is not symplectic. The field is
the Boozer field of a chartmap (``orbit_coord = 1``), evaluated at the
particle position by inverting the chartmap; :math:`\mathbf{B} = \nabla \times
\mathbf{A}` is built from the Boozer vector potential, so that the field is
exactly tangent to the flux surfaces. Full-orbit tracing resolves the
gyration and needs 100 to 1000 times more steps than guiding-centre tracing.
It is a reference for finite-gyroradius effects and does not yet support
walls or collisions. ``DOC/full-orbit-integration.md``
(:doc:`/reference/design-notes`) discusses the design choices.

SPECTRE equilibria
------------------

Multi-region relaxed MHD (MRxMHD) equilibria from SPECTRE or SPEC consist of
nested volumes with a smooth field inside each volume and current sheets on
the interfaces between volumes, where the tangential field jumps. SIMPLE
represents each volume in its own chart with the stacked radial label
:math:`\rho_g \in [0, M_\text{vol}]` (``integ_coords = 6``).

Inside a volume the guiding-centre equations are integrated in per-volume
Meiss canonical coordinates with the symplectic integrators. At a current
sheet the guiding-centre ordering breaks down, because the field changes on a
scale smaller than the gyroradius. When a guiding centre comes closer than
``spectre_fo_entry_larmor`` radial gyroradii to an interface or to the outer
boundary, SIMPLE switches to the full Lorentz orbit in the sharp per-volume
field, integrated by an adaptive Runge-Kutta method with localised interface
crossings. It switches back to the guiding centre one gyroradius farther
away. Guiding-centre maps across the interfaces with a conserved magnetic
moment were tested and rejected because they produce spurious confinement;
the record and the current status are in ``DOC/spectre-interface-crossing.md``
and ``DOC/spectre-rejected-interface-models.md``
(:doc:`/reference/design-notes`).

GPU loss objective
------------------

The CUDA backend (:doc:`/user-guide/gpu`) integrates large ensembles of
guiding-centre orbits in Boozer coordinates with adaptive Runge-Kutta methods
on the GPU. With ``gpu_landreman = .true.`` it evaluates a segmented loss
objective: the tracing is divided into segments of length ``gpu_t_block``,
losses are weighted by :math:`\exp(-t/\tau_\text{loss})` with
``gpu_loss_tau``, and the run stops early when the weighted loss fraction
exceeds ``gpu_maxloss``. This follows the protocol used for alpha-particle
optimisation by Landreman et al. (2026, arXiv:2606.19523).
