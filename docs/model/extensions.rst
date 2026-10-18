Full orbit, SPECTRE and GPU objective
=====================================

This page states the models that replace or supplement the guiding-centre
system: full-orbit tracing, the per-volume treatment of SPECTRE equilibria,
and the segmented loss objective of the CUDA backend.
Usage: :doc:`/user-guide/running`.

Full-orbit tracing
------------------

``orbit_model = 7`` (with ``orbit_coord = 1``) integrates the Lorentz motion
in Cartesian position and velocity,

.. math::

   \frac{d\mathbf x}{d\bar t} = \bar{\mathbf v}, \qquad
   \frac{d\bar{\mathbf v}}{d\bar t} = \frac{1}{\bar\rho_0}\,\bar{\mathbf v}\times\mathbf B(\mathbf x),

in the normalisation of :doc:`guiding-centre` (``src/orbit_fo_boris.f90``).
:math:`\mathbf B(\mathbf x)` comes from a Boozer chartmap: the chartmap is
inverted to :math:`u = (\rho,\vartheta_B,\varphi_B)` by Newton iteration and
:math:`\mathbf B = \nabla\times\mathbf A` is assembled from the Boozer
potentials :math:`A_\vartheta(s)`, :math:`A_\varphi(s)`, so that
:math:`B^s = 0` exactly; a pseudo-Cartesian chart
:math:`(\rho\cos\vartheta,\rho\sin\vartheta,\varphi)` keeps the inversion
regular at the axis. The step is the Boris drift–rotate–drift scheme with
:math:`\bar{\mathbf t} = \mathbf B\,h/(2\bar\rho_0)`:

.. math::

   \mathbf x^{n+1/2} = \mathbf x^n + \tfrac h2\bar{\mathbf v}^n,\quad
   \bar{\mathbf v}' = \bar{\mathbf v}^n + \bar{\mathbf v}^n\times\bar{\mathbf t},\quad
   \bar{\mathbf v}^{n+1} = \bar{\mathbf v}^n
     + \bar{\mathbf v}'\times\frac{2\bar{\mathbf t}}{1+|\bar{\mathbf t}|^2},\quad
   \mathbf x^{n+1} = \mathbf x^{n+1/2} + \tfrac h2\bar{\mathbf v}^{n+1},

with :math:`\mathbf B` at :math:`\mathbf x^{n+1/2}`. It rotates the velocity
exactly for constant :math:`\mathbf B` and preserves phase-space volume.
The step is the guiding-centre microstep, :math:`h = \Delta\tau/\sqrt2`, so
the gyro-angle per step is :math:`\Delta\tau\,B/\rho_0`: ``npoiper2`` must
make :math:`2\pi R_0 B/(\texttt{npoiper2}\,\rho_0) \ll 1`. The particle
starts one gyroradius from the sampled guiding centre (``init_fo`` in
``src/simple.f90``). A marker is lost when its reconstructed guiding centre reaches
:math:`\rho \ge 1` and the particle is within 0.05 of the edge. Full-orbit
tracing does not support walls, collisions or classification.

SPECTRE equilibria
------------------

A multi-region relaxed MHD equilibrium has :math:`M_\text{vol}` nested
volumes with a smooth field inside each and current sheets on the
interfaces, across which the tangential field jumps. SIMPLE labels position
by :math:`\rho_g \in [0, M_\text{vol}]`, volume :math:`k` covering
:math:`[k-1,k]` (``integ_coords = 6``), and builds a Meiss canonical chart per
volume (:doc:`coordinates`). Inside a volume the symplectic guiding-centre
integrators run unchanged.

At a sheet the guiding-centre ordering fails, because :math:`\mathbf b`
changes over a distance smaller than the gyroradius. Multi-volume tracing
therefore requires ``integmode > 0`` and ``crossing_level = 3``: when the
guiding centre comes within ``spectre_fo_entry_larmor`` (default 1.5) radial
gyroradii :math:`\rho_r = \rho_L|\nabla\rho_g|` of an interface or of the
outer boundary, SIMPLE switches to the full Lorentz orbit in the sharp
per-volume field (RK45 with located sheet crossings), and returns beyond
``spectre_fo_entry_larmor`` + 1 radial gyroradii. The algorithm is in
``DOC/spectre-interface-crossing.md``; ``DOC/spectre-rejected-interface-models.md``
derives why fixed-:math:`\mu` guiding-centre maps across a sheet are invalid
(:doc:`/development/design-notes`).

GPU loss objective
------------------

The CUDA backend (:doc:`/user-guide/gpu`) traces guiding-centre orbits in
Boozer coordinates in segments of ``gpu_t_block`` and after each segment
evaluates

.. math::

   W = \frac1N\sum_{i\ \text{lost}} e^{-t_i/\tau_\text{loss}}, \qquad
   \tau_\text{loss} = \texttt{gpu\_loss\_tau},

stopping when :math:`W > \texttt{gpu\_maxloss}`
(``src/simple_cuda_native_rk54.cu``). This follows the objective of Landreman
et al. (2026).

References
----------

* M. Landreman, M. Czekanski, A. Giuliani, B. Jang, R. Conlin, Bayesian
  optimization of stellarator alpha-particle confinement using data-informed
  parameter spaces and dimensionality reduction, arXiv:2606.19523 (2026),
  https://arxiv.org/abs/2606.19523
