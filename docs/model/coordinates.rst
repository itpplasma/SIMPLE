Coordinates and field representation
====================================

The canonical form of the guiding-centre Lagrangian
(:doc:`guiding-centre`, Eq. :eq:`canonical`) needs coordinates
:math:`(r, \vartheta, \varphi)` in which the radial covariant components of
the vector potential and of the field direction vanish,

.. math::

   A_r = 0, \qquad h_r = 0 .

The integrators then need only five functions of position and their first and
second derivatives: :math:`A_\vartheta`, :math:`A_\varphi`,
:math:`h_\vartheta`, :math:`h_\varphi` and :math:`B`. SIMPLE stores them as
three-dimensional splines (``field_can_t`` in
``src/field/field_can_base.f90``). This page describes how the coordinates
are obtained. ``DOC/coordinates-and-fields.md`` (:doc:`/reference/design-notes`)
documents the implementation in detail.

Reference and integration coordinates
-------------------------------------

SIMPLE distinguishes two coordinate systems.

* **Reference coordinates** describe the input and the output: start
  positions, final positions and orbit files. They are VMEC coordinates
  :math:`(s, \vartheta, \varphi)`, with :math:`s` the normalised toroidal flux
  and :math:`\varphi` the geometric toroidal angle, or chartmap coordinates
  :math:`(\rho, \vartheta, \zeta)` defined by a table of Cartesian positions.
* **Integration coordinates** are the canonical coordinates in which the
  orbits are integrated. They are selected by ``integ_coords``.

Transformations between the two are applied at the start and at the end of
each orbit, and at every macrostep for orbit output.

Flux coordinates
----------------

In straight-field-line flux coordinates :math:`(r, \vartheta, \varphi)` the
vector potential has the form

.. math::

   \mathbf{A} = A_\vartheta(r) \nabla\vartheta + A_\varphi(r) \nabla\varphi,
   \qquad
   A_\vartheta = \psi_\text{tor}(r), \quad A_\varphi = -\psi_\text{pol}(r),

with the toroidal and poloidal fluxes divided by :math:`2\pi`. Thus
:math:`A_r = 0` holds already. The radial covariant component :math:`B_r` of
the field is non-zero in general.

**Canonical flux coordinates** (``integ_coords = 0``) remove :math:`B_r` by a
shift of the angles along the field line,

.. math::

   \vartheta_c = \vartheta + \iota(r) G(r, \vartheta, \varphi), \qquad
   \varphi_c = \varphi + G(r, \vartheta, \varphi),

with a function :math:`G` that is obtained, in its inverse form
:math:`G_c(r, \vartheta_c, \varphi_c)`, from a linear first-order
differential equation in :math:`r` at fixed canonical angles, so that
:math:`B_r^c = 0`. The transformation keeps the field lines straight and
:math:`A_k` unchanged. It is constructed directly from VMEC data
(``src/get_canonical_coordinates.F90``). The derivation is in
``DOC/canonical_and_boozer_flux_coords_via_VMEC.tex``.

**Boozer coordinates** (``integ_coords = 2``) are flux coordinates in which
the covariant angular components of :math:`\mathbf{B}` are flux functions,

.. math::

   B_\vartheta = I(r), \qquad B_\varphi = G(r).

SIMPLE computes them from VMEC by a generating function of the angle shift
(``src/boozer_converter.F90`` of libneo and ``src/field/field_can_boozer.f90``)
or reads them from a Boozer chartmap. The remaining radial component
:math:`B_r^B` is not used: ``field_can_t`` carries only the angular
components, and the log reports ``B_r is not computed``. The Boozer
Lagrangian therefore neglects :math:`h_r`, the standard approximation of
Boozer guiding-centre codes. In return, :math:`A_\vartheta` and
:math:`A_\varphi` are one-dimensional and :math:`h_\vartheta`,
:math:`h_\varphi` depend on angles only through :math:`B`, which makes the
evaluation fast. This is the default for production runs.

Boozer chartmaps store :math:`B`, :math:`B_\vartheta`, :math:`B_\varphi`
and :math:`A_\varphi` over a uniform grid together with the geometry; their
format and the left-handed sign convention are specified in
:doc:`/boozer-chartmap-schema`.

Meiss coordinates
-----------------

Flux coordinates need nested flux surfaces. Fields with islands, stochastic
regions or a region outside the last closed flux surface need another
construction. Following Meiss and Hazeltine (1990), SIMPLE starts from any
reference coordinates :math:`(r, \vartheta, \varphi)` (VMEC, chartmap or
cylindrical) and changes the toroidal angle and the gauge,

.. math::

   \varphi = \varphi_c + \lambda(r, \vartheta, \varphi_c), \qquad
   \mathbf{A}' = \mathbf{A} - \nabla\chi(r, \vartheta, \varphi_c),

with :math:`\lambda` and :math:`\chi` chosen so that :math:`h_r` and
:math:`A_r` vanish in the new coordinates
:math:`(r, \vartheta, \varphi_c)`. This gives the ordinary differential
equations in :math:`r`

.. math::

   \frac{\partial \lambda}{\partial r} = -\frac{h_r}{h_\varphi}, \qquad
   \frac{\partial \chi}{\partial r} = A_r + A_\varphi \frac{\partial \lambda}{\partial r},

integrated outward from the axis with :math:`\lambda = \chi = 0` there, on
each :math:`(\vartheta, \varphi_c)` node of a grid with
``canonical_grid_nr`` :math:`\times` ``canonical_grid_ntheta`` :math:`\times`
``canonical_grid_nphi`` points (relative tolerance ``canonical_ode_relerr``).
The new components are

.. math::

   A'_\vartheta = A_\vartheta + A_\varphi \partial_\vartheta\lambda - \partial_\vartheta\chi,
   \quad
   A'_\varphi = A_\varphi (1 + \partial_{\varphi}\lambda) - \partial_{\varphi}\chi,
   \quad
   h'_\vartheta = h_\vartheta + h_\varphi \partial_\vartheta\lambda,
   \quad
   h'_\varphi = h_\varphi (1 + \partial_{\varphi}\lambda),

which are splined on the same grid (``src/field/field_can_meiss.f90``). The
only requirement is a non-vanishing toroidal covariant field component
:math:`h_\varphi`, which holds in stellarators and tokamaks. Meiss
coordinates (``integ_coords = 3``) are used for coil fields and for tracing
beyond the last closed flux surface.

Albert coordinates
------------------

The integrators solve for the radial variable with a Newton iteration on
:math:`p_\vartheta(r, \ldots)`. **Albert coordinates** (``integ_coords = 4``)
replace :math:`r` of the Meiss coordinates by

.. math::

   \psi = A'_\vartheta(r, \vartheta, \varphi_c) / A_{\vartheta,\text{norm}},

so that :math:`A_\vartheta = A_{\vartheta,\text{norm}}\, \psi` is linear in
the radial variable and :math:`\partial_\psi A_\vartheta` is constant. The
other field components are regridded from :math:`r` to :math:`\psi`
(``src/field/field_can_albert.f90``). This simplifies the radial Newton
solve.

Field input and splines
-----------------------

VMEC
   :math:`R`, :math:`Z` and :math:`\lambda` of VMEC are splined in :math:`s`
   (order ``ns_s``) and in the angles (order ``ns_tp``, grid ``multharm``
   times the largest mode number). Near the axis the harmonics are continued
   by ``axis_healing``. The canonical and Boozer transformations are computed
   from these splines.
Coils
   The Biot-Savart field of filamentary coils is evaluated on the grid of the
   reference coordinates and splined there (``src/field/field_coils.f90``,
   ``field_splined.f90``) before the Meiss or Albert construction.
Chartmaps
   A chartmap gives the Cartesian position :math:`\mathbf{x}(\rho,
   \vartheta, \zeta)` on a uniform grid. Its metric defines the reference
   coordinates. Chartmaps extend reference coordinates beyond the last closed
   flux surface for wall studies, and allow equilibria from GVEC and other
   codes.
SPECTRE
   Multi-region relaxed MHD equilibria are represented per volume, with
   sharp interfaces between volumes (:doc:`extensions`).

References
----------

* J. D. Meiss and R. D. Hazeltine, Canonical coordinates for guiding center
  particles, Phys. Fluids B 2, 2563 (1990), https://doi.org/10.1063/1.859380
* C. G. Albert, S. V. Kasilov, W. Kernbichler, Accelerated methods for direct
  computation of fusion alpha particle losses within stellarator
  optimization, J. Plasma Phys. 86, 815860201 (2020),
  https://doi.org/10.1017/S0022377820000203
