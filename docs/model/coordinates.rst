Coordinates and field representation
====================================

This page derives the coordinate constructions that make the guiding-centre
Lagrangian canonical (:doc:`guiding-centre`, Eq. :eq:`canonical`) and states
what is splined for each ``integ_coords`` choice.
Usage: :doc:`/user-guide/field-sources`.

All constructions produce coordinates :math:`(r,\vartheta,\varphi)` with
:math:`A_r = h_r = 0` and splines of the five functions :math:`A_\vartheta, A_\varphi, h_\vartheta,
h_\varphi, B` with first and second derivatives (``field_can_t`` in
``src/field/field_can_base.f90``). The integrators use nothing else.
Implementation details and the libneo interfaces are in
``DOC/coordinates-and-fields.md`` (:doc:`/development/design-notes`).

Reference and integration coordinates
-------------------------------------

*Reference coordinates* :math:`u` describe inputs and outputs: VMEC
:math:`(s,\vartheta,\varphi)` with the geometric toroidal angle
:math:`\varphi`, or chartmap :math:`(\rho,\vartheta,\zeta)` defined by a table
of Cartesian positions :math:`\mathbf x(u)`. *Integration coordinates* are
the canonical coordinates selected by ``integ_coords``. ``ref_to_integ`` and
``integ_to_ref`` (``src/field_can.f90``) convert start points, end points,
orbit samples and wall-hit points.

Canonical flux coordinates (``integ_coords = 0``)
--------------------------------------------------

VMEC provides straight-field-line flux coordinates :math:`(s,\vartheta,
\varphi)` with :math:`\vartheta = \theta + \lambda_\text{V}` and

.. math::

   A_s = 0, \qquad A_\vartheta = \psi_\text{tor}(s), \qquad
   A_\varphi = -\psi_\text{pol}(s), \qquad \psi_\text{pol}' = \iota\,\psi_\text{tor}',

but :math:`B_s \neq 0`. SIMPLE shifts both angles along the field line by a
generating function :math:`G_c(s,\vartheta_c,\varphi_c)`:

.. math::

   \vartheta = \vartheta_c + \iota(s)\,G_c, \qquad \varphi = \varphi_c + G_c .

The shift keeps field lines straight. It changes the one-form
:math:`A_\vartheta d\vartheta + A_\varphi d\varphi` only by
:math:`d\bigl((\iota\psi_\text{tor} - \psi_\text{pol})G_c\bigr)`, because
:math:`(\iota\psi_\text{tor}-\psi_\text{pol})' = \iota'\psi_\text{tor}`; this
gauge term is dropped, so :math:`A_\vartheta` and :math:`A_\varphi` are
unchanged. The radial covariant field component becomes
:math:`B_s^c = B_s + \iota' G_c B_\vartheta + (\iota B_\vartheta + B_\varphi)
\partial_s G_c`, and :math:`B_s^c = 0` gives a linear ODE along :math:`s` at
fixed :math:`(\vartheta_c,\varphi_c)`:

.. math::

   \frac{\partial G_c}{\partial s} =
   -\frac{B_s + \iota' B_\vartheta G_c}{\iota B_\vartheta + B_\varphi}.

It is integrated outward from the axis, starting from :math:`G_c = 10^{-8}`,
on a grid uniform in :math:`\sqrt s` (``rhs_cancoord`` in
``src/get_canonical_coordinates.F90``), with the VMEC fields evaluated at
:math:`(\vartheta,\varphi)` from the equations above. The new components are

.. math::

   B_\vartheta^c = B_\vartheta + (\iota B_\vartheta + B_\varphi)\,\partial_{\vartheta_c}G_c,
   \quad
   B_\varphi^c = B_\varphi + (\iota B_\vartheta + B_\varphi)\,\partial_{\varphi_c}G_c,
   \quad
   \sqrt{g_c} = \sqrt g\,(1 + \iota\,\partial_{\vartheta_c}G_c + \partial_{\varphi_c}G_c),

with angular derivatives from finite-difference stencils. The integration
radius is :math:`r = s`, with :math:`A_\vartheta = \psi_\text{tor}^\text{edge}\,s`
(``field_can_flux.f90``). The derivation from VMEC data is in
``DOC/canonical_and_boozer_flux_coords_via_VMEC.tex``.

Boozer coordinates (``integ_coords = 2``)
-----------------------------------------

Boozer coordinates are flux coordinates with

.. math::

   B_\vartheta = I(s), \qquad B_\varphi = G(s), \qquad
   B = B(s,\vartheta_B,\varphi_B).

SIMPLE builds them from VMEC by an angle shift computed in libneo
(``boozer_converter.F90``) or reads them from a Boozer chartmap
(:doc:`/boozer-chartmap-schema`). The radial covariant component
:math:`B_s` is not computed and not used: the Boozer Lagrangian sets
:math:`h_s = 0`. :math:`A_\vartheta` and :math:`A_\varphi` are then
one-dimensional, and :math:`h_k = B_k/B` vary on a surface only through
:math:`B`; a field evaluation needs a three-dimensional spline of :math:`B`
and one-dimensional splines of :math:`A_\varphi`, :math:`B_\vartheta` and
:math:`B_\varphi` (:math:`A_\vartheta` is linear in :math:`s`). The integration radius is
:math:`r = s`; a chartmap stores :math:`\rho = \sqrt s`
(``field_can_boozer.f90``).

Meiss coordinates (``integ_coords = 3``)
----------------------------------------

Meiss and Hazeltine (1990) remove :math:`h_r` and :math:`A_r` in any
coordinates with :math:`h_\varphi \neq 0`, without nested flux surfaces.
Starting from reference coordinates :math:`(r,\vartheta,\varphi)` (VMEC with
:math:`r=\sqrt s`, a chartmap with :math:`r=\rho`, a coil field splined on
such a grid, or a SPECTRE volume), SIMPLE changes the toroidal angle and the
gauge,

.. math::

   \varphi = \varphi_c + \lambda(r,\vartheta,\varphi_c), \qquad
   \mathbf A' = \mathbf A - \nabla\chi(r,\vartheta,\varphi_c).

The radial components in :math:`(r,\vartheta,\varphi_c)` are
:math:`h_r + h_\varphi\partial_r\lambda` and
:math:`A_r + A_\varphi\partial_r\lambda - \partial_r\chi`. Setting both to
zero gives two ODEs in :math:`r` at fixed :math:`(\vartheta,\varphi_c)`,

.. math::

   \frac{\partial\lambda}{\partial r} = -\frac{h_r}{h_\varphi}, \qquad
   \frac{\partial\chi}{\partial r} = A_r + A_\varphi\frac{\partial\lambda}{\partial r},

evaluated at :math:`\varphi = \varphi_c + \lambda` and integrated outward
from :math:`r_\text{min} = 10^{-6}` with :math:`\lambda = \chi = 0`
(``rh_can`` in ``src/field/field_can_meiss.f90``; adaptive RK with relative
tolerance ``canonical_ode_relerr``; identity where :math:`|h_r| < 10^{-9}`
along the slice). The new angular components,

.. math::

   A'_\vartheta = A_\vartheta + A_\varphi\partial_\vartheta\lambda - \partial_\vartheta\chi,
   \quad
   A'_\varphi = A_\varphi(1 + \partial_{\varphi_c}\lambda) - \partial_{\varphi_c}\chi,
   \quad
   h'_\vartheta = h_\vartheta + h_\varphi\partial_\vartheta\lambda,
   \quad
   h'_\varphi = h_\varphi(1 + \partial_{\varphi_c}\lambda),

and :math:`B` are tabulated on a grid of ``canonical_grid_nr`` :math:`\times`
``canonical_grid_ntheta`` :math:`\times` ``canonical_grid_nphi`` points over one
field period and splined in three dimensions (``init_canonical_field_components``).
All five functions depend on all three coordinates, so a field evaluation
costs one batched three-dimensional spline of five components.

Albert coordinates (``integ_coords = 4``)
-----------------------------------------

The Newton solves of the integrators find :math:`r` from :math:`p_\vartheta`
(:doc:`integrators`). Albert coordinates replace the Meiss radius by

.. math::

   \psi = A'_\vartheta(r,\vartheta,\varphi_c)/A_{\vartheta,\text{norm}}, \qquad
   A_{\vartheta,\text{norm}} = \pm\max_\text{grid}|A'_\vartheta|,

with the sign of :math:`A'_\vartheta` at the outer boundary, so that
:math:`A_\vartheta = A_{\vartheta,\text{norm}}\,\psi` is linear and
:math:`\partial_\psi A_\vartheta` constant. The other four functions are
regridded from :math:`r` to a uniform :math:`\psi` grid covered by every
:math:`(\vartheta,\varphi_c)` line (``src/field/field_can_albert.f90``).

Field sources
-------------

VMEC :math:`R`, :math:`Z`, :math:`\lambda_\text{V}` are splined in :math:`s`
(order ``ns_s``) and the angles (order ``ns_tp``, ``multharm`` points per
highest harmonic), with near-axis harmonics continued by ``axis_healing``.
Coil fields (Biot–Savart, ``src/field/field_coils.f90``) are splined on a
reference grid before the Meiss or Albert construction. A chartmap
tabulates :math:`\mathbf x(\rho,\vartheta,\zeta)` (and :math:`B`,
:math:`B_\vartheta`, :math:`B_\varphi`, :math:`A_\varphi` for Boozer
chartmaps) and may extend beyond the LCFS. Each SPECTRE volume has its own
Meiss chart on :math:`\rho_g \in [k-1, k]` (:doc:`extensions`).

References
----------

* J. D. Meiss and R. D. Hazeltine, Canonical coordinates for guiding center
  particles, Phys. Fluids B 2, 2563 (1990), https://doi.org/10.1063/1.859380
* C. G. Albert, G. S. Grassler, S. V. Kasilov, M. Markl, J. Schatzlmayr,
  Symplectic integration of guiding-center equations in canonical
  coordinates for general toroidal fields, arXiv:2503.15278 (2025),
  https://arxiv.org/abs/2503.15278
