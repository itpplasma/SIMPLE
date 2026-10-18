Guiding-centre model
====================

This page derives the four-dimensional canonical system that the symplectic
integrators advance, its normalisation in the code, and the invariants and
the trapping parameter reported in the outputs.
Usage: :doc:`/user-guide/running`.

Phase-space Lagrangian
----------------------

A particle of mass :math:`m = n_d m_p` and charge :math:`e = n_e e_0`
(``n_d``, ``n_e``) moves in a static field
:math:`\mathbf B = \nabla\times\mathbf A` without electric field. In Gaussian
units and arbitrary spatial coordinates :math:`x^k`, the guiding-centre
phase-space Lagrangian to first order in :math:`\rho_L/L` is (Littlejohn 1983)

.. math::
   :label: lagrangian

   L = \Bigl(m v_\parallel h_k + \frac{e}{c} A_k\Bigr)\dot x^k
       + \frac{mc}{e}\,\mu\,\dot\gamma - H,
   \qquad
   H = \frac{m v_\parallel^2}{2} + \mu B ,

with covariant components :math:`h_k = B_k/B`, gyrophase :math:`\gamma` and
magnetic moment :math:`\mu = m v_\perp^2/(2B)`. :math:`\gamma` is cyclic, so
:math:`\mu` is constant; :math:`H` has no explicit time dependence, so the
energy is constant. Position and :math:`v_\parallel` remain: four
phase-space variables.

Canonical form
--------------

If the coordinates :math:`(r,\vartheta,\varphi)` satisfy
:math:`A_r = h_r = 0` (construction in :doc:`coordinates`), the
:math:`\dot r` term of :eq:`lagrangian` vanishes and

.. math::
   :label: canonical

   L = p_\vartheta\dot\vartheta + p_\varphi\dot\varphi - H, \qquad
   p_\vartheta = m v_\parallel h_\vartheta + \frac{e}{c}A_\vartheta, \qquad
   p_\varphi = m v_\parallel h_\varphi + \frac{e}{c}A_\varphi ,

so :math:`(\vartheta,p_\vartheta)` and :math:`(\varphi,p_\varphi)` are
canonical pairs. The integrators carry

.. math::

   z = (r,\vartheta,\varphi,p_\varphi)

and evaluate the field at :math:`(r,\vartheta,\varphi)`. Everything else is a
function of :math:`z`:

.. math::
   :label: vpar

   v_\parallel(z) = \frac{p_\varphi - (e/c)A_\varphi}{m h_\varphi}, \qquad
   p_\vartheta(z) = m v_\parallel h_\vartheta + \frac{e}{c}A_\vartheta, \qquad
   H(z) = \frac{m v_\parallel^2}{2} + \mu B .

``get_val`` in ``src/field_can.f90`` implements :eq:`vpar`;
``get_derivatives`` and ``get_derivatives2`` add first and second
:math:`z`-derivatives by the chain rule. Since :math:`r` replaces
:math:`p_\vartheta`, Hamilton's equations follow by the implicit-function
theorem at fixed :math:`p_\vartheta(z)`:

.. math::
   :label: hamilton

   \dot\vartheta = \frac{\partial_r H}{\partial_r p_\vartheta}, \quad
   \dot p_\vartheta = -\partial_\vartheta H
      + \frac{\partial_r H}{\partial_r p_\vartheta}\,\partial_\vartheta p_\vartheta, \quad
   \dot\varphi = \partial_{p_\varphi} H
      - \frac{\partial_r H}{\partial_r p_\vartheta}\,\partial_{p_\varphi} p_\vartheta, \quad
   \dot p_\varphi = -\partial_\varphi H
      + \frac{\partial_r H}{\partial_r p_\vartheta}\,\partial_\varphi p_\vartheta .

With :math:`\partial_{p_\varphi} H = v_\parallel/h_\varphi` and
:math:`\partial_{p_\varphi} p_\vartheta = h_\vartheta/h_\varphi` (normalised
units) the third equation becomes
:math:`\dot\varphi = (v_\parallel - h_\vartheta\,\partial_r H/\partial_r
p_\vartheta)/h_\varphi`, as coded in ``sympl_euler1_advance_angles`` and
``gauss_canfield_residual``. Albert, Kasilov and Kernbichler (2020, JCP)
introduced this non-canonical quadrature. It requires
:math:`\partial_r p_\vartheta \approx (e/c)\,\psi_\text{tor}'(r) \neq 0`,
which fails at the axis for :math:`r=\rho`; for :math:`r=s` the polar
singularity at :math:`s=0` remains (treatment in :doc:`integrators`).

Normalisation
-------------

Inside the code all lengths are in cm, fields in G, and the reference speed
is the birth speed

.. math::

   v_0 = \sqrt{2E_0/m}, \qquad E_0 = 3.5\,\text{MeV}/\texttt{facE\_al},
   \qquad \rho_0 = \frac{m c v_0}{e}\ \ (\texttt{ro0},\ \text{G cm}),

set in ``src/params.f90`` (``params_init``). Time is the path length
:math:`\tau = v_0 t`. The public state of a marker is
:math:`z_\text{out} = (r,\vartheta,\varphi,v/v_0,\lambda)`.

The symplectic integrators use a second set of normalised quantities, marked
by a bar, in which the Hamiltonian has the form
:math:`\bar H = \bar v_\parallel^2/2 + J_\perp B`:

.. csv-table::
   :header: "Quantity", "Definition", "In the code"
   :widths: 24 46 30

   "parallel velocity", ":math:`\bar v_\parallel = \sqrt2\,v_\parallel/v_0 = \sqrt2\,(v/v_0)\lambda`", "``f%vpar``"
   "perpendicular invariant", ":math:`J_\perp = 2\mu/(m v_0^2) = (v/v_0)^2(1-\lambda^2)/B`", "``f%mu`` (G\ :sup:`-1`)"
   "Hamiltonian", ":math:`\bar H = 2H/(m v_0^2) = (v/v_0)^2`", "``f%H``"
   "gyroradius parameter", ":math:`\bar\rho_0 = \rho_0/\sqrt2`", "``f%ro0``"
   "canonical momenta", ":math:`\bar p_k = p_k/(m v_0/\sqrt2) = \bar v_\parallel h_k + A_k/\bar\rho_0`", "``si%z(4)`` :math:`=\bar p_\varphi`, ``f%pth``"
   "time", ":math:`\bar t = \tau/\sqrt2`", "``si%dt`` :math:`=\Delta\tau/\sqrt2`"

With these definitions Hamilton's equations keep their form,
:math:`d\bar q/d\bar t = \partial\bar H/\partial\bar p`, because
:math:`\partial H/\partial p = (v_0/\sqrt2)\,\partial\bar H/\partial\bar p`.
``src/simple.f90`` (``encode_symplectic_state``) is the single conversion
from :math:`z_\text{out}` to :math:`(r,\vartheta,\varphi,\bar p_\varphi)`;
the reverse map is
:math:`v/v_0 = \sqrt{\bar H} = \sqrt{J_\perp B + \bar v_\parallel^2/2}` and
:math:`\lambda = \bar v_\parallel/(\sqrt2\,v/v_0)`. After a collision
``update_momentum`` re-encodes the state in the same way (:doc:`collisions`).

Invariants and trapping parameter
---------------------------------

For each marker the outputs report, at the start point
(``compute_pitch_angle_params`` in ``src/simple_main.f90``),

.. math::

   J_\perp = \frac{(v/v_0)^2(1-\lambda^2)}{B}, \qquad
   \kappa = \Bigl((1-\lambda^2)\frac{B_\text{max}}{B} - 1\Bigr)
            \frac{B_\text{min}}{B_\text{max}-B_\text{min}},

where :math:`B_\text{min}` and :math:`B_\text{max}` are the extrema of
:math:`B` on the start surface (``get_bminmax``). :math:`\kappa < 0` is
equivalent to :math:`\lambda^2 > 1 - B/B_\text{max}`, i.e. a passing orbit;
:math:`\kappa = 0` is the trapped-passing boundary and :math:`\kappa = 1` a
deeply trapped particle at :math:`B_\text{min}`; ``contr_pp`` skips passing
markers by :math:`\kappa`. Without collisions :math:`\bar H`,
:math:`J_\perp` and, in axisymmetry, :math:`\bar p_\varphi` are invariants;
``src/orbit_invariants.f90`` converts between states and invariants.

References
----------

* R. G. Littlejohn, Variational principles of guiding centre motion,
  J. Plasma Phys. 29, 111 (1983), https://doi.org/10.1017/S002237780000060X
* C. G. Albert, S. V. Kasilov, W. Kernbichler, Symplectic integration with
  non-canonical quadrature for guiding-center orbits in magnetic confinement
  devices, J. Comput. Phys. 403, 109065 (2020),
  https://doi.org/10.1016/j.jcp.2019.109065
