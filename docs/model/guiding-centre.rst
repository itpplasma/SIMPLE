Guiding-centre model
====================

Phase-space Lagrangian
----------------------

SIMPLE follows the guiding centre of a charged particle of mass :math:`m` and
charge :math:`e` (here :math:`e = Z e_0` with charge number ``n_e``) in a
static magnetic field :math:`\mathbf{B} = \nabla \times \mathbf{A}` without
electric field. In Gaussian units, the guiding-centre phase-space Lagrangian
in arbitrary spatial coordinates :math:`x^k` is (Littlejohn 1983)

.. math::
   :label: lagrangian

   L = \left( m v_\parallel h_k + \frac{e}{c} A_k \right) \dot x^k
       - J_\perp \dot\phi - H ,
   \qquad
   H = \frac{m v_\parallel^2}{2} + J_\perp \omega_c ,

with the covariant components :math:`h_k = B_k / B` of the unit vector along
the field, the gyrophase :math:`\phi`, the perpendicular adiabatic invariant
:math:`J_\perp = m v_\perp^2 / (2 \omega_c)` and the cyclotron frequency
:math:`\omega_c = e B / (m c)`. The magnetic moment
:math:`\mu = m v_\perp^2 / (2B)` is related by :math:`J_\perp \omega_c = \mu B`.
The guiding-centre position is accurate to first order in the ratio of the
gyroradius :math:`\rho` to the scale length of the field.

:math:`\phi` is cyclic, so :math:`J_\perp` and :math:`\mu` are constants of
motion. :math:`H` does not depend on time, so the energy is conserved. Four
phase-space variables remain: three for the position and one for the parallel
motion.

Canonical form
--------------

In coordinates :math:`(r, \vartheta, \varphi)` in which the radial covariant
components of both :math:`\mathbf{A}` and :math:`\mathbf{h}` vanish,
:math:`A_r = h_r = 0`, the Lagrangian :eq:`lagrangian` becomes canonical:

.. math::
   :label: canonical

   L = p_\vartheta \dot\vartheta + p_\varphi \dot\varphi - H,
   \qquad
   p_\vartheta = m v_\parallel h_\vartheta + \frac{e}{c} A_\vartheta ,
   \qquad
   p_\varphi = m v_\parallel h_\varphi + \frac{e}{c} A_\varphi .

The angles :math:`(\vartheta, \varphi)` and the momenta
:math:`(p_\vartheta, p_\varphi)` form two canonical pairs. How such coordinates
are constructed for VMEC equilibria, Boozer chartmaps and general fields is
described in :doc:`coordinates`.

The integrators do not use :math:`p_\vartheta` as an independent variable.
They use the four variables

.. math::

   z = (r, \vartheta, \varphi, p_\varphi),

in which all field quantities are evaluated. Given :math:`z`, the parallel
velocity and the other momentum follow from :eq:`canonical`:

.. math::
   :label: vpar

   v_\parallel = \frac{p_\varphi - (e/c) A_\varphi(r, \vartheta, \varphi)}
                      {m h_\varphi(r, \vartheta, \varphi)},
   \qquad
   p_\vartheta(z) = m v_\parallel h_\vartheta + \frac{e}{c} A_\vartheta ,
   \qquad
   H(z) = \frac{m v_\parallel^2}{2} + \mu B(r, \vartheta, \varphi).

The radial coordinate :math:`r` is a non-canonical variable, defined
implicitly by :math:`p_\vartheta(z)`. Hamilton's equations for the canonical
pairs read

.. math::
   :label: hamilton

   \dot p_\vartheta = -\left.\frac{\partial H}{\partial \vartheta}\right|_{p},
   \quad
   \dot \vartheta = \left.\frac{\partial H}{\partial p_\vartheta}\right|_{q},
   \quad
   \dot p_\varphi = -\left.\frac{\partial H}{\partial \varphi}\right|_{p},
   \quad
   \dot \varphi = \left.\frac{\partial H}{\partial p_\varphi}\right|_{q},

where the derivatives at constant canonical variables are obtained from
derivatives in :math:`z` by the implicit function theorem, for example

.. math::

   \left.\frac{\partial H}{\partial p_\vartheta}\right|_{q}
     = \frac{\partial_r H}{\partial_r p_\vartheta},
   \qquad
   \left.\frac{\partial H}{\partial \vartheta}\right|_{p}
     = \partial_\vartheta H
       - \frac{\partial_r H}{\partial_r p_\vartheta}\,
         \partial_\vartheta p_\vartheta .

This formulation, symplectic integration with non-canonical quadrature
points, is the subject of Albert, Kasilov and Kernbichler (2020, J. Comput.
Phys.). It requires :math:`\partial_r p_\vartheta \neq 0`. In flux coordinates
:math:`\partial_r p_\vartheta \approx (e/c)\, \partial_r A_\vartheta =
(e/c)\, \psi_\text{tor}'`, which is non-zero away from the magnetic axis. The
code is in ``src/field_can.f90`` (``get_val`` and its derivatives) and
``src/orbit_symplectic*.f90``.

Normalisation
-------------

The code uses Gaussian units with lengths in cm and fields in G. The particle
energy is :math:`E_0 = 3.5\,\text{MeV}/\text{facE\_al}` and the reference
speed is the birth speed

.. math::

   v_0 = \sqrt{2 E_0 / m}.

The physical path-length time is :math:`\tau = v_0 t`. The symplectic
integrators use :math:`\bar\tau = \tau/\sqrt{2}`, with velocities normalised
to :math:`v_0/\sqrt{2}`, energies to
:math:`E_0` and momenta divided by :math:`m v_0/\sqrt{2}`, which introduces the
reference gyroradius :math:`\rho_0 = m c v_0 / (\sqrt{2} e)` (unit G cm) in
the normalised momenta, :math:`\bar p_k = \bar v_\parallel h_k + A_k/\rho_0`.
The output state uses :math:`v/v_0` and the pitch
:math:`\lambda = v_\parallel / v`.

The time step is set by ``npoiper2``: the step in :math:`\tau` is
:math:`2\pi R_0 / \text{npoiper2}`, with the major radius :math:`R_0` of the
equilibrium. A macrostep, after which the results are recorded, contains an
integer number of such microsteps.

Invariants and the trapping parameter
-------------------------------------

The perpendicular invariant reported in the outputs is
:math:`J_\perp \propto (1 - \lambda^2) v^2 / B` at the start point. The
trapping parameter

.. math::

   \kappa = \left((1 - \lambda^2) \frac{B_\text{max}}{B} - 1\right)
            \frac{B_\text{min}}{B_\text{max} - B_\text{min}}

uses the extrema :math:`B_\text{min}`, :math:`B_\text{max}` of the field
strength on the start surface. :math:`\kappa < 0` is passing,
:math:`0 \le \kappa \le 1` trapped.

Collisions
----------

With ``swcoll = .true.``, Coulomb collisions with a Maxwellian background of
electrons and two ion species are added by operator splitting: after each
integration step of length :math:`\Delta\tau`, the speed and the pitch are
changed by a Monte Carlo step (``src/collis_alphas.f90``, ``stost``):

.. math::

   \Delta\lambda = \sqrt{2 D_{\lambda} (1 - \lambda^2) \Delta\tau}\, \xi_1
                   - 2 \lambda D_{\lambda} \Delta\tau ,
   \qquad
   \Delta p = \sqrt{2 D_{p} \Delta\tau}\, \xi_2 + F_p \Delta\tau ,

with independent random numbers :math:`\xi_{1,2}` of zero mean and unit
variance (:math:`\pm 1` for the pitch, uniformly distributed for the
speed), the pitch-angle scattering coefficient :math:`D_\lambda`, the energy
diffusion coefficient :math:`D_p` and the drag :math:`F_p` for the normalised
momentum :math:`p = v/v_0`. The coefficients are computed for the local
background densities and temperatures from the ``&profiles`` namelist. Drag
on electrons and ions slows the alphas down; pitch-angle scattering moves them
between trapped and passing orbits.

References
----------

* R. G. Littlejohn, Variational principles of guiding centre motion,
  J. Plasma Phys. 29, 111 (1983), https://doi.org/10.1017/S002237780000060X
* C. G. Albert, S. V. Kasilov, W. Kernbichler, Symplectic integration with
  non-canonical quadrature for guiding-center orbits in magnetic confinement
  devices, J. Comput. Phys. 403, 109065 (2020),
  https://doi.org/10.1016/j.jcp.2019.109065
