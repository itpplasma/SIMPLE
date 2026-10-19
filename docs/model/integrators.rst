Integrators
===========

This page states the discrete equations each ``integmode`` solves, the step
size, how a loss crossing is located inside a step, and what SIMPLE does
when an implicit step fails. Usage: :doc:`/user-guide/running`.

Step size
---------

With :math:`\Delta\tau_\text{macro} = v_0\,\texttt{trace\_time}/(\texttt{ntimstep}-1)`
the microstep is

.. math::

   \Delta\tau = \frac{\Delta\tau_\text{macro}}{N_\tau}, \qquad
   N_\tau = \Bigl\lceil \frac{\Delta\tau_\text{macro}\,\texttt{npoiper2}}{2\pi R_0} \Bigr\rceil ,

so :math:`\Delta\tau \le 2\pi R_0/\texttt{npoiper2}` and every macrostep holds
an integer number of microsteps (``params_init`` in ``src/params.f90``;
:math:`R_0` from the VMEC or chartmap geometry). The symplectic schemes take
:math:`h = \Delta\tau/\sqrt2` in the time :math:`\bar t` of
:doc:`guiding-centre`. Losses are resolved within a microstep; confined
fractions are recorded per macrostep.

Symplectic schemes
------------------

All schemes advance :math:`z = (r,\vartheta,\varphi,\bar p_\varphi)` with the
canonical equations :eq:`hamilton`, written with
:math:`H' \equiv \partial_r\bar H/\partial_r\bar p_\vartheta` and the
normalised quantities of :doc:`guiding-centre` (bars omitted below). The
implicit equations are multiplied by :math:`\partial_r p_\vartheta` to avoid
the division and are solved by Newton's method with analytic Jacobians from
the second field derivatives.

.. csv-table::
   :header: "Scheme (``integmode``)", "Order", "Newton unknowns", "Field evaluated at"
   :widths: 22 10 32 36

   "Euler 1 (``1``)", "1", ":math:`r^{n+1}, p_\varphi^{n+1}`", ":math:`(r^{n+1},\vartheta^n,\varphi^n,p_\varphi^{n+1})`"
   "Euler 2 (``2``)", "1", ":math:`r^{n+1},\vartheta^{n+1},\varphi^{n+1}`", ":math:`(r^{n+1},\vartheta^{n+1},\varphi^{n+1},p_\varphi^n)`"
   "Midpoint (``3``)", "2", ":math:`z^{n+1}` and :math:`r^{n+1/2}` (5)", "averaged angles and :math:`p_\varphi`, radius :math:`r^{n+1/2}`"
   "Gauss–Legendre (``4``–``7``)", ":math:`2s`", "all :math:`s` stage states (:math:`4s`)", "the :math:`s = 1..4` stages"
   "Lobatto IIIA–IIIB (``15``)", "4", ":math:`4s-2 = 10`", "the :math:`s = 3` stages"

**Euler 1** (default) is explicit in the angles and implicit in the momenta.
``sympl_euler1_residual`` in ``src/orbit_symplectic_euler1.f90`` solves

.. math::

   \partial_r p_\vartheta\,(p_\vartheta - p_\vartheta^n)
     + h\,(\partial_\vartheta H\,\partial_r p_\vartheta - \partial_r H\,\partial_\vartheta p_\vartheta) = 0,
   \qquad
   \partial_r p_\vartheta\,(p_\varphi^{n+1} - p_\varphi^n)
     + h\,(\partial_\varphi H\,\partial_r p_\vartheta - \partial_r H\,\partial_\varphi p_\vartheta) = 0,

with all field quantities at :math:`(r^{n+1},\vartheta^n,\varphi^n,p_\varphi^{n+1})`,
then updates the angles explicitly,
:math:`\vartheta^{n+1} = \vartheta^n + hH'`,
:math:`\varphi^{n+1} = \varphi^n + h(v_\parallel - H' h_\vartheta)/h_\varphi`.
After convergence the field values are extrapolated to first order in the
last Newton correction instead of being re-evaluated. **Euler 2** is the
adjoint: :math:`p_\vartheta(r^{n+1},\vartheta^{n+1},\varphi^{n+1},p_\varphi^n)
= p_\vartheta^n` fixes :math:`r^{n+1}`, the angle equations are implicit, and
:math:`p_\varphi` is updated explicitly. **Midpoint** evaluates the field
once at :math:`(r^{n+1/2},\bar\vartheta,\bar\varphi,\bar p_\varphi)` with
averaged angles and momentum; the extra unknown :math:`r^{n+1/2}` satisfies
:math:`p_\vartheta(z^{n+1/2}) = p_\vartheta^n + \tfrac h2\dot p_\vartheta`,
which together with the full-step equation makes :math:`p_\vartheta^{n+1/2}`
the average of both ends (``f_midpoint_part1``/``part2`` in
``src/orbit_symplectic.f90``). **Gauss–Legendre** solves the stage equations
:math:`z_k = z^n + h\sum_l a_{kl}\dot z_l` for
:math:`(p_\vartheta,\vartheta,\varphi,p_\varphi)` and updates them with the
weights :math:`b_l`; the stored radius after the step is the last stage
radius, and :math:`r^{n+1}` is fixed implicitly by :math:`p_\vartheta^{n+1}`
at the next step (``gauss_canfield_residual`` in ``src/orbit_rk_core.f90``).

Newton stops when all residuals are below :math:`10^{-15}` or all
corrections below ``relerr`` times a reference scale (Euler 1: 1 for
:math:`r`, :math:`10\,\psi_\text{tor}^\text{edge}/\bar\rho_0` for
:math:`p_\varphi`), after at most 32 iterations. An outward correction past
the loss boundary :math:`r_\text{max}` is cut to 0.8 of the remaining
distance; a step ending at :math:`r < 0` is mapped through the axis,
:math:`(r,\vartheta)\to(-r,\vartheta+\pi)`. The schemes conserve the
canonical two-form to the Newton tolerance, so backward-error analysis
bounds the long-time energy error for small steps in smooth regions (Hairer,
Lubich and Wanner 2006); with early classification, loss computations cost
about one order of magnitude less than with adaptive Runge–Kutta (Albert,
Kasilov and Kernbichler 2020, JPP).

Cartesian canonical chart (``canonical_chart = 'cartesian'``)
-------------------------------------------------------------

The polar schemes above are singular at the magnetic axis, where
:math:`(r,\vartheta)` is not a chart: steps that pass the axis closer than
about one step length have no solution in :math:`(r,\vartheta)`, and the
axis flip and the recovery path take over. With ``canonical_chart =
'cartesian'`` midpoint and Gauss (``integmode`` 3–7) step instead in the
global Darboux chart

.. math::

   X = \sqrt{2\kappa}\,x, \qquad Y = \sigma\sqrt{2\kappa}\,y, \qquad
   \kappa = \frac{\sigma p_\vartheta}{s}, \qquad
   (x, y) = \sqrt s\,(\cos\vartheta, \sin\vartheta),

with :math:`\sigma = \operatorname{sign}(p_\vartheta/s)`. Then
:math:`X^2 + Y^2 = 2\sigma p_\vartheta` and :math:`dX\wedge dY = dp_\vartheta
\wedge d\vartheta`, so :math:`(X, Y, \varphi, p_\varphi)` are canonical and
one chart covers the axis :math:`X = Y = 0`. For a field that is regular in
:math:`(x, y)` (the libneo regular Boozer field of :doc:`coordinates`, or the
analytic test field) :math:`\kappa` is regular and positive, and so is the
Hamiltonian :math:`K(X, Y, \varphi, p_\varphi)`. Midpoint is Gauss-1; for
:math:`s` stages the unknowns are the stage points :math:`W_i` and their
field points :math:`u_i = (x_i, y_i)`:

.. math::

   W_i = w^n + h\sum_j a_{ij}\,\mathbb J\nabla K(W_j), \qquad
   C(u_i, \varphi_i, p_{\varphi,i}) = (X_i, Y_i),

where :math:`C` is the chart map above and :math:`\nabla K` follows from
:math:`\nabla H` at :math:`u_i` by implicit differentiation of :math:`C`.
This is the Gauss method for :math:`\dot w = \mathbb J\nabla K(w)` and
therefore symplectic, with no radius floor, axis flip, chart handoff or
recovery: a step either converges or ends the marker with exit code 102 or
103. Newton's first iteration uses the second field derivatives; the chord
iterations reuse that Jacobian with first derivatives. The end point
:math:`w^{n+1} = w^n + \sum_i d_i (W_i - w^n)`, :math:`d = b^\top A^{-1}`, is
mapped back to :math:`(s, \vartheta)` for output, loss and wall checks.
A collision kick changes :math:`p_\varphi` and :math:`\mu` at fixed
:math:`(s, \vartheta, \varphi)`; the next step re-encodes :math:`(X, Y)`.
The pointwise algebra (:math:`\kappa`, :math:`H`, :math:`C`, the vector field
and the stage Jacobian blocks) is generated with FortSym by
``gen/gen_darboux_chart.f90`` into ``src/generated/``; the step is
``src/orbit_symplectic_darboux.f90``. SPECTRE, canonical flux, Meiss and
Albert fields and the GPU paths are not supported.

Adaptive Runge–Kutta (``integmode = 0``)
----------------------------------------

The Cash–Karp RK5(4) integrator of libneo (``odeint_allroutines``) solves the
non-canonical guiding-centre equations in any coordinates, with relative
tolerance ``relerr`` and :math:`w = v/v_0` (``velo_can`` in
``src/sub_alpha_lifetime_can.f90``):

.. math::

   \frac{dx^i}{d\tau} = \frac{w\lambda\,h^{*i}
      + \dfrac{\rho_0 w^2(1-\lambda^2)}{2B\sqrt g}\,\varepsilon^{ijk}h_j\partial_k\ln B}
      {h^*_\parallel},
   \qquad
   \frac{d\lambda}{d\tau} = -\frac{(1-\lambda^2)\,w}{2h^*_\parallel}\,
      \mathbf h^*\!\cdot\nabla\ln B,
   \qquad \frac{dw}{d\tau} = 0,

with :math:`\mathbf h^* = \mathbf h + w\lambda(\rho_0/B)\nabla\times\mathbf h`
and :math:`h^*_\parallel = 1 + w\lambda(\rho_0/B)\,\mathbf h\cdot\nabla\times
\mathbf h`. It is not symplectic. Below :math:`s = 0.01` the step switches to
:math:`(\sqrt s\cos\vartheta, \sqrt s\sin\vartheta, \varphi)`, which is regular
at the axis (``orbit_timestep_axis``). The same equations serve as the first
recovery method below.

Loss crossing inside a step
---------------------------

If a symplectic trial step ends beyond :math:`r_\text{max}` (LCFS, outer
chartmap boundary, or :math:`M_\text{vol}` for SPECTRE),
``locate_symplectic_boundary`` repeats it from the accepted state with a
fraction :math:`f` of :math:`h`, chosen by a secant estimate of
:math:`r(f) = r_\text{max}` damped to 0.8 of the way, else by bisection;
trials outside the domain become the upper bracket. It converges when the
bracket width is at most ``boundary_event_fraction_tolerance`` and
:math:`|r - r_\text{max}|` at most ``boundary_event_radial_tolerance``
(defaults :math:`\max(10^{-12}, 10\,\texttt{relerr})`,
:math:`\max(10^{-10}, 10\,\texttt{relerr})`), within 64 trials. The loss
time is :math:`(n+f)\Delta\tau`; residual and bracket width go to
``orbit_exit_code.dat``.

Failure recovery
----------------

With ``symplectic_newton_warning_mode = .true.`` (default) a failed
microstep, typically near the axis, goes through these stages
(``src/orbit_symplectic.f90``, ``src/simple_main.f90``):

1. At the iteration limit, accept a final correction of at most 100
   (SPECTRE: 10) times the convergence threshold.
2. Retry as two half steps, recursively to depth 8.
3. RK equations above from the last accepted state, at most 10,000 steps;
   then, at reference radius :math:`< 0.1` with :math:`|B^*_\parallel/B| \ge
   0.2`, with axis floor :math:`10^{-3}`; then a toroidally regularised RK
   step (at most 100,000 steps); then, below reference radius 0.02 with a
   chartmap field, a Boris full-orbit bridge. Each must keep :math:`v/v_0`
   within 0.5 %.
4. Resume the symplectic scheme (``reseed_sympl``) after 256 consecutive RK
   steps with :math:`|B^*_\parallel/B| \ge 0.3` and reference radius
   :math:`\ge 0.12` (recovery started below 0.1) or :math:`\le 0.98`
   (started above 0.98), otherwise after :math:`2^{7+k}` such steps for the
   :math:`k`-th recovery (:math:`k\le9`).
5. Otherwise roll back and skip the microstep; a second consecutive failure
   ends the marker with exit code 101–105, or 4 if collisionless at
   :math:`s < 0.01` (:doc:`/user-guide/outputs`).

With the warning mode off, the first failure ends the marker. Stages 2–5
are not a fixed symplectic map.

References
----------

* C. G. Albert, S. V. Kasilov, W. Kernbichler, Symplectic integration with
  non-canonical quadrature for guiding-center orbits in magnetic confinement
  devices, J. Comput. Phys. 403, 109065 (2020),
  https://doi.org/10.1016/j.jcp.2019.109065
* C. G. Albert, S. V. Kasilov, W. Kernbichler, Accelerated methods for direct
  computation of fusion alpha particle losses within stellarator
  optimization, J. Plasma Phys. 86, 815860201 (2020),
  https://doi.org/10.1017/S0022377820000203
* E. Hairer, C. Lubich, G. Wanner, Geometric Numerical Integration,
  2nd ed., Springer (2006), https://doi.org/10.1007/3-540-30666-8
