Symplectic integrators
======================

Why symplectic
--------------

Collisionless alpha losses in stellarators are dominated by trapped orbits
that drift slowly over thousands of bounce periods, and by the boundary
between regular and chaotic orbits. A numerical method must therefore keep
the long-time structure of the guiding-centre motion. Symplectic integrators
conserve the symplectic form of the Hamiltonian system exactly (up to the
tolerance of the implicit solve). As a consequence, phase-space volume is
conserved, the energy error stays bounded without secular drift, and
invariant tori of regular orbits are preserved at step sizes much larger than
those of standard Runge-Kutta methods of the same accuracy. Albert, Kasilov
and Kernbichler (2020, J. Comput. Phys.) introduced these integrators for
guiding-centre orbits; together with early orbit classification they reduce
the cost of loss computations by about one order of magnitude compared with
adaptive Runge-Kutta integration (Albert, Kasilov and Kernbichler 2020,
J. Plasma Phys.).

Formulation
-----------

All integrators advance the four variables :math:`z = (r, \vartheta,
\varphi, p_\varphi)` of :doc:`guiding-centre`. Let
:math:`q = (\vartheta, \varphi)` and :math:`p = (p_\vartheta, p_\varphi)` be
the canonical pairs and :math:`h` the step in normalised time. The schemes
are standard symplectic methods for :math:`(q, p)`. Their implicit equations
are written in :math:`z`, where :math:`p_\vartheta(z)` and :math:`H(z)` are
evaluated from the splined field (non-canonical quadrature points). The
derivatives at constant canonical variables follow from the implicit function
theorem (:doc:`guiding-centre`, Eq. :eq:`hamilton`). The resulting nonlinear
equations are solved by Newton's method with analytic Jacobians from the
second derivatives of the field splines, to the relative tolerance
``relerr``.

Symplectic Euler (``integmode = 1``, ``2``)
   Variant 1 (explicit in :math:`q`, implicit in :math:`p`):

   .. math::

      p^{n+1} = p^n - h\, \partial_q H(q^n, p^{n+1}), \qquad
      q^{n+1} = q^n + h\, \partial_p H(q^n, p^{n+1}).

   The unknowns of the Newton iteration are :math:`r^{n+1}` and
   :math:`p_\varphi^{n+1}` (``src/orbit_symplectic_euler1.f90``). Variant 2
   is the adjoint method, implicit in :math:`q` and explicit in :math:`p`.
   Both are first order. Variant 1 is the default: it needs a
   two-dimensional Newton solve and one field evaluation per iteration.

Implicit midpoint (``integmode = 3``)
   .. math::

      p^{n+1} = p^n - h\, \partial_q H(z^{n+1/2}), \qquad
      q^{n+1} = q^n + h\, \partial_p H(z^{n+1/2}),

   with the field evaluated at the midpoint of the angles and of
   :math:`p_\varphi`. The midpoint radius is an additional unknown. Second
   order and symmetric.

Gauss-Legendre collocation (``integmode = 4`` to ``7``)
   Implicit Runge-Kutta methods with 1 to 4 stages at the Gauss-Legendre
   nodes, of order 2, 4, 6 and 8. They are symplectic and symmetric. All
   stages are solved together in one Newton iteration; the stage field values
   are extrapolated to first order in the Newton correction
   (``DOC/coordinates-and-fields.md``, section 8.3).

Lobatto IIIA-IIIB (``integmode = 15``)
   Partitioned Runge-Kutta method with 3 stages, order 4, symplectic for the
   partitioned system.

Adaptive Runge-Kutta (``integmode = 0``)
   The Cash-Karp RK5(4) method with adaptive step size (libneo
   ``odeint_allroutines``) on the non-canonical guiding-centre equations of
   motion, in any coordinates. It is not symplectic. It serves as reference,
   as fallback, and for fields without canonical coordinates
   (``integ_coords`` 1 and 5).

Step size and output times
--------------------------

The microstep is :math:`h = 2\pi R_0 / \text{npoiper2}` in path length, 256
steps per :math:`2\pi R_0` by default. The ``ntimstep`` macrosteps
between 0 and ``trace_time`` contain an integer number of microsteps each;
losses are recorded with the resolution of a microstep.

Boundary events
---------------

When a step crosses the last closed flux surface (:math:`s = 1` or the outer
boundary of a chartmap), SIMPLE locates the crossing within the step by
bracketing the step fraction and bisecting the bracket.
``boundary_event_fraction_tolerance`` and ``boundary_event_radial_tolerance``
control the accuracy; by default they are
derived from ``relerr``. The loss time and position are those of the crossing.
``orbit_exit_code.dat`` reports the radial residual and the time width of the
final bracket.

Failure handling
----------------

Close to the magnetic axis, the radial coordinate of flux coordinates is
singular and :math:`\partial_r p_\vartheta` becomes small, so the Newton
iteration can fail. In the default warning mode
(``symplectic_newton_warning_mode = .true.``) SIMPLE

1. accepts a bounded final Newton correction when the iteration limit is
   reached,
2. retries a rejected step with dyadically smaller symplectic steps,
3. falls back to the adaptive Runge-Kutta method from the last accepted state
   for at most 10,000 internal steps,
4. as a last resort, skips one isolated microstep.

Orbits for which all methods fail end with a numerical exit code (101 to 105)
and are reported as unresolved, not as lost. ``DOC/config.md``
(:doc:`/reference/design-notes`) gives the exact rules, including the
treatment of steps that cross the polar axis.

GPU integrators
---------------

The CUDA backend integrates the guiding-centre equations in Boozer
coordinates with adaptive Dormand-Prince or Cash-Karp Runge-Kutta methods on
the GPU (:doc:`/user-guide/gpu`). It is not symplectic and is designed for
large ensembles over short to moderate tracing times, in particular for the
segmented loss objective.

References
----------

* C. G. Albert, S. V. Kasilov, W. Kernbichler, Symplectic integration with
  non-canonical quadrature for guiding-center orbits in magnetic confinement
  devices, J. Comput. Phys. 403, 109065 (2020),
  https://doi.org/10.1016/j.jcp.2019.109065
* E. Hairer, C. Lubich, G. Wanner, Geometric Numerical Integration,
  2nd ed., Springer (2006), https://doi.org/10.1007/3-540-30666-8
