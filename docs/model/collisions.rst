Collisions
==========

This page derives the Monte Carlo Coulomb collision operator applied with
``swcoll = .true.``: coefficients, discrete kicks, substep bounds and random
streams. On this page only, :math:`p = v/v_0`.
Usage: :doc:`/user-guide/running`.

Scope and splitting
-------------------

A test alpha particle scatters on a Maxwellian background of electrons and
two ion species with densities :math:`n_a` and temperatures :math:`T_a` from
the ``&profiles`` namelist. The operator always uses the alpha mass
:math:`m_\alpha = 4m_p` and charge :math:`Z_\alpha = 2`
(``src/collis_alphas.f90:107,126``): ``n_e`` does not enter, and a run with
``n_d`` :math:`\neq 4` stops in ``init_collisions`` with ``v0_coll != v0``.
Collisions cannot be combined with classification.

Each orbit microstep :math:`\Delta\tau` is followed by a collision step of
the same length at fixed position (``collide`` in ``src/simple_main.f90``)
that changes :math:`p` and :math:`\lambda`; ``update_momentum`` re-encodes
:math:`J_\perp`, :math:`\bar v_\parallel`, :math:`\bar p_\varphi`
(:doc:`guiding-centre`). A failed collision step rejects the microstep as a
numerical failure (:doc:`integrators`).

Coefficients
------------

For species :math:`a \in \{1, 2, e\}` let :math:`v_{T,a} = \sqrt{2T_a/m_a}`,
:math:`x_a = p\,v_0/v_{T,a}`, :math:`n_e = Z_1n_1 + Z_2n_2`, :math:`Z_e = 1`,
:math:`E_\alpha = m_\alpha v_0^2/2 = 3.5\,\text{MeV}/\texttt{facE\_al}`, and

.. math::

   \nu_a = \frac{2\pi Z_\alpha^2 Z_a^2 e^4 n_a \ln\Lambda_a}{m_\alpha^2 v_0^3\, v_{T,a}}
   \quad [\text{cm}^{-1}].

With the Chandrasekhar function
:math:`G(x) = [\operatorname{erf}x - x\operatorname{erf}'x]/(2x^2)`,
``onseff`` evaluates (series for :math:`x < 0.01`, asymptotes for
:math:`x > 6`)

.. math::

   d_p = \frac{2G(x)}{x}, \qquad
   d_\lambda = \frac{\operatorname{erf}x - G(x)}{x}, \qquad
   d_F = \frac1x\frac{d(x^2 d_p)}{dx} = \frac{4e^{-x^2}}{\sqrt\pi} - d_p ,

and ``coleff`` sums the species:

.. math::

   D_p = \sum_a \nu_a d_p(x_a), \qquad
   D_\lambda = \frac{1}{p^2}\sum_a \nu_a d_\lambda(x_a), \qquad
   F_p = \sum_a \nu_a\Bigl[\frac{d_F(x_a)}{p} - 2p\frac{E_\alpha}{T_a}d_p(x_a)\Bigr]
       = \frac{1}{p^2}\frac{\partial(p^2 D_p)}{\partial p}
         - 2p\sum_a\frac{E_\alpha}{T_a}\nu_a d_p(x_a).

:math:`F_p` and :math:`D_p` are the Itô drift and diffusion of :math:`p`, and
:math:`D_\lambda` the pitch diffusion, of the Fokker–Planck equation for
:math:`f(p,\lambda)` per unit velocity volume

.. math::

   \frac{\partial f}{\partial\tau} = \sum_a \frac{1}{p^2}\frac{\partial}{\partial p}
     \Bigl[p^2\nu_a d_p(x_a)\Bigl(\frac{\partial f}{\partial p}
       + 2p\frac{E_\alpha}{T_a}f\Bigr)\Bigr]
   + D_\lambda\frac{\partial}{\partial\lambda}\Bigl[(1-\lambda^2)\frac{\partial f}{\partial\lambda}\Bigr],

whose speed flux vanishes for each species on the Maxwellian
:math:`f\propto\exp(-m_\alpha v^2/2T_a)`.

The Coulomb logarithms (``compute_collision_coeffs``; :math:`n` in
cm\ :sup:`-3`, :math:`T`, :math:`E_\alpha` in eV, :math:`A_i` = ``am1``,
``am2``) are the NRL electron–ion form and the NRL mixed ion–ion form with
:math:`m_\alpha = 4m_p`, :math:`T_\alpha \to E_\alpha` and only the
background-ion density under the root:

.. math::

   \ln\Lambda_e = 24 - \ln\frac{\sqrt{n_e}}{T_e}, \qquad
   \ln\Lambda_i = 23 - \ln\Bigl[\frac{Z_\alpha Z_i(4 + A_i)}{4T_i + A_iE_\alpha}
                  \sqrt{\frac{n_iZ_i^2}{T_i}}\Bigr].

``init_collision_profiles`` tabulates :math:`\nu_a`, :math:`v_0/v_{T,a}` and
:math:`E_\alpha/T_a` on 101 equidistant points in :math:`[0,1]`, with floors
:math:`T_a \ge 1` eV and :math:`n_a \ge 10^6` m\ :sup:`-3` (both ion
densities zero is vacuum). Each step interpolates linearly at the clamped
integration radius :math:`r`, which is the normalised toroidal flux
:math:`s` for ``integ_coords`` 0 and 2 only.

Discrete step
-------------

``stost`` (mode 1) advances :math:`(p,\lambda)` over :math:`\Delta\tau` in
substeps :math:`\delta`, with coefficients at the current :math:`p`:

.. math::

   \lambda' = \lambda\cos\alpha + \sqrt{1-\lambda^2}\,\sin\alpha\cos\psi,
   \quad \alpha = \sqrt{4D_\lambda\delta},\ \psi\sim\mathcal U[0,2\pi);
   \qquad
   p' = p + F_p\delta + \sqrt{2D_p\delta}\,\xi,\quad \xi\sim\mathcal U[-\sqrt3,\sqrt3].

The pitch rotation has :math:`\langle\Delta\lambda\rangle =
-2\lambda D_\lambda\delta` and :math:`\langle\Delta\lambda^2\rangle =
2(1-\lambda^2)D_\lambda\delta` to leading order, the generator of the
Lorentz operator, and keeps :math:`|\lambda| \le 1` exactly; :math:`\xi` has
zero mean and unit variance. Before the random draws,

.. math::

   \delta = \min\Bigl(\text{remaining},\ \frac{f^2}{D_\lambda},\
            \frac{f\,p}{|F_p|},\ \frac{f^2p^2}{6D_p}\Bigr), \qquad f = 0.1,

which bounds the rotation angle by :math:`2f` and the drift and the largest
noise amplitude :math:`\sqrt{6D_p\delta}` by :math:`fp`. ``stost`` accepts
:math:`0 < f \le 0.25` (``max_substep_fraction``); ``simple.x`` uses 0.1. A
non-finite or negative coefficient, :math:`p \le 0`, a substep that no longer
advances time, or more than 100,000 substeps restores the initial state and
returns an error. In vacuum neither the state nor the generator changes.

Random streams
--------------

With ``deterministic = .true.``, ``reset_seed_if_deterministic``
(``src/params.f90``) reseeds the generator before each marker from
``ran_seed`` and :math:`k` = ``particle_index_offset`` + source row:

.. math::

   s_0 = (\texttt{ran\_seed} + 104729\,k) \bmod (2^{31}-2) + 1, \qquad
   s_i = (48271\,s_{i-1} + 7919\,i) \bmod (2^{31}-2) + 1 .

A marker thus draws the same kicks for any thread count or contiguous
shard; distinct seeds do not prove statistical independence.

References
----------

* J. D. Huba, NRL Plasma Formulary, NRL/PU/6790--19-640, Naval Research
  Laboratory (2019)
